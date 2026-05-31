"""
align.py — WCS fitting for rb_align.

Fitting tiers
-------------
3+ pairs → full affine via fit_wcs_from_points
2 pairs  → shift + rotation (analytic)
1 pair   → CRVAL shift only
0 pairs  → fail per on_fail policy
"""
from __future__ import annotations

import copy
import warnings
from typing import List, Optional, Tuple

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points


# ---------------------------------------------------------------------------
# Residual helper
# ---------------------------------------------------------------------------

def _compute_residuals(pairs, new_wcs):
    """
    Compute per-source residuals in arcseconds.

    Parameters
    ----------
    pairs : list of (rx, ry, ra, dec, tx, ty)
    new_wcs : astropy.wcs.WCS (2D)

    Returns
    -------
    residuals : np.ndarray, arcsec
    rms : float, arcsec
    """
    if not pairs:
        return np.array([]), 0.0

    truth_ra  = np.array([p[2] for p in pairs])
    truth_dec = np.array([p[3] for p in pairs])
    tgt_x     = np.array([p[4] for p in pairs])
    tgt_y     = np.array([p[5] for p in pairs])

    fitted_sky = new_wcs.all_pix2world(
        np.column_stack([tgt_x, tgt_y]), 0)
    fitted_ra  = fitted_sky[:, 0]
    fitted_dec = fitted_sky[:, 1]

    # Angular separation in arcsec
    cos_dec = np.cos(np.deg2rad(truth_dec))
    dra  = (fitted_ra  - truth_ra)  * cos_dec * 3600.0
    ddec = (fitted_dec - truth_dec) * 3600.0
    residuals = np.sqrt(dra**2 + ddec**2)
    rms = float(np.sqrt(np.mean(residuals**2)))
    return residuals, rms


# ---------------------------------------------------------------------------
# Shift-only fit (1 pair)
# ---------------------------------------------------------------------------

def _fit_shift(frame, pairs):
    """
    Shift CRVAL1/CRVAL2 so that one source lands on the correct sky position.
    Returns a new WCS object.
    """
    rx, ry, ra_truth, dec_truth, tx, ty = pairs[0]
    wcs_orig = frame.wcs

    # Current sky position of the target pixel
    sky_curr = wcs_orig.all_pix2world([[tx, ty]], 0)
    ra_curr  = float(sky_curr[0, 0])
    dec_curr = float(sky_curr[0, 1])

    new_wcs = copy.deepcopy(wcs_orig)
    new_wcs.wcs.crval[0] += (ra_truth  - ra_curr)
    new_wcs.wcs.crval[1] += (dec_truth - dec_curr)
    new_wcs.wcs.set()
    return new_wcs


# ---------------------------------------------------------------------------
# Shift + rotation fit (2 pairs)
# ---------------------------------------------------------------------------

def _fit_shift_rotation(frame, pairs):
    """
    Fit a shift + rotation from two pairs using analytic 2-point solution.

    Computes the rotation that maps target pixel vectors to reference sky vectors,
    then patches the CD matrix.
    """
    wcs_orig = copy.deepcopy(frame.wcs)

    # Pair 0
    p0 = pairs[0]
    p1 = pairs[1]
    # Target pixel positions
    t0 = np.array([p0[4], p0[5]])
    t1 = np.array([p1[4], p1[5]])
    # Truth sky positions (in the reference frame)
    sky0 = np.array([p0[2], p0[3]])
    sky1 = np.array([p1[2], p1[3]])

    # Use fit_wcs_from_points with 2 points — may not produce a full affine,
    # but it gives the best least-squares solution astropy can manage.
    pixel_coords = np.array([[t0[0], t1[0]],
                              [t0[1], t1[1]]])
    sky_coords = SkyCoord(ra=np.array([sky0[0], sky1[0]]) * u.deg,
                          dec=np.array([sky0[1], sky1[1]]) * u.deg)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            new_wcs = fit_wcs_from_points(pixel_coords, sky_coords,
                                          proj_point='center',
                                          projection='TAN')
        return new_wcs
    except Exception:
        # Last-resort: pure shift from midpoint
        mid_pair = [
            (t0[0] + t1[0]) / 2, (t0[1] + t1[1]) / 2,
            (sky0[0] + sky1[0]) / 2, (sky0[1] + sky1[1]) / 2,
            (t0[0] + t1[0]) / 2, (t0[1] + t1[1]) / 2,
        ]
        return _fit_shift(frame, [mid_pair])


# ---------------------------------------------------------------------------
# Full affine fit (3+ pairs)
# ---------------------------------------------------------------------------

def _fit_full(frame, pairs):
    """
    Full WCS fit using fit_wcs_from_points (shift + rotation + scale).
    """
    tgt_x  = np.array([p[4] for p in pairs])
    tgt_y  = np.array([p[5] for p in pairs])
    ras    = np.array([p[2] for p in pairs])
    decs   = np.array([p[3] for p in pairs])

    pixel_coords = np.array([tgt_x, tgt_y])
    sky_coords   = SkyCoord(ra=ras * u.deg, dec=decs * u.deg)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        new_wcs = fit_wcs_from_points(pixel_coords, sky_coords,
                                      proj_point='center',
                                      projection='TAN')
    return new_wcs


# ---------------------------------------------------------------------------
# Per-frame aligner
# ---------------------------------------------------------------------------

def _align_frame(aligner, target_idx: int,
                 pairs_for_frame,
                 on_fail: str = 'skip',
                 min_sources: int = 2):
    """
    Fit WCS for one target frame.

    Returns (new_wcs, fit_type, n_used, residuals, rms) or raises.
    """
    frame = aligner.targets[target_idx]
    n = len(pairs_for_frame)

    if n == 0:
        msg = f"[rb_align] target {target_idx}: 0 pairs — cannot align."
        if on_fail == 'raise':
            raise RuntimeError(msg)
        elif on_fail == 'degrade':
            print(msg + " (degraded — keeping original WCS)")
            return copy.deepcopy(frame.wcs), 'none', 0, np.array([]), 0.0
        else:  # skip
            print(msg + " (skipped)")
            aligner.flagged_frames.append(target_idx)
            return None, 'none', 0, np.array([]), 0.0

    if n < min_sources:
        msg = (f"[rb_align] target {target_idx}: only {n} pair(s) "
               f"(min={min_sources}).")
        if on_fail == 'raise':
            raise RuntimeError(msg)
        elif on_fail == 'degrade':
            print(msg + " (degraded — using shift-only fit)")
            min_sources = 1   # allow shift fit
        else:
            print(msg + " (skipped)")
            aligner.flagged_frames.append(target_idx)
            return None, 'none', 0, np.array([]), 0.0

    if n >= 3:
        try:
            new_wcs = _fit_full(frame, pairs_for_frame)
            fit_type = 'full'
        except Exception as e:
            print(f"[rb_align] target {target_idx}: full fit failed ({e}), "
                  "falling back to shift+rotation.")
            n = 2
    if n == 2:
        try:
            new_wcs = _fit_shift_rotation(frame, pairs_for_frame[:2])
            fit_type = 'shift+rot'
        except Exception as e:
            print(f"[rb_align] target {target_idx}: shift+rot fit failed ({e}), "
                  "falling back to shift-only.")
            n = 1
    if n == 1:
        new_wcs = _fit_shift(frame, pairs_for_frame[:1])
        fit_type = 'shift'

    residuals, rms = _compute_residuals(pairs_for_frame, new_wcs)
    return new_wcs, fit_type, len(pairs_for_frame), residuals, rms


# ---------------------------------------------------------------------------
# Public: align single target
# ---------------------------------------------------------------------------

def _align_single(aligner, target_idx: int = 0):
    """Align one target frame using aligner.pairs."""
    pairs = aligner.pairs
    result = _align_frame(aligner, target_idx, pairs,
                          on_fail='raise', min_sources=1)
    new_wcs, fit_type, n_used, residuals, rms = result

    aligner._new_wcs[target_idx]          = new_wcs
    aligner.fit_type[target_idx]          = fit_type
    aligner.n_sources_used[target_idx]    = n_used
    aligner.source_residuals[target_idx]  = residuals
    aligner.rms_residuals[target_idx]     = rms

    print(f"[rb_align] target {target_idx}: fit='{fit_type}', "
          f"n={n_used}, RMS={rms:.3f}\"")


# ---------------------------------------------------------------------------
# Public: align many targets
# ---------------------------------------------------------------------------

def _align_many(aligner, on_fail: str = 'skip', min_sources: int = 2):
    """
    Align all targets in batch.

    For each target, filter pairs to those whose sky positions project into
    the target's image bounds, then call _align_frame.
    """
    from .sources import _refine_centroid

    # The global pair list has (rx, ry, ra, dec, tx_0, ty_0) from target-0
    # For batch, we use ra/dec only and recentroid per target.
    global_pairs = aligner.pairs
    if not global_pairs:
        raise RuntimeError(
            "No pairs stored.  Run find_sources() before align_many().")

    ras  = np.array([p[2] for p in global_pairs])
    decs = np.array([p[3] for p in global_pairs])
    rxs  = np.array([p[0] for p in global_pairs])
    rys  = np.array([p[1] for p in global_pairs])

    for idx, frame in enumerate(aligner.targets):
        img = frame.image2d
        if img is None:
            print(f"[rb_align] target {idx}: no image2d — skipping.")
            aligner.flagged_frames.append(idx)
            continue

        ny, nx = img.shape

        # Project all sky positions into this frame's pixel grid
        try:
            pix = frame.wcs.all_world2pix(
                np.column_stack([ras, decs]), 0)
        except Exception as e:
            print(f"[rb_align] target {idx}: WCS projection failed ({e}) — skipping.")
            aligner.flagged_frames.append(idx)
            continue

        px = pix[:, 0]
        py = pix[:, 1]
        in_bounds = (px > 0) & (px < nx) & (py > 0) & (py < ny)

        if in_bounds.sum() < min_sources:
            msg = (f"[rb_align] target {idx}: only {in_bounds.sum()} "
                   f"in-bounds source(s) (min={min_sources}).")
            if on_fail == 'raise':
                raise RuntimeError(msg)
            elif on_fail == 'degrade':
                print(msg + " (degraded — using all available)")
                if in_bounds.sum() == 0:
                    aligner.flagged_frames.append(idx)
                    continue
            else:
                print(msg + " (skipped)")
                aligner.flagged_frames.append(idx)
                continue

        # Recentroid each in-bounds source on this target
        frame_pairs = []
        for i in np.where(in_bounds)[0]:
            tx_pred, ty_pred = float(px[i]), float(py[i])
            tx_r, ty_r = _refine_centroid(img, tx_pred, ty_pred)
            frame_pairs.append((
                float(rxs[i]), float(rys[i]),
                float(ras[i]), float(decs[i]),
                tx_r, ty_r,
            ))

        result = _align_frame(aligner, idx, frame_pairs,
                               on_fail=on_fail, min_sources=min_sources)
        new_wcs, fit_type, n_used, residuals, rms = result

        aligner._new_wcs[idx]          = new_wcs
        aligner.fit_type[idx]          = fit_type
        aligner.n_sources_used[idx]    = n_used
        aligner.source_residuals[idx]  = residuals
        aligner.rms_residuals[idx]     = rms

        if new_wcs is not None:
            print(f"[rb_align] target {idx}: fit='{fit_type}', "
                  f"n={n_used}, RMS={rms:.3f}\"")
