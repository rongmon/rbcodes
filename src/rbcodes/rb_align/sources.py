"""
sources.py — source detection and matching strategies for rb_align.

Strategies
----------
interactive  — two-panel matplotlib state machine (the primary method)
gaia         — Gaia DR3 catalog (astroquery required)
dao          — DAOStarFinder (photutils optional, fallback to scipy)
knots        — segmentation map → bright knot centroids
cross_corr   — cross-correlation shift (scipy only)
batch        — reproject saved catalog onto each target and recentroid
auto         — cascade: gaia → dao → knots → cross_corr → interactive
"""
from __future__ import annotations

import warnings
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Centroiding — photutils optional
# ---------------------------------------------------------------------------

try:
    from photutils.centroids import centroid_com, centroid_2dg
    HAS_PHOTUTILS = True
except ImportError:
    from scipy.ndimage import center_of_mass as _com

    def centroid_com(cutout):          # same behaviour as photutils version
        cy, cx = _com(np.nan_to_num(np.abs(cutout), nan=0.0))
        return cx, cy

    def centroid_2dg(cutout):
        return centroid_com(cutout)

    HAS_PHOTUTILS = False


# ---------------------------------------------------------------------------
# Display normalisation helper
# ---------------------------------------------------------------------------

def _compute_norm(data2d: np.ndarray, stretch='zscale'):
    """
    Build a matplotlib Normalize object for the given 2D array.

    Parameters
    ----------
    data2d : np.ndarray
    stretch : str or (vmin, vmax)
        'zscale' | '99.5%' | '99%' | '98%' | '97%' | '95%' | 'minmax'
        or a (vmin, vmax) tuple/list for manual limits.
    """
    import matplotlib.colors as mcolors
    from astropy.visualization import ZScaleInterval

    finite = data2d[np.isfinite(data2d)]
    if len(finite) == 0:
        return mcolors.Normalize(vmin=0, vmax=1)

    if isinstance(stretch, (tuple, list)):
        vmin, vmax = float(stretch[0]), float(stretch[1])
    elif stretch == 'zscale':
        try:
            vmin, vmax = ZScaleInterval().get_limits(finite)
        except Exception:
            vmin, vmax = np.nanpercentile(finite, [1, 99])
    elif stretch == 'minmax':
        vmin, vmax = float(finite.min()), float(finite.max())
    else:
        pct = float(str(stretch).rstrip('%'))
        lo = 100.0 - pct
        vmin, vmax = np.nanpercentile(finite, [lo, pct])

    return mcolors.Normalize(vmin=vmin, vmax=vmax)


# ---------------------------------------------------------------------------
# Centroid refinement
# ---------------------------------------------------------------------------

def _refine_centroid(image2d: np.ndarray,
                     cx: float, cy: float,
                     box: int = 15):
    """
    Refine centroid within a box around (cx, cy).

    Returns (refined_x, refined_y) in image pixel coordinates.
    Falls back to input position on any error.
    """
    ny, nx = image2d.shape
    x0 = int(round(cx))
    y0 = int(round(cy))
    x1 = max(0, x0 - box)
    x2 = min(nx, x0 + box + 1)
    y1 = max(0, y0 - box)
    y2 = min(ny, y0 + box + 1)

    cutout = image2d[y1:y2, x1:x2].copy()
    if cutout.size == 0:
        return cx, cy

    # Replace NaNs/infs with local minimum
    bad = ~np.isfinite(cutout)
    if bad.all():
        return cx, cy
    cutout[bad] = np.nanmin(cutout[~bad])

    # Clip negatives
    cutout = np.clip(cutout - np.nanmin(cutout), 0, None)

    try:
        rcx, rcy = centroid_2dg(cutout)
    except Exception:
        try:
            rcx, rcy = centroid_com(cutout)
        except Exception:
            return cx, cy

    # Sanity check — stay inside cutout
    if not (0 <= rcx < cutout.shape[1] and 0 <= rcy < cutout.shape[0]):
        return cx, cy

    return x1 + rcx, y1 + rcy


# ---------------------------------------------------------------------------
# Interactive strategy
# ---------------------------------------------------------------------------

def _interactive(aligner, stretch='zscale', box=15, save_catalog=None):
    """
    Open a two-panel matplotlib window for interactive source matching.

    State machine:
        IDLE    — left-click reference panel selects a source
        PENDING — left-click target panel confirms it (right-click cancels)
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    ref_frame = aligner.reference
    tgt_frame = aligner.targets[0]   # interactive always acts on first target

    ref_img = ref_frame.image2d
    tgt_img = tgt_frame.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError(
            "Call preprocess() before find_sources(strategy='interactive').")

    # ---- Figure setup -------------------------------------------------------
    fig, (ax_ref, ax_tgt) = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle(
        'Left-click reference to select source  |  u: undo  |  Enter: done',
        fontsize=10)

    ax_ref.imshow(ref_img, origin='lower', cmap='gray',
                  norm=_compute_norm(ref_img, stretch), interpolation='nearest')
    ax_tgt.imshow(tgt_img, origin='lower', cmap='gray',
                  norm=_compute_norm(tgt_img, stretch), interpolation='nearest')
    ax_ref.set_title('Reference')
    ax_tgt.set_title('Target')

    # ---- State ---------------------------------------------------------------
    state = {'mode': 'IDLE'}           # 'IDLE' | 'PENDING'
    pending = {
        'ref_x': None, 'ref_y': None,
        'ra': None,    'dec': None,
        'pred_x': None,'pred_y': None,
        'marker_ref': [],              # matplotlib artists on ax_ref
        'marker_tgt': [],              # matplotlib artists on ax_tgt
    }
    confirmed_markers_ref = []         # list of lists (one per stored pair)
    confirmed_markers_tgt = []

    # ---- Helper: draw confirmed pair markers --------------------------------
    def _draw_confirmed(ax, x, y, n, container):
        m1, = ax.plot(x, y, 'r+', ms=14, mew=2)
        m2 = ax.text(x + 3, y + 3, str(n), color='red', fontsize=9)
        container.append([m1, m2])

    def _remove_markers(mlist):
        for artist in mlist:
            try:
                artist.remove()
            except Exception:
                pass

    def _remove_all_confirmed():
        for grp in confirmed_markers_ref:
            _remove_markers(grp)
        for grp in confirmed_markers_tgt:
            _remove_markers(grp)
        confirmed_markers_ref.clear()
        confirmed_markers_tgt.clear()

    def _redraw_all_confirmed():
        _remove_all_confirmed()
        for i, pair in enumerate(aligner.pairs):
            rx, ry, ra, dec, tx, ty = pair
            grp_r = []
            m1, = ax_ref.plot(rx, ry, 'r+', ms=14, mew=2)
            m2 = ax_ref.text(rx + 3, ry + 3, str(i + 1),
                              color='red', fontsize=9)
            grp_r.extend([m1, m2])
            confirmed_markers_ref.append(grp_r)

            grp_t = []
            m3, = ax_tgt.plot(tx, ty, 'r+', ms=14, mew=2)
            m4 = ax_tgt.text(tx + 3, ty + 3, str(i + 1),
                              color='red', fontsize=9)
            grp_t.extend([m3, m4])
            confirmed_markers_tgt.append(grp_t)

    def _set_title():
        if state['mode'] == 'IDLE':
            fig.suptitle(
                'Left-click reference to select source  |  u: undo  |  Enter: done',
                fontsize=10)
        else:
            n = len(aligner.pairs) + 1
            fig.suptitle(
                f'Source {n} pending — click target to confirm  '
                '|  Space: accept prediction  |  Right-click: cancel',
                fontsize=10)

    # ---- Scroll zoom ---------------------------------------------------------
    def _on_scroll(event):
        ax = None
        if event.inaxes is ax_ref:
            ax = ax_ref
        elif event.inaxes is ax_tgt:
            ax = ax_tgt
        if ax is None or event.xdata is None or event.ydata is None:
            return
        factor = 0.85 if event.button == 'up' else 1.15
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        cx, cy = event.xdata, event.ydata
        ax.set_xlim([cx + (x - cx) * factor for x in xlim])
        ax.set_ylim([cy + (y - cy) * factor for y in ylim])
        fig.canvas.draw_idle()

    # ---- Accept pair (shared between click and spacebar) --------------------
    def _accept_pair(tgt_x_refined, tgt_y_refined):
        rx = pending['ref_x']
        ry = pending['ref_y']
        ra = pending['ra']
        dec = pending['dec']
        aligner.pairs.append((rx, ry, ra, dec, tgt_x_refined, tgt_y_refined))
        if aligner.on_pair_callback is not None:
            aligner.on_pair_callback(rx, ry, tgt_x_refined, tgt_y_refined, ra, dec)
        # Remove pending markers
        for m in pending['marker_ref'] + pending['marker_tgt']:
            _remove_markers(m if isinstance(m, list) else [m])
        pending['marker_ref'].clear()
        pending['marker_tgt'].clear()
        # Redraw all confirmed
        _redraw_all_confirmed()
        state['mode'] = 'IDLE'
        _set_title()
        fig.canvas.draw_idle()

    # ---- Click handler -------------------------------------------------------
    def _on_click(event):
        # Ignore clicks during zoom/pan
        try:
            if fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != '':
                return
        except Exception:
            pass
        if event.xdata is None or event.ydata is None:
            return
        if event.inaxes not in (ax_ref, ax_tgt):
            return

        if state['mode'] == 'IDLE':
            if event.inaxes is not ax_ref:
                return
            if event.button != 1:
                return
            # Centroid on reference
            cx, cy = _refine_centroid(ref_img, event.xdata, event.ydata, box)
            # Convert to sky
            sky = ref_frame.wcs.all_pix2world([[cx, cy]], 0)
            ra, dec = float(sky[0, 0]), float(sky[0, 1])
            # Project to target
            try:
                tgt_pred = tgt_frame.wcs.all_world2pix([[ra, dec]], 0)
                px, py = float(tgt_pred[0, 0]), float(tgt_pred[0, 1])
            except Exception:
                px, py = tgt_img.shape[1] / 2, tgt_img.shape[0] / 2

            pending['ref_x'] = cx
            pending['ref_y'] = cy
            pending['ra'] = ra
            pending['dec'] = dec
            pending['pred_x'] = px
            pending['pred_y'] = py

            # Draw solid red + on reference
            m1, = ax_ref.plot(cx, cy, 'r+', ms=14, mew=2)
            pending['marker_ref'].clear()
            pending['marker_ref'].append(m1)

            # Draw dashed cyan circle on target at predicted position
            circ = mpatches.Circle((px, py), radius=box,
                                   edgecolor='cyan', facecolor='none',
                                   linestyle='--', linewidth=1.5)
            ax_tgt.add_patch(circ)
            pending['marker_tgt'].clear()
            pending['marker_tgt'].append(circ)

            state['mode'] = 'PENDING'
            _set_title()
            fig.canvas.draw_idle()

        elif state['mode'] == 'PENDING':
            if event.button == 3:          # right-click → cancel
                for m in pending['marker_ref']:
                    _remove_markers(m if isinstance(m, list) else [m])
                for m in pending['marker_tgt']:
                    _remove_markers(m if isinstance(m, list) else [m])
                pending['marker_ref'].clear()
                pending['marker_tgt'].clear()
                state['mode'] = 'IDLE'
                _set_title()
                fig.canvas.draw_idle()
                return

            if event.button != 1:
                return
            if event.inaxes is not ax_tgt:
                return
            # Centroid on target
            tx, ty = _refine_centroid(tgt_img, event.xdata, event.ydata, box)
            _accept_pair(tx, ty)

    # ---- Key handler ---------------------------------------------------------
    def _on_key(event):
        if event.key == 'u':
            # Undo last confirmed pair
            if aligner.pairs:
                aligner.pairs.pop()
                _redraw_all_confirmed()
                fig.canvas.draw_idle()
        elif event.key == ' ':
            # Spacebar — auto-accept around predicted position
            if state['mode'] == 'PENDING':
                px = pending['pred_x']
                py = pending['pred_y']
                tx, ty = _refine_centroid(tgt_img, px, py, box)
                _accept_pair(tx, ty)
        elif event.key == 'enter':
            plt.close(fig)

    # ---- Window close --------------------------------------------------------
    def _on_close(event):
        pass  # pairs already stored in aligner.pairs

    # ---- Connect events -------------------------------------------------------
    fig.canvas.mpl_connect('button_press_event', _on_click)
    fig.canvas.mpl_connect('key_press_event', _on_key)
    fig.canvas.mpl_connect('close_event', _on_close)
    fig.canvas.mpl_connect('scroll_event', _on_scroll)

    plt.tight_layout()
    plt.show()

    # ---- Save catalog if requested ------------------------------------------
    if save_catalog is not None:
        _save_pairs_catalog(aligner.pairs, save_catalog)


# ---------------------------------------------------------------------------
# Batch strategy
# ---------------------------------------------------------------------------

def _batch(aligner, catalog: str, box: int = 15):
    """
    Load a saved pairs catalog and recentroid on each target frame.

    The catalog provides (ra, dec) sky positions (truth from reference WCS).
    For each target, these are projected to pixels and recentroided.
    """
    from astropy.table import Table

    tbl = Table.read(catalog)
    ras  = np.asarray(tbl['ra'])
    decs = np.asarray(tbl['dec'])
    rxs  = np.asarray(tbl['x_ref'])
    rys  = np.asarray(tbl['y_ref'])

    # We operate on target 0 by default for single-target, but for
    # align_many the per-frame recentroiding happens in align.py.
    # Here we populate aligner.pairs for target 0.
    if not aligner.targets:
        raise RuntimeError("No targets loaded.")
    tgt = aligner.targets[0]
    img = tgt.image2d
    if img is None:
        raise RuntimeError("Call preprocess() before find_sources(strategy='batch').")
    ny, nx = img.shape

    aligner.pairs.clear()
    for rx, ry, ra, dec in zip(rxs, rys, ras, decs):
        try:
            px_arr = tgt.wcs.all_world2pix([[ra, dec]], 0)
            px, py = float(px_arr[0, 0]), float(px_arr[0, 1])
        except Exception:
            continue
        if not (0 <= px < nx and 0 <= py < ny):
            continue
        tx, ty = _refine_centroid(img, px, py, box)
        aligner.pairs.append((rx, ry, ra, dec, tx, ty))


# ---------------------------------------------------------------------------
# Gaia strategy
# ---------------------------------------------------------------------------

def _gaia(aligner, box: int = 15):
    """
    Match Gaia DR3 catalog to reference image, then project to target.

    Requires astroquery.
    """
    try:
        from astroquery.gaia import Gaia
    except ImportError:
        raise ImportError(
            "astroquery is required for strategy='gaia'. "
            "Install with: pip install astroquery"
        )
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    ref = aligner.reference
    img = ref.image2d
    if img is None:
        raise RuntimeError("Call preprocess() before find_sources(strategy='gaia').")

    ny, nx = img.shape
    # Compute centre and radius from WCS
    try:
        centre_sky = ref.wcs.all_pix2world([[nx / 2, ny / 2]], 0)
        ra_c, dec_c = float(centre_sky[0, 0]), float(centre_sky[0, 1])
        corner_sky = ref.wcs.all_pix2world([[0, 0]], 0)
        dra = abs(float(corner_sky[0, 0]) - ra_c)
        ddec = abs(float(corner_sky[0, 1]) - dec_c)
        radius_deg = min(max(dra, ddec) * 1.5, 0.5)   # cap at 0.5 deg
    except Exception:
        raise RuntimeError("Reference WCS is not valid — cannot query Gaia.")

    coord = SkyCoord(ra=ra_c * u.deg, dec=dec_c * u.deg)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            j = Gaia.cone_search_async(coord, radius=radius_deg * u.deg)
            result = j.get_results()
    except Exception as e:
        raise RuntimeError(f"Gaia query failed: {e}")

    if len(result) == 0:
        raise RuntimeError("No Gaia sources found in field.")

    # Filter bright, well-measured stars
    mag_col = 'phot_g_mean_mag' if 'phot_g_mean_mag' in result.colnames else result.colnames[0]
    result.sort(mag_col)
    result = result[:50]   # use up to 50 brightest

    tgt = aligner.targets[0]
    tgt_img = tgt.image2d
    if tgt_img is None:
        raise RuntimeError("Call preprocess() first.")
    t_ny, t_nx = tgt_img.shape

    aligner.pairs.clear()
    for row in result:
        ra = float(row['ra'])
        dec = float(row['dec'])
        # Project to reference pixels
        try:
            r_pix = ref.wcs.all_world2pix([[ra, dec]], 0)
            rx, ry = float(r_pix[0, 0]), float(r_pix[0, 1])
        except Exception:
            continue
        if not (0 <= rx < nx and 0 <= ry < ny):
            continue
        # Refine on reference
        rx_r, ry_r = _refine_centroid(img, rx, ry, box)
        # Project to target
        try:
            t_pix = tgt.wcs.all_world2pix([[ra, dec]], 0)
            tx, ty = float(t_pix[0, 0]), float(t_pix[0, 1])
        except Exception:
            continue
        if not (0 <= tx < t_nx and 0 <= ty < t_ny):
            continue
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, box)
        aligner.pairs.append((rx_r, ry_r, ra, dec, tx_r, ty_r))

    if not aligner.pairs:
        raise RuntimeError("Gaia sources found but none projected into both images.")


# ---------------------------------------------------------------------------
# DAO strategy
# ---------------------------------------------------------------------------

def _dao(aligner, box: int = 15):
    """
    Detect point sources with DAOStarFinder, then match between images.

    Requires photutils; falls back to SExtractor-style peak finding if absent.
    """
    ref = aligner.reference
    tgt = aligner.targets[0]
    ref_img = ref.image2d
    tgt_img = tgt.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError("Call preprocess() before find_sources(strategy='dao').")

    ref_sources = _detect_dao_sources(ref_img, box)
    if len(ref_sources) == 0:
        raise RuntimeError("DAOStarFinder found no sources in reference image.")

    tgt_ny, tgt_nx = tgt_img.shape
    aligner.pairs.clear()
    for (rx, ry) in ref_sources:
        sky = ref.wcs.all_pix2world([[rx, ry]], 0)
        ra, dec = float(sky[0, 0]), float(sky[0, 1])
        try:
            t_pix = tgt.wcs.all_world2pix([[ra, dec]], 0)
            tx, ty = float(t_pix[0, 0]), float(t_pix[0, 1])
        except Exception:
            continue
        if not (0 <= tx < tgt_nx and 0 <= ty < tgt_ny):
            continue
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, box)
        aligner.pairs.append((rx, ry, ra, dec, tx_r, ty_r))


def _detect_dao_sources(img: np.ndarray, box: int):
    """Return list of (x, y) source positions from DAOStarFinder or fallback."""
    if HAS_PHOTUTILS:
        try:
            from photutils.detection import DAOStarFinder
            from astropy.stats import sigma_clipped_stats
            _, median, std = sigma_clipped_stats(img[np.isfinite(img)])
            finder = DAOStarFinder(fwhm=2.5, threshold=5.0 * std)
            sources = finder(img - median)
            if sources is None or len(sources) == 0:
                return []
            return list(zip(
                np.asarray(sources['xcentroid']),
                np.asarray(sources['ycentroid'])
            ))
        except Exception:
            pass
    # Fallback: simple peak finding via scipy
    return _simple_peak_find(img, box)


def _simple_peak_find(img: np.ndarray, box: int, n_max: int = 30):
    """Basic peak finding: find local maxima above 3-sigma threshold."""
    from scipy.ndimage import maximum_filter, label
    clean = np.nan_to_num(img, nan=0.0)
    flat = clean.ravel()
    flat_finite = flat[np.isfinite(flat)]
    threshold = np.median(flat_finite) + 3.0 * np.std(flat_finite)
    local_max = maximum_filter(clean, size=box) == clean
    above_thresh = clean > threshold
    peaks = local_max & above_thresh
    labeled, n = label(peaks)
    sources = []
    for i in range(1, min(n + 1, n_max + 1)):
        ys, xs = np.where(labeled == i)
        sources.append((float(xs.mean()), float(ys.mean())))
    return sources


# ---------------------------------------------------------------------------
# Knots strategy
# ---------------------------------------------------------------------------

def _knots(aligner, box: int = 15):
    """
    Detect bright knots in reference via segmentation, project to target.

    Designed for extended sources (arcs, nebulae, etc.).
    """
    ref = aligner.reference
    tgt = aligner.targets[0]
    ref_img = ref.image2d
    tgt_img = tgt.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError("Call preprocess() before find_sources(strategy='knots').")

    knot_positions = _detect_knots(ref_img, box)
    if len(knot_positions) == 0:
        raise RuntimeError("No knots detected in reference image.")

    tgt_ny, tgt_nx = tgt_img.shape
    aligner.pairs.clear()
    for (rx, ry) in knot_positions:
        sky = ref.wcs.all_pix2world([[rx, ry]], 0)
        ra, dec = float(sky[0, 0]), float(sky[0, 1])
        try:
            t_pix = tgt.wcs.all_world2pix([[ra, dec]], 0)
            tx, ty = float(t_pix[0, 0]), float(t_pix[0, 1])
        except Exception:
            continue
        if not (0 <= tx < tgt_nx and 0 <= ty < tgt_ny):
            continue
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, box)
        aligner.pairs.append((rx, ry, ra, dec, tx_r, ty_r))


def _detect_knots(img: np.ndarray, box: int, n_knots: int = 10):
    """Segment image and return centroids of brightest segments."""
    if HAS_PHOTUTILS:
        try:
            from photutils.segmentation import detect_sources, SourceCatalog
            from astropy.stats import sigma_clipped_stats
            _, _, std = sigma_clipped_stats(img[np.isfinite(img)])
            seg = detect_sources(np.nan_to_num(img, nan=0.0),
                                  threshold=2.0 * std, npixels=5)
            if seg is None:
                return _simple_peak_find(img, box, n_max=n_knots)
            cat = SourceCatalog(np.nan_to_num(img, nan=0.0), seg)
            tbl = cat.to_table()
            tbl.sort('segment_flux', reverse=True)
            tbl = tbl[:n_knots]
            return list(zip(
                np.asarray(tbl['xcentroid']),
                np.asarray(tbl['ycentroid'])
            ))
        except Exception:
            pass
    return _simple_peak_find(img, box, n_max=n_knots)


# ---------------------------------------------------------------------------
# Cross-correlation strategy
# ---------------------------------------------------------------------------

def _cross_corr(aligner):
    """
    Find the shift between reference and target via cross-correlation.

    Works best for same-instrument pairs with similar PSFs.
    Returns one pair representing the centre of the reference after shift.
    """
    from scipy.signal import fftconvolve

    ref = aligner.reference
    tgt = aligner.targets[0]
    ref_img = ref.image2d
    tgt_img = tgt.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError("Call preprocess() before find_sources(strategy='cross_corr').")

    # Pad/crop to same size for correlation
    r_ny, r_nx = ref_img.shape
    t_ny, t_nx = tgt_img.shape
    ny = max(r_ny, t_ny)
    nx = max(r_nx, t_nx)

    def _pad(arr, out_ny, out_nx):
        padded = np.zeros((out_ny, out_nx), dtype=np.float64)
        src_ny, src_nx = arr.shape
        padded[:src_ny, :src_nx] = np.nan_to_num(arr, nan=0.0)
        return padded

    ref_pad = _pad(ref_img, ny, nx)
    tgt_pad = _pad(tgt_img, ny, nx)

    # Normalise
    ref_norm = ref_pad - ref_pad.mean()
    tgt_norm = tgt_pad - tgt_pad.mean()

    # Cross-correlate: convolve with flipped target
    xcorr = fftconvolve(ref_norm, tgt_norm[::-1, ::-1], mode='full')
    peak = np.unravel_index(np.argmax(xcorr), xcorr.shape)
    dy = peak[0] - (ny - 1)
    dx = peak[1] - (nx - 1)

    # Build one synthetic pair: centre of reference shifted by (dx, dy)
    ref_cx = r_nx / 2.0
    ref_cy = r_ny / 2.0
    tgt_cx = ref_cx - dx
    tgt_cy = ref_cy - dy

    sky = ref.wcs.all_pix2world([[ref_cx, ref_cy]], 0)
    ra, dec = float(sky[0, 0]), float(sky[0, 1])

    aligner.pairs.clear()
    aligner.pairs.append((ref_cx, ref_cy, ra, dec, tgt_cx, tgt_cy))


# ---------------------------------------------------------------------------
# Auto cascade
# ---------------------------------------------------------------------------

def _auto(aligner, stretch='zscale', box=15):
    """
    Try strategies in cascade order, fall back gracefully.

    Order: gaia → dao → knots → cross_corr → interactive
    """
    strategies = [
        ('gaia', lambda: _gaia(aligner, box=box)),
        ('dao',  lambda: _dao(aligner, box=box)),
        ('knots', lambda: _knots(aligner, box=box)),
        ('cross_corr', lambda: _cross_corr(aligner)),
        ('interactive', lambda: _interactive(aligner, stretch=stretch, box=box)),
    ]
    for name, fn in strategies:
        try:
            fn()
            if aligner.pairs:
                print(f"[rb_align] auto: used strategy '{name}', "
                      f"found {len(aligner.pairs)} pair(s).")
                return
        except (ImportError, RuntimeError, Exception) as e:
            print(f"[rb_align] auto: strategy '{name}' failed ({e}). Trying next…")
    print("[rb_align] auto: all automated strategies exhausted; "
          "falling through to interactive.")


# ---------------------------------------------------------------------------
# Catalog save / load
# ---------------------------------------------------------------------------

def _save_pairs_catalog(pairs, path: str):
    """Save matched pairs to a FITS table catalog."""
    from astropy.table import Table

    if not pairs:
        warnings.warn("[rb_align] No pairs to save.")
        return

    rows = {
        'x_ref': [], 'y_ref': [], 'ra': [], 'dec': [],
        'x_target': [], 'y_target': [], 'residual_arcsec': [],
    }
    for p in pairs:
        rx, ry, ra, dec, tx, ty = p
        rows['x_ref'].append(float(rx))
        rows['y_ref'].append(float(ry))
        rows['ra'].append(float(ra))
        rows['dec'].append(float(dec))
        rows['x_target'].append(float(tx))
        rows['y_target'].append(float(ty))
        rows['residual_arcsec'].append(np.nan)   # filled by align()

    tbl = Table(rows)
    tbl.write(path, overwrite=True)
    print(f"[rb_align] saved source catalog → {path}")


# ---------------------------------------------------------------------------
# Public dispatcher
# ---------------------------------------------------------------------------

def _run_strategy(aligner,
                  strategy: str = 'interactive',
                  stretch='zscale',
                  box: int = 15,
                  save_catalog: Optional[str] = None,
                  catalog: Optional[str] = None):
    """Dispatch to the selected source detection strategy."""

    strategy = strategy.lower()

    if strategy == 'interactive':
        _interactive(aligner, stretch=stretch, box=box, save_catalog=save_catalog)

    elif strategy == 'gaia':
        _gaia(aligner, box=box)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'dao':
        _dao(aligner, box=box)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'knots':
        _knots(aligner, box=box)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'cross_corr':
        _cross_corr(aligner)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'batch':
        if catalog is None:
            raise ValueError("strategy='batch' requires catalog='path/to/catalog.fits'.")
        _batch(aligner, catalog=catalog, box=box)

    elif strategy == 'auto':
        _auto(aligner, stretch=stretch, box=box)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    else:
        raise ValueError(
            f"Unknown strategy '{strategy}'. Choose from: "
            "interactive, gaia, dao, knots, cross_corr, batch, auto."
        )

    print(f"[rb_align] find_sources('{strategy}'): "
          f"{len(aligner.pairs)} pair(s) stored.")
