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
            if not aligner.pairs:
                return
            # Delete the pair whose reference marker is nearest to the cursor.
            # Falls back to the last pair if the cursor is not over an axes.
            if (event.inaxes in (ax_ref, ax_tgt)
                    and event.xdata is not None and event.ydata is not None):
                cx, cy = event.xdata, event.ydata
                if event.inaxes is ax_ref:
                    dists = [np.sqrt((p[0]-cx)**2 + (p[1]-cy)**2)
                             for p in aligner.pairs]
                else:  # ax_tgt
                    dists = [np.sqrt((p[4]-cx)**2 + (p[5]-cy)**2)
                             for p in aligner.pairs]
                del aligner.pairs[int(np.argmin(dists))]
            else:
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

def _gaia(aligner, box: int = 5, radius_deg: float = 0.5,
          max_stars: int = 50):
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
        radius_deg = min(max(dra, ddec) * 1.5, radius_deg)
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
    result = result[:max_stars]

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

def _dao(aligner, box: int = 5, fwhm: float = 0.0,
         threshold_sigma: float = 3.0, max_sources: int = 50):
    """
    Detect point sources with DAOStarFinder, then match between images.

    Requires photutils; falls back to peak finding if absent.
    """
    ref = aligner.reference
    tgt = aligner.targets[0]
    ref_img = ref.image2d
    tgt_img = tgt.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError("Call preprocess() before find_sources(strategy='dao').")

    ref_sources = _detect_dao_sources(ref_img, box, fwhm=fwhm,
                                       threshold_sigma=threshold_sigma)
    if len(ref_sources) == 0:
        raise RuntimeError("DAOStarFinder found no sources in reference image.")

    # Sort by brightness and cap
    ref_sources = sorted(
        ref_sources,
        key=lambda xy: -float(ref_img[min(int(round(xy[1])), ref_img.shape[0]-1),
                                       min(int(round(xy[0])), ref_img.shape[1]-1)])
    )[:max_sources]

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


def _detect_dao_sources(img: np.ndarray, box: int,
                        fwhm: float = 0.0,
                        threshold_sigma: float = 3.0):
    """
    Return list of (x, y) source positions from DAOStarFinder or fallback.

    Parameters
    ----------
    fwhm : float
        Expected source FWHM in pixels.  0 → auto as ``box * 0.5``.
    threshold_sigma : float
        Detection threshold in units of image σ.
    """
    if fwhm <= 0:
        fwhm = max(2.0, box * 0.5)

    if HAS_PHOTUTILS:
        try:
            from photutils.detection import DAOStarFinder
            from astropy.stats import sigma_clipped_stats
            finite = img[np.isfinite(img)]
            _, median, std = sigma_clipped_stats(finite)
            finder = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * std)
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
    return _simple_peak_find(img, box, threshold_sigma=threshold_sigma)


def _simple_peak_find(img: np.ndarray, box: int, n_max: int = 30,
                      threshold_sigma: float = 3.0):
    """
    Basic peak finding: local maxima above threshold_sigma × σ_background.

    Uses sigma-clipped stats to estimate background and noise, matching
    the behaviour of the DAOStarFinder path.
    """
    from scipy.ndimage import maximum_filter, label
    clean = np.nan_to_num(img, nan=0.0)
    finite = clean.ravel()
    finite = finite[np.isfinite(finite)]

    # Sigma-clipped background estimate (2 iterations, 3-sigma clip)
    bg = finite.copy()
    for _ in range(2):
        med = np.median(bg)
        std = np.std(bg)
        bg = bg[np.abs(bg - med) < 3.0 * std]
    bg_median = float(np.median(bg))
    bg_std    = float(np.std(bg))

    threshold = bg_median + threshold_sigma * bg_std
    local_max = maximum_filter(clean, size=max(3, box)) == clean
    peaks = local_max & (clean > threshold)
    labeled, n = label(peaks)
    sources = []
    for i in range(1, min(n + 1, n_max + 1)):
        ys, xs = np.where(labeled == i)
        sources.append((float(xs.mean()), float(ys.mean())))
    return sources


# ---------------------------------------------------------------------------
# Knots strategy
# ---------------------------------------------------------------------------

def _knots(aligner, box: int = 5, threshold_sigma: float = 2.0,
           max_knots: int = 10):
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

    knot_positions = _detect_knots(ref_img, box, n_knots=max_knots,
                                    threshold_sigma=threshold_sigma)
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


def _detect_knots(img: np.ndarray, box: int, n_knots: int = 10,
                  threshold_sigma: float = 2.0):
    """Segment image and return centroids of brightest segments."""
    if HAS_PHOTUTILS:
        try:
            from photutils.segmentation import detect_sources, SourceCatalog
            from astropy.stats import sigma_clipped_stats
            _, _, std = sigma_clipped_stats(img[np.isfinite(img)])
            seg = detect_sources(np.nan_to_num(img, nan=0.0),
                                  threshold=threshold_sigma * std, npixels=5)
            if seg is None:
                return _simple_peak_find(img, box, n_max=n_knots,
                                          threshold_sigma=threshold_sigma)
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
    return _simple_peak_find(img, box, n_max=n_knots,
                              threshold_sigma=threshold_sigma)


# ---------------------------------------------------------------------------
# Cross-correlation strategy
# ---------------------------------------------------------------------------

def _estimate_pixel_scale_arcsec(frame, img: np.ndarray) -> float:
    """Estimate pixel scale in arcsec/px from the frame's WCS."""
    try:
        ny, nx = img.shape
        sky0 = frame.wcs.all_pix2world([[nx / 2,     ny / 2    ]], 0)[0]
        sky1 = frame.wcs.all_pix2world([[nx / 2 + 1, ny / 2    ]], 0)[0]
        sky2 = frame.wcs.all_pix2world([[nx / 2,     ny / 2 + 1]], 0)[0]
        cos_dec = np.cos(np.deg2rad(sky0[1]))
        dx = abs(sky1[0] - sky0[0]) * cos_dec * 3600.0
        dy = abs(sky2[1] - sky0[1]) * 3600.0
        return float((dx + dy) / 2.0)
    except Exception:
        return 0.3   # fallback: KCWI-like scale


def _cross_corr_image_fallback(aligner):
    """
    Single-pair shift via FFT image cross-correlation.

    Fallback used when source detection yields nothing in either image.
    Only reliable for same pixel-scale images. Returns 1 pair at the
    reference centre displaced by the measured shift.
    """
    from scipy.signal import fftconvolve

    ref = aligner.reference
    tgt = aligner.targets[0]
    ref_img = ref.image2d
    tgt_img = tgt.image2d

    r_ny, r_nx = ref_img.shape
    t_ny, t_nx = tgt_img.shape
    ny = max(r_ny, t_ny)
    nx = max(r_nx, t_nx)

    def _pad(arr, out_ny, out_nx):
        padded = np.zeros((out_ny, out_nx), dtype=np.float64)
        padded[:arr.shape[0], :arr.shape[1]] = np.nan_to_num(arr, nan=0.0)
        return padded

    ref_norm = _pad(ref_img, ny, nx)
    tgt_norm = _pad(tgt_img, ny, nx)
    ref_norm -= ref_norm.mean()
    tgt_norm -= tgt_norm.mean()

    xcorr = fftconvolve(ref_norm, tgt_norm[::-1, ::-1], mode='full')
    peak  = np.unravel_index(np.argmax(xcorr), xcorr.shape)
    dy    = peak[0] - (ny - 1)
    dx    = peak[1] - (nx - 1)

    ref_cx = r_nx / 2.0
    ref_cy = r_ny / 2.0
    tgt_cx = ref_cx - dx
    tgt_cy = ref_cy - dy

    sky = ref.wcs.all_pix2world([[ref_cx, ref_cy]], 0)
    ra, dec = float(sky[0, 0]), float(sky[0, 1])

    aligner.pairs.clear()
    aligner.pairs.append((ref_cx, ref_cy, ra, dec, tgt_cx, tgt_cy))


def _cross_corr(aligner, box: int = 5, fwhm: float = 0.0,
                threshold_sigma: float = 3.0, max_sources: int = 50):
    """
    Multi-source cross-match via Hough-style catalog voting.

    Algorithm
    ---------
    1. Detect sources in each image independently (in their own pixel space).
    2. Convert pixel positions to sky (RA/Dec) via each frame's own WCS.
    3. Build all pairwise (ref_i, tgt_j) sky-shift votes → 2-D histogram.
    4. The peak bin gives the best (ΔRA, ΔDec) offset between the catalogs.
    5. Match sources within that offset and recentroid each → N pairs.

    Works for different pixel scales and FoVs. Falls back to single-pair
    image FFT cross-correlation if either image yields no detectable sources.

    Parameters
    ----------
    box : int
        Centroid box half-width (px); also sets the Hough bin size as
        ``clip(box × pixel_scale, 1.5, 10)`` arcsec.
    """
    ref = aligner.reference
    tgt = aligner.targets[0]
    ref_img = ref.image2d
    tgt_img = tgt.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError(
            "Call preprocess() before find_sources(strategy='cross_corr').")

    r_ny, r_nx = ref_img.shape
    t_ny, t_nx = tgt_img.shape

    # ------------------------------------------------------------------
    # Step 1: source detection — keep up to 50 brightest per image
    # ------------------------------------------------------------------
    ref_sources = _detect_dao_sources(ref_img, box, fwhm=fwhm,
                                       threshold_sigma=threshold_sigma)
    tgt_sources = _detect_dao_sources(tgt_img, box, fwhm=fwhm,
                                       threshold_sigma=threshold_sigma)

    if not ref_sources or not tgt_sources:
        print("[rb_align] cross_corr: source detection failed; "
              "falling back to image FFT cross-correlation (1 pair).")
        _cross_corr_image_fallback(aligner)
        return

    def _brightest(img, sources, n_max=50):
        ny_, nx_ = img.shape
        vals = []
        for x, y in sources:
            xi, yi = int(round(x)), int(round(y))
            xi = max(0, min(nx_ - 1, xi))
            yi = max(0, min(ny_ - 1, yi))
            vals.append(float(img[yi, xi]))
        order = np.argsort(vals)[::-1][:n_max]
        return [sources[i] for i in order]

    ref_sources = _brightest(ref_img, ref_sources, n_max=max_sources)
    tgt_sources = _brightest(tgt_img, tgt_sources, n_max=max_sources)

    # ------------------------------------------------------------------
    # Step 2: pixel → sky for each catalog independently
    # ------------------------------------------------------------------
    try:
        ref_sky = ref.wcs.all_pix2world(np.array(ref_sources, dtype=np.float64), 0)
        tgt_sky = tgt.wcs.all_pix2world(np.array(tgt_sources, dtype=np.float64), 0)
    except Exception as e:
        print(f"[rb_align] cross_corr: WCS projection failed ({e}); "
              "falling back to image FFT cross-correlation.")
        _cross_corr_image_fallback(aligner)
        return

    ref_ra  = ref_sky[:, 0]
    ref_dec = ref_sky[:, 1]
    tgt_ra  = tgt_sky[:, 0]
    tgt_dec = tgt_sky[:, 1]

    mean_dec = float(np.mean(np.concatenate([ref_dec, tgt_dec])))
    cos_dec  = np.cos(np.deg2rad(mean_dec))

    # ------------------------------------------------------------------
    # Step 3: Hough vote — all pairwise (tgt_j − ref_i) shifts in arcsec
    # ------------------------------------------------------------------
    dra_all  = ((tgt_ra[:, None]  - ref_ra[None, :])  * cos_dec * 3600.0).ravel()
    ddec_all = ((tgt_dec[:, None] - ref_dec[None, :]) *           3600.0).ravel()

    px_scale   = _estimate_pixel_scale_arcsec(ref, ref_img)
    bin_arcsec = float(np.clip(box * px_scale, 1.5, 10.0))

    dra_lo,  dra_hi  = dra_all.min()  - bin_arcsec, dra_all.max()  + bin_arcsec
    ddec_lo, ddec_hi = ddec_all.min() - bin_arcsec, ddec_all.max() + bin_arcsec

    n_bins_ra  = max(3, int((dra_hi  - dra_lo)  / bin_arcsec) + 1)
    n_bins_dec = max(3, int((ddec_hi - ddec_lo) / bin_arcsec) + 1)

    hist, dra_edges, ddec_edges = np.histogram2d(
        dra_all, ddec_all,
        bins=[n_bins_ra, n_bins_dec],
        range=[[dra_lo, dra_hi], [ddec_lo, ddec_hi]],
    )

    peak_idx   = np.unravel_index(np.argmax(hist), hist.shape)
    peak_votes = int(hist[peak_idx])
    best_dra   = (dra_edges[peak_idx[0]]     + dra_edges[peak_idx[0] + 1])  / 2.0
    best_ddec  = (ddec_edges[peak_idx[1]]    + ddec_edges[peak_idx[1] + 1]) / 2.0

    print(f"[rb_align] cross_corr: Hough peak = {peak_votes} vote(s), "
          f"ΔRA={best_dra:.2f}\", ΔDec={best_ddec:.2f}\"  "
          f"(bin={bin_arcsec:.1f}\", "
          f"N_ref={len(ref_sources)}, N_tgt={len(tgt_sources)})")

    if peak_votes < 2 and len(ref_sources) > 3 and len(tgt_sources) > 3:
        print("[rb_align] cross_corr: weak Hough peak — shift may be unreliable.")

    # ------------------------------------------------------------------
    # Step 4: match sources under the best shift; recentroid each pair
    # ------------------------------------------------------------------
    match_tol_arcsec = bin_arcsec * 1.5

    aligner.pairs.clear()
    used_tgt = set()

    for i in range(len(ref_sources)):
        rx, ry   = ref_sources[i]
        ra_r     = float(ref_ra[i])
        dec_r    = float(ref_dec[i])

        # Expected sky position of this source in the target frame
        exp_ra  = ra_r  + best_dra  / (cos_dec * 3600.0)
        exp_dec = dec_r + best_ddec / 3600.0

        # Distance from every target source to the expected position
        dra_t  = (tgt_ra  - exp_ra)  * cos_dec * 3600.0
        ddec_t = (tgt_dec - exp_dec) * 3600.0
        sep    = np.sqrt(dra_t**2 + ddec_t**2)

        for j_used in used_tgt:
            sep[j_used] = np.inf

        j = int(np.argmin(sep))
        if sep[j] > match_tol_arcsec:
            continue

        used_tgt.add(j)
        tx, ty = tgt_sources[j]

        if not (0 <= rx < r_nx and 0 <= ry < r_ny):
            continue
        if not (0 <= tx < t_nx and 0 <= ty < t_ny):
            continue

        rx_r, ry_r = _refine_centroid(ref_img, rx, ry, box)
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, box)

        sky = ref.wcs.all_pix2world([[rx_r, ry_r]], 0)
        ra, dec = float(sky[0, 0]), float(sky[0, 1])

        aligner.pairs.append((rx_r, ry_r, ra, dec, tx_r, ty_r))

    if not aligner.pairs:
        print("[rb_align] cross_corr: no source matches after Hough vote; "
              "falling back to image FFT cross-correlation (1 pair).")
        _cross_corr_image_fallback(aligner)


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
        ('cross_corr', lambda: _cross_corr(aligner, box=box)),
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
                  box: int = 5,
                  save_catalog: Optional[str] = None,
                  catalog: Optional[str] = None,
                  # cross_corr / dao options
                  fwhm: float = 0.0,
                  threshold_sigma: float = 3.0,
                  max_sources: int = 50,
                  # knots options
                  max_knots: int = 10,
                  # gaia options
                  radius_deg: float = 0.5,
                  max_stars: int = 50):
    """Dispatch to the selected source detection strategy."""

    strategy = strategy.lower()

    if strategy == 'interactive':
        _interactive(aligner, stretch=stretch, box=box, save_catalog=save_catalog)

    elif strategy == 'gaia':
        _gaia(aligner, box=box, radius_deg=radius_deg, max_stars=max_stars)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'dao':
        _dao(aligner, box=box, fwhm=fwhm,
             threshold_sigma=threshold_sigma, max_sources=max_sources)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'knots':
        _knots(aligner, box=box, threshold_sigma=threshold_sigma,
               max_knots=max_knots)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'cross_corr':
        _cross_corr(aligner, box=box, fwhm=fwhm,
                    threshold_sigma=threshold_sigma, max_sources=max_sources)
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
