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

def _compute_norm(data2d: np.ndarray, stretch='zscale',
                  scale='linear', vmin=None, vmax=None):
    """
    Build a matplotlib Normalize object for the given 2D array.

    Parameters
    ----------
    data2d : np.ndarray
    stretch : str or (vmin, vmax)
        'zscale' | percentile string | 'minmax' | 'manual'
        or a (vmin, vmax) tuple/list for manual limits.
    scale : str
        'linear' (default) | 'log' | 'sqrt'
        Transfer function applied after the stretch limits are determined.
    vmin, vmax : float or None
        Explicit limits; used when stretch == 'manual' or a tuple is passed.
    """
    import matplotlib.colors as mcolors
    from astropy.visualization import ZScaleInterval

    finite = data2d[np.isfinite(data2d)]
    if len(finite) == 0:
        return mcolors.Normalize(vmin=0, vmax=1)

    # Step 1 — determine lo/hi from stretch
    if stretch == 'manual' and vmin is not None and vmax is not None:
        lo, hi = float(vmin), float(vmax)
    elif isinstance(stretch, (tuple, list)):
        lo, hi = float(stretch[0]), float(stretch[1])
    elif stretch == 'zscale':
        try:
            lo, hi = ZScaleInterval().get_limits(finite)
        except Exception:
            lo, hi = np.nanpercentile(finite, [1, 99])
    elif stretch == 'minmax':
        lo, hi = float(finite.min()), float(finite.max())
    else:  # percentile string like '99%'
        pct = float(str(stretch).rstrip('%'))
        lo, hi = np.nanpercentile(finite, [100.0 - pct, pct])

    # Step 2 — build norm for the requested scale type
    if scale == 'log':
        pos = finite[finite > 0]
        if lo <= 0:
            lo = float(pos.min()) if len(pos) else 1e-10
        if hi <= lo:
            hi = lo * 1000
        return mcolors.LogNorm(vmin=lo, vmax=hi)
    elif scale == 'sqrt':
        lo = max(lo, 0.0)   # sqrt of negative is undefined
        return mcolors.PowerNorm(gamma=0.5, vmin=lo, vmax=hi)
    else:  # linear
        return mcolors.Normalize(vmin=lo, vmax=hi)


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

    # Replace NaNs/infs with local median
    bad = ~np.isfinite(cutout)
    if bad.all():
        return cx, cy
    finite = cutout[~bad]
    bg_median = float(np.median(finite))
    bg_std    = float(1.4826 * np.median(np.abs(finite - bg_median)))
    cutout[bad] = bg_median

    # Subtract background and keep only signal above 1σ
    cutout = np.clip(cutout - (bg_median + bg_std), 0, None)
    if cutout.sum() == 0:
        return cx, cy

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

def _interactive(aligner, stretch='zscale', box=0.1, save_catalog=None):
    """
    Open a two-panel matplotlib window for interactive source matching.

    Parameters
    ----------
    box : float
        Centroid box half-width in arcsec (converted to pixels per frame via WCS).
        Default 1.0".

    State machine:
        IDLE    — left-click reference panel selects a new source
                  left-click target panel re-places an edit-mode circle
        PENDING — left-click target panel confirms pending source

    Keyboard / mouse:
        double-click — toggle edit mode for nearest pair
                       (target marker becomes a numbered cyan dashed circle;
                        left-click on target re-places it)
        u            — hover over any marker → delete nearest pair
        right-click  — IDLE: delete nearest pair  |  PENDING: cancel
        Space        — PENDING: auto-accept centroid at predicted position
        Enter        — close window, keep all confirmed pairs

    If existing pairs are in aligner.pairs when called, they are drawn
    immediately — useful after strategy='batch' or remove_pair().
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    ref_frame = aligner.reference
    tgt_frame = aligner.targets[0]

    ref_img = ref_frame.image2d
    tgt_img = tgt_frame.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError(
            "Call preprocess() before find_sources(strategy='interactive').")

    # Convert box from arcsec to pixels for each frame independently
    ref_box_px = _box_to_pixels(box, ref_frame, ref_img)
    tgt_box_px = _box_to_pixels(box, tgt_frame, tgt_img)

    # ---- Figure setup -------------------------------------------------------
    fig, (ax_ref, ax_tgt) = plt.subplots(1, 2, figsize=(12, 6))

    ax_ref.imshow(ref_img, origin='lower', cmap='gray',
                  norm=_compute_norm(ref_img, stretch), interpolation='nearest')
    ax_tgt.imshow(tgt_img, origin='lower', cmap='gray',
                  norm=_compute_norm(tgt_img, stretch), interpolation='nearest')
    ax_ref.set_title('Reference')
    ax_tgt.set_title('Target')

    # ---- State ---------------------------------------------------------------
    state = {'mode': 'IDLE'}
    pending = {
        'ref_x': None, 'ref_y': None,
        'ra': None,    'dec': None,
        'pred_x': None, 'pred_y': None,
        'marker_ref': [],
        'marker_tgt': [],
    }

    # confirmed_markers_ref/tgt: list of artist-lists, one per pair (parallel to aligner.pairs)
    confirmed_markers_ref = []
    confirmed_markers_tgt = []

    # edit_mode: {pair_idx: {'circle': patch, 'label': text}}
    # pairs in edit_mode have their target red-cross hidden and a cyan circle shown instead
    edit_mode = {}

    # ---- Marker helpers ------------------------------------------------------
    def _remove_artists(alist):
        for a in alist:
            try:
                a.remove()
            except Exception:
                pass

    def _clear_edit_mode():
        for info in edit_mode.values():
            try:
                info['circle'].remove()
                info['label'].remove()
            except Exception:
                pass
        edit_mode.clear()

    def _redraw_all_confirmed():
        """Rebuild all confirmed markers from aligner.pairs.
        Pairs in edit_mode get a cyan circle on the target instead of a red cross."""
        for grp in confirmed_markers_ref:
            _remove_artists(grp)
        for grp in confirmed_markers_tgt:
            _remove_artists(grp)
        confirmed_markers_ref.clear()
        confirmed_markers_tgt.clear()

        # Also remove old edit circles — will be re-created below
        for info in edit_mode.values():
            try:
                info['circle'].remove()
                info['label'].remove()
            except Exception:
                pass

        for i, (rx, ry, ra, dec, tx, ty) in enumerate(aligner.pairs):
            n = i + 1
            m1, = ax_ref.plot(rx, ry, 'r+', ms=14, mew=2)
            m2 = ax_ref.text(rx + 3, ry + 3, str(n), color='red', fontsize=9)
            confirmed_markers_ref.append([m1, m2])

            if i in edit_mode:
                # Re-draw cyan circle (artists were removed above)
                circ = mpatches.Circle((tx, ty), radius=tgt_box_px,
                                       edgecolor='cyan', facecolor='none',
                                       linestyle='--', linewidth=1.5, zorder=5)
                ax_tgt.add_patch(circ)
                lbl = ax_tgt.text(tx + tgt_box_px + 2, ty + tgt_box_px + 2, str(n),
                                  color='cyan', fontsize=9, zorder=6)
                edit_mode[i]['circle'] = circ
                edit_mode[i]['label']  = lbl
                confirmed_markers_tgt.append([])   # placeholder — no red cross
            else:
                m3, = ax_tgt.plot(tx, ty, 'r+', ms=14, mew=2)
                m4 = ax_tgt.text(tx + 3, ty + 3, str(n), color='red', fontsize=9)
                confirmed_markers_tgt.append([m3, m4])

    def _enter_edit_mode(idx):
        """Toggle edit mode for pair idx."""
        if idx in edit_mode:
            # Toggle off: remove circle, restore red cross
            try:
                edit_mode[idx]['circle'].remove()
                edit_mode[idx]['label'].remove()
            except Exception:
                pass
            del edit_mode[idx]
            _redraw_all_confirmed()
            fig.canvas.draw_idle()
            return

        # Enter edit mode: hide red target cross, show cyan circle
        tx, ty = aligner.pairs[idx][4], aligner.pairs[idx][5]
        if idx < len(confirmed_markers_tgt):
            _remove_artists(confirmed_markers_tgt[idx])
            confirmed_markers_tgt[idx] = []

        circ = mpatches.Circle((tx, ty), radius=tgt_box_px,
                               edgecolor='cyan', facecolor='none',
                               linestyle='--', linewidth=1.5, zorder=5)
        ax_tgt.add_patch(circ)
        lbl = ax_tgt.text(tx + tgt_box_px + 2, ty + tgt_box_px + 2, str(idx + 1),
                          color='cyan', fontsize=9, zorder=6)
        edit_mode[idx] = {'circle': circ, 'label': lbl}
        fig.canvas.draw_idle()

    def _delete_pair(idx):
        """Delete pair at idx, clear all edit mode (indices shift)."""
        del aligner.pairs[idx]
        _clear_edit_mode()
        _redraw_all_confirmed()
        fig.canvas.draw_idle()

    def _set_title():
        if state['mode'] == 'IDLE':
            n_edit = len(edit_mode)
            edit_hint = f'  |  {n_edit} pair(s) in edit mode' if n_edit else ''
            fig.suptitle(
                'Left-click ref: add pair  |  double-click: edit/un-edit pair  |  '
                'u / right-click: delete  |  Enter: done' + edit_hint,
                fontsize=9)
        else:
            n = len(aligner.pairs) + 1
            fig.suptitle(
                f'Source {n} pending — click target to confirm  '
                '|  Space: accept  |  Right-click: cancel',
                fontsize=10)

    # ---- Scroll zoom ---------------------------------------------------------
    def _on_scroll(event):
        ax = event.inaxes
        if ax not in (ax_ref, ax_tgt) or event.xdata is None:
            return
        factor = 0.85 if event.button == 'up' else 1.15
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        cx, cy = event.xdata, event.ydata
        ax.set_xlim([cx + (x - cx) * factor for x in xlim])
        ax.set_ylim([cy + (y - cy) * factor for y in ylim])
        fig.canvas.draw_idle()

    # ---- Accept new pair (PENDING → IDLE) ------------------------------------
    def _accept_pair(tx, ty):
        rx, ry = pending['ref_x'], pending['ref_y']
        ra, dec = pending['ra'], pending['dec']
        aligner.pairs.append((rx, ry, ra, dec, tx, ty))
        if aligner.on_pair_callback is not None:
            aligner.on_pair_callback(rx, ry, tx, ty, ra, dec)
        for m in pending['marker_ref'] + pending['marker_tgt']:
            _remove_artists(m if isinstance(m, list) else [m])
        pending['marker_ref'].clear()
        pending['marker_tgt'].clear()
        _redraw_all_confirmed()
        state['mode'] = 'IDLE'
        _set_title()
        fig.canvas.draw_idle()

    # ---- Click handler -------------------------------------------------------
    def _on_click(event):
        try:
            if fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != '':
                return
        except Exception:
            pass
        if event.xdata is None or event.ydata is None:
            return
        if event.inaxes not in (ax_ref, ax_tgt):
            return

        cx, cy = event.xdata, event.ydata

        # ---- Double-click: toggle edit mode for nearest pair ----------------
        if event.dblclick and event.button == 1:
            if aligner.pairs and state['mode'] == 'IDLE':
                if event.inaxes is ax_ref:
                    dists = [np.sqrt((p[0]-cx)**2 + (p[1]-cy)**2)
                             for p in aligner.pairs]
                else:
                    dists = [np.sqrt((p[4]-cx)**2 + (p[5]-cy)**2)
                             for p in aligner.pairs]
                _enter_edit_mode(int(np.argmin(dists)))
                _set_title()
            return

        # ---- Right-click ------------------------------------------------
        if event.button == 3:
            if state['mode'] == 'PENDING':
                # Cancel pending
                for m in pending['marker_ref'] + pending['marker_tgt']:
                    _remove_artists(m if isinstance(m, list) else [m])
                pending['marker_ref'].clear()
                pending['marker_tgt'].clear()
                state['mode'] = 'IDLE'
                _set_title()
                fig.canvas.draw_idle()
                return
            # IDLE: delete nearest pair
            if not aligner.pairs:
                return
            if event.inaxes is ax_ref:
                dists = [np.sqrt((p[0]-cx)**2 + (p[1]-cy)**2) for p in aligner.pairs]
            else:
                dists = [np.sqrt((p[4]-cx)**2 + (p[5]-cy)**2) for p in aligner.pairs]
            _delete_pair(int(np.argmin(dists)))
            _set_title()
            return

        if event.button != 1:
            return

        # ---- IDLE left-click on target: re-place nearest edit circle --------
        if state['mode'] == 'IDLE' and event.inaxes is ax_tgt and edit_mode:
            dists = {idx: np.sqrt((info['circle'].center[0] - cx)**2 +
                                  (info['circle'].center[1] - cy)**2)
                     for idx, info in edit_mode.items()}
            nearest_idx = min(dists, key=dists.get)
            tx, ty = _refine_centroid(tgt_img, cx, cy, tgt_box_px)
            rx, ry, ra, dec, _, _ = aligner.pairs[nearest_idx]
            aligner.pairs[nearest_idx] = (rx, ry, ra, dec, tx, ty)
            # Remove artist BEFORE deleting from dict so _redraw_all_confirmed
            # doesn't encounter a missing key and leave a ghost circle on canvas
            info = edit_mode.pop(nearest_idx)
            try:
                info['circle'].remove()
                info['label'].remove()
            except Exception:
                pass
            _redraw_all_confirmed()
            _set_title()
            fig.canvas.draw_idle()
            return

        # ---- IDLE left-click on reference: start new pair -------------------
        if state['mode'] == 'IDLE' and event.inaxes is ax_ref:
            cx, cy = _refine_centroid(ref_img, cx, cy, ref_box_px)
            sky = ref_frame.wcs.all_pix2world([[cx, cy]], 0)
            ra, dec = float(sky[0, 0]), float(sky[0, 1])
            try:
                tgt_pred = tgt_frame.wcs.all_world2pix([[ra, dec]], 0)
                px, py = float(tgt_pred[0, 0]), float(tgt_pred[0, 1])
            except Exception:
                px, py = tgt_img.shape[1] / 2, tgt_img.shape[0] / 2

            pending['ref_x'], pending['ref_y'] = cx, cy
            pending['ra'], pending['dec'] = ra, dec
            pending['pred_x'], pending['pred_y'] = px, py

            m1, = ax_ref.plot(cx, cy, 'r+', ms=14, mew=2)
            pending['marker_ref'].clear()
            pending['marker_ref'].append(m1)

            circ = mpatches.Circle((px, py), radius=tgt_box_px,
                                   edgecolor='cyan', facecolor='none',
                                   linestyle='--', linewidth=1.5)
            ax_tgt.add_patch(circ)
            pending['marker_tgt'].clear()
            pending['marker_tgt'].append(circ)

            state['mode'] = 'PENDING'
            _set_title()
            fig.canvas.draw_idle()

        # ---- PENDING left-click on target: confirm pair ---------------------
        elif state['mode'] == 'PENDING' and event.inaxes is ax_tgt:
            tx, ty = _refine_centroid(tgt_img, cx, cy, tgt_box_px)
            _accept_pair(tx, ty)

    # ---- Key handler ---------------------------------------------------------
    def _on_key(event):
        if event.key == 'e':
            if not aligner.pairs or state['mode'] != 'IDLE':
                return
            if event.xdata is None or event.ydata is None:
                return
            cx, cy = event.xdata, event.ydata
            if event.inaxes is ax_ref:
                dists = [np.sqrt((p[0]-cx)**2 + (p[1]-cy)**2) for p in aligner.pairs]
            elif event.inaxes is ax_tgt:
                dists = [np.sqrt((p[4]-cx)**2 + (p[5]-cy)**2) for p in aligner.pairs]
            else:
                return
            _enter_edit_mode(int(np.argmin(dists)))
            _set_title()

        elif event.key == 'u':
            if not aligner.pairs:
                return
            if (event.inaxes in (ax_ref, ax_tgt)
                    and event.xdata is not None and event.ydata is not None):
                cx, cy = event.xdata, event.ydata
                if event.inaxes is ax_ref:
                    dists = [np.sqrt((p[0]-cx)**2 + (p[1]-cy)**2)
                             for p in aligner.pairs]
                else:
                    dists = [np.sqrt((p[4]-cx)**2 + (p[5]-cy)**2)
                             for p in aligner.pairs]
                _delete_pair(int(np.argmin(dists)))
            else:
                _delete_pair(len(aligner.pairs) - 1)
            _set_title()

        elif event.key == ' ':
            if state['mode'] == 'PENDING':
                tx, ty = _refine_centroid(tgt_img,
                                          pending['pred_x'], pending['pred_y'], tgt_box_px)
                _accept_pair(tx, ty)

        elif event.key == 'enter':
            _clear_edit_mode()
            plt.close(fig)

    # ---- Window close --------------------------------------------------------
    def _on_close(event):
        pass

    # ---- Connect events -------------------------------------------------------
    fig.canvas.mpl_connect('button_press_event', _on_click)
    fig.canvas.mpl_connect('key_press_event', _on_key)
    fig.canvas.mpl_connect('close_event', _on_close)
    fig.canvas.mpl_connect('scroll_event', _on_scroll)

    # Draw any pairs already in aligner.pairs (e.g. loaded by batch)
    if aligner.pairs:
        _redraw_all_confirmed()

    _set_title()
    plt.tight_layout()
    plt.show()

    # ---- Save catalog if requested ------------------------------------------
    if save_catalog is not None:
        _save_pairs_catalog(aligner.pairs, save_catalog)


# ---------------------------------------------------------------------------
# Batch strategy
# ---------------------------------------------------------------------------

def _batch(aligner, catalog: str, box: float = 0.1,
           stretch: str = 'zscale', save_catalog=None, headless: bool = False):
    """
    Load a saved pairs catalog, recentroid on current reference+target, then
    open the interactive window so the user can inspect and edit any pairs.

    The catalog's (ra, dec) sky positions are re-projected onto whatever reference
    and target are currently loaded — so the catalog is reusable across different
    reference files (e.g. HST vs KCWI, or two different HST pointings).
    Stored x_ref/y_ref values are ignored; positions are always recomputed from WCS.

    After loading, the interactive window opens with all pairs shown as numbered
    red markers.  Press 'e' near a target marker to toggle edit mode (cyan circle);
    left-click to re-place it.  Right-click or 'u' to delete.  Enter to finish.
    """
    from astropy.table import Table

    tbl = Table.read(catalog)
    ras  = np.asarray(tbl['ra'])
    decs = np.asarray(tbl['dec'])

    if not aligner.targets:
        raise RuntimeError("No targets loaded.")

    ref = aligner.reference
    tgt = aligner.targets[0]

    ref_img = ref.image2d
    tgt_img = tgt.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError("Call preprocess() before find_sources(strategy='batch').")

    ref_ny, ref_nx = ref_img.shape
    tgt_ny, tgt_nx = tgt_img.shape

    # Convert box from arcsec to pixels for each frame independently
    ref_box_px = _box_to_pixels(box, ref, ref_img)
    tgt_box_px = _box_to_pixels(box, tgt, tgt_img)

    aligner.pairs.clear()
    for ra, dec in zip(ras, decs):
        # Re-project onto current reference WCS
        try:
            ref_pix = ref.wcs.all_world2pix([[ra, dec]], 0)
            ref_px, ref_py = float(ref_pix[0, 0]), float(ref_pix[0, 1])
        except Exception:
            continue
        if not (0 <= ref_px < ref_nx and 0 <= ref_py < ref_ny):
            continue
        rx, ry = _refine_centroid(ref_img, ref_px, ref_py, ref_box_px)

        # Re-project onto current target WCS
        try:
            tgt_pix = tgt.wcs.all_world2pix([[ra, dec]], 0)
            tgt_px, tgt_py = float(tgt_pix[0, 0]), float(tgt_pix[0, 1])
        except Exception:
            continue
        if not (0 <= tgt_px < tgt_nx and 0 <= tgt_py < tgt_ny):
            continue
        tx, ty = _refine_centroid(tgt_img, tgt_px, tgt_py, tgt_box_px)

        aligner.pairs.append((rx, ry, ra, dec, tx, ty))

    print(f"[rb_align] batch: loaded {len(aligner.pairs)} pair(s) from '{catalog}'.")

    # Open interactive window for inspection / editing (skipped when called from GUI)
    if not headless:
        _interactive(aligner, stretch=stretch, box=box, save_catalog=save_catalog)


# ---------------------------------------------------------------------------
# Gaia strategy
# ---------------------------------------------------------------------------

def _gaia(aligner, box: float = 0.1, radius_deg: float = 0.5,
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

    # Convert box from arcsec to pixels for each frame
    ref_box_px = _box_to_pixels(box, ref, img)
    tgt_box_px = _box_to_pixels(box, tgt, tgt_img)

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
        rx_r, ry_r = _refine_centroid(img, rx, ry, ref_box_px)
        # Project to target
        try:
            t_pix = tgt.wcs.all_world2pix([[ra, dec]], 0)
            tx, ty = float(t_pix[0, 0]), float(t_pix[0, 1])
        except Exception:
            continue
        if not (0 <= tx < t_nx and 0 <= ty < t_ny):
            continue
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, tgt_box_px)
        aligner.pairs.append((rx_r, ry_r, ra, dec, tx_r, ty_r))

    if not aligner.pairs:
        raise RuntimeError("Gaia sources found but none projected into both images.")


# ---------------------------------------------------------------------------
# DAO strategy
# ---------------------------------------------------------------------------

def _dao(aligner, box: float = 0.1, fwhm: float = 0.0,
         threshold_sigma: float = 3.0, max_sources: int = 50,
         bg_method: str = 'sigma_clip'):
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

    # Convert box from arcsec to pixels for each frame
    ref_box_px = _box_to_pixels(box, ref, ref_img)
    tgt_box_px = _box_to_pixels(box, tgt, tgt_img)

    ref_sources = _detect_dao_sources(ref_img, ref_box_px, fwhm=fwhm,
                                       threshold_sigma=threshold_sigma,
                                       bg_method=bg_method)
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
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, tgt_box_px)
        aligner.pairs.append((rx, ry, ra, dec, tx_r, ty_r))


def _bg_stats(finite: np.ndarray, bg_method: str = 'sigma_clip'):
    """
    Estimate background median and noise std from a 1-D array of finite values.

    Parameters
    ----------
    bg_method : 'sigma_clip' | 'mad'
        'sigma_clip' — iterative sigma-clipped mean/std (astropy).
        'mad'        — median + 1.4826 × MAD (more robust for IFU data with
                       correlated noise, non-Gaussian backgrounds, or many sources).
    Returns
    -------
    bg_median, bg_std : float
    """
    if bg_method == 'mad':
        from astropy.stats import mad_std
        bg_median = float(np.median(finite))
        bg_std    = float(mad_std(finite))
    else:
        from astropy.stats import sigma_clipped_stats
        _, bg_median, bg_std = sigma_clipped_stats(finite)
        bg_median = float(bg_median)
        bg_std    = float(bg_std)
    return bg_median, bg_std


def _detect_dao_sources(img: np.ndarray, box: int,
                        fwhm: float = 0.0,
                        threshold_sigma: float = 3.0,
                        bg_method: str = 'sigma_clip'):
    """
    Return list of (x, y) source positions from DAOStarFinder or fallback.

    Parameters
    ----------
    fwhm : float
        Expected source FWHM in pixels.  0 → auto as ``box * 0.5``.
    threshold_sigma : float
        Detection threshold in units of background σ.
    bg_method : 'sigma_clip' | 'mad'
        Background noise estimator.  'mad' is more robust for IFU whitelight
        images with correlated noise or a high source fraction.
    """
    if fwhm <= 0:
        fwhm = max(2.0, box * 0.5)

    finite = img[np.isfinite(img)]
    bg_median, bg_std = _bg_stats(finite, bg_method)

    if HAS_PHOTUTILS:
        try:
            from photutils.detection import DAOStarFinder
            finder = DAOStarFinder(fwhm=fwhm, threshold=threshold_sigma * bg_std)
            sources = finder(img - bg_median)
            if sources is None or len(sources) == 0:
                return []
            return list(zip(
                np.asarray(sources['xcentroid']),
                np.asarray(sources['ycentroid'])
            ))
        except Exception:
            pass
    # Fallback: simple peak finding via scipy
    return _simple_peak_find(img, box, threshold_sigma=threshold_sigma,
                              bg_median=bg_median, bg_std=bg_std)


def _simple_peak_find(img: np.ndarray, box: int, n_max: int = 30,
                      threshold_sigma: float = 3.0,
                      bg_method: str = 'sigma_clip',
                      bg_median: float = None,
                      bg_std: float = None):
    """
    Basic peak finding: local maxima above threshold_sigma × σ_background.

    bg_median / bg_std may be passed in pre-computed (from _detect_dao_sources)
    to avoid re-estimating.  If None, estimated here using bg_method.
    """
    from scipy.ndimage import maximum_filter, label
    clean = np.nan_to_num(img, nan=0.0)

    if bg_median is None or bg_std is None:
        finite = clean.ravel()
        finite = finite[np.isfinite(finite)]
        bg_median, bg_std = _bg_stats(finite, bg_method)

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

def _knots(aligner, box: float = 0.1, threshold_sigma: float = 2.0,
           max_knots: int = 10, bg_method: str = 'sigma_clip'):
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

    # Convert box from arcsec to pixels for each frame
    ref_box_px = _box_to_pixels(box, ref, ref_img)
    tgt_box_px = _box_to_pixels(box, tgt, tgt_img)

    knot_positions = _detect_knots(ref_img, ref_box_px, n_knots=max_knots,
                                    threshold_sigma=threshold_sigma,
                                    bg_method=bg_method)
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
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, tgt_box_px)
        aligner.pairs.append((rx, ry, ra, dec, tx_r, ty_r))


def _detect_knots(img: np.ndarray, box: int, n_knots: int = 10,
                  threshold_sigma: float = 2.0,
                  bg_method: str = 'sigma_clip'):
    """Segment image and return centroids of brightest segments."""
    finite = img[np.isfinite(img)]
    bg_median, bg_std = _bg_stats(finite, bg_method)

    if HAS_PHOTUTILS:
        try:
            from photutils.segmentation import detect_sources, SourceCatalog
            seg = detect_sources(np.nan_to_num(img, nan=0.0),
                                  threshold=bg_median + threshold_sigma * bg_std,
                                  npixels=5)
            if seg is None:
                return _simple_peak_find(img, box, n_max=n_knots,
                                          threshold_sigma=threshold_sigma,
                                          bg_median=bg_median, bg_std=bg_std)
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
                              threshold_sigma=threshold_sigma,
                              bg_median=bg_median, bg_std=bg_std)


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


def _box_to_pixels(box_arcsec: float, frame, img: np.ndarray,
                   minimum: int = 3) -> int:
    """
    Convert a centroid box size from arcsec to pixels for a given frame.

    Parameters
    ----------
    box_arcsec : float
        Half-width of the centroid box in arcsec.
    frame : Frame
        The frame whose WCS pixel scale is used.
    img : np.ndarray
        The 2D image (used for pixel scale estimation).
    minimum : int
        Minimum box size in pixels (default 3).

    Returns
    -------
    int
        Box half-width in pixels, at least ``minimum``.
    """
    scale = _estimate_pixel_scale_arcsec(frame, img)   # arcsec/px
    if scale <= 0:
        scale = 0.3
    px = max(minimum, int(round(box_arcsec / scale)))
    return px


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


def _cross_corr(aligner, box: float = 0.1, fwhm: float = 0.0,
                threshold_sigma: float = 3.0, max_sources: int = 50,
                bg_method: str = 'sigma_clip'):
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
        Centroid box half-width in arcsec; also used directly as the Hough
        bin size (clipped to 1.5–10 arcsec).
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

    # Convert box from arcsec to pixels for each frame independently
    ref_box_px = _box_to_pixels(box, ref, ref_img)
    tgt_box_px = _box_to_pixels(box, tgt, tgt_img)

    # ------------------------------------------------------------------
    # Step 1: source detection — keep up to 50 brightest per image
    # ------------------------------------------------------------------
    ref_sources = _detect_dao_sources(ref_img, ref_box_px, fwhm=fwhm,
                                       threshold_sigma=threshold_sigma,
                                       bg_method=bg_method)
    tgt_sources = _detect_dao_sources(tgt_img, tgt_box_px, fwhm=fwhm,
                                       threshold_sigma=threshold_sigma,
                                       bg_method=bg_method)

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

    # box is already in arcsec — use directly as the Hough bin size
    bin_arcsec = float(np.clip(box, 1.5, 10.0))

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

        rx_r, ry_r = _refine_centroid(ref_img, rx, ry, ref_box_px)
        tx_r, ty_r = _refine_centroid(tgt_img, tx, ty, tgt_box_px)

        sky = ref.wcs.all_pix2world([[rx_r, ry_r]], 0)
        ra, dec = float(sky[0, 0]), float(sky[0, 1])

        aligner.pairs.append((rx_r, ry_r, ra, dec, tx_r, ty_r))

    if not aligner.pairs:
        print("[rb_align] cross_corr: no source matches after Hough vote; "
              "falling back to image FFT cross-correlation (1 pair).")
        _cross_corr_image_fallback(aligner)


# ---------------------------------------------------------------------------
# Arc NCC strategy — NCC template matching on a user-selected region
# ---------------------------------------------------------------------------

def _apply_log(img: np.ndarray, sigma: float) -> np.ndarray:
    """
    Negated Laplacian-of-Gaussian filter.

    Suppresses smooth/diffuse emission (seeing halos, continuum gradients)
    and enhances compact structures (arc knots, emission peaks).
    The negation makes brightness peaks positive.

    Parameters
    ----------
    img : np.ndarray
        2D image array.
    sigma : float
        Gaussian smoothing scale in pixels.  Match to expected knot size.

    Returns
    -------
    np.ndarray
        LoG-filtered image, same shape as input.
    """
    from scipy.ndimage import gaussian_laplace
    return -gaussian_laplace(np.nan_to_num(img, nan=0.0), sigma=sigma)


def _ncc_subpixel_peak(surface: np.ndarray):
    """
    Locate the peak of a 2D NCC surface with sub-pixel precision.

    Uses independent 1-D parabolic interpolation along each axis through
    the integer-pixel maximum.

    Parameters
    ----------
    surface : np.ndarray
        2D NCC surface.

    Returns
    -------
    (sub_row, sub_col) : (float, float)
        Sub-pixel peak position in array coordinates.
    """
    peak_idx = np.unravel_index(np.argmax(surface), surface.shape)
    py, px = int(peak_idx[0]), int(peak_idx[1])
    ny, nx = surface.shape

    # Parabolic sub-pixel in y
    if 0 < py < ny - 1:
        vy = (surface[py - 1, px], surface[py, px], surface[py + 1, px])
        denom = vy[0] - 2.0 * vy[1] + vy[2]
        sub_y = py - 0.5 * (vy[0] - vy[2]) / denom if abs(denom) > 1e-10 else float(py)
    else:
        sub_y = float(py)

    # Parabolic sub-pixel in x
    if 0 < px < nx - 1:
        vx = (surface[py, px - 1], surface[py, px], surface[py, px + 1])
        denom = vx[0] - 2.0 * vx[1] + vx[2]
        sub_x = px - 0.5 * (vx[0] - vx[2]) / denom if abs(denom) > 1e-10 else float(px)
    else:
        sub_x = float(px)

    return sub_y, sub_x


def _arc_ncc_interactive(aligner,
                         search_radius: float = 5.0,
                         use_LoG: bool = False,
                         LoG_sigma: float = 0.0,
                         bg_method: str = 'mad',
                         stretch: str = 'zscale'):
    """
    Interactive version of arc_ncc: click-drag on the reference panel to draw
    a bounding box, then NCC runs automatically on mouse release.

    Controls
    --------
    Click-drag on reference  — draw box, run NCC immediately on release
    r                        — clear result and draw a new box
    Enter                    — accept current result and close
    Escape / close window    — close without saving any pair

    The NCC result is shown as:
      - a cyan rectangle on the reference panel (selected box)
      - a red cross on the target panel (matched position)
      - a third panel showing the NCC surface (peak = best shift)
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.gridspec as gridspec

    ref_frame = aligner.reference
    tgt_frame = aligner.targets[0]
    ref_img   = ref_frame.image2d
    tgt_img   = tgt_frame.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError(
            "Call preprocess() before find_sources(strategy='arc_ncc').")

    ref_norm = _compute_norm(ref_img, stretch)
    tgt_norm = _compute_norm(tgt_img, stretch)

    # ---- Figure: reference | target | NCC surface ---------------------------
    fig = plt.figure(figsize=(15, 5))
    gs  = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.35)
    ax_ref = fig.add_subplot(gs[0])
    ax_tgt = fig.add_subplot(gs[1])
    ax_ncc = fig.add_subplot(gs[2])

    ax_ref.imshow(ref_img, origin='lower', cmap='gray',
                  norm=ref_norm, interpolation='nearest')
    ax_tgt.imshow(tgt_img, origin='lower', cmap='gray',
                  norm=tgt_norm, interpolation='nearest')
    ax_ref.set_title('Reference — click-drag to select region')
    ax_tgt.set_title('Target — NCC match')
    ax_ncc.set_title('NCC surface')
    ax_ncc.set_visible(False)

    state = {
        'press_xy':            None,   # (x, y) where mouse button went down on ref
        'rect':                None,   # cyan box on reference panel
        'tgt_rect':            None,   # cyan box on target panel (same size, at match)
        'tgt_mark':            None,   # red cross on target panel
        'ncc_im':              None,   # imshow on ax_ncc
        'result':              None,   # (bx0, by0, bx1, by1) if a box has been drawn
        'edit_mode':           False,  # True after double-click on target
        'seed_xy':             None,   # (x, y) manual search seed in target (or None → WCS)
        'ignore_next_release': False,  # swallow the release that follows a double-click
    }

    def _clear_artists():
        for key in ('rect', 'tgt_rect', 'tgt_mark', 'ncc_im'):
            if state[key] is not None:
                try:
                    state[key].remove()
                except Exception:
                    pass
                state[key] = None
        ax_ncc.cla()
        ax_ncc.set_visible(False)
        state['result']    = None
        state['edit_mode'] = False
        state['seed_xy']   = None

    def _set_title():
        if state['edit_mode']:
            fig.suptitle(
                'Edit mode — left-click on target to set new search centre',
                fontsize=10, color='cyan')
        elif state['result'] is not None:
            bx0, by0, bx1, by1 = state['result']
            if aligner.pairs:
                rx, ry, ra, dec, tgt_cx, tgt_cy = aligner.pairs[0]
                tgt_scale = _estimate_pixel_scale_arcsec(tgt_frame, tgt_img)
                try:
                    pred = tgt_frame.wcs.all_world2pix([[ra, dec]], 0)
                    shift_dx = tgt_cx - float(pred[0, 0])
                    shift_dy = tgt_cy - float(pred[0, 1])
                except Exception:
                    shift_dx = shift_dy = 0.0
                shift_as = np.sqrt((shift_dx * tgt_scale)**2 +
                                   (shift_dy * tgt_scale)**2)
                seed_note = '  [manual seed]' if state['seed_xy'] else ''
                fig.suptitle(
                    f'NCC shift = ({shift_dx:+.2f}, {shift_dy:+.2f}) px'
                    f' = {shift_as:.2f}"{seed_note}'
                    '  |  double-click target: re-seed  |  r: redo  |  Enter: accept',
                    fontsize=9)
            else:
                fig.suptitle('NCC failed.  r: redo', fontsize=10, color='red')
        else:
            fig.suptitle(
                'Click-drag on reference to select a region'
                '  |  r: redo  |  Enter: accept  |  Esc: cancel',
                fontsize=10)

    def _run_ncc(bx0, by0, bx1, by1, seed_xy=None):
        """Run NCC for the given reference box, optionally with a manual seed."""
        # Remove old NCC result artists but keep the ref box if already drawn
        for key in ('tgt_rect', 'tgt_mark', 'ncc_im'):
            if state[key] is not None:
                try:
                    state[key].remove()
                except Exception:
                    pass
                state[key] = None
        ax_ncc.cla()
        ax_ncc.set_visible(False)
        state['edit_mode'] = False
        aligner.pairs.clear()

        try:
            _arc_ncc(aligner, box_ref=(bx0, by0, bx1, by1),
                     search_radius=search_radius,
                     use_LoG=use_LoG, LoG_sigma=LoG_sigma, bg_method=bg_method,
                     _seed_xy=seed_xy)
        except Exception as e:
            fig.suptitle(f'NCC failed: {e}', fontsize=9, color='red')
            fig.canvas.draw_idle()
            return

        if not aligner.pairs:
            fig.suptitle('NCC produced no result.', fontsize=9, color='red')
            fig.canvas.draw_idle()
            return

        _, _, _, _, tgt_cx, tgt_cy = aligner.pairs[0]

        # Cyan box in target — same pixel size as the reference box
        box_w = bx1 - bx0
        box_h = by1 - by0
        tgt_rect = mpatches.Rectangle(
            (tgt_cx - box_w / 2.0, tgt_cy - box_h / 2.0), box_w, box_h,
            edgecolor='cyan', facecolor='none', linewidth=1.5, linestyle='--')
        ax_tgt.add_patch(tgt_rect)
        state['tgt_rect'] = tgt_rect

        mark, = ax_tgt.plot(tgt_cx, tgt_cy, 'r+', ms=18, mew=2.5)
        state['tgt_mark'] = mark

        _show_ncc_surface(ax_ncc, ref_img, tgt_img, ref_frame, tgt_frame,
                          bx0, by0, bx1, by1,
                          search_radius, use_LoG, LoG_sigma, bg_method,
                          seed_xy=seed_xy)
        ax_ncc.set_visible(True)
        _set_title()
        fig.canvas.draw_idle()

    # ---- Mouse handlers -----------------------------------------------------
    def _on_press(event):
        if fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != '':
            return
        if event.xdata is None or event.ydata is None:
            return

        # Start a drag on the reference panel (single left-click press)
        if event.inaxes is ax_ref and event.button == 1 and not event.dblclick:
            state['press_xy'] = (event.xdata, event.ydata)

    def _on_release(event):
        if fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != '':
            return
        if event.xdata is None or event.ydata is None:
            return

        # Swallow the release that immediately follows a double-click press
        if state['ignore_next_release']:
            state['ignore_next_release'] = False
            return

        # ---- Edit-mode: left-click on target sets new search seed -----------
        if (state['edit_mode'] and event.inaxes is ax_tgt
                and event.button == 1 and not event.dblclick):
            seed_x, seed_y = event.xdata, event.ydata
            state['seed_xy'] = (seed_x, seed_y)
            state['edit_mode'] = False
            bx0, by0, bx1, by1 = state['result']
            _run_ncc(bx0, by0, bx1, by1, seed_xy=(seed_x, seed_y))
            return

        # ---- Normal drag on reference: draw box and run NCC -----------------
        if event.inaxes is not ax_ref or event.button != 1:
            state['press_xy'] = None
            return
        if state['press_xy'] is None:
            return

        x0, y0 = state['press_xy']
        x1, y1 = event.xdata, event.ydata
        state['press_xy'] = None

        bx0, bx1 = sorted([int(round(x0)), int(round(x1))])
        by0, by1 = sorted([int(round(y0)), int(round(y1))])
        if bx1 - bx0 < 3 or by1 - by0 < 3:
            fig.suptitle('Box too small — try again.', fontsize=10, color='red')
            fig.canvas.draw_idle()
            return

        # Remove previous artists, draw new reference box
        _clear_artists()
        rect = mpatches.Rectangle(
            (bx0, by0), bx1 - bx0, by1 - by0,
            edgecolor='cyan', facecolor='none', linewidth=1.5, linestyle='--')
        ax_ref.add_patch(rect)
        state['rect']   = rect
        state['result'] = (bx0, by0, bx1, by1)

        _run_ncc(bx0, by0, bx1, by1)

    def _on_click(event):
        """Separate handler for double-click on target → enter edit mode."""
        if fig.canvas.toolbar is not None and fig.canvas.toolbar.mode != '':
            return
        if not event.dblclick or event.button != 1:
            return
        if event.inaxes is not ax_tgt:
            return
        if state['result'] is None or not aligner.pairs:
            return

        # Always swallow the release that pairs with this double-click press
        state['ignore_next_release'] = True

        # Toggle edit mode — no visual decoration, just change the title
        state['edit_mode'] = not state['edit_mode']
        _set_title()
        fig.canvas.draw_idle()

    def _on_key(event):
        if event.key == 'r':
            _clear_artists()
            aligner.pairs.clear()
            ax_ref.set_title('Reference — click-drag to select region')
            _set_title()
            fig.canvas.draw_idle()
        elif event.key == 'enter':
            plt.close(fig)
        elif event.key == 'escape':
            aligner.pairs.clear()
            plt.close(fig)

    def _on_scroll(event):
        ax = event.inaxes
        if ax not in (ax_ref, ax_tgt) or event.xdata is None:
            return
        factor = 0.85 if event.button == 'up' else 1.15
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        cx, cy = event.xdata, event.ydata
        ax.set_xlim([cx + (x - cx) * factor for x in xlim])
        ax.set_ylim([cy + (y - cy) * factor for y in ylim])
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('button_press_event',   _on_press)
    fig.canvas.mpl_connect('button_release_event', _on_release)
    fig.canvas.mpl_connect('button_press_event',   _on_click)
    fig.canvas.mpl_connect('key_press_event',      _on_key)
    fig.canvas.mpl_connect('scroll_event',         _on_scroll)

    _set_title()
    plt.tight_layout()
    plt.show()


def _show_ncc_surface(ax, ref_img, tgt_img, ref_frame, tgt_frame,
                      x0, y0, x1, y1,
                      search_radius, use_LoG, LoG_sigma, bg_method,
                      seed_xy=None):
    """
    Recompute the NCC surface for a given box and display it in ax.
    Used only for visualisation inside _arc_ncc_interactive.
    """
    from scipy.signal import fftconvolve

    template = np.nan_to_num(ref_img[y0:y1, x0:x1], nan=0.0).copy()
    tmpl_h, tmpl_w = template.shape
    ref_cx = (x0 + x1) / 2.0
    ref_cy = (y0 + y1) / 2.0

    if seed_xy is not None:
        pred_tx, pred_ty = float(seed_xy[0]), float(seed_xy[1])
    else:
        try:
            sky = ref_frame.wcs.all_pix2world([[ref_cx, ref_cy]], 0)
            ra, dec = float(sky[0, 0]), float(sky[0, 1])
            pred = tgt_frame.wcs.all_world2pix([[ra, dec]], 0)
            pred_tx, pred_ty = float(pred[0, 0]), float(pred[0, 1])
        except Exception:
            t_ny, t_nx = tgt_img.shape
            pred_tx, pred_ty = t_nx / 2.0, t_ny / 2.0

    ref_scale = _estimate_pixel_scale_arcsec(ref_frame, ref_img)
    tgt_scale = _estimate_pixel_scale_arcsec(tgt_frame, tgt_img)
    zoom_factor = tgt_scale / ref_scale if ref_scale > 0 else 1.0
    search_px   = max(5, int(np.ceil(search_radius / tgt_scale)))

    t_ny, t_nx = tgt_img.shape
    half_w = (tmpl_w / zoom_factor) / 2.0 + search_px
    half_h = (tmpl_h / zoom_factor) / 2.0 + search_px
    sx0 = max(0, int(np.floor(pred_tx - half_w)))
    sx1 = min(t_nx, int(np.ceil(pred_tx + half_w)))
    sy0 = max(0, int(np.floor(pred_ty - half_h)))
    sy1 = min(t_ny, int(np.ceil(pred_ty + half_h)))

    search_raw = np.nan_to_num(tgt_img[sy0:sy1, sx0:sx1], nan=0.0).copy()

    if abs(zoom_factor - 1.0) > 0.05:
        from scipy.ndimage import zoom as _zoom
        search_img = _zoom(search_raw, zoom_factor, order=3)
    else:
        search_img = search_raw
        zoom_factor = 1.0

    if use_LoG:
        sigma = LoG_sigma if LoG_sigma > 0 else 1.5
        template   = _apply_log(template,   sigma)
        search_img = _apply_log(search_img, sigma)

    def _sub_bg(arr):
        fin = arr[np.isfinite(arr)].ravel()
        if len(fin) == 0:
            return arr
        bg, _ = _bg_stats(fin, bg_method)
        return arr - bg

    template   = _sub_bg(template)
    search_img = _sub_bg(search_img)

    tmpl_std = template.std()
    if tmpl_std < 1e-10:
        return
    tmpl_norm = template / (tmpl_std * template.size)
    ncc = fftconvolve(search_img, tmpl_norm[::-1, ::-1], mode='valid')
    search_std = search_img.std()
    if search_std > 1e-10:
        ncc = ncc / search_std

    if ncc.size == 0:
        return

    peak_row, peak_col = _ncc_subpixel_peak(ncc)

    ax.cla()
    im = ax.imshow(ncc, origin='lower', cmap='hot', interpolation='nearest')
    ax.plot(peak_col, peak_row, 'c+', ms=14, mew=2)
    ax.set_title(f'NCC surface  (peak={ncc.max():.3f})')
    ax.set_xlabel('Δx (px in search)')
    ax.set_ylabel('Δy (px in search)')
    try:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    except Exception:
        pass


def _arc_ncc(aligner,
             box_ref,
             search_radius: float = 5.0,
             use_LoG: bool = False,
             LoG_sigma: float = 0.0,
             bg_method: str = 'mad',
             _seed_xy=None):
    """
    Align two IFU images by NCC template matching on a user-selected region.

    The user specifies a bounding box in the reference image that contains a
    compact feature common to both observations (a bright arc knot, an emission
    peak, a foreground star).  The function extracts that cutout as a template,
    builds a search region in the target centred on the WCS-predicted position,
    optionally reprojects to a common pixel scale, and finds the shift that
    maximises the Normalised Cross-Correlation (NCC) between the two cutouts.
    A single pair is stored in ``aligner.pairs``.

    Unlike the ``cross_corr`` strategy, this requires *no* point-source
    detection — it works directly on extended emission.

    Parameters
    ----------
    box_ref : (x0, y0, x1, y1)
        Bounding box in **reference** pixel coordinates (floats are rounded).
        Can be given in any corner order; will be sorted internally.
    search_radius : float
        Half-width of the search region around the WCS-predicted target
        position, in arcsec.  Default 5.0".  Increase if the WCS offset is
        known to be large.
    use_LoG : bool
        Apply a Laplacian-of-Gaussian (LoG) filter to both cutouts before
        computing the NCC.  Suppresses seeing halos and broad backgrounds;
        useful when the two datasets have noticeably different seeing.
        Default False.
    LoG_sigma : float
        LoG smoothing scale in pixels.  0 (default) → auto set to 1.5 px.
    bg_method : str
        Background estimator: ``'mad'`` (default, robust) or
        ``'sigma_clip'``.

    Notes
    -----
    * If the two frames have pixel scales that differ by more than 5 %, the
      target search region is resampled to the reference pixel scale before
      NCC so that the template and search image are in the same pixel units.
    * The NCC peak is refined to sub-pixel precision by 1-D parabolic
      interpolation along each axis.
    * Only one pair is produced, representing the centre of the selected box
      matched to its best-fit position in the target.  Pass this to
      ``aligner.compute_offset()`` (translation) or combine with additional
      pairs from other strategies for a full WCS solution.

    Examples
    --------
    >>> aligner.find_sources(strategy='arc_ncc',
    ...                      box_ref=(45, 30, 75, 60),
    ...                      search_radius=8.0,
    ...                      use_LoG=True)
    """
    from scipy.signal import fftconvolve

    ref = aligner.reference
    tgt = aligner.targets[0]
    ref_img = ref.image2d
    tgt_img = tgt.image2d
    if ref_img is None or tgt_img is None:
        raise RuntimeError(
            "Call preprocess() before find_sources(strategy='arc_ncc').")

    # --- Parse and clamp bounding box ---
    x0, y0, x1, y1 = [int(round(v)) for v in box_ref]
    x0, x1 = sorted([x0, x1])
    y0, y1 = sorted([y0, y1])
    r_ny, r_nx = ref_img.shape
    x0 = max(0, x0);  x1 = min(r_nx, x1)
    y0 = max(0, y0);  y1 = min(r_ny, y1)
    if x1 - x0 < 3 or y1 - y0 < 3:
        raise ValueError(
            f"box_ref ({x0},{y0})→({x1},{y1}) is too small (< 3 px on a side).")

    template = np.nan_to_num(ref_img[y0:y1, x0:x1], nan=0.0).copy()
    tmpl_h, tmpl_w = template.shape

    # --- Reference box centre → sky ---
    ref_cx = (x0 + x1) / 2.0
    ref_cy = (y0 + y1) / 2.0
    sky = ref.wcs.all_pix2world([[ref_cx, ref_cy]], 0)
    ra, dec = float(sky[0, 0]), float(sky[0, 1])

    # --- Search centre in target: manual seed overrides WCS prediction ---
    if _seed_xy is not None:
        pred_tx, pred_ty = float(_seed_xy[0]), float(_seed_xy[1])
        print(f"[rb_align] arc_ncc: using manual seed "
              f"({pred_tx:.1f}, {pred_ty:.1f}) as search centre.")
    else:
        try:
            pred = tgt.wcs.all_world2pix([[ra, dec]], 0)
            pred_tx, pred_ty = float(pred[0, 0]), float(pred[0, 1])
        except Exception:
            t_ny0, t_nx0 = tgt_img.shape
            pred_tx, pred_ty = t_nx0 / 2.0, t_ny0 / 2.0
            print("[rb_align] arc_ncc: WCS projection failed; "
                  "centring search on target image centre.")

    # --- Pixel scales (arcsec/px) ---
    ref_scale = _estimate_pixel_scale_arcsec(ref, ref_img)
    tgt_scale = _estimate_pixel_scale_arcsec(tgt, tgt_img)

    # zoom_factor: rescale target → reference pixel scale
    # zoom > 1 if target pixels are larger (coarser) than reference
    zoom_factor = tgt_scale / ref_scale if ref_scale > 0 else 1.0

    # search_radius in original target pixels
    search_px = max(5, int(np.ceil(search_radius / tgt_scale)))

    t_ny, t_nx = tgt_img.shape

    # --- Extract search region from target (original pixel coords) ---
    # After zoom_factor rescaling, the template occupies (tmpl_h, tmpl_w) pixels.
    # We need the rescaled search image to be at least (tmpl + 2*search_px) on
    # each side, which in original target pixels is (tmpl/zoom + 2*search_px).
    half_w_tgt = (tmpl_w / zoom_factor) / 2.0 + search_px
    half_h_tgt = (tmpl_h / zoom_factor) / 2.0 + search_px

    sx0 = max(0, int(np.floor(pred_tx - half_w_tgt)))
    sx1 = min(t_nx, int(np.ceil(pred_tx  + half_w_tgt)))
    sy0 = max(0, int(np.floor(pred_ty - half_h_tgt)))
    sy1 = min(t_ny, int(np.ceil(pred_ty  + half_h_tgt)))

    if sx1 - sx0 < 3 or sy1 - sy0 < 3:
        raise RuntimeError(
            "Target search region is entirely outside the image. "
            "Check WCS headers or increase search_radius.")

    search_raw = np.nan_to_num(tgt_img[sy0:sy1, sx0:sx1], nan=0.0).copy()

    # --- Resample target cutout to reference pixel scale if needed ---
    if abs(zoom_factor - 1.0) > 0.05:
        from scipy.ndimage import zoom as _zoom
        search_img = _zoom(search_raw, zoom_factor, order=3)
        _rescaled = True
        print(f"[rb_align] arc_ncc: resampling target cutout "
              f"(tgt {tgt_scale:.3f}\"/px → ref {ref_scale:.3f}\"/px, "
              f"zoom={zoom_factor:.3f})")
    else:
        search_img = search_raw
        zoom_factor = 1.0
        _rescaled = False

    sh, sw = search_img.shape
    if sh < tmpl_h or sw < tmpl_w:
        raise RuntimeError(
            f"Rescaled search region ({sw}×{sh} px) is smaller than the "
            f"template ({tmpl_w}×{tmpl_h} px). Increase search_radius.")

    # --- Optional LoG preprocessing ---
    if use_LoG:
        sigma = LoG_sigma if LoG_sigma > 0 else 1.5
        template   = _apply_log(template,   sigma)
        search_img = _apply_log(search_img, sigma)
        print(f"[rb_align] arc_ncc: LoG preprocessing (sigma={sigma:.1f} px)")

    # --- Background subtraction ---
    def _sub_bg(arr):
        fin = arr[np.isfinite(arr)].ravel()
        if len(fin) == 0:
            return arr
        bg, _ = _bg_stats(fin, bg_method)
        return arr - bg

    template   = _sub_bg(template)
    search_img = _sub_bg(search_img)

    # --- Normalised Cross-Correlation via FFT (valid mode = template matching) ---
    # valid mode: only positions where template fully overlaps the search image
    # NCC output shape = (sh - tmpl_h + 1, sw - tmpl_w + 1)
    tmpl_std = template.std()
    if tmpl_std < 1e-10:
        raise RuntimeError(
            "Reference cutout has no variance — "
            "choose a box that contains signal.")
    tmpl_norm = template / (tmpl_std * template.size)

    ncc = fftconvolve(search_img, tmpl_norm[::-1, ::-1], mode='valid')

    if ncc.size == 0:
        raise RuntimeError(
            "NCC surface is empty — search region is too small for the template. "
            "Reduce box_ref size or increase search_radius.")

    search_std = search_img.std()
    if search_std > 1e-10:
        ncc = ncc / search_std

    # --- Sub-pixel peak location ---
    peak_row, peak_col = _ncc_subpixel_peak(ncc)
    peak_val = float(ncc[int(round(peak_row)), int(round(peak_col))])

    # --- Convert NCC peak → original target pixel coordinates ---
    # Position (r, c) in the valid NCC output means the template top-left corner
    # is at (r, c) in the (rescaled) search_img.
    # Template centre in rescaled search coords:
    cx_rescaled = peak_col + tmpl_w / 2.0
    cy_rescaled = peak_row + tmpl_h / 2.0

    # Map from rescaled search coords back to original target pixels:
    # rescaled pixel j → original target pixel sx0 + j / zoom_factor
    tgt_cx = sx0 + cx_rescaled / zoom_factor
    tgt_cy = sy0 + cy_rescaled / zoom_factor
    tgt_cx = float(np.clip(tgt_cx, 0, t_nx - 1))
    tgt_cy = float(np.clip(tgt_cy, 0, t_ny - 1))

    shift_dx    = tgt_cx - pred_tx
    shift_dy    = tgt_cy - pred_ty
    shift_arcsec = np.sqrt((shift_dx * tgt_scale) ** 2 +
                           (shift_dy * tgt_scale) ** 2)

    print(f"[rb_align] arc_ncc: NCC peak = {peak_val:.4f} | "
          f"WCS pred = ({pred_tx:.1f}, {pred_ty:.1f}) px | "
          f"NCC found = ({tgt_cx:.2f}, {tgt_cy:.2f}) px | "
          f"shift = ({shift_dx:+.2f}, {shift_dy:+.2f}) px "
          f"= {shift_arcsec:.2f}\"")

    aligner.pairs.clear()
    aligner.pairs.append((ref_cx, ref_cy, ra, dec, tgt_cx, tgt_cy))


# ---------------------------------------------------------------------------
# Auto cascade
# ---------------------------------------------------------------------------

def _auto(aligner, stretch='zscale', box=0.1):
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
                  box: float = 0.1,
                  save_catalog: Optional[str] = None,
                  catalog: Optional[str] = None,
                  # cross_corr / dao / knots options
                  fwhm: float = 0.0,
                  threshold_sigma: float = 3.0,
                  max_sources: int = 50,
                  bg_method: str = 'sigma_clip',
                  # knots options
                  max_knots: int = 10,
                  # gaia options
                  radius_deg: float = 0.5,
                  max_stars: int = 50,
                  # batch options
                  headless: bool = False,
                  # arc_ncc options
                  box_ref=None,
                  search_radius: float = 5.0,
                  use_LoG: bool = False,
                  LoG_sigma: float = 0.0):
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
             threshold_sigma=threshold_sigma, max_sources=max_sources,
             bg_method=bg_method)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'knots':
        _knots(aligner, box=box, threshold_sigma=threshold_sigma,
               max_knots=max_knots, bg_method=bg_method)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'cross_corr':
        _cross_corr(aligner, box=box, fwhm=fwhm,
                    threshold_sigma=threshold_sigma, max_sources=max_sources,
                    bg_method=bg_method)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'arc_ncc':
        if box_ref is None:
            # No box given → open interactive draw-and-match window
            _arc_ncc_interactive(aligner, search_radius=search_radius,
                                 use_LoG=use_LoG, LoG_sigma=LoG_sigma,
                                 bg_method=bg_method, stretch=stretch)
        else:
            _arc_ncc(aligner, box_ref=box_ref, search_radius=search_radius,
                     use_LoG=use_LoG, LoG_sigma=LoG_sigma, bg_method=bg_method)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    elif strategy == 'batch':
        if catalog is None:
            raise ValueError("strategy='batch' requires catalog='path/to/catalog.fits'.")
        _batch(aligner, catalog=catalog, box=box,
               stretch=stretch, save_catalog=save_catalog, headless=headless)

    elif strategy == 'auto':
        _auto(aligner, stretch=stretch, box=box)
        if save_catalog is not None:
            _save_pairs_catalog(aligner.pairs, save_catalog)

    else:
        raise ValueError(
            f"Unknown strategy '{strategy}'. Choose from: "
            "interactive, gaia, dao, knots, cross_corr, arc_ncc, batch, auto."
        )

    print(f"[rb_align] find_sources('{strategy}'): "
          f"{len(aligner.pairs)} pair(s) stored.")
