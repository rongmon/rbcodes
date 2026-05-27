"""
Profile-weighted (optimal) spectral extraction for IFU cubes.

Implements the static spatial-profile approach from extract_weighted_spectrum
(kcwitools / Horne 1986 in spirit), adapted for arbitrary boolean aperture masks
and without any kcwitools / linetools dependency.

Two profile methods
-------------------
'data'
    Spatial weights = collapsed whitelight image within the aperture.
    Fast, no fitting, works even for irregular/extended sources.

'gaussian'
    Spatial weights = 2-D Gaussian fitted to the whitelight image.
    Smoother profile; more robust when individual spaxels are noisy.
    Requires astropy.modeling.
"""
import numpy as np
import warnings


def extract_optimal_weighted(flux, var, mask, method='data'):
    """
    Profile-weighted optimal extraction over a spatial boolean mask.

    Parameters
    ----------
    flux   : ndarray, shape (n_wave, ny, nx)
    var    : ndarray, same shape, or None
        If None the returned *sig* is None.
    mask   : ndarray bool, shape (ny, nx)
        Aperture mask — any shape (circular, rectangular, freehand).
    method : {'data', 'gaussian'}
        How the spatial weight profile is built (see module docstring).

    Returns
    -------
    fl  : ndarray (n_wave,)
        Extracted, renormalized flux spectrum.
    sig : ndarray (n_wave,) or None
        Propagated 1-sigma error spectrum (None when var is None).

    Notes
    -----
    The flux is renormalized after weighted extraction so that its total
    integral matches the simple boxcar sum, preserving physical units.
    Error is scaled by the same factor.
    """
    if not mask.any():
        raise ValueError("Aperture mask is empty — no spaxels selected.")

    # ------------------------------------------------------------------ #
    # Build whitelight image within mask
    # ------------------------------------------------------------------ #
    wl = np.nansum(flux, axis=0).copy()   # (ny, nx)
    wl[~mask] = 0.0
    wl[~np.isfinite(wl)] = 0.0
    wl[wl < 0] = 0.0                      # negative flux → no weight

    # ------------------------------------------------------------------ #
    # Build spatial weight map
    # ------------------------------------------------------------------ #
    if method == 'gaussian':
        weights = _gaussian_profile(wl, mask)
    else:
        weights = wl.copy()

    # Fall back to uniform if the profile is flat/zero (e.g. all-negative cube)
    w_sum = weights[mask].sum()
    if w_sum <= 0:
        weights[mask] = 1.0
        w_sum = float(mask.sum())

    # Normalize within aperture
    w_norm = np.zeros_like(weights)
    w_norm[mask] = weights[mask] / w_sum   # sums to 1 over aperture

    # ------------------------------------------------------------------ #
    # Vectorized extraction
    # ------------------------------------------------------------------ #
    w_flat = w_norm[mask]                  # (n_pix,)
    f_flat = flux[:, mask]                 # (n_wave, n_pix)

    fl = np.nansum(f_flat * w_flat, axis=1)   # (n_wave,)

    if var is not None:
        v = var.copy()
        v[~np.isfinite(v)] = 0.0
        v[v < 0] = 0.0
        v_flat = v[:, mask]                # (n_wave, n_pix)
        sig = np.sqrt(np.nansum(v_flat * w_flat ** 2, axis=1))
    else:
        sig = None

    # ------------------------------------------------------------------ #
    # Renormalize to boxcar-sum flux scale
    # ------------------------------------------------------------------ #
    fl_box = np.nansum(f_flat, axis=1)
    total_box = np.nansum(fl_box)
    total_opt = np.nansum(fl)
    if total_opt != 0 and np.isfinite(total_opt):
        norm = total_box / total_opt
    else:
        norm = 1.0

    fl  = fl  * norm
    if sig is not None:
        sig = sig * norm

    return fl, sig


# ---------------------------------------------------------------------------
# Gaussian profile helper
# ---------------------------------------------------------------------------

def _gaussian_profile(wl, mask):
    """
    Fit a 2-D Gaussian to *wl* (whitelight image, zeros outside mask)
    and return the evaluated model — zeroed outside *mask*.

    Falls back to the raw data if the fit fails.
    """
    try:
        from astropy.modeling import models, fitting
    except ImportError:
        warnings.warn("astropy.modeling not available; falling back to Data profile.")
        return wl.copy()

    ydim, xdim = wl.shape
    yi, xi = np.indices(wl.shape)

    # Initial guesses from the data
    amp = float(wl.max()) if wl.max() > 0 else 1.0
    # centroid from intensity-weighted mean within mask
    total = wl[mask].sum()
    if total > 0:
        cx = float((xi[mask] * wl[mask]).sum() / total)
        cy = float((yi[mask] * wl[mask]).sum() / total)
    else:
        cy, cx = ydim / 2.0, xdim / 2.0

    stdev = max(float(mask.sum()) ** 0.5 / 2.0, 1.0)

    g_init = models.Gaussian2D(
        amplitude=amp,
        x_mean=cx, y_mean=cy,
        x_stddev=stdev, y_stddev=stdev,
        theta=0,
    )

    fitter = fitting.LevMarLSQFitter()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            g_fit = fitter(g_init, xi, yi, wl)
        profile = g_fit(xi, yi)
    except Exception:
        warnings.warn("Gaussian profile fit failed; falling back to Data profile.")
        return wl.copy()

    profile[profile < 0] = 0.0
    profile[~mask] = 0.0
    return profile
