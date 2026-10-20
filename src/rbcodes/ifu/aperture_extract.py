"""
Aperture spectral extraction helpers.

Phase 9.
"""
import numpy as np


# ---------------------------------------------------------------------------
# Mask builders
# ---------------------------------------------------------------------------

def make_circular_mask(ny, nx, cx, cy, radius):
    """Boolean mask: pixels within *radius* of (cx, cy)."""
    y, x = np.mgrid[0:ny, 0:nx]
    return (x - cx) ** 2 + (y - cy) ** 2 <= radius ** 2


def make_square_mask(ny, nx, cx, cy, half_width):
    """Boolean mask: pixels within a square of half-width centred on (cx, cy)."""
    y, x = np.mgrid[0:ny, 0:nx]
    return (np.abs(x - cx) <= half_width) & (np.abs(y - cy) <= half_width)


def make_annulus_mask(ny, nx, cx, cy, inner_r, outer_r):
    """Boolean mask: pixels in the annulus [inner_r, outer_r]."""
    y, x = np.mgrid[0:ny, 0:nx]
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    return (r2 >= inner_r ** 2) & (r2 <= outer_r ** 2)


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def extract_with_method(flux, mask, method='sum'):
    """
    Extract spectrum over *mask* using the given collapse method.

    Parameters
    ----------
    flux   : ndarray (n_wave, ny, nx)
    mask   : ndarray bool (ny, nx)
    method : 'sum' | 'mean' | 'median'

    Returns
    -------
    spec : ndarray (n_wave,)
    """
    f = flux[:, mask]   # (n_wave, n_pix)
    if method == 'sum':
        return np.nansum(f, axis=1)
    elif method == 'mean':
        return np.nanmean(f, axis=1)
    elif method == 'median':
        return np.nanmedian(f, axis=1)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'sum', 'mean', or 'median'.")


def extract_aperture(flux, var, mask):
    """
    Simple summed extraction over a spatial boolean mask.

    Parameters
    ----------
    flux : ndarray, shape (n_wave, ny, nx)
    var  : ndarray same shape, or None
    mask : ndarray bool, shape (ny, nx)

    Returns
    -------
    spec : ndarray (n_wave,)
    err  : ndarray (n_wave,) or None
    """
    spec = np.nansum(flux[:, mask], axis=1)
    if var is not None:
        with np.errstate(invalid='ignore'):
            err = np.sqrt(np.nansum(var[:, mask], axis=1))
    else:
        err = None
    return spec, err


def extract_variance_weighted(flux, var, mask):
    """
    Variance-weighted (quasi-optimal) extraction over a spatial mask.

    Each spaxel is weighted by 1/variance at every wavelength channel.
    This is the standard optimal extraction approach (Horne 1986 in spirit).

    Parameters
    ----------
    flux : ndarray, shape (n_wave, ny, nx)
    var  : ndarray same shape — required (raises ValueError if None)
    mask : ndarray bool, shape (ny, nx)

    Returns
    -------
    spec : ndarray (n_wave,)
    err  : ndarray (n_wave,)
    """
    if var is None:
        raise ValueError("Variance array required for variance-weighted extraction")

    f = flux[:, mask]   # (n_wave, n_spaxels)
    v = var[:, mask]

    w = np.where(v > 0, 1.0 / v, 0.0)
    w_sum = np.nansum(w, axis=1)   # (n_wave,)

    with np.errstate(invalid='ignore', divide='ignore'):
        spec = np.where(w_sum > 0,
                        np.nansum(f * w, axis=1) / w_sum,
                        np.nan)
        err  = np.where(w_sum > 0,
                        1.0 / np.sqrt(w_sum),
                        np.nan)

    return spec, err


# ---------------------------------------------------------------------------
# Background subtraction
# ---------------------------------------------------------------------------

def subtract_background(spec, flux, bg_mask, method='mean'):
    """
    Subtract a per-pixel background estimate from *spec*.

    The background level per wavelength is estimated from flux[:, bg_mask]
    using *method*, then subtracted as a constant per-pixel level.

    Parameters
    ----------
    spec    : ndarray (n_wave,)  — source spectrum already extracted
    flux    : ndarray (n_wave, ny, nx)
    bg_mask : ndarray bool (ny, nx)
    method  : 'mean' | 'median'

    Returns
    -------
    spec_sub : ndarray (n_wave,)
    """
    f_bg = flux[:, bg_mask]   # (n_wave, n_bg_pix)
    if method == 'median':
        bg_per_pix = np.nanmedian(f_bg, axis=1)
    else:
        bg_per_pix = np.nanmean(f_bg, axis=1)
    return spec - bg_per_pix
