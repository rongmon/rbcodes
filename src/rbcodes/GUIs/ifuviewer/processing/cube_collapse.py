"""
Cube collapse functions — whitelight, narrowband, continuum-subtracted.

Phase 3 implements build_whitelight.
Phases 8+ add build_narrowband and build_continuum_sub.
"""
import warnings
import numpy as np


def build_whitelight(flux, wave, wmin=None, wmax=None, method='mean'):
    """
    Collapse a 3-D flux cube along the wavelength axis.

    Parameters
    ----------
    flux : np.ndarray, shape (n_wave, ny, nx)
    wave : np.ndarray, shape (n_wave,)
    wmin, wmax : float or None
        Wavelength limits in Angstroms.  None = no limit.
    method : {'mean', 'sum', 'median'}

    Returns
    -------
    np.ndarray, shape (ny, nx)
    """
    mask = np.ones(len(wave), dtype=bool)
    if wmin is not None:
        mask &= wave >= wmin
    if wmax is not None:
        mask &= wave <= wmax

    if not mask.any():
        raise ValueError(f"No wavelength channels in range [{wmin}, {wmax}].")

    data = flux[mask, :, :]
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('ignore', RuntimeWarning)
        if method == 'mean':
            return np.nanmean(data, axis=0)
        elif method == 'sum':
            return np.nansum(data, axis=0)
        elif method == 'median':
            return np.nanmedian(data, axis=0)
    raise ValueError(f"Unknown method '{method}'. Use 'mean', 'sum', or 'median'.")


def build_narrowband(flux, wave, wmin, wmax, method='mean'):
    """Alias for build_whitelight with explicit wavelength bounds (Phase 8)."""
    return build_whitelight(flux, wave, wmin=wmin, wmax=wmax, method=method)


def build_continuum_sub(flux, wave, wmin, wmax,
                        c1min, c1max, c2min=None, c2max=None, method='mean'):
    """
    Continuum-subtracted narrowband image (Phase 8).

    Subtracts the mean of one or two continuum windows from the on-band image.
    """
    nb   = build_narrowband(flux, wave, wmin, wmax, method)
    cont = build_narrowband(flux, wave, c1min, c1max, method)
    if c2min is not None and c2max is not None:
        cont2 = build_narrowband(flux, wave, c2min, c2max, method)
        cont  = (cont + cont2) / 2.0
    return nb - cont
