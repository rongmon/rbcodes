"""
Moment map computation for IFU cubes.

Phase 10.
"""
import numpy as np

C_KMS = 2.998e5   # speed of light in km/s


def velocity_array(wave, lambda_rest):
    """
    Wavelength array → velocity in km/s relative to lambda_rest.

    lambda_rest is the line center in the *observed* frame (the wavelength
    that should map to v=0).  It must be in the same units as wave (Å).
    """
    return C_KMS * (wave - lambda_rest) / lambda_rest


def moment0(flux, wave, wmin, wmax):
    """
    Zeroth moment — integrated flux (flux · Å).

    Parameters
    ----------
    flux  : ndarray (n_wave, ny, nx)
    wave  : ndarray (n_wave,)
    wmin, wmax : float  wavelength window in Å

    Returns
    -------
    ndarray (ny, nx)
    """
    mask = (wave >= wmin) & (wave <= wmax)
    if not mask.any():
        raise ValueError(f"No channels in wavelength window [{wmin:.1f}, {wmax:.1f}] Å.")
    dw = np.gradient(wave[mask])
    with np.errstate(all='ignore'):
        return np.nansum(flux[mask] * dw[:, None, None], axis=0)


def moment1(flux, wave, wmin, wmax, lambda_rest):
    """
    First moment — flux-weighted velocity centroid (km/s).

    Spaxels where M0 ≤ 0 are masked to NaN.
    """
    mask = (wave >= wmin) & (wave <= wmax)
    if not mask.any():
        raise ValueError(f"No channels in wavelength window [{wmin:.1f}, {wmax:.1f}] Å.")
    vel = velocity_array(wave[mask], lambda_rest)
    dw  = np.gradient(wave[mask])
    with np.errstate(all='ignore'):
        m0 = np.nansum(flux[mask] * dw[:, None, None], axis=0)
        m1 = np.nansum(flux[mask] * vel[:, None, None] * dw[:, None, None], axis=0)
        return np.where(m0 > 0, m1 / m0, np.nan)


def moment2(flux, wave, wmin, wmax, lambda_rest):
    """
    Second moment — flux-weighted velocity dispersion (km/s).

    Spaxels where M0 ≤ 0 are masked to NaN.
    """
    mask = (wave >= wmin) & (wave <= wmax)
    if not mask.any():
        raise ValueError(f"No channels in wavelength window [{wmin:.1f}, {wmax:.1f}] Å.")
    vel = velocity_array(wave[mask], lambda_rest)
    dw  = np.gradient(wave[mask])
    with np.errstate(all='ignore'):
        m0  = np.nansum(flux[mask] * dw[:, None, None], axis=0)
        m1  = moment1(flux, wave, wmin, wmax, lambda_rest)
        dv2 = (vel[:, None, None] - m1[None]) ** 2
        m2  = np.nansum(flux[mask] * dv2 * dw[:, None, None], axis=0)
        return np.where(m0 > 0, np.sqrt(np.abs(m2 / m0)), np.nan)


def moment_map(flux, wave, wmin, wmax, order, lambda_rest=None):
    """
    Compute moment map of given order.

    Parameters
    ----------
    order       : int  0, 1, or 2
    lambda_rest : float or None  rest wavelength in Å (required for order 1 and 2)

    Returns
    -------
    ndarray (ny, nx)
    """
    if order == 0:
        return moment0(flux, wave, wmin, wmax)
    if order in (1, 2):
        if lambda_rest is None:
            raise ValueError("lambda_rest is required for moment order 1 and 2.")
        if order == 1:
            return moment1(flux, wave, wmin, wmax, lambda_rest)
        return moment2(flux, wave, wmin, wmax, lambda_rest)
    raise ValueError(f"Unsupported moment order {order}. Use 0, 1, or 2.")
