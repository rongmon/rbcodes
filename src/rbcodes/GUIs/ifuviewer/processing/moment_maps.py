"""
Moment map computation for IFU cubes.

Phase 10.
"""
import warnings
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


def compute_snr_map(m0, flux, wave, wmin, wmax,
                    var=None, sky_mask=None, cont1=None, cont2=None):
    """
    Compute a per-spaxel SNR map for use as a moment-map mask.

    Parameters
    ----------
    m0       : ndarray (ny, nx)  — moment-0 map (integrated flux)
    flux     : ndarray (n_wave, ny, nx)
    wave     : ndarray (n_wave,)
    wmin, wmax : float  — line window boundaries (Å)
    var      : ndarray (n_wave, ny, nx) or None  — variance cube
    sky_mask : ndarray (ny, nx) bool or None     — sky region mask
    cont1, cont2 : (float, float) or None        — continuum windows (Å)

    Returns
    -------
    snr : ndarray (ny, nx) or None if no noise estimate is possible
    """
    line_mask = (wave >= wmin) & (wave <= wmax)
    if not line_mask.any():
        return None
    dw = np.gradient(wave[line_mask])
    N  = int(line_mask.sum())

    if var is not None:
        # Per-spaxel: σ_M0 = sqrt(∫ var · dw²)
        with np.errstate(all='ignore'):
            sigma = np.sqrt(np.nansum(
                var[line_mask] * dw[:, None, None] ** 2, axis=0))

    elif sky_mask is not None and sky_mask.any():
        # Sky region: RMS per channel across sky spaxels, then propagate to M0
        sky_flux = flux[line_mask][:, sky_mask]   # (N_line, N_sky)
        with np.errstate(all='ignore'):
            rms_per_ch = np.nanstd(sky_flux, axis=1)          # (N_line,)
            sigma = np.sqrt(np.nansum(
                (rms_per_ch[:, None, None] * dw[:, None, None]) ** 2,
                axis=0))                                        # (ny, nx) — same value everywhere

    elif cont1 is not None:
        # Continuum window(s): RMS per spaxel → scale to line window
        cmask = (wave >= cont1[0]) & (wave <= cont1[1])
        if cont2 is not None:
            cmask |= (wave >= cont2[0]) & (wave <= cont2[1])
        if not cmask.any():
            return None
        with np.errstate(all='ignore'):
            rms = np.nanstd(flux[cmask], axis=0)               # (ny, nx)
            sigma = rms * np.sqrt(N) * float(np.mean(np.abs(dw)))

    else:
        return None

    with np.errstate(all='ignore'):
        snr = np.where(sigma > 0, m0 / sigma, np.nan)
    return snr


def subtract_linear_continuum(flux, wave, bcont_min, bcont_max, rcont_min, rcont_max):
    """
    Subtract a per-spaxel linear continuum from a flux cube.

    A linear baseline is anchored at the mean flux in the blue and red
    continuum windows, evaluated at the midpoint wavelength of each window.
    The fit is extrapolated / interpolated across the full wavelength range.

    Parameters
    ----------
    flux       : ndarray (n_wave, ny, nx)
    wave       : ndarray (n_wave,)
    bcont_min, bcont_max : float  blue continuum window (Å)
    rcont_min, rcont_max : float  red  continuum window (Å)

    Returns
    -------
    flux_sub : ndarray (n_wave, ny, nx)  continuum-subtracted cube
    """
    bmask = (wave >= bcont_min) & (wave <= bcont_max)
    rmask = (wave >= rcont_min) & (wave <= rcont_max)
    if not bmask.any():
        raise ValueError(f"Blue continuum window [{bcont_min:.1f}, {bcont_max:.1f}] Å has no channels.")
    if not rmask.any():
        raise ValueError(f"Red continuum window [{rcont_min:.1f}, {rcont_max:.1f}] Å has no channels.")

    lambda_b = float(wave[bmask].mean())
    lambda_r = float(wave[rmask].mean())
    if lambda_b >= lambda_r:
        raise ValueError("Blue continuum window must be at shorter wavelength than red window.")

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        flux_b = np.nanmean(flux[bmask], axis=0)   # (ny, nx)
        flux_r = np.nanmean(flux[rmask], axis=0)   # (ny, nx)

    slope     = (flux_r - flux_b) / (lambda_r - lambda_b)            # (ny, nx)
    intercept = flux_b - slope * lambda_b                             # (ny, nx)
    continuum = (slope[np.newaxis] * wave[:, np.newaxis, np.newaxis]
                 + intercept[np.newaxis])                             # (n_wave, ny, nx)
    return flux - continuum


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
