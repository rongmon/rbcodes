"""
Noise cubes for IFU data.

    from rbcodes.ifu.noise import noise_cube
    sigma, how = noise_cube(cube.flux, var=cube.var)                     # sqrt(variance)
    sigma, how = noise_cube(cube.flux, line_free=lf, sky_mask=sky)       # no variance: estimate

Without a variance cube the noise is modelled as sigma(lambda, y, x) = sigma_sky(lambda) * s(y, x):
sigma_sky is the robust spread across sky spaxels in each channel (follows sky-line residuals),
s is the robust spread of flux / sigma_sky over line-free channels in each spaxel (follows the
exposure map).
"""
import numpy as np

__all__ = ['robust_std', 'noise_cube']


def robust_std(x, axis=None):
    """Standard deviation from the median absolute deviation (1.4826 MAD), ignoring NaN."""
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x, axis=axis, keepdims=True)
    return 1.4826 * np.nanmedian(np.abs(x - med), axis=axis)


def noise_cube(flux, var=None, line_free=None, sky_mask=None):
    """
    Per-channel, per-spaxel 1-sigma noise cube.

    Parameters
    ----------
    flux : ndarray (n_wave, ny, nx)
    var : ndarray like flux, optional
        Variance cube. If given, the noise is sqrt(var) (non-positive or non-finite -> NaN).
    line_free : bool ndarray (n_wave,), optional
        Channels without line emission; required when ``var`` is None.
    sky_mask : bool ndarray (ny, nx), optional
        Sky spaxels. Default: the 10% of spaxels with the lowest |summed flux|.

    Returns
    -------
    sigma : ndarray like flux
    description : str
    """
    if var is not None:
        v = np.where(np.isfinite(var) & (var > 0), var, np.nan)
        return np.sqrt(v), 'variance cube'
    if line_free is None:
        raise ValueError('Without a variance cube, give line_free channels to estimate the noise.')
    line_free = np.asarray(line_free, dtype=bool)
    if line_free.sum() < 5:
        raise ValueError(f'Need at least 5 line-free channels, got {line_free.sum()}.')
    if sky_mask is None:
        tot = np.abs(np.nansum(flux, axis=0))
        sky_mask = tot <= np.nanpercentile(tot, 10)
    sky_mask = np.asarray(sky_mask, dtype=bool)
    if sky_mask.sum() < 5:
        raise ValueError(f'Need at least 5 sky spaxels, got {sky_mask.sum()}.')

    sig_sky = robust_std(flux[:, sky_mask], axis=1)
    sig_sky = np.where(np.isfinite(sig_sky) & (sig_sky > 0), sig_sky, np.nanmedian(sig_sky))
    scale = robust_std(flux[line_free] / sig_sky[line_free, None, None], axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, np.nanmedian(scale))
    return sig_sky[:, None, None] * scale[None], 'sky spectrum x per-spaxel scale'
