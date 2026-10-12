"""
Noise cubes for IFU data.

    from rbcodes.ifu.noise import noise_cube
    sigma, how = noise_cube(cube.flux, var=cube.var)                     # sqrt(variance)
    sky = sky_mask_from_region('sky.reg', cube.wcs, cube.flux.shape[1:])  # or any bool (ny, nx)
    sigma, how = noise_cube(cube.flux, line_free=lf, sky_mask=sky)       # no variance: estimate

Without a variance cube the noise is modelled as sigma(lambda, y, x) = sigma_sky(lambda) * s(y, x):
sigma_sky is the robust spread across the sky spaxels in each channel (follows sky-line
residuals), s is the robust spread of flux / sigma_sky over line-free channels in each spaxel
(follows the exposure map). The sky region must be chosen by the user: blank sky, no sources,
away from the cube edges.
"""
import numpy as np

__all__ = ['robust_std', 'noise_cube', 'sky_mask_from_region']


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
    sky_mask : bool ndarray (ny, nx)
        Sky spaxels; required when ``var`` is None (see ``sky_mask_from_region``).

    Returns
    -------
    sigma : ndarray like flux
    description : str
    """
    if var is not None:
        v = np.where(np.isfinite(var) & (var > 0), var, np.nan)
        return np.sqrt(v), 'variance cube'
    if line_free is None or sky_mask is None:
        raise ValueError('Without a variance cube, give line_free channels and a sky_mask '
                         '(blank-sky spaxels) to estimate the noise.')
    line_free = np.asarray(line_free, dtype=bool)
    if line_free.sum() < 5:
        raise ValueError(f'Need at least 5 line-free channels, got {line_free.sum()}.')
    sky_mask = np.asarray(sky_mask, dtype=bool)
    if sky_mask.shape != flux.shape[1:]:
        raise ValueError(f'sky_mask shape {sky_mask.shape} does not match the cube {flux.shape[1:]}.')
    if sky_mask.sum() < 5:
        raise ValueError(f'Need at least 5 sky spaxels, got {sky_mask.sum()}.')

    sig_sky = robust_std(flux[:, sky_mask], axis=1)
    sig_sky = np.where(np.isfinite(sig_sky) & (sig_sky > 0), sig_sky, np.nanmedian(sig_sky))
    scale = robust_std(flux[line_free] / sig_sky[line_free, None, None], axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, np.nanmedian(scale))
    return sig_sky[:, None, None] * scale[None], 'sky spectrum x per-spaxel scale'


def sky_mask_from_region(region, wcs, shape):
    """
    Boolean (ny, nx) sky mask from ds9 regions (the union of all extraction shapes).

    Parameters
    ----------
    region : str or path
        A ds9 region file, or region text such as ``'image\nbox(5,3,8,5,0)'``.
    wcs : astropy.wcs.WCS or None
        2D WCS of the cube (needed for sky-coordinate regions).
    shape : (ny, nx)
    """
    import os
    from .spatial_mask import parse_ds9_regions, region_to_mask, EXTRACTION_SHAPES
    text = region = os.fspath(region)
    if '\n' not in region and os.path.isfile(region):
        with open(region) as fh:
            text = fh.read()
    regs = [r for r in parse_ds9_regions(text) if r.get('shape') in EXTRACTION_SHAPES]
    if not regs:
        raise ValueError('No region shapes found for the sky mask.')
    mask = np.zeros(shape, dtype=bool)
    for r in regs:
        mask |= region_to_mask(r, wcs, *shape)
    return mask
