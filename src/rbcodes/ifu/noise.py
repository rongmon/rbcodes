"""
Noise model for IFU cubes.

    from rbcodes.ifu.noise import NoiseModel
    nm = NoiseModel.from_cube(cube)                                   # variance extension present
    nm = NoiseModel.from_cube(cube, linefree=lm, sky_box=(0, 5, 0, 8)) # no variance: hybrid estimate
    nm.sigma            # 1-sigma noise cube (n_wave, ny, nx)
    nm.correlation(cube.flux, linefree=lm)   # adjacent-spaxel / adjacent-channel correlation

Without a variance cube the noise is sigma(lambda, y, x) = sigma_sky(lambda) * s(y, x):
sigma_sky(lambda) is the robust spread across the sky spaxels in each channel (captures sky-line
residuals); s(y, x) is the robust spread over line-free channels of f(lambda)/sigma_sky(lambda)
(captures exposure-map variations). ``spaxel_scale=False`` drops s(y, x). The sky region is
chosen by the user (``sky_box``, ``sky_mask`` or a ds9 ``sky_region``).
"""
import os

import numpy as np

__all__ = ['robust_std', 'NoiseModel', 'sky_mask_from_region', 'narrowband', 'NarrowBand']


def robust_std(x, axis=None):
    """1.4826 x median absolute deviation, ignoring NaN."""
    med = np.nanmedian(x, axis=axis, keepdims=True)
    return 1.4826 * np.nanmedian(np.abs(x - med), axis=axis)


def sky_mask_from_region(region, wcs, shape):
    """Boolean (ny, nx) mask from ds9 regions (file or text; union of all extraction shapes)."""
    from .spatial_mask import parse_ds9_regions, region_to_mask, EXTRACTION_SHAPES
    text = region = os.fspath(region)
    if '\n' not in region and os.path.isfile(region):
        with open(region) as fh:
            text = fh.read()
    regs = [r for r in parse_ds9_regions(text) if r.get('shape') in EXTRACTION_SHAPES]
    if not regs:
        raise ValueError('No region shapes found.')
    mask = np.zeros(shape, dtype=bool)
    for r in regs:
        mask |= region_to_mask(r, wcs, *shape)
    return mask


def _sky(shape, sky_box=None, sky_mask=None, sky_region=None, wcs=None):
    given = [x is not None for x in (sky_box, sky_mask, sky_region)]
    if sum(given) != 1:
        raise ValueError('Give exactly one of sky_box=(y0, y1, x0, x1), sky_mask (bool ny x nx) '
                         'or sky_region (ds9 file or text).')
    if sky_box is not None:
        y0, y1, x0, x1 = sky_box
        m = np.zeros(shape, bool)
        m[y0:y1, x0:x1] = True
    elif sky_mask is not None:
        m = np.asarray(sky_mask, dtype=bool)
        if m.shape != tuple(shape):
            raise ValueError(f'sky_mask shape {m.shape} does not match the cube {tuple(shape)}.')
    else:
        m = sky_mask_from_region(sky_region, wcs, shape)
    if m.sum() < 5:
        raise ValueError(f'The sky region has {m.sum()} spaxels; need at least 5.')
    return m


def _linefree(linefree, wave, n):
    """bool (n_wave,) from a bool array, a list of (lo, hi) Angstrom ranges, or a LineModel."""
    if linefree is None:
        return None
    if hasattr(linefree, 'masks'):
        return linefree.masks(wave)[2]
    arr = np.asarray(linefree)
    if arr.dtype == bool:
        if arr.shape != (n,):
            raise ValueError(f'linefree mask has shape {arr.shape}, expected ({n},).')
        return arr
    if wave is None:
        raise ValueError('Wavelength ranges for linefree need the wavelength array.')
    m = np.zeros(n, bool)
    for lo, hi in linefree:
        m |= (wave > lo) & (wave < hi)
    return m


class NoiseModel:
    """
    Per-channel, per-spaxel 1-sigma noise.

    Attributes
    ----------
    sigma : ndarray (n_wave, ny, nx)
    description : str
    sky_mask : bool (ny, nx) or None
    wave : ndarray (n_wave,) or None
    sig_sky : ndarray (n_wave,) or None      (hybrid model)
    scale : ndarray (ny, nx) or None         (hybrid model, s(y, x))
    """

    def __init__(self, sigma, description='', sky_mask=None, sig_sky=None, scale=None, wave=None):
        self.sigma, self.description = sigma, description
        self.sky_mask, self.sig_sky, self.scale = sky_mask, sig_sky, scale
        self.wave = None if wave is None else np.asarray(wave, dtype=float)

    def __repr__(self):
        return f'NoiseModel({self.description}, shape={self.sigma.shape})'

    @classmethod
    def from_cube(cls, cube, linefree=None, sky_box=None, sky_mask=None, sky_region=None,
                  spaxel_scale=True, use_variance=True, wave=None, var=None, wcs=None):
        """
        Parameters
        ----------
        cube : IFUCube (``rbcodes.ifu.load_fits``) or ndarray (n_wave, ny, nx)
        linefree : bool (n_wave,), list of (lo, hi) Angstrom, or a LineModel
            Line-free channels (needed for the hybrid model).
        sky_box : (y0, y1, x0, x1) | sky_mask : bool (ny, nx) | sky_region : ds9 file/text
            Blank-sky spaxels (exactly one; needed for the hybrid model).
        spaxel_scale : bool
            Include the per-spaxel scale s(y, x) in the hybrid model.
        use_variance : bool
            Use the variance cube when present (False forces the hybrid model).
        wave, var, wcs : for an array cube.
        """
        if hasattr(cube, 'flux'):
            flux, wave, wcs = cube.flux, cube.wave, cube.wcs if wcs is None else wcs
            var = cube.var if var is None else var
        else:
            flux = np.asarray(cube, dtype=float)
        if var is not None and use_variance:
            v = np.where(np.isfinite(var) & (var > 0), var, np.nan)
            sky = (_sky(flux.shape[1:], sky_box, sky_mask, sky_region, wcs)
                   if any(x is not None for x in (sky_box, sky_mask, sky_region)) else None)
            return cls(np.sqrt(v), 'variance extension', sky_mask=sky, wave=wave)

        lf = _linefree(linefree, wave, flux.shape[0])
        if lf is None or lf.sum() < 5:
            raise ValueError('No usable variance: give linefree channels (>= 5) and a sky region '
                             'for the hybrid noise model.')
        sky = _sky(flux.shape[1:], sky_box, sky_mask, sky_region, wcs)
        sig_sky = robust_std(flux[:, sky], axis=1)
        sig_sky = np.where(np.isfinite(sig_sky) & (sig_sky > 0), sig_sky, np.nanmedian(sig_sky))
        if spaxel_scale:
            norm = flux[lf] / sig_sky[lf, None, None]
            s = robust_std(norm, axis=0)
            s = np.where(np.isfinite(s) & (s > 0), s, np.nanmedian(s))
            desc = 'hybrid sky-spectrum x per-spaxel scale'
        else:
            s = np.ones(flux.shape[1:])
            desc = 'sky spectrum only (no per-spaxel scale)'
        return cls(sig_sky[:, None, None] * s[None], desc, sky_mask=sky, sig_sky=sig_sky, scale=s, wave=wave)

    def correlation(self, flux, linefree, sky_mask=None, wave=None):
        """
        Correlation of the normalised noise between neighbouring spaxels and channels.

        Uses r = (f - median over line-free channels) / sigma in sky spaxels and line-free channels.

        Returns
        -------
        dict(rho_x, rho_y, rho_lambda, n_x, n_y, n_lambda): Pearson coefficients for pairs
        of sky spaxels adjacent along x, along y, and for adjacent line-free channels.
        """
        sky = self.sky_mask if sky_mask is None else np.asarray(sky_mask, bool)
        if sky is None:
            raise ValueError('Give sky_mask (the model has none).')
        lf = _linefree(linefree, self.wave if wave is None else wave, flux.shape[0])
        if lf is None:
            raise ValueError('Give linefree channels.')
        r = (flux - np.nanmedian(flux[lf], axis=0)) / self.sigma

        def pearson(a, b):
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < 10:
                return np.nan, int(ok.sum())
            return float(np.corrcoef(a[ok], b[ok])[0, 1]), int(ok.sum())

        rl = r[lf]
        out = {}
        pairs = dict(x=(rl[:, :, :-1], rl[:, :, 1:], sky[:, :-1] & sky[:, 1:]),
                     y=(rl[:, :-1, :], rl[:, 1:, :], sky[:-1, :] & sky[1:, :]))
        for key, (a, b, m) in pairs.items():
            out[f'rho_{key}'], out[f'n_{key}'] = pearson(a[:, m], b[:, m])
        both = lf[:-1] & lf[1:]
        out['rho_lambda'], out['n_lambda'] = pearson(r[:-1][both][:, sky], r[1:][both][:, sky])
        return out


class NarrowBand:
    """
    Continuum-subtracted narrow-band image with its noise (from ``narrowband``).

    Attributes: ``image`` (flux x Angstrom), ``sigma_stat`` (propagated from the variance cube),
    ``k`` (sky calibration), ``sigma`` = k sigma_stat, ``snr`` = image / sigma, ``continuum``
    (mean continuum level per spaxel), ``sky_mask``, ``band``, ``windows``.
    """

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def __repr__(self):
        return (f'NarrowBand(band={self.band}, k={self.k:.3f}, sky spaxels={int(self.sky_mask.sum())}, '
                f'S/N>3: {int(np.nansum(self.snr > 3))})')


def _auto_sky(cont, finite, nb_over_sig, clip, edge=3, grow=2):
    """Blank sky: interior, off continuum sources, then without strong line emitters."""
    from scipy.ndimage import binary_dilation, binary_erosion
    interior = binary_erosion(finite, iterations=edge)
    c = cont[interior]
    med, rs = np.nanmedian(c), robust_std(c)
    src = np.zeros_like(finite)
    src[interior] = cont[interior] > med + 3 * rs
    src = binary_dilation(src, iterations=grow)
    sky = interior & ~src & np.isfinite(nb_over_sig)
    k = robust_std(nb_over_sig[sky])
    sky &= np.abs(nb_over_sig - np.nanmedian(nb_over_sig[sky])) < clip * k
    return sky


def narrowband(cube, band, windows, var=None, wave=None, wcs=None, continuum='linear', sky_box=None,
               sky_mask=None, sky_region=None, clip=4.):
    """
    Continuum-subtracted narrow-band image and its calibrated noise.

    NB = sum over the band of (f - continuum) dlambda, with the continuum a constant or a straight
    line fitted per spaxel to the channels in ``windows``. sigma_stat is propagated exactly from the
    variance (band channels and continuum channels with their weights). k = robust scatter of
    NB / sigma_stat over blank sky, and sigma = k sigma_stat: the variance cube sets the spatial
    pattern, the sky sets the level (MUSE variance cubes miss the correlated noise, so k > 1).

    Parameters
    ----------
    cube : IFUCube or ndarray (n_wave, ny, nx)
    band : (lo, hi) Angstrom, or a list of (lo, hi) whose channels are summed together (e.g. one
        window per component of a doublet, skipping an absorption trough between them)
    windows : list of (lo, hi) Angstrom, continuum (line-free) windows
    var : variance cube (default ``cube.var``); without one sigma_stat = 1 and k is the sky noise
    wave, wcs : for an array cube
    continuum : 'linear' or 'constant'
    sky_box / sky_mask / sky_region : blank sky; default automatic (see ``NarrowBand``)
    clip : line emitters with |NB / sigma_stat| > clip k are removed from the automatic sky

    Returns
    -------
    NarrowBand
    """
    if hasattr(cube, 'flux'):
        flux, wave = cube.flux, cube.wave
        var = cube.var if var is None else var
        wcs = cube.wcs if wcs is None else wcs
    else:
        flux = np.asarray(cube, dtype=float)
        if wave is None:
            raise ValueError('With an array cube, give wave=.')
    wave = np.asarray(wave, dtype=float)
    if continuum not in ('linear', 'constant'):
        raise ValueError("continuum must be 'linear' or 'constant'.")
    bands = [tuple(band)] if np.ndim(band) == 1 else [tuple(x) for x in band]
    b = np.zeros(len(wave), bool)
    for lo, hi in bands:
        b |= (wave >= lo) & (wave <= hi)
    c = np.zeros(len(wave), bool)
    for lo, hi in windows:
        c |= (wave >= lo) & (wave <= hi)
    c &= ~b
    if b.sum() == 0 or c.sum() < 2:
        raise ValueError(f'The band has {b.sum()} channels and the continuum windows {c.sum()}; need >= 1 and >= 2.')
    dl = np.gradient(wave)
    lam0 = np.mean(wave[b])
    if continuum == 'linear':
        X = np.c_[np.ones(c.sum()), wave[c] - lam0]
        P = np.linalg.pinv(X)                                     # (2, Nc)
        W = np.c_[np.ones(b.sum()), wave[b] - lam0] @ P           # (Nb, Nc): continuum at band channels
    else:
        W = np.full((b.sum(), c.sum()), 1. / c.sum())
    w = dl[b] @ W                                                 # weight of each continuum channel
    fb, fc = flux[b], flux[c]
    image = np.tensordot(dl[b], fb, 1) - np.tensordot(w, fc, 1)
    if var is not None:
        vb, vc = var[b], var[c]
        sigma_stat = np.sqrt(np.tensordot(dl[b] ** 2, vb, 1) + np.tensordot(w ** 2, vc, 1))
    else:
        sigma_stat = np.ones(image.shape)
    cont = np.nanmean(fc, axis=0)
    finite = np.isfinite(image) & np.isfinite(sigma_stat) & (sigma_stat > 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        r = np.where(finite, image / sigma_stat, np.nan)
    if any(x is not None for x in (sky_box, sky_mask, sky_region)):
        sky = _sky(image.shape, sky_box, sky_mask, sky_region, wcs) & finite
    else:
        sky = _auto_sky(cont, finite, r, clip)
    if sky.sum() < 20:
        raise ValueError(f'Only {sky.sum()} blank-sky spaxels; give sky_box / sky_mask / sky_region.')
    k = float(robust_std(r[sky]))
    sigma = k * sigma_stat
    with np.errstate(invalid='ignore', divide='ignore'):
        snr = np.where(finite, image / sigma, np.nan)
    return NarrowBand(image=image, sigma_stat=sigma_stat, k=k, sigma=sigma, snr=snr, continuum=cont,
                      sky_mask=sky, band=bands[0] if len(bands) == 1 else bands, windows=list(windows),
                      continuum_kind=continuum)
