"""
Spatial binning of IFU cubes.

    from rbcodes.ifu.binning import bin_cube
    cube3 = bin_cube(cube, 3)          # 3x3 blocks: summed flux and variance, rescaled WCS

Non-overlapping blocks keep the errors of different bins independent (what the debiased
structure function assumes); the result is an ordinary cube, so narrow-band images, line fits,
VelocityMap, the scattering transform and ray tracing work on it unchanged.
"""
import numpy as np

__all__ = ['bin_cube']


def _block_sum(a, f):
    nw, ny, nx = a.shape
    a = a[:, :ny // f * f, :nx // f * f]
    return a.reshape(nw, ny // f, f, nx // f, f).sum(axis=(2, 4))


def bin_cube(cube, factor, require_all=True):
    """
    Sum a cube over non-overlapping factor x factor blocks of spaxels.

    Parameters
    ----------
    cube : IFUCube (``rbcodes.ifu.load_fits``)
    factor : int
        Block size in spaxels. Edge rows/columns that do not fill a whole block are dropped.
    require_all : bool
        A block with any non-finite spaxel (in a channel) is NaN in that channel (default), so
        partial sums are never mixed with full ones; False sums the finite spaxels.

    Returns
    -------
    IFUCube with flux and var (sums), the same wavelengths, a WCS whose pixels are the blocks
    (centre of new pixel = centre of its block), header NAXIS/CRPIX/CDELT (or CD) updated, and
    ``binning`` = factor.
    """
    from .io.base import IFUCube
    f = int(factor)
    if f < 1:
        raise ValueError('factor must be a positive integer.')
    flux = np.asarray(cube.flux, dtype=float)
    if flux.ndim != 3:
        raise ValueError('bin_cube needs a 3-D cube.')

    def bsum(a):
        s = _block_sum(np.nan_to_num(a, nan=0.), f)
        if require_all:
            bad = _block_sum((~np.isfinite(a)).astype(float), f) > 0
            s[bad] = np.nan
        return s

    out = IFUCube()
    out.flux = bsum(flux)
    out.var = None if cube.var is None else bsum(np.asarray(cube.var, dtype=float))
    out.wave = None if cube.wave is None else np.array(cube.wave, copy=True)
    out.airvac = getattr(cube, 'airvac', 'unknown')
    out.binning = f * getattr(cube, 'binning', 1)
    out.name = f'{cube.name}_bin{f}' if getattr(cube, 'name', '') else f'bin{f}'
    out.path = ''

    out.wcs = None
    if cube.wcs is not None:
        w = cube.wcs.deepcopy()
        w.wcs.crpix = (np.asarray(w.wcs.crpix) - 0.5) / f + 0.5
        if w.wcs.has_cd():
            w.wcs.cd = np.asarray(w.wcs.cd) * f
        else:
            w.wcs.cdelt = np.asarray(w.wcs.cdelt) * f
        w.wcs.set()
        out.wcs = w

    out.header = None
    if getattr(cube, 'header', None) is not None:
        h = cube.header.copy()
        ny, nx = out.flux.shape[1:]
        h['NAXIS1'], h['NAXIS2'] = nx, ny
        for ax in (1, 2):
            if f'CRPIX{ax}' in h:
                h[f'CRPIX{ax}'] = (h[f'CRPIX{ax}'] - 0.5) / f + 0.5
            if f'CDELT{ax}' in h and not any(f'CD{ax}_{j}' in h for j in (1, 2)):
                h[f'CDELT{ax}'] = h[f'CDELT{ax}'] * f          # (a PC matrix stays as it is)
            for j in (1, 2):
                if f'CD{ax}_{j}' in h:
                    h[f'CD{ax}_{j}'] = h[f'CD{ax}_{j}'] * f
        h.add_history(f'rbcodes: summed over {f}x{f} spaxel blocks')
        out.header = h
    return out
