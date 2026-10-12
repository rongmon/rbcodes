"""
Synthetic data with a known answer, for tests and tutorials.

    from rbcodes.turbulence.demo import demo_velocity_field, demo_cube
    flux, vel, vturb, gal = demo_velocity_field()
    cube, truth = demo_cube('OII', z=0.785273)      # IFUCube-like object + the true maps

The demo nebula is an elongated exponential flux distribution with a velocity gradient plus a
turbulent field (alpha = 2/3, L = 10 kpc, rms 20 km/s after the flux-weighted PSF average);
the cube is continuum-subtracted with Gaussian noise and a variance cube.
Ported from vsf_pipeline (``demo_velocity_field``, ``demo_cube``).
"""
import numpy as np
from scipy.ndimage import gaussian_filter

from .forward import grf, psf_weight

__all__ = ['DEMO_TRUTH', 'demo_velocity_field', 'demo_cube']

DEMO_TRUTH = dict(alpha=2 / 3, L_kpc=10., sig_turb=20., grad=(-7.0, -5.8), centre=(24., 21.), seed=7)


def demo_velocity_field(shape=(40, 50), kpp=1.49, psf_sig=1.486, truth=DEMO_TRUTH):
    """
    Returns (flux, vel, vturb, gal): flux map, total velocity (gradient + turbulence), the
    turbulent part, and a small 'galaxy' mask at the centre. ``truth['seed']`` sets the field.
    """
    rng = np.random.default_rng(truth['seed'])
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    cx, cy = truth['centre']
    th = np.deg2rad(35.)
    xr = (xx - cx) * np.cos(th) + (yy - cy) * np.sin(th)
    yr = -(xx - cx) * np.sin(th) + (yy - cy) * np.cos(th)
    flux = gaussian_filter(150 * np.exp(-np.hypot(xr / 11., yr / 6.)), psf_sig)
    g = grf(truth['alpha'], truth['L_kpc'], shape, kpp, rng)
    vt = psf_weight(g, flux, psf_sig)
    vt = vt - vt[flux > 3].mean()
    vt *= truth['sig_turb'] / vt[flux > 3].std()
    vel = truth['grad'][0] * (xx - cx) + truth['grad'][1] * (yy - cy) + vt
    gal = np.hypot(xr / 2.5, yr / 1.8) < 1
    return flux, vel, vt, gal


def demo_cube(line, z, shape=(40, 50), kpp=1.49, psf_sig=1.486, air=True, flux_scale=1.0, seed=11, sig_int=60.,
              lsf=None, pix_arcsec=0.2, truth=DEMO_TRUTH):
    """
    Synthetic continuum-subtracted cube of ``line`` (a key of ``rbcodes.ifu.lines.LINES``) with the
    demo velocity field. Doublets with a free ratio get R = 1.35; intrinsic dispersion ``sig_int``;
    LSF default MUSE. Noise sigma 0.6 x (1 + 0.3 U(0,1)) per spaxel; the variance cube is exact.

    Returns
    -------
    cube : IFUCube with flux, var, wave, wcs (TAN, ``pix_arcsec``), name
    truth : dict(vturb, vel, gal, flux, sig_int, window) — window is the (lo, hi) range of the cube
    """
    from astropy.wcs import WCS
    from rbcodes.ifu.io.base import IFUCube
    from rbcodes.ifu.lines import LineModel, lsf_muse

    rng = np.random.default_rng(seed + sum(map(ord, line)))
    flux, vel, vt, gal = demo_velocity_field(shape, kpp, psf_sig, truth)
    rest = LineModel(line, z, dict(fit=(0., np.inf)), air=air).obs
    lo, hi = rest.min() * 0.985, rest.max() * 1.015
    m = LineModel(line, z, dict(fit=(lo, hi)), air=air, lsf=lsf_muse if lsf is None else lsf)
    wave = np.arange(lo, hi, 1.25)
    pars = [0., 1., 1.35, 0., sig_int] if m.kind == 'free' else [0., 1., 0., sig_int]
    data = np.zeros((len(wave),) + tuple(shape))
    for iy in range(shape[0]):
        for ix in range(shape[1]):
            p = list(pars)
            p[1] = flux[iy, ix] * flux_scale
            p[-2] = vel[iy, ix]
            data[:, iy, ix] = m(wave, *p)
    nmap = 0.6 * (1 + 0.3 * rng.random(shape))
    data += rng.normal(size=data.shape) * nmap
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [150., 2.]
    w.wcs.crpix = [shape[1] / 2 + 0.5, shape[0] / 2 + 0.5]
    w.wcs.cdelt = [-pix_arcsec / 3600., pix_arcsec / 3600.]
    cube = IFUCube()
    cube.flux, cube.var, cube.wave, cube.wcs = data, np.broadcast_to(nmap ** 2, data.shape).copy(), wave, w
    cube.name, cube.path, cube.header = f'demo_{line}', '', None
    return cube, dict(vturb=vt, vel=vel, gal=gal, flux=flux, sig_int=sig_int, window=(lo, hi))
