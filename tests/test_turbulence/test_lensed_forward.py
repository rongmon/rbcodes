"""Tests for the source-plane (lensed) forward model in rbcodes.turbulence.forward."""
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from rbcodes.ifu.vmap import VelocityMap
from rbcodes.lensing.lens_model import LensModel, match_images
from rbcodes.turbulence.forward import Observation, TurbulenceField
from rbcodes.turbulence.vsf import StructureFunction

N = 40


class Ramp(TurbulenceField):
    """Deterministic 'field': the source-grid column index."""
    def draw(self, shape, kpp, rng):
        return dict(g=np.tile(np.arange(shape[1], dtype=float), (shape[0], 1)))


def _grid():
    yy, xx = np.mgrid[0:N, 0:N].astype(float)
    return yy, xx


def _obs(beta, mask=None, **kw):
    yy, xx = _grid()
    mask = np.ones((N, N), bool) if mask is None else mask
    return Observation(mask, np.ones((N, N)), 1.0, 1e-3, np.full((N, N), 2.), beta=beta, **kw)


def test_sampling_at_beta_is_exact():
    yy, xx = _grid()
    bx, by = 0.37 * xx + 0.1 * yy, 0.5 * yy                      # a sheared, magnified mapping
    obs = _obs((bx, by), src_kpp=0.2)
    m = obs(Ramp(10., 1., None), np.random.default_rng(0), noise=False)
    v, b = m.v[m.mask], bx[m.mask]
    assert np.corrcoef(v, b)[0, 1] == pytest.approx(1., abs=1e-10)
    assert m.has_source_positions and m.meta['mock']['lensed']
    np.testing.assert_array_equal(m.beta_x, bx)


def test_multiple_images_share_the_field():
    yy, xx = _grid()
    xc = xx - (N - 1) / 2
    bx, by = np.abs(xc) - 5., yy * 1.0                            # x<0 and x>0 halves map to the same beta
    img = np.where(xc > 0, 1, 2)
    keep = np.abs(xc) > 0.6
    obs = _obs((bx, by), mask=keep, image=img)
    m = obs(TurbulenceField(10., 1.0, None), np.random.default_rng(1), noise=False)
    left, right = m.v[:, :N // 2][:, ::-1], m.v[:, N // 2:]
    np.testing.assert_allclose(left, right, atol=1e-9)
    assert set(np.unique(m.image)) == {1, 2}


def test_source_plane_statistics():
    # beta = image positions / 2: the source-plane S2 of lensed mocks must match that of the same
    # field made directly (unlensed) on a 0.5 kpc grid
    yy, xx = _grid()
    lensed = _obs((0.5 * xx, 0.5 * yy))
    assert lensed.src_kpp == pytest.approx(0.5)
    direct = Observation(np.ones((N, N), bool), np.ones((N, N)), 0.5, 1e-3, np.full((N, N), 2.))
    slopes = {}
    for name, obs, pos in (('lensed', lensed, 'source'), ('direct', direct, 'image')):
        S = []
        for seed in range(12):
            m = obs(TurbulenceField(10., 1.0, None), np.random.default_rng(seed), noise=False)
            sf = StructureFunction(m, n_bins=12, positions=pos)
            S.append(sf.S_raw)
        S = np.nanmean(S, 0)
        k = sf.ok & (sf.rc > 1.0) & (sf.rc < 5.)
        slopes[name] = np.polyfit(np.log(sf.rc[k]), np.log(S[k]), 1)[0]
    assert slopes['lensed'] == pytest.approx(slopes['direct'], abs=0.1)
    assert slopes['lensed'] == pytest.approx(1.0, abs=0.2)


def test_from_raytraced_map_and_matching(tmp_path):
    # SIS lens (theta_E = 2") on 0.1" maps; image grid 0.1" with velocities from a lensed mock
    n, pix = 201, 0.1
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [150., 2.]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-pix / 3600., pix / 3600.]
    y, x = np.mgrid[0:n, 0:n] - (n - 1) / 2
    r = np.hypot(x, y)
    r[r == 0] = 1
    for nm, a in (('dplx', 2.0 * x / r), ('dply', 2.0 * y / r)):
        fits.PrimaryHDU(a, header=w.to_header()).writeto(tmp_path / f'{nm}_1.fits')
    lens = LensModel.from_directory(str(tmp_path), zl=0.43, zs=2.76)
    m_img, pimg = 61, 0.1
    wi = WCS(naxis=2)
    wi.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    wi.wcs.crval = [150., 2.]
    wi.wcs.crpix = [(m_img + 1) / 2, (m_img + 1) / 2]
    wi.wcs.cdelt = [-pimg / 3600., pimg / 3600.]
    yy, xx = np.mgrid[0:m_img, 0:m_img] - (m_img - 1) / 2
    slit = np.abs(yy) < 0.5
    vm = VelocityMap(np.zeros((m_img, m_img)), np.ones((m_img, m_img)), slit, wcs=wi, kpc_per_pix=0.8,
                     psf_fwhm_kpc=1e-3, flux=np.ones((m_img, m_img)))
    traced = lens.raytrace(vm, 2.76, images={1: xx * pimg > 0.3, 2: xx * pimg < -0.3}, interp='bilinear')
    obs = Observation.from_vmap(traced, positions='source')
    mock = obs(TurbulenceField(10., 1.0, None), np.random.default_rng(2), noise=False)
    pairs = match_images(mock, tol=0.05 * lens.kpc_per_arcsec(2.76))
    assert len(pairs['a']) > 5 and np.all(np.abs(pairs['dv']) < 1.0)   # same beta, same velocity
    with pytest.raises(ValueError):
        Observation.from_vmap(vm, positions='source')


def test_3d_and_kpp_override():
    yy, xx = _grid()
    obs = _obs((0.5 * xx, 0.5 * yy))
    m = obs(TurbulenceField(20., 0.8, 20., mode='3D', depth_kpc=5.), np.random.default_rng(3), noise=False)
    assert np.all(np.isfinite(m.v[m.mask])) and np.all(m.width[m.mask] >= 2.3)
    with pytest.raises(ValueError):
        obs(TurbulenceField(10., 1., None), np.random.default_rng(4), kpp=2.)


def test_grid_guards():
    yy, xx = _grid()
    with pytest.raises(ValueError):                               # would need a 4e8-pixel source grid
        _obs((xx, yy), src_kpp=1e-3)(TurbulenceField(10., 1., None), np.random.default_rng(5))
    with pytest.raises(ValueError):
        _obs((np.zeros((N, N)), np.zeros((N, N))))                # no two distinct positions
