"""Tests for rbcodes.lensing.lens_model against a singular isothermal sphere (SIS)."""
import numpy as np
import pytest

pytestmark = pytest.mark.filterwarnings("ignore:divide by zero:RuntimeWarning", "ignore:invalid value:RuntimeWarning")
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS

from rbcodes.ifu.vmap import VelocityMap
from rbcodes.lensing.lens_model import LensModel, match_images

ZL, ZS = 0.43, 2.76
PIX = 0.1               # arcsec / pixel of the deflection maps
N = 201
CENTRE = (N - 1) / 2
RA0, DEC0 = 150.0, 2.0


def _wcs(pix, n):
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [RA0, DEC0]
    w.wcs.crpix = [(n - 1) / 2 + 1, (n - 1) / 2 + 1]
    w.wcs.cdelt = [-pix / 3600., pix / 3600.]
    return w


def _write_sis(path, label, theta_e):
    w = _wcs(PIX, N)
    y, x = np.mgrid[0:N, 0:N]
    dx, dy = x - CENTRE, y - CENTRE
    r = np.hypot(dx, dy)
    r[r == 0] = 1.
    fits.PrimaryHDU(theta_e * dx / r, header=w.to_header()).writeto(path / f'dplx_{label}.fits')
    fits.PrimaryHDU(theta_e * dy / r, header=w.to_header()).writeto(path / f'dply_{label}.fits')


@pytest.fixture(scope='module')
def lensdir(tmp_path_factory):
    d = tmp_path_factory.mktemp('lens')
    _write_sis(d, '1', 2.0)          # best fit
    _write_sis(d, '2', 2.2)
    _write_sis(d, '10', 1.8)
    return d


def _vm(n=41, pix=0.2):
    """Velocity map on a 0.2" grid centred on the lens; v = x offset (arcsec)."""
    w = _wcs(pix, n)
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    v = (xx - (n - 1) / 2) * pix
    return VelocityMap(v, np.ones_like(v), np.ones_like(v, bool), wcs=w, kpc_per_pix=1.0, line='MgII')


def _kpc(z):
    return FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0.).kpc_proper_per_arcmin(z).value / 60.


def test_from_directory_orders_best_first(lensdir):
    lens = LensModel.from_directory(str(lensdir), ZL, ZS)
    assert lens.n_realizations == 3
    assert lens.files[0][0].endswith('dplx_1.fits') and lens.files[1][0].endswith('dplx_2.fits')
    assert lens.files[2][0].endswith('dplx_10.fits')
    with pytest.raises(FileNotFoundError):
        LensModel.from_directory(str(lensdir), ZL, ZS, best='7')


def test_raytrace_matches_sis(lensdir):
    lens = LensModel(str(lensdir / 'dplx_1.fits'), str(lensdir / 'dply_1.fits'), ZL, ZS)
    vm = _vm()
    out = lens.raytrace(vm, ZS, centre=(RA0, DEC0), interp='bilinear')
    # image-plane offsets (arcsec, x towards the west = +pixel x)
    tx, ty = (vm.x - 20) * 0.2, (vm.y - 20) * 0.2
    r = np.hypot(tx, ty)
    far = r > 0.5
    k = _kpc(ZS)
    assert lens.kpc_per_arcsec(ZS) == pytest.approx(k, rel=1e-10)
    np.testing.assert_allclose(out.beta_x[far], (tx * (1 - 2.0 / r))[far] * k, atol=0.02 * k)
    np.testing.assert_allclose(out.beta_y[far], (ty * (1 - 2.0 / r))[far] * k, atol=0.02 * k)
    np.testing.assert_array_equal(out.v, vm.v)                     # velocities stay with their spaxels
    assert out.meta['lens']['realization'] == 0 and 'ray traced' in out.provenance[-1]
    # nearest-pixel lookup agrees to within the deflection-map pixel
    near = lens.raytrace(vm, ZS, centre=(RA0, DEC0))
    assert np.nanmax(np.abs(near.beta_x - out.beta_x)[r > 1]) < 0.1 * k


def test_lower_plane_scales_deflection(lensdir):
    lens = LensModel(str(lensdir / 'dplx_1.fits'), str(lensdir / 'dply_1.fits'), ZL, ZS)
    z = 1.70
    s = lens.model().deflection_scale(z)
    assert 0 < s < 1
    out = lens.raytrace(_vm(), z, centre=(RA0, DEC0), interp='bilinear')
    tx = (out.x - 20) * 0.2
    ty = (out.y - 20) * 0.2
    r = np.hypot(tx, ty)
    far = r > 0.5
    np.testing.assert_allclose(out.beta_x[far], (tx * (1 - s * 2.0 / r))[far] * _kpc(z), atol=0.02 * _kpc(z))


def test_realizations(lensdir):
    lens = LensModel.from_directory(str(lensdir), ZL, ZS)
    vm = _vm()
    outs = list(lens.raytrace_realizations(vm, ZS, interp='bilinear'))
    assert [k for k, _ in outs] == [0, 1, 2]
    centres = {tuple(o.meta['lens']['centre']) for _, o in outs}
    assert len(centres) == 1                                         # same origin for all
    b0 = outs[0][1]                                                  # default origin: mean source position (best fit)
    assert abs(np.nanmean(b0.beta_x[b0.mask])) < 1e-6 and abs(np.nanmean(b0.beta_y[b0.mask])) < 1e-6
    # a spaxel at +3" on the x axis maps to 3 - theta_E
    iy, ix = 20, 35
    k = _kpc(ZS)
    c = outs[0][1].meta['lens']['centre']
    for (_, o), te in zip(outs, (2.0, 2.2, 1.8)):
        assert o.beta_x[iy, ix] == pytest.approx((3.0 - te) * k + outs[0][1].beta_x[20, 20] * 0, abs=0.02 * k)


def test_magnification_matches_sis(lensdir):
    lens = LensModel(str(lensdir / 'dplx_1.fits'), str(lensdir / 'dply_1.fits'), ZL, ZS)
    w = _wcs(PIX, N)
    for r in (1.2, 3.0, 5.0):
        for ang in (0.3, 2.0):
            px = CENTRE + r / PIX * np.cos(ang)
            py = CENTRE + r / PIX * np.sin(ang)
            ra, dec = w.pixel_to_world_values(px, py)
            for z in (ZS, 1.7):
                s = lens.model().deflection_scale(z)
                m = lens.magnification(ra, dec, z, interp='bilinear')
                # central differences on the 0.1" map: ~(0.1/r)^2 errors, amplified near the critical curve
                assert float(m['mu']) == pytest.approx(r / (r - s * 2.0), rel=0.02), (r, ang, z)
                assert float(m['kappa']) == pytest.approx(s * 2.0 / (2 * r), rel=0.02)
                # nearest map pixel: r off by <= 0.071", so mu off by <= |dmu/dr| x 0.071
                near = lens.magnification(ra, dec, z)
                dmu = 2.0 * s / (r - s * 2.0) ** 2 * 0.071
                assert float(near['mu']) == pytest.approx(r / (r - s * 2.0), abs=dmu + 0.02 * abs(r / (r - s * 2.0)))
    assert np.isnan(lens.magnification(RA0 + 1, DEC0, ZS)['mu'])        # outside the maps


def test_images_and_matching(lensdir):
    lens = LensModel(str(lensdir / 'dplx_1.fits'), str(lensdir / 'dply_1.fits'), ZL, ZS)
    n, pix = 81, 0.1                                                # +-4": both images of sources up to 1.7"
    w = _wcs(pix, n)
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    tx, ty = (xx - 40) * pix, (yy - 40) * pix
    r = np.hypot(tx, ty)
    r[r == 0] = 1
    bx = tx * (1 - 2.0 / r)                                          # true source position (arcsec)
    v = 50. * bx                                                     # velocity field fixed in the source plane
    keep = np.abs(ty) < 0.05                                         # a slit along x
    vm = VelocityMap(v, np.ones_like(v), keep, wcs=w)
    img = {1: (tx > 0.3), 2: (tx < -0.3)}
    out = lens.raytrace(vm, ZS, images=img, centre=(RA0, DEC0), interp='bilinear')
    assert set(np.unique(out.image)) == {-1, 1, 2}
    m = match_images(out, tol=0.01 * _kpc(ZS))
    assert len(m['a']) >= 15 and np.all(m['image_a'] != m['image_b'])
    assert np.all(np.abs(m['pull']) < 2)                             # same source point, same velocity
    with pytest.raises(ValueError):
        lens.raytrace(vm, ZS, images={1: tx > 0, 2: tx > -1})        # overlapping images
    with pytest.raises(ValueError):
        match_images(vm, tol=1.)                                     # no source positions
