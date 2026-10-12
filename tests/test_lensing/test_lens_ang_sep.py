"""
Tests for rbcodes.lensing.lens_ang_sep.

Distances are checked against astropy; ray tracing against a singular
isothermal sphere (SIS) whose deflection field is known analytically.
"""
import numpy as np
import pytest
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS

from rbcodes._selftest import angular_diameter_distance_12 as d_a12
from rbcodes.lensing import lens_ang_sep as las

COSMO = FlatLambdaCDM(H0=70., Om0=0.3, Tcmb0=0.)
MPC_M = 3.0856775814913673e22
# lens_ang_sep uses c/H0 = 9.26e25/h m, 0.1% above the exact value
RTOL_D = 2e-3

ZL, ZS = 0.564, 1.701
THETA_E = 2.0         # arcsec
PIX = 0.1             # arcsec / pixel of the deflection maps
N = 201
CENTRE = (N - 1) / 2  # 0-based pixel of the lens centre


@pytest.mark.parametrize('z', [0.1, 0.564, 1.701, 3.0])
def test_cosmic_d_matches_astropy(z):
    out = las.cosmic_D(0.3, 0.7, z)
    assert out['D_C'] / MPC_M == pytest.approx(COSMO.comoving_distance(z).value, rel=RTOL_D)
    assert out['D_A'] / MPC_M == pytest.approx(COSMO.angular_diameter_distance(z).value, rel=RTOL_D)
    assert out['D_L'] / MPC_M == pytest.approx(COSMO.luminosity_distance(z).value, rel=RTOL_D)


@pytest.mark.parametrize('z1,z2', [(0.392, 2.76), (0.564, 1.701), (0.564, 0.8)])
def test_ang_d12_matches_astropy(z1, z2):
    got = las.ang_D12(0.3, 0.7, z1, z2) / MPC_M
    assert got == pytest.approx(d_a12(COSMO, z1, z2).value, rel=RTOL_D)


@pytest.mark.parametrize('sep_arcsec', [1e-3, 0.1, 1., 100., 3600. * 30, 3600. * 170])
def test_ang_sep_d_matches_astropy(sep_arcsec):
    from astropy.coordinates import SkyCoord
    rng = np.random.default_rng(1)
    ra1, dec1 = rng.uniform(0, 360, 500), rng.uniform(-89, 89, 500)
    ra2 = ra1 + rng.normal(0, sep_arcsec / 3600., 500)
    dec2 = np.clip(dec1 + rng.normal(0, sep_arcsec / 3600., 500), -90, 90)
    ref = SkyCoord(ra1, dec1, unit='deg').separation(SkyCoord(ra2, dec2, unit='deg')).rad
    np.testing.assert_allclose(las.ang_sep_D(ra1, dec1, ra2, dec2), ref, rtol=1e-12)


@pytest.fixture(scope='module')
def sis_maps(tmp_path_factory):
    """SIS deflection maps in arcsec along the pixel axes, with a TAN WCS."""
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [150.0, 2.0]
    w.wcs.crpix = [CENTRE + 1, CENTRE + 1]
    w.wcs.cdelt = [-PIX / 3600., PIX / 3600.]
    y, x = np.mgrid[0:N, 0:N]
    dx, dy = x - CENTRE, y - CENTRE
    r = np.hypot(dx, dy)
    r[r == 0] = 1.
    ax = THETA_E * dx / r
    ay = THETA_E * dy / r
    d = tmp_path_factory.mktemp('lens')
    fx, fy = d / 'dplx.fits', d / 'dply.fits'
    fits.PrimaryHDU(ax, header=w.to_header()).writeto(fx)
    fits.PrimaryHDU(ay, header=w.to_header()).writeto(fy)
    return str(fx), str(fy), w


def _trace(sis_maps, px, py, zs_n):
    fx, fy, w = sis_maps
    ra, dec = w.pixel_to_world_values(px, py)
    s = las.model_delens(fx, fy, zl=ZL, zs_o=ZS)
    s.raytrace_new_z(ra, dec, zs_n=zs_n)
    return w.world_to_pixel_values(s.src_ra, s.src_dec)


def test_sis_raytrace_at_model_redshift(sis_maps):
    # integer pixels, so the nearest-pixel lookup is exact
    px = np.array([CENTRE + 40., CENTRE - 25., CENTRE + 10.])
    py = np.array([CENTRE, CENTRE + 30., CENTRE - 60.])
    bx, by = _trace(sis_maps, px, py, zs_n=ZS)
    dx, dy = px - CENTRE, py - CENTRE
    r = np.hypot(dx, dy)
    shrink = 1. - (THETA_E / PIX) / r
    np.testing.assert_allclose(bx - CENTRE, dx * shrink, atol=1e-6)
    np.testing.assert_allclose(by - CENTRE, dy * shrink, atol=1e-6)


def test_einstein_ring_maps_to_centre(sis_maps):
    rE = THETA_E / PIX
    bx, by = _trace(sis_maps, np.array([CENTRE + rE]), np.array([CENTRE]), zs_n=ZS)
    assert bx[0] == pytest.approx(CENTRE, abs=1e-6)
    assert by[0] == pytest.approx(CENTRE, abs=1e-6)


def test_rescaling_to_intermediate_plane(sis_maps):
    z_new = 0.9
    scale = ((d_a12(COSMO, ZL, z_new) / COSMO.angular_diameter_distance(z_new)) /
             (d_a12(COSMO, ZL, ZS) / COSMO.angular_diameter_distance(ZS))).value
    px, py = np.array([CENTRE + 50.]), np.array([CENTRE])
    bx, _ = _trace(sis_maps, px, py, zs_n=z_new)
    expected = 50. - scale * THETA_E / PIX
    assert bx[0] - CENTRE == pytest.approx(expected, rel=RTOL_D)


def test_plane_at_lens_has_no_deflection(sis_maps):
    px, py = np.array([CENTRE + 50.]), np.array([CENTRE + 20.])
    bx, by = _trace(sis_maps, px, py, zs_n=ZL)
    assert bx[0] == pytest.approx(px[0], abs=1e-6)
    assert by[0] == pytest.approx(py[0], abs=1e-6)


def _write_maps(tmp_path, ax, ay, header):
    fx, fy = tmp_path / 'dplx.fits', tmp_path / 'dply.fits'
    fits.PrimaryHDU(ax, header=header).writeto(fx)
    fits.PrimaryHDU(ay, header=header).writeto(fy)
    return str(fx), str(fy)


def test_raytrace_returns_and_stores(sis_maps):
    fx, fy, w = sis_maps
    ra, dec = w.pixel_to_world_values(np.array([CENTRE + 40.]), np.array([CENTRE]))
    s = las.model_delens(fx, fy, zl=ZL, zs_o=ZS)
    out = s.raytrace_new_z(ra, dec, zs_n=ZS)
    np.testing.assert_array_equal(out[0], s.src_ra)
    np.testing.assert_array_equal(out[1], s.src_dec)


def test_outside_map_is_nan(sis_maps):
    fx, fy, w = sis_maps
    px = np.array([-5., CENTRE + 40., N + 3., CENTRE])
    py = np.array([CENTRE, CENTRE, CENTRE, -0.6])
    ra, dec = w.pixel_to_world_values(px, py)
    s = las.model_delens(fx, fy, zl=ZL, zs_o=ZS)
    for interp in ('nearest', 'bilinear'):
        sra, sdec = s.raytrace_new_z(ra, dec, zs_n=ZS, interp=interp)
        assert np.isnan(sra[[0, 2, 3]]).all() and np.isnan(sdec[[0, 2, 3]]).all()
        assert np.isfinite(sra[1]) and np.isfinite(sdec[1])


def test_nearest_uses_containing_pixel(sis_maps):
    # +-0.4 px offsets stay in the same pixel: identical deflection
    fx, fy, w = sis_maps
    s = las.model_delens(fx, fy, zl=ZL, zs_o=ZS)
    ax0, ay0 = s.deflection_at(np.array([CENTRE + 40.]), np.array([CENTRE + 7.]))
    for d in (-0.4, 0.4):
        ax, ay = s.deflection_at(np.array([CENTRE + 40. + d]), np.array([CENTRE + 7. - d]))
        assert ax[0] == ax0[0] and ay[0] == ay0[0]


def test_bilinear_exact_for_linear_field(tmp_path):
    # alpha linear in pixel coordinates -> bilinear interpolation is exact
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [150.0, 2.0]
    w.wcs.crpix = [26., 26.]
    w.wcs.cdelt = [-PIX / 3600., PIX / 3600.]
    y, x = np.mgrid[0:51, 0:51].astype(float)
    fx, fy = _write_maps(tmp_path, 0.02 * x - 0.01 * y + 0.3, 0.015 * y + 0.005 * x, w.to_header())
    s = las.model_delens(fx, fy, zl=ZL, zs_o=ZS)
    rng = np.random.default_rng(3)
    px, py = rng.uniform(2, 48, 50), rng.uniform(2, 48, 50)
    ax, ay = s.deflection_at(px, py, interp='bilinear')
    np.testing.assert_allclose(ax * PIX, 0.02 * px - 0.01 * py + 0.3, rtol=1e-10)
    np.testing.assert_allclose(ay * PIX, 0.015 * py + 0.005 * px, rtol=1e-10)


def test_cd_matrix_header_matches_cdelt(sis_maps, tmp_path):
    fx, fy, w = sis_maps
    hdr = w.to_header()
    cd = fits.Header()
    for k in ('CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2'):
        cd[k] = hdr[k]
    cd['CD1_1'], cd['CD1_2'] = -PIX / 3600., 0.
    cd['CD2_1'], cd['CD2_2'] = 0., PIX / 3600.
    gx, gy = _write_maps(tmp_path, fits.getdata(fx), fits.getdata(fy), cd)
    px, py = np.array([CENTRE + 40., CENTRE - 25.]), np.array([CENTRE, CENTRE + 30.])
    ra, dec = w.pixel_to_world_values(px, py)
    a = las.model_delens(fx, fy, zl=ZL, zs_o=ZS).raytrace_new_z(ra, dec, zs_n=ZS)
    b = las.model_delens(gx, gy, zl=ZL, zs_o=ZS).raytrace_new_z(ra, dec, zs_n=ZS)
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)


def test_invalid_interp_raises(sis_maps):
    fx, fy, _ = sis_maps
    s = las.model_delens(fx, fy, zl=ZL, zs_o=ZS)
    with pytest.raises(ValueError):
        s.deflection_at(np.array([1.]), np.array([1.]), interp='cubic')
