"""Tests for rbcodes.ifu (cube I/O, collapse, moments, extraction, region masks)."""
import numpy as np
import pytest
from astropy.io import fits

from rbcodes.ifu import load_fits, load_image
from rbcodes.ifu.io.auto_cube import GenericCube, KCWICube, MUSECube
from rbcodes.ifu.io.image2d import FITSImage
from rbcodes.ifu.cube_collapse import build_whitelight, build_narrowband
from rbcodes.ifu.moment_maps import C_KMS, moment0, moment1, moment2, moment_map
from rbcodes.ifu.aperture_extract import (
    make_circular_mask, extract_aperture, extract_variance_weighted)
from rbcodes.ifu.spatial_mask import parse_ds9_regions, region_to_mask

NW, NY, NX = 200, 12, 15
CRVAL3, CDELT3 = 6500., 0.5
RA0, DEC0, PIX = 150.0, 2.0, 0.2 / 3600.    # 0.2" pixels


def _cube_header(instrume=None):
    h = fits.Header()
    h['CTYPE1'], h['CTYPE2'], h['CTYPE3'] = 'RA---TAN', 'DEC--TAN', 'AWAV'
    h['CRVAL1'], h['CRVAL2'], h['CRVAL3'] = RA0, DEC0, CRVAL3
    h['CRPIX1'], h['CRPIX2'], h['CRPIX3'] = 8., 6., 1.
    h['CDELT1'], h['CDELT2'], h['CDELT3'] = -PIX, PIX, CDELT3
    h['CUNIT1'], h['CUNIT2'], h['CUNIT3'] = 'deg', 'deg', 'Angstrom'
    if instrume:
        h['INSTRUME'] = instrume
    return h


def _wave():
    return CRVAL3 + CDELT3 * np.arange(NW)


@pytest.fixture
def line_cube():
    """Gaussian emission line whose centre shifts linearly across x; constant width."""
    wave = _wave()
    lam0, sig_a = 6550., 2.0
    shift = np.linspace(-3., 3., NX)                 # Angstrom
    cen = lam0 + shift[None, None, :] * np.ones((1, NY, 1))
    flux = 1.0 + 10. * np.exp(-0.5 * ((wave[:, None, None] - cen) / sig_a) ** 2)
    var = np.full_like(flux, 0.04)
    return dict(flux=flux, var=var, wave=wave, lam0=lam0, sig_a=sig_a, shift=shift)


def _write_muse(path, flux, var):
    hdus = [fits.PrimaryHDU(header=fits.Header({'INSTRUME': 'MUSE'})),
            fits.ImageHDU(flux, header=_cube_header(), name='DATA'),
            fits.ImageHDU(var, header=_cube_header(), name='STAT')]
    fits.HDUList(hdus).writeto(path)


class TestLoad:

    def test_muse_layout(self, tmp_path, line_cube):
        p = tmp_path / 'muse.fits'
        _write_muse(p, line_cube['flux'], line_cube['var'])
        cube = load_fits(p)
        assert type(cube) is MUSECube
        np.testing.assert_allclose(cube.flux, line_cube['flux'])
        np.testing.assert_allclose(cube.var, line_cube['var'])
        np.testing.assert_allclose(cube.wave, _wave())
        assert cube.wcs.naxis == 2
        ra, dec = cube.wcs.all_pix2world(7., 5., 0)      # CRPIX is 1-based
        assert ra == pytest.approx(RA0) and dec == pytest.approx(DEC0)
        assert cube.name == 'muse'

    def test_kcwi_primary_and_generic(self, tmp_path, line_cube):
        for inst, cls in (('KCWI', KCWICube), (None, GenericCube)):
            p = tmp_path / f'{inst}.fits'
            fits.HDUList([fits.PrimaryHDU(line_cube['flux'], header=_cube_header(inst)),
                          fits.ImageHDU(line_cube['var'], name='VAR')]).writeto(p)
            cube = load_fits(p)
            assert type(cube) is cls
            np.testing.assert_allclose(cube.var, line_cube['var'])
            np.testing.assert_allclose(cube.wave, _wave())

    def test_external_variance_overrides(self, tmp_path, line_cube):
        p, pv = tmp_path / 'c.fits', tmp_path / 'v.fits'
        _write_muse(p, line_cube['flux'], line_cube['var'])
        fits.PrimaryHDU(2 * line_cube['var']).writeto(pv)
        np.testing.assert_allclose(load_fits(p, var=pv).var, 2 * line_cube['var'])

    def test_2d_image(self, tmp_path):
        h = _cube_header()
        for k in ('CTYPE3', 'CRVAL3', 'CRPIX3', 'CDELT3', 'CUNIT3'):
            del h[k]
        img = np.arange(NY * NX, dtype=float).reshape(NY, NX)
        p = tmp_path / 'img.fits'
        fits.PrimaryHDU(img, header=h).writeto(p)
        for loaded in (load_fits(p), load_image(p)):
            assert isinstance(loaded, FITSImage)
            assert loaded.wave is None
            np.testing.assert_allclose(loaded.flux, img)


class TestMoments:

    def test_velocity_and_dispersion(self, line_cube):
        c = line_cube
        win = (c['lam0'] - 20., c['lam0'] + 20.)
        m0 = moment0(c['flux'] - 1., c['wave'], *win)
        m1 = moment1(c['flux'] - 1., c['wave'], *win, lambda_rest=c['lam0'])
        m2 = moment2(c['flux'] - 1., c['wave'], *win, lambda_rest=c['lam0'])
        np.testing.assert_allclose(m0, 10. * c['sig_a'] * np.sqrt(2 * np.pi), rtol=1e-3)
        np.testing.assert_allclose(m1[0], C_KMS * c['shift'] / c['lam0'], atol=0.5)
        np.testing.assert_allclose(m2, C_KMS * c['sig_a'] / c['lam0'], rtol=2e-3)
        np.testing.assert_array_equal(moment_map(c['flux'] - 1., c['wave'], *win, 1,
                                                 lambda_rest=c['lam0']), m1)

    def test_empty_window_raises(self, line_cube):
        with pytest.raises(ValueError):
            moment0(line_cube['flux'], line_cube['wave'], 1000., 1001.)


class TestCollapseAndExtract:

    def test_whitelight_and_narrowband(self, line_cube):
        c = line_cube
        np.testing.assert_allclose(build_whitelight(c['flux'], c['wave']), c['flux'].mean(0))
        sel = (c['wave'] >= 6540.) & (c['wave'] <= 6560.)
        np.testing.assert_allclose(build_narrowband(c['flux'], c['wave'], 6540., 6560.),
                                   c['flux'][sel].mean(0))

    def test_aperture_extraction(self, line_cube):
        c = line_cube
        mask = make_circular_mask(NY, NX, 7, 5, 2.5)
        n = mask.sum()
        spec, err = extract_aperture(c['flux'], c['var'], mask)
        np.testing.assert_allclose(spec, c['flux'][:, mask].sum(1))
        np.testing.assert_allclose(err, np.sqrt(0.04 * n))
        wspec, werr = extract_variance_weighted(c['flux'], c['var'], mask)
        np.testing.assert_allclose(wspec, c['flux'][:, mask].mean(1))
        np.testing.assert_allclose(werr, np.sqrt(0.04 / n))
        with pytest.raises(ValueError):
            extract_variance_weighted(c['flux'], None, mask)


class TestRegions:

    def test_image_and_sky_circles_agree(self, tmp_path, line_cube):
        p = tmp_path / 'muse.fits'
        _write_muse(p, line_cube['flux'], line_cube['var'])
        wcs = load_fits(p).wcs
        # same circle: centre pixel (8,6) 1-based = CRVAL, radius 2.5 pix = 0.5"
        # (no pixel centre within 0.25 pix of the edge; icrs = the header's frame)
        img = parse_ds9_regions('image\ncircle(8,6,2.5)')
        sky = parse_ds9_regions(f'icrs\ncircle({RA0},{DEC0},0.5")')
        m_img = region_to_mask(img[0], wcs, NY, NX)
        m_sky = region_to_mask(sky[0], wcs, NY, NX)
        expected = make_circular_mask(NY, NX, 7, 5, 2.5)
        assert m_img.sum() > 0
        np.testing.assert_array_equal(m_img, m_sky)
        np.testing.assert_array_equal(m_img, expected)

    def test_region_outside_fov_is_empty(self, line_cube):
        reg = parse_ds9_regions('image\ncircle(500,500,3)')
        assert not region_to_mask(reg[0], None, NY, NX).any()
