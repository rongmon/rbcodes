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


class TestPartialAndVacuum:

    @staticmethod
    def _write(tmp_path, ctype='AWAV', name='c.fits'):
        rng = np.random.default_rng(0)
        flux, var = rng.normal(size=(NW, NY, NX)), rng.random((NW, NY, NX)) + 0.5
        h = _cube_header()
        h['CTYPE3'] = ctype
        p = tmp_path / name
        fits.HDUList([fits.PrimaryHDU(header=fits.Header({'INSTRUME': 'MUSE'})),
                      fits.ImageHDU(flux, header=h, name='DATA'),
                      fits.ImageHDU(var, header=h, name='STAT')]).writeto(p)
        return p, flux, var

    def test_wave_range_reads_a_slab(self, tmp_path):
        from rbcodes.ifu.io.auto_cube import _build_wave
        p, flux, var = self._write(tmp_path)
        full = load_fits(p)
        part = load_fits(p, wave_range=(6540., 6560.))
        sel = (full.wave >= 6540.) & (full.wave <= 6560.)
        np.testing.assert_array_equal(part.wave, full.wave[sel])
        np.testing.assert_array_equal(part.flux, flux[sel])
        np.testing.assert_array_equal(part.var, var[sel])
        np.testing.assert_allclose(_build_wave(part.header), part.wave)       # header matches the slab
        assert part.header['NAXIS3'] == sel.sum() and full.header['NAXIS3'] == NW
        np.testing.assert_array_equal(full.flux, flux)                         # default: everything
        with pytest.raises(ValueError):
            load_fits(p, wave_range=(9000., 9100.))

    def test_wave_range_with_external_variance(self, tmp_path):
        p, flux, var = self._write(tmp_path)
        pv = tmp_path / 'v.fits'
        fits.PrimaryHDU(2 * var).writeto(pv)
        part = load_fits(p, var=pv, wave_range=(6510., 6520.))
        sel = (_wave() >= 6510.) & (_wave() <= 6520.)
        np.testing.assert_array_equal(part.var, 2 * var[sel])

    def test_air_to_vacuum(self, tmp_path):
        from rbcodes.utils.rb_spectrum import air_to_vacuum
        p, _, _ = self._write(tmp_path)
        air = load_fits(p)
        vac = load_fits(p, air_to_vacuum=True, wave_range=(6540., 6560.))
        assert air.airvac == 'air' and vac.airvac == 'vac'
        sel = (_wave() >= 6540.) & (_wave() <= 6560.)
        np.testing.assert_allclose(vac.wave, air_to_vacuum(_wave()[sel]))
        assert np.all(vac.wave - _wave()[sel] > 1.7)                           # ~1.8 A at 6550 A
        pw, _, _ = self._write(tmp_path, ctype='WAVE', name='vac.fits')
        with pytest.warns(UserWarning, match='already in vacuum'):
            already = load_fits(pw, air_to_vacuum=True)
        np.testing.assert_array_equal(already.wave, _wave()) and already.airvac == 'vac'


def test_rb_spectrum_air_vacuum_functions():
    from astropy import units as u
    from rbcodes.utils.rb_spectrum import air_to_vacuum, rb_spectrum, vacuum_to_air
    w = np.array([1500., 3726.032, 5006.843, 6562.80, 8000.])
    vac = air_to_vacuum(w)
    assert vac[0] == 1500. and vac[3] == pytest.approx(6564.61, abs=0.01)
    np.testing.assert_allclose(vacuum_to_air(vac), w, atol=2e-3)
    s = rb_spectrum.from_tuple((w * u.AA, np.ones(5), np.ones(5)))
    s.meta['airvac'] = 'air'                                                  # from_tuple assumes vacuum
    s.air2vac()
    np.testing.assert_array_equal(s.wavelength.to(u.AA).value, vac)          # the method uses the function
    s.vac2air()
    np.testing.assert_array_equal(s.wavelength.to(u.AA).value, vacuum_to_air(vac))


class TestBinCube:

    def _cube(self, tmp_path, cd=False):
        rng = np.random.default_rng(3)
        flux, var = rng.normal(size=(NW, NY, NX)), rng.random((NW, NY, NX)) + 0.5
        h = _cube_header()
        if cd:
            for k in ('CDELT1', 'CDELT2'):
                del h[k]
            h['CD1_1'], h['CD2_2'], h['CD1_2'], h['CD2_1'] = -PIX, PIX, 0., 0.
        p = tmp_path / ('cd.fits' if cd else 'c.fits')
        fits.HDUList([fits.PrimaryHDU(header=fits.Header({'INSTRUME': 'MUSE'})), fits.ImageHDU(flux, h, name='DATA'),
                      fits.ImageHDU(var, h, name='STAT')]).writeto(p)
        return load_fits(p), flux, var

    @pytest.mark.parametrize('cd', [False, True])
    def test_sums_and_wcs(self, tmp_path, cd):
        from rbcodes.ifu.binning import bin_cube
        cube, flux, var = self._cube(tmp_path, cd)
        b = bin_cube(cube, 3)
        assert b.flux.shape == (NW, NY // 3, NX // 3) and b.binning == 3
        np.testing.assert_allclose(b.flux[:, 1, 2], flux[:, 3:6, 6:9].sum((1, 2)))
        np.testing.assert_allclose(b.var[:, 1, 2], var[:, 3:6, 6:9].sum((1, 2)))
        np.testing.assert_array_equal(b.wave, cube.wave)
        # the centre of a new pixel is the centre of its block
        for iy, ix in ((0, 0), (1, 2), (3, 4)):
            got = b.wcs.pixel_to_world_values(ix, iy)
            ref = cube.wcs.pixel_to_world_values(3 * ix + 1, 3 * iy + 1)
            np.testing.assert_allclose(got, ref, atol=1e-10)
        from rbcodes.ifu.io.auto_cube import _make_2d_wcs
        hw = _make_2d_wcs(b.header)                                  # header agrees with b.wcs
        np.testing.assert_allclose(hw.pixel_to_world_values(2, 1), b.wcs.pixel_to_world_values(2, 1), atol=1e-10)

    def test_nan_blocks_and_noise(self, tmp_path):
        from rbcodes.ifu.binning import bin_cube
        cube, flux, var = self._cube(tmp_path)
        cube.flux[:, 0, 0] = np.nan
        b = bin_cube(cube, 2)
        assert np.all(np.isnan(b.flux[:, 0, 0])) and np.all(np.isfinite(b.flux[:, 0, 1]))
        part = bin_cube(cube, 2, require_all=False)
        np.testing.assert_allclose(part.flux[:, 0, 0], np.nansum(cube.flux[:, 0:2, 0:2], (1, 2)))
        # independent noise: binned variance = sum, so sum / sqrt(var) keeps unit scatter
        rng = np.random.default_rng(4)
        cube.var = np.full((NW, NY, NX), 4.)
        cube.flux = rng.normal(0, 2., (NW, NY, NX))
        b = bin_cube(cube, 2)
        assert np.std(b.flux / np.sqrt(b.var)) == pytest.approx(1., abs=0.05)
        with pytest.raises(ValueError):
            bin_cube(cube, 0)
