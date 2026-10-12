"""Tests for rbcodes.ifu.linefit and rbcodes.ifu.noise against known answers."""
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from rbcodes.ifu import load_fits
from rbcodes.ifu.linefit import (
    C_KMS, EmissionLine, LineFitResult, fit_cube, fit_spectrum, lsf_constant, lsf_muse,
    lsf_resolution, plot_spaxel, quality_mask, spaxel_report, vac_to_air, REJECT)
from rbcodes.ifu.noise import noise_cube, robust_std, sky_mask_from_region

Z = 0.785273
NY, NX = 14, 16


def _wave():
    return np.arange(6600., 6720., 1.25)            # MUSE-like sampling around [O II] at z=0.785


def _oii_cube(rng, noise_sigma=0.05, F=3.0, sig=45., R=1.25, with_line=None):
    """[O II] doublet with a velocity gradient; flux in arbitrary units. Returns cube, noise, truth."""
    line = EmissionLine('OII', z=Z, air=True, lsf='muse')
    wave = _wave()
    yy, xx = np.mgrid[0:NY, 0:NX]
    v_true = 8. * (xx - NX / 2) - 5. * (yy - NY / 2)
    model = np.empty((len(wave), NY, NX))
    for iy in range(NY):
        for ix in range(NX):
            model[:, iy, ix] = line(wave, 0.02, F, R, v_true[iy, ix], sig)
    if with_line is not None:
        model[:, ~with_line] = 0.            # continuum-subtracted sky
    noise = np.full_like(model, noise_sigma)
    flux = model + rng.normal(0, noise_sigma, model.shape)
    return line, wave, flux, noise, dict(v=v_true, sig=sig, R=R, F=F, C=0.02)


class TestWavelengthsAndModel:

    @pytest.mark.parametrize('vac, air', [(3727.092, 3726.032), (3729.875, 3728.815),
                                          (4862.683, 4861.333), (5008.240, 5006.843),
                                          (6564.614, 6562.819), (2796.354, 2795.528)])
    def test_vac_to_air(self, vac, air):
        # published air values (NIST); conversion formulae agree to ~0.02 A
        assert vac_to_air(vac) == pytest.approx(air, abs=0.02)

    def test_presets_use_rbcodes_line_lists(self):
        oii = EmissionLine('OII')
        np.testing.assert_allclose(oii.rest, [3727.092, 3729.874])
        assert oii.kind == 'free' and oii.npar == 5
        assert EmissionLine('NII').kind == 'fixed'
        np.testing.assert_allclose(EmissionLine('OII', air=True).rest, vac_to_air(oii.rest))
        with pytest.raises(KeyError):
            EmissionLine('NotALine')

    def test_custom_line_and_bad_input(self):
        ln = EmissionLine('CIII', z=2., rest=(1906.68, 1908.73), ratio_bounds=(0.5, 2.))
        assert ln.kind == 'free'
        np.testing.assert_allclose(ln.obs, np.array([1906.68, 1908.73]) * 3.)
        with pytest.raises(ValueError):
            EmissionLine('x', rest=(1., 2., 3.), kind='free')

    @pytest.mark.parametrize('name', ['OII', 'OIII_doublet', 'Halpha'])
    def test_flux_is_integral_and_centroid(self, name):
        ln = EmissionLine(name, z=0.3, lsf='muse')
        p = (0., 7.5, 1.3, 120., 50.) if ln.kind == 'free' else (0., 7.5, 120., 50.)
        lam = np.linspace(ln.obs.min() - 100, ln.obs.max() + 100, 40001)
        prof = ln(lam, *p)
        assert np.trapezoid(prof, lam) == pytest.approx(7.5, rel=1e-6)
        # each component sits at lambda_rest (1+z)(1+v/c)
        for comp, lk in zip(ln.components(lam, *p), ln.obs * (1 + 120. / C_KMS)):
            assert np.trapezoid(comp * lam, lam) / np.trapezoid(comp, lam) == pytest.approx(lk, abs=1e-6)

    def test_lsf_width_adds_in_quadrature(self):
        lam = np.linspace(6500, 6620, 20001)
        for lsf, s_lsf in ((lsf_constant(1.2), 1.2), (lsf_resolution(3000.), 6562.8 / 3000. / 2.3548),
                           ('muse', lsf_muse(6562.8))):
            ln = EmissionLine('Ha', rest=(6562.8,), lsf=lsf)
            prof = ln(lam, 0., 1., 0., 40.)
            var = np.trapezoid(prof * (lam - 6562.8) ** 2, lam)
            assert np.sqrt(var) == pytest.approx(np.hypot(6562.8 * 40. / C_KMS, s_lsf), rel=1e-4)


class TestFitSpectrum:

    def test_noiseless_recovery_any_flux_units(self):
        ln = EmissionLine('OII', z=Z, air=True, lsf='muse')
        wave = _wave()
        for scale in (1., 1e-20, 1e5):
            true = (0.1 * scale, 4. * scale, 0.8, -75., 35.)
            r = fit_spectrum(ln, wave, ln(wave, *true), np.full(wave.size, 0.01 * scale))
            assert r['status'] == 0
            np.testing.assert_allclose(r['p'], true, rtol=1e-4, atol=1e-4 * scale)
            assert r['chi2'] < 1e-6

    def test_bad_data_status(self):
        ln = EmissionLine('OIII', z=0.1)
        wave = np.linspace(5490, 5530, 40)
        f, e = np.ones(40), np.ones(40)
        f[3] = np.nan
        assert fit_spectrum(ln, wave, f, e)['status'] == 1
        assert fit_spectrum(ln, wave[:4], f[:4], e[:4])['status'] == 1


class TestFitCube:

    def test_recovers_truth_and_errors_are_calibrated(self):
        rng = np.random.default_rng(1)
        line, wave, flux, noise, t = _oii_cube(rng)
        res = fit_cube(flux, line, noise=noise, wave=wave, verbose=False)
        ok = np.isfinite(res['v'])
        assert ok.all()
        for key, err in (('v', 'v_err'), ('sig', 'sig_err'), ('R', 'R_err'), ('F', 'F_err')):
            truth = t[key] if np.isscalar(t[key]) else t[key][ok]
            pull = (res[key][ok] - truth) / res[err][ok]
            assert abs(np.mean(pull)) < 0.25, key
            assert 0.8 < np.std(pull) < 1.2, key          # quoted errors match the scatter
        assert 0.8 < np.median(res['chi2']) < 1.2

    def test_detection_and_overrides(self):
        rng = np.random.default_rng(2)
        has = np.zeros((NY, NX), bool)
        has[3:11, 4:12] = True
        line, wave, flux, noise, _ = _oii_cube(rng, with_line=has)
        res = fit_cube(flux, line, noise=noise, wave=wave, overrides={(0, 0): dict(p0_v=0.)},
                       verbose=False)
        fitted = np.isfinite(res['status'])
        assert fitted[has].all()
        assert fitted.sum() <= has.sum() + 3            # only a chance noise detection or two
        assert fitted[0, 0]                             # override forces a fit

    def test_ifucube_input_and_save_load(self, tmp_path):
        rng = np.random.default_rng(3)
        line, wave, flux, noise, _ = _oii_cube(rng)
        h = fits.Header()
        h['CTYPE1'], h['CTYPE2'], h['CTYPE3'] = 'RA---TAN', 'DEC--TAN', 'AWAV'
        h['CRVAL1'], h['CRVAL2'], h['CRVAL3'] = 231.94, 6.87, wave[0]
        h['CRPIX1'], h['CRPIX2'], h['CRPIX3'] = 8., 7., 1.
        h['CDELT1'], h['CDELT2'], h['CDELT3'] = -0.2 / 3600, 0.2 / 3600, 1.25
        h['CUNIT3'], h['INSTRUME'] = 'Angstrom', 'MUSE'
        p = tmp_path / 'cube.fits'
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(flux, h, name='DATA'),
                      fits.ImageHDU(noise ** 2, h, name='STAT')]).writeto(p)
        cube = load_fits(p)
        res = fit_cube(cube, line, exclude=[(6680., 6685.)], verbose=False)
        direct = fit_cube(flux, line, noise=noise, wave=cube.wave, exclude=[(6680., 6685.)],
                          verbose=False)
        np.testing.assert_allclose(res['v'], direct['v'])

        out = tmp_path / 'fit.fits'
        res.save(out)
        back = LineFitResult.load(out)
        assert set(back.keys()) == set(res.keys())
        for k in res.keys():
            np.testing.assert_array_equal(back[k], res[k])
        for k in ('line', 'z', 'kind', 'vmax', 'window', 'exclude', 'sig_bounds'):
            assert back.meta[k] == pytest.approx(res.meta[k]) if k in ('z', 'vmax') else \
                back.meta[k] == res.meta[k] or np.allclose(back.meta[k], res.meta[k]), k
        np.testing.assert_allclose(back.meta['rest'], line.rest)
        assert back.wcs.wcs.ctype[0] == 'RA---TAN'
        np.testing.assert_allclose(back.wcs.wcs.crval, [231.94, 6.87])

    def test_no_variance_needs_noise(self):
        rng = np.random.default_rng(4)
        line, wave, flux, noise, _ = _oii_cube(rng)
        with pytest.raises(ValueError):
            fit_cube(flux, line, wave=wave)


class TestQualityMask:

    def _res(self):
        n = (10, 10)
        maps = dict(v=np.zeros(n), v_err=np.full(n, 5.), sig=np.full(n, 40.), sig_err=np.full(n, 3.), F=np.full(n, 10.), F_err=np.ones(n),
                    chi2=np.ones(n), status=np.zeros(n), R=np.full(n, np.nan))
        maps['F'][2, 2] = 2.                    # low S/N
        maps['v_err'][2, 6] = 100.              # large error
        maps['chi2'][6, 2] = 50.                # bad fit
        maps['v'][6, 6] = 599.                  # at the bound
        for k in ('v', 'v_err', 'sig', 'sig_err', 'F', 'F_err', 'chi2'):
            maps[k][9, :] = np.nan              # bottom row not fitted
        maps['status'][9, :] = np.nan
        maps['status'][9, 0] = 2                # fit failed
        return LineFitResult(maps, dict(vmax=600.))

    def test_codes(self):
        res = self._res()
        keep, code = quality_mask(res, min_neigh=3, min_island=5)
        assert code[2, 2] == 4 and code[2, 6] == 5 and code[6, 2] == 6 and code[6, 6] == 7
        assert code[9, 0] == 3 and code[9, 5] == 1
        assert keep.sum() == 90 - 4
        assert set(np.unique(code)) <= set(REJECT)
        assert 'low S/N' in spaxel_report(res, code, 2, 2)

    def test_isolated_island_and_forcing(self):
        res = self._res()
        res.maps['F'][:, :] = 1.                # everything low S/N ...
        res.maps['F'][0, 0] = 10.               # ... except one isolated spaxel
        res.maps['F'][3:6, 3:6] = 10.           # and a 3x3 island
        keep, code = quality_mask(res, min_neigh=3, min_island=15)
        assert code[0, 0] == 8 and (code[3:6, 3:6] == 9).all() and not keep.any()
        keep, code = quality_mask(res, force_include=[(0, 0)], force_exclude=[(4, 4)])
        assert keep[0, 0] and code[0, 0] == 10 and code[4, 4] == 11


class TestNoise:

    def test_variance_cube(self):
        var = np.full((5, 3, 3), 4.)
        var[0, 0, 0] = -1.
        s, how = noise_cube(np.zeros_like(var), var=var)
        assert np.isnan(s[0, 0, 0]) and np.all(s.ravel()[1:] == 2.)

    def test_estimate_without_variance(self):
        rng = np.random.default_rng(5)
        nw, ny, nx = 300, 20, 20
        sig_sky = 1. + 3. * (np.arange(nw) % 37 == 0)          # sky-line channels
        scale = np.linspace(0.8, 1.6, nx)[None, :] * np.ones((ny, 1))
        true = sig_sky[:, None, None] * scale[None]
        flux = rng.normal(0, 1, (nw, ny, nx)) * true
        flux[130:170, 8:12, 8:12] += 5.                         # a source
        line_free = np.ones(nw, bool)
        line_free[120:180] = False
        sky = np.zeros((ny, nx), bool)
        sky[:5, :] = True
        est, _ = noise_cube(flux, line_free=line_free, sky_mask=sky)
        ratio = est / true
        assert abs(np.median(ratio) - 1) < 0.05
        assert robust_std(ratio) < 0.15
        with pytest.raises(ValueError):
            noise_cube(flux)                                    # no sky region given
        with pytest.raises(ValueError):
            noise_cube(flux, line_free=line_free)
        with pytest.raises(ValueError):
            noise_cube(flux, line_free=line_free, sky_mask=sky[:4])

    def test_sky_mask_from_region(self, tmp_path):
        # ds9 image coords are 1-based: box centre (8.5, 3) size 16x5 covers rows 0-4, all 16 columns
        text = 'image\nbox(8.5,3,16,5,0)\ncircle(14,12,1.2)'
        mask = sky_mask_from_region(text, None, (14, 16))
        expected = np.zeros((14, 16), bool)
        expected[0:5, :] = True
        expected[11, 13] = expected[10, 13] = expected[12, 13] = True
        expected[11, 12] = expected[11, 14] = True
        np.testing.assert_array_equal(mask, expected)
        f = tmp_path / 'sky.reg'
        f.write_text(text)
        np.testing.assert_array_equal(sky_mask_from_region(f, None, (14, 16)), mask)


def test_plot_spaxel_runs():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(6)
    line, wave, flux, noise, _ = _oii_cube(rng)
    res = fit_cube(flux, line, noise=noise, wave=wave, exclude=[(6680., 6685.)], verbose=False)
    alt = fit_spectrum(line, wave, flux[:, 3, 3], noise[:, 3, 3])
    fig, axes = plt.subplots(2, 1)
    plot_spaxel(axes, line, wave, flux, noise, 3, 3, res, alt=alt)
    assert len(axes[0].lines) >= 5 and len(axes[1].lines) >= 2
    plt.close(fig)
