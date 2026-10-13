"""Known-answer tests for rbcodes.ifu.lines, rbcodes.ifu.measure and rbcodes.ifu.noise."""
import numpy as np
import pytest
from astropy.io import fits
from scipy.ndimage import gaussian_filter

from rbcodes.ifu import load_fits
from rbcodes.ifu.lines import (C_KMS, LINES, LineModel, lsf_constant, lsf_muse, lsf_resolution,
                               window_masks)
from rbcodes.ifu.measure import (REJECT, FitResult, GaussianFitter, MCErrors, MomentMeasurer,
                                 QualityMask, SpaxelDecisions, plot_spaxel, spaxel_report)
from rbcodes.ifu.noise import NoiseModel, robust_std, sky_mask_from_region

Z = 0.785273
NY, NX = 12, 14
OII_WIN = dict(flux=(6640., 6680.), fit=(6640., 6680.))


def _wave():
    return np.arange(6600., 6720., 1.25)


def _oii_cube(seed, noise_sigma=0.05, F=3.0, sig=45., R=1.25, has_line=None):
    lm = LineModel('OII', z=Z, windows=OII_WIN, lsf=lsf_muse)
    wave = _wave()
    yy, xx = np.mgrid[0:NY, 0:NX]
    v_true = 8. * (xx - NX / 2) - 5. * (yy - NY / 2)
    model = np.zeros((len(wave), NY, NX))
    dom = (wave > 6640) & (wave < 6680)
    for iy in range(NY):
        for ix in range(NX):
            if has_line is None or has_line[iy, ix]:
                model[dom, iy, ix] = lm(wave[dom], 0., F, R, v_true[iy, ix], sig)
    sigma = np.full_like(model, noise_sigma)
    flux = model + np.random.default_rng(seed).normal(0, noise_sigma, model.shape)
    return lm, wave, flux, sigma, dict(v=v_true, sig=sig, R=R, F=F)


# ---------------------------------------------------------------------------
# lines
# ---------------------------------------------------------------------------

class TestLineModel:

    def test_air_vacuum_is_an_option(self):
        vac = LineModel('OII', Z, OII_WIN, air=False)
        air = LineModel('OII', Z, OII_WIN, air=True)
        np.testing.assert_allclose(vac.groups[0]['rest'], LINES['OII']['vac'])
        np.testing.assert_allclose(air.groups[0]['rest'], LINES['OII']['air'])
        assert vac.param_names == ('C', 'F', 'R', 'v', 'sig')

    def test_windows_required_and_no_overlap(self):
        with pytest.raises(ValueError):
            LineModel('OII', Z, dict(flux=(6640, 6680)))
        with pytest.raises(ValueError):
            LineModel(['OII', 'Hbeta'], Z, {'OII': dict(fit=(6640, 6680))})
        with pytest.raises(ValueError):
            LineModel(['OII', 'Hbeta'], Z, {'OII': dict(fit=(6640, 6680)), 'Hbeta': dict(fit=(6670, 6700))})

    def test_ratio_bounds_and_ratios_are_options(self):
        assert LineModel('MgII', 0.5, dict(fit=(4180, 4220))).groups[0]['ratio_bounds'] == LINES['MgII']['ratio_bounds']
        assert LineModel('MgII', 0.5, dict(fit=(4180, 4220)), ratio_bounds=(0.5, 1.0)).groups[0]['ratio_bounds'] == (0.5, 1.0)
        nii = LineModel('NII', 0.1, dict(fit=(7190, 7260)), ratios=(1.0, 0.5))
        np.testing.assert_allclose(nii.groups[0]['w'], [2 / 3, 1 / 3])

    def test_custom_catalogue(self):
        cat = {**LINES, 'CIII': dict(label='C III]', air=(1906.68, 1908.73), vac=(1906.68, 1908.73),
                                     ratio='free', ratio_bounds=(0.5, 2.))}
        lm = LineModel('CIII', 2.0, dict(fit=(5700, 5750)), catalogue=cat)
        np.testing.assert_allclose(lm.obs, np.array([1906.68, 1908.73]) * 3)

    @pytest.mark.parametrize('name, p', [('OII', (0., 7.5, 1.3, 120., 50.)), ('OIII_doublet', (0., 7.5, 120., 50.)),
                                         ('Halpha', (0., 7.5, 120., 50.))])
    def test_flux_and_centroid(self, name, p):
        obs = np.array(LINES[name]['vac']) * 1.3
        lm = LineModel(name, 0.3, dict(fit=(obs.min() - 100, obs.max() + 100)), lsf=lsf_muse)
        lam = np.linspace(obs.min() - 90, obs.max() + 90, 40001)
        assert np.trapezoid(lm(lam, *p), lam) == pytest.approx(7.5, rel=1e-6)
        for comp, lk in zip(lm.components(lam, *p), lm.obs * (1 + 120. / C_KMS)):
            assert np.trapezoid(comp * lam, lam) / np.trapezoid(comp, lam) == pytest.approx(lk, abs=1e-6)

    @pytest.mark.parametrize('lsf, s_lsf', [(lsf_constant(1.2), 1.2), (lsf_resolution(3000.), 6564.614 / 3000. / 2.3548),
                                            (lsf_muse, lsf_muse(6564.614)), (None, 0.)])
    def test_lsf_adds_in_quadrature(self, lsf, s_lsf):
        lm = LineModel('Halpha', 0., dict(fit=(6400, 6700)), lsf=lsf)
        lam = np.linspace(6500, 6630, 20001)
        prof = lm(lam, 0., 1., 0., 40.)
        width = np.sqrt(np.trapezoid(prof * (lam - 6564.614) ** 2, lam))
        assert width == pytest.approx(np.hypot(6564.614 * 40. / C_KMS, s_lsf), rel=1e-4)

    def test_window_masks(self):
        w = _wave()
        fl, ft, lf = window_masks(w, dict(flux=(6650, 6660), fit=(6640, 6680), exclude=[(6670, 6672)]))
        assert fl.sum() == 7 and not ft[(w > 6670) & (w < 6672)].any()
        assert not lf[(w > 6635) & (w < 6685)].any() and lf[w < 6634].all()


# ---------------------------------------------------------------------------
# Gaussian fitting
# ---------------------------------------------------------------------------

class TestGaussianFitter:

    def test_noiseless_recovery_in_any_units(self):
        lm = LineModel('OII', Z, OII_WIN, lsf=lsf_muse)
        wave = _wave()
        dom = lm.fit_mask(wave)
        for scale in (1., 1e-20, 1e5):
            true = (0.1 * scale, 4. * scale, 0.8, -75., 35.)
            r = GaussianFitter(lm).fit_spectrum(wave[dom], lm(wave[dom], *true), np.full(dom.sum(), 0.05 * scale))
            assert r['status'] == 0 and r['chi2'] < 1e-6
            np.testing.assert_allclose(r['p'], true, rtol=1e-4, atol=1e-4 * scale)

    def test_bad_data(self):
        lm = LineModel('OIII', 0.1, dict(fit=(5490, 5530)))
        wave = np.linspace(5490, 5530, 40)
        f = np.ones(40)
        f[3] = np.nan
        assert GaussianFitter(lm).fit_spectrum(wave, f, np.ones(40))['status'] == 1
        assert GaussianFitter(lm).fit_spectrum(wave[:4], f[:4], np.ones(4))['status'] == 1

    def test_tied_lines(self):
        z = 0.3
        oii, hb = np.array(LINES['OII']['vac']) * 1.3, LINES['Hbeta']['vac'][0] * 1.3
        lm = LineModel(['OII', 'Hbeta'], z, {'OII': dict(fit=(oii[0] - 20, oii[1] + 20)),
                                            'Hbeta': dict(fit=(hb - 20, hb + 20))}, lsf=lsf_resolution(3000))
        assert lm.param_names == ('C_OII', 'F_OII', 'R_OII', 'C_Hbeta', 'F_Hbeta', 'v', 'sig')
        wave = np.arange(4800., 6400., 0.5)
        true = (0., 5., 1.4, 0.05, 2., 33., 55.)      # continua within the default +-5 sigma bound
        dom = lm.fit_mask(wave)
        r = GaussianFitter(lm).fit_spectrum(wave[dom], lm(wave[dom], *true), np.full(dom.sum(), 0.02))
        assert r['status'] == 0
        np.testing.assert_allclose(r['p'], true, rtol=1e-4, atol=1e-4)
        assert r['F'] == pytest.approx(5.) and r['F_Hbeta'] == pytest.approx(2.)

    def test_cube_recovery_and_error_calibration(self):
        lm, wave, flux, sigma, t = _oii_cube(seed=1)
        fit = GaussianFitter(lm).fit_cube(flux, sigma, wave=wave, verbose=False)
        ok = np.isfinite(fit['v'])
        assert ok.all()
        for key in ('v', 'sig', 'R', 'F'):
            truth = t[key][ok] if key == 'v' else t[key]
            pull = (fit[key][ok] - truth) / fit[key + '_err'][ok]
            assert abs(pull.mean()) < 0.3, key
            assert 0.8 < pull.std() < 1.2, key
        assert 0.8 < np.median(fit['chi2']) < 1.2

    def test_detection_overrides_and_ifucube_input(self, tmp_path):
        has = np.zeros((NY, NX), bool)
        has[3:10, 4:11] = True
        lm, wave, flux, sigma, _ = _oii_cube(seed=2, has_line=has)
        h = fits.Header()
        for k, v in dict(CTYPE1='RA---TAN', CTYPE2='DEC--TAN', CTYPE3='WAVE', CUNIT3='Angstrom', CRVAL1=231.9,
                         CRVAL2=6.87, CRVAL3=wave[0], CRPIX1=1., CRPIX2=1., CRPIX3=1., CDELT1=-0.2 / 3600,
                         CDELT2=0.2 / 3600, CDELT3=1.25, INSTRUME='MUSE').items():
            h[k] = v
        p = tmp_path / 'c.fits'
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(flux, h, name='DATA'),
                      fits.ImageHDU(sigma ** 2, h, name='STAT')]).writeto(p)
        cube = load_fits(p)
        nm = NoiseModel.from_cube(cube)
        fit = GaussianFitter(lm).fit_cube(cube, nm, overrides={(0, 0): dict(p0_v=0., vmax=300.)}, verbose=False)
        attempted = np.isfinite(fit['status'])
        assert attempted[has].all() and attempted[0, 0]
        assert attempted.sum() <= has.sum() + 3
        assert fit.detect[has].all()

        out = tmp_path / 'fit.fits'
        fit.save(out)
        back = FitResult.load(out)
        assert set(back.keys()) == set(fit.keys())
        for k in fit.keys():
            np.testing.assert_array_equal(back[k], fit[k])
        np.testing.assert_array_equal(back.detect, fit.detect)
        assert back.meta['line'] == 'OII' and back.meta['air'] is False and back.meta['vmax'] == 600.
        assert back.meta['overrides'] == {'0,0': {'p0_v': 0.0, 'vmax': 300.0}}
        np.testing.assert_allclose(back.wcs.wcs.crval, [231.9, 6.87])


# ---------------------------------------------------------------------------
# Quality mask and decisions
# ---------------------------------------------------------------------------

def _maps():
    n = (10, 10)
    maps = dict(v=np.zeros(n), v_err=np.full(n, 5.), sig=np.full(n, 40.), sig_err=np.full(n, 3.),
                F=np.full(n, 10.), F_err=np.ones(n), chi2=np.ones(n), status=np.zeros(n),
                R=np.full(n, np.nan), R_err=np.full(n, np.nan), C=np.zeros(n))
    maps['F'][2, 2] = 2.                    # low S/N
    maps['v_err'][2, 6] = 100.              # large error
    maps['chi2'][6, 2] = 50.                # bad fit
    maps['v'][6, 6] = 599.                  # at the bound
    for k in ('v', 'v_err', 'sig', 'sig_err', 'F', 'F_err', 'chi2'):
        maps[k][9, :] = np.nan
    maps['status'][9, :] = np.nan           # not attempted ...
    maps['status'][9, 0] = 2                # ... except one failed fit
    detect = np.ones(n, bool)
    detect[9, 1:] = False
    return FitResult(maps, dict(line='OII', vmax=600.), detect)


class TestQualityMask:

    def test_codes_follow_reject_table(self):
        fit = _maps()
        q = QualityMask(min_neigh=3, min_island=5).apply(fit)
        assert (q.code[2, 2], q.code[2, 6], q.code[6, 2], q.code[6, 6]) == (3, 4, 5, 6)
        assert q.code[9, 0] == 2 and q.code[9, 5] == 1
        assert q.mask.sum() == 90 - 4
        assert set(np.unique(q.code)) <= set(REJECT) and max(REJECT) == 10
        assert 'low S/N' in spaxel_report(fit, q.code, 2, 2)

    def test_isolated_island_and_forcing(self):
        fit = _maps()
        fit['F'][:, :] = 1.
        fit['F'][0, 0] = 10.
        fit['F'][3:6, 3:6] = 10.
        q = QualityMask(min_island=15).apply(fit)
        assert q.code[0, 0] == 7 and (q.code[3:6, 3:6] == 8).all() and not q.mask.any()
        q = QualityMask(force_include=[(0, 0)], force_exclude=[(4, 4)]).apply(fit)
        assert q.mask[0, 0] and q.code[0, 0] == 9 and q.code[4, 4] == 10

    def test_decisions_update_maps_without_refitting(self, tmp_path):
        lm, wave, flux, sigma, _ = _oii_cube(seed=3)
        fitter = GaussianFitter(lm)
        fit = fitter.fit_cube(flux, sigma, wave=wave, verbose=False)
        path = tmp_path / 'decisions.json'
        dec = SpaxelDecisions(path)
        r = fitter.refit(flux, sigma, 4, 5, wave=wave, fit=(6645., 6675.), sig_bounds=(5., 200.))
        dec.accept(fit, 4, 5, r)
        assert fit['v'][4, 5] == r['v']
        dec.reject('OII', 6, 6)
        dec.accept(fit, 7, 7)
        q = QualityMask().apply(fit, decisions=dec)
        assert q.code[4, 5] == 9 and q.code[7, 7] == 9 and q.code[6, 6] == 10

        again = SpaxelDecisions(path)                              # persisted
        assert again.excludes('OII') == [(6, 6)] and sorted(again.includes('OII')) == [(4, 5), (7, 7)]
        refit_all = fitter.fit_cube(flux, sigma, wave=wave, decisions=again, verbose=False)
        assert refit_all['v'][4, 5] == pytest.approx(r['v'])       # stored override reproduces the refit
        again.forget('OII', 6, 6)
        assert SpaxelDecisions(path).excludes('OII') == []


# ---------------------------------------------------------------------------
# Moments and Monte Carlo errors
# ---------------------------------------------------------------------------

class TestMoments:

    def test_single_line_moments_and_errors(self):
        rng = np.random.default_rng(4)
        lm = LineModel('Halpha', 0.1, dict(fit=(7180, 7260)), lsf=lsf_constant(1.0))
        wave = np.arange(7150., 7290., 0.8)
        n, v_true, s_int = 400, 37., 60.
        clean = lm(wave, 0., 20., v_true, s_int)
        flux = clean[:, None, None] + rng.normal(0, 0.3, (len(wave), 1, n))
        sigma = np.full_like(flux, 0.3)
        mm = MomentMeasurer(lm, vwin=(-250, 250))          # ~ +-3.4 x the observed width
        res = mm.measure_cube(flux, sigma, wave=wave, detect=np.ones((1, n), bool))
        v, ve = res['v'][0], res['v_err'][0]
        assert abs(np.mean(v) - v_true) < 3 * np.std(v) / np.sqrt(n)
        assert np.std(v) / np.median(ve) == pytest.approx(1., abs=0.12)
        s_lsf = 1.0 / (6564.614 * 1.1) * C_KMS
        assert np.median(res['sig'][0]) == pytest.approx(np.hypot(s_int, s_lsf), rel=0.03)
        assert np.median(res['F'][0]) == pytest.approx(20., rel=0.02)

    def test_stacking_tied_lines_reduces_error(self):
        # inverse-variance stacking weights by noise, so it helps for lines of similar strength
        rng = np.random.default_rng(5)
        lm = LineModel('OIII_doublet', 0.2, dict(fit=(5920, 6030)), lsf=lsf_constant(1.0), ratios=(1., 1.))
        wave = np.arange(5900., 6060., 0.8)
        clean = lm(wave, 0., 20., -40., 50.)
        n = 200
        flux = clean[:, None, None] + rng.normal(0, 0.3, (len(wave), 1, n))
        sigma = np.full_like(flux, 0.3)
        both = MomentMeasurer(lm, vwin=(-400, 400)).measure_cube(flux, sigma, wave=wave, detect=np.ones((1, n), bool))
        one = MomentMeasurer(lm, vwin=(-400, 400), components=[LINES['OIII']['vac'][0]]).measure_cube(
            flux, sigma, wave=wave, detect=np.ones((1, n), bool))
        assert abs(np.mean(both['v']) + 40.) < 1.5
        assert np.median(both['v_err']) == pytest.approx(np.median(one['v_err']) / np.sqrt(2), rel=0.1)


class TestMCErrors:

    def test_white_noise_errors_match_and_correlated_exceed(self):
        lm, wave, flux, sigma, _ = _oii_cube(seed=6)
        fitter = GaussianFitter(lm)
        fit = fitter.fit_cube(flux, sigma, wave=wave, verbose=False)
        white = MCErrors(fitter, n_real=60, seed=1).run(fit, flux, sigma, wave=wave, n_spaxels=15)
        assert white['median_ratio'] == pytest.approx(1., abs=0.2)
        corr = MCErrors(fitter, n_real=60, seed=1, corr_spectral=2.).run(fit, flux, sigma, wave=wave,
                                                                         spaxels=white['spaxels'])
        assert corr['median_ratio'] > 1.4
        assert corr['spaxels'] == white['spaxels']


def test_plot_spaxel_runs():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lm, wave, flux, sigma, _ = _oii_cube(seed=7)
    fitter = GaussianFitter(lm)
    r = fitter.refit(flux, sigma, 3, 3, wave=wave)
    alt = fitter.refit(flux, sigma, 3, 3, wave=wave, fit=(6645., 6675.))
    fig, axes = plt.subplots(2, 1)
    plot_spaxel(axes, lm, wave, flux, sigma, 3, 3, r, alt)
    assert len(axes[0].lines) >= 6 and len(axes[1].lines) >= 2
    plt.close(fig)


# ---------------------------------------------------------------------------
# Noise
# ---------------------------------------------------------------------------

class TestNoiseModel:

    @staticmethod
    def _cube(seed, smooth=(0, 0, 0)):
        rng = np.random.default_rng(seed)
        nw, ny, nx = 300, 20, 20
        sig_sky = 1. + 3. * (np.arange(nw) % 37 == 0)
        scale = np.linspace(0.8, 1.6, nx)[None, :] * np.ones((ny, 1))
        true = sig_sky[:, None, None] * scale[None]
        z = rng.normal(0, 1, (nw, ny, nx))
        if any(smooth):
            z = gaussian_filter(z, smooth, mode='wrap')
            z /= z.std()
        flux = z * true
        flux[130:170, 8:12, 8:12] += 5.
        lf = np.ones(nw, bool)
        lf[120:180] = False
        return flux, true, lf

    def test_variance(self):
        var = np.full((8, 4, 4), 4.)
        var[0, 0, 0] = -1.
        nm = NoiseModel.from_cube(np.zeros_like(var), var=var)
        assert np.isnan(nm.sigma[0, 0, 0]) and np.all(nm.sigma.ravel()[1:] == 2.)

    def test_hybrid_recovers_truth(self):
        flux, true, lf = self._cube(8)
        nm = NoiseModel.from_cube(flux, linefree=lf, sky_box=(0, 5, 0, 20))
        ratio = nm.sigma / true
        assert abs(np.median(ratio) - 1) < 0.05 and robust_std(ratio) < 0.15
        flat = NoiseModel.from_cube(flux, linefree=lf, sky_box=(0, 5, 0, 20), spaxel_scale=False)
        assert np.all(flat.scale == 1) and flat.sigma.shape == flux.shape

    def test_sky_region_is_required_and_unique(self):
        flux, _, lf = self._cube(9)
        with pytest.raises(ValueError):
            NoiseModel.from_cube(flux, linefree=lf)
        with pytest.raises(ValueError):
            NoiseModel.from_cube(flux, linefree=lf, sky_box=(0, 5, 0, 20), sky_mask=np.ones((20, 20), bool))
        reg = NoiseModel.from_cube(flux, linefree=lf, sky_region='image\nbox(10.5,3,20,5,0)')
        box = NoiseModel.from_cube(flux, linefree=lf, sky_box=(0, 5, 0, 20))
        np.testing.assert_array_equal(reg.sky_mask, box.sky_mask)
        np.testing.assert_allclose(reg.sigma, box.sigma)

    def test_linefree_as_ranges_or_linemodel(self):
        lm, wave, flux, sigma, _ = _oii_cube(seed=10)
        a = NoiseModel.from_cube(flux, wave=wave, linefree=lm, sky_box=(0, 3, 0, NX))
        b = NoiseModel.from_cube(flux, wave=wave, linefree=[(6590, 6634.9), (6685.1, 6730)], sky_box=(0, 3, 0, NX))
        np.testing.assert_allclose(a.sigma, b.sigma)
        assert a.correlation(flux, lm) == a.correlation(flux, [(6590, 6634.9), (6685.1, 6730)])

    def test_correlation(self):
        flux, _, lf = self._cube(11)
        nm = NoiseModel.from_cube(flux, linefree=lf, sky_box=(0, 8, 0, 20))
        c = nm.correlation(flux, lf)
        assert all(abs(c[k]) < 0.05 for k in ('rho_x', 'rho_y', 'rho_lambda'))
        flux, _, lf = self._cube(11, smooth=(1.5, 0, 0))
        c = NoiseModel.from_cube(flux, linefree=lf, sky_box=(0, 8, 0, 20)).correlation(flux, lf)
        assert c['rho_lambda'] > 0.6 and abs(c['rho_x']) < 0.05
        flux, _, lf = self._cube(11, smooth=(0, 1.0, 1.0))
        c = NoiseModel.from_cube(flux, linefree=lf, sky_box=(0, 8, 0, 20)).correlation(flux, lf)
        assert c['rho_x'] > 0.4 and c['rho_y'] > 0.4 and abs(c['rho_lambda']) < 0.05


def test_sky_mask_from_region_file(tmp_path):
    f = tmp_path / 'sky.reg'
    f.write_text('image\nbox(8.5,3,16,5,0)')
    m = sky_mask_from_region(f, None, (14, 16))
    assert m[:5].all() and not m[5:].any()


# ---------------------------------------------------------------------------
# linear continuum and C III]
# ---------------------------------------------------------------------------

class TestLinearContinuum:

    @staticmethod
    def _spec(scale=1.):
        lm = LineModel('CIII', 1.701, dict(fit=(5130., 5180.)), lsf=lsf_muse, continuum='linear')
        wave = np.arange(5130., 5180., 1.25)
        true = (2.0 * scale, 0.03 * scale, 6. * scale, 0.9, 45., 40.)       # C, S, F, R, v, sig
        return lm, wave, true

    def test_ciii_entry_and_parameters(self):
        lm, _, _ = self._spec()
        np.testing.assert_allclose(lm.groups[0]['rest'], [1906.683, 1908.734])
        assert lm.param_names == ('C', 'S', 'F', 'R', 'v', 'sig')
        assert LineModel('CIII', 1.701, dict(fit=(5130., 5180.)), air=True).groups[0]['rest'][0] == 1906.683
        with pytest.raises(ValueError):
            LineModel('CIII', 1.7, dict(fit=(5130., 5180.)), continuum='quadratic')

    def test_noiseless_recovery_any_units(self):
        for scale in (1., 1e-20, 1e4):
            lm, wave, true = self._spec(scale)
            r = GaussianFitter(lm).fit_spectrum(wave, lm(wave, *true), np.full(wave.size, 0.05 * scale))
            assert r['status'] == 0 and r['chi2'] < 1e-6
            np.testing.assert_allclose(r['p'], true, rtol=1e-4, atol=1e-6 * scale)
            assert r['S'] == pytest.approx(true[1]) and r['C'] == pytest.approx(true[0])

    def test_constant_model_is_biased_by_a_sloped_continuum(self):
        lm, wave, true = self._spec()
        f = lm(wave, *true)
        e = np.full(wave.size, 0.05)
        const = LineModel('CIII', 1.701, dict(fit=(5130., 5180.)), lsf=lsf_muse)
        bad = GaussianFitter(const).fit_spectrum(wave, f, e)
        assert bad['status'] != 0 or bad['chi2'] > 100                   # continuum outside +-5 sigma, slope ignored
        good = GaussianFitter(lm).fit_spectrum(wave, f, e)
        assert good['v'] == pytest.approx(45., abs=1e-3)

    def test_cube_with_continuum(self):
        lm, wave, true = self._spec()
        rng = np.random.default_rng(12)
        ny, nx = 6, 7
        v = rng.uniform(-80, 80, (ny, nx))
        cube = np.empty((wave.size, ny, nx))
        for iy in range(ny):
            for ix in range(nx):
                cube[:, iy, ix] = lm(wave, 2.0 + 0.2 * ix, 0.02 * iy, 6., 0.9, v[iy, ix], 40.)
        sig = np.full_like(cube, 0.05)
        cube += rng.normal(0, 0.05, cube.shape)
        fit = GaussianFitter(lm).fit_cube(cube, sig, wave=wave, detect=np.ones((ny, nx), bool), verbose=False)
        pull = (fit['v'] - v) / fit['v_err']
        assert np.all(np.isfinite(pull)) and np.abs(pull).max() < 4
        assert 'S' in fit.keys() and fit['S'][3, 0] == pytest.approx(0.06, abs=0.01)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1)
        plot_spaxel(axes, lm, wave, cube, sig, 2, 3, GaussianFitter(lm).refit(cube, sig, 2, 3, wave=wave))
        plt.close(fig)
