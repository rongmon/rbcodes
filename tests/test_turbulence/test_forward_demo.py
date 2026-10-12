"""Known-answer tests for rbcodes.turbulence.forward and rbcodes.turbulence.demo."""
import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from rbcodes.ifu.lines import LineModel, lsf_muse
from rbcodes.ifu.measure import GaussianFitter, QualityMask
from rbcodes.ifu.noise import NoiseModel
from rbcodes.ifu.vmap import VelocityMap
from rbcodes.turbulence.demo import DEMO_TRUTH, demo_cube, demo_velocity_field
from rbcodes.turbulence.forward import Observation, TurbulenceField, grf, psf_weight
from rbcodes.turbulence.vsf import StructureFunction


def _vm(f):
    return VelocityMap(f, np.full(f.shape, 1e-3), np.ones(f.shape, bool), kpc_per_pix=1.)


class TestGRF:

    def test_reproducible_and_shapes(self):
        a = grf(0.67, 10., (20, 30), 1.0, np.random.default_rng(1))
        b = grf(0.67, 10., (20, 30), 1.0, np.random.default_rng(1))
        np.testing.assert_array_equal(a, b)
        assert a.shape == (20, 30) and grf(1., 10., (6, 8, 9), 1., np.random.default_rng(2)).shape == (6, 8, 9)

    @pytest.mark.parametrize('alpha', [0.67, 1.2])
    def test_structure_function_slope(self, alpha):
        S = []
        for seed in range(8):
            sf = StructureFunction(_vm(grf(alpha, None, (48, 48), 1.0, np.random.default_rng(seed))), n_bins=14)
            S.append(sf.S)
        S = np.nanmean(S, 0)
        k = sf.ok & (sf.rc > 2.5) & (sf.rc < 12)
        slope = np.polyfit(np.log(sf.rc[k]), np.log(S[k]), 1)[0]
        assert slope == pytest.approx(alpha, abs=0.12)

    def test_outer_scale_flattens(self):
        S = {}
        for L in (None, 3.):
            acc = []
            for seed in range(6):
                sf = StructureFunction(_vm(grf(1.0, L, (40, 40), 1.0, np.random.default_rng(seed))), n_bins=12)
                acc.append(sf.S / np.nanmean(sf.S[sf.ok][-3:]))
            S[L] = np.nanmean(acc, 0)
        k = sf.ok & (sf.rc > 6) & (sf.rc < 20)
        slope = {L: np.polyfit(np.log(sf.rc[k]), np.log(S[L][k]), 1)[0] for L in S}
        assert slope[3.] < slope[None] - 0.3                          # saturates beyond L


class TestPSFWeight:

    def test_constant_and_uniform_flux(self):
        rng = np.random.default_rng(3)
        x = rng.normal(size=(30, 30))
        flux = np.ones((30, 30))
        np.testing.assert_allclose(psf_weight(np.full((30, 30), 7.), rng.random((30, 30)) + 0.1, 1.5), 7.)
        np.testing.assert_allclose(psf_weight(x, flux, 1.5), gaussian_filter(x, 1.5) / gaussian_filter(flux, 1.5))
        assert np.all(psf_weight(x, np.zeros((30, 30)), 1.5) == 0)


def _obs(**kw):
    yy, xx = np.mgrid[0:30, 0:36]
    flux = np.exp(-np.hypot(xx - 18, (yy - 15) * 1.5) / 8.)
    mask = flux > 0.2
    return Observation(mask, flux, 1.5, 1.2, np.full(mask.shape, 4.), **kw)


class TestObservation:

    def test_2d_normalisation_and_noise(self):
        obs = _obs()
        f = TurbulenceField(12., 0.67, 15.)
        clean = obs(f, np.random.default_rng(4), noise=False)
        v = clean.v[obs.mask]
        assert v.std() == pytest.approx(12., rel=1e-10) and abs(v.mean()) < 1e-10
        assert np.all(np.isnan(clean.v[~obs.mask])) and clean.mask.sum() == obs.mask.sum()
        noisy = obs(f, np.random.default_rng(4))
        assert np.std(noisy.v[obs.mask] - v) == pytest.approx(4., rel=0.1)
        assert isinstance(noisy, VelocityMap) and noisy.meta['mock']['alpha'] == 0.67
        assert noisy.kpc_per_pix == 1.5 and noisy.psf_fwhm_kpc == pytest.approx(1.2 * 2.3548 * 1.5)

    def test_bulk_and_theta(self):
        bulk = np.mgrid[0:30, 0:36][1] * 3.
        a = _obs(bulk_model=bulk).simulate((10., 1.0, 1.0, 0.), np.random.default_rng(5), noise=False)
        b = _obs().simulate((10., 1.0, 1.0, 0.), np.random.default_rng(5), noise=False)
        np.testing.assert_allclose((a.v - b.v)[a.mask], bulk[a.mask])

    def test_3d_slab(self):
        obs = _obs()
        f = TurbulenceField(20., 0.67, 30., mode='3D', depth_kpc=15.)
        m = obs(f, np.random.default_rng(6), noise=False)
        assert abs(np.nanmean(m.v[m.mask])) < 1e-10
        assert np.all(m.width[m.mask] >= 2.3)
        # line-of-sight averaging and PSF reduce the centroid rms below the 3-D rms
        assert np.nanstd(m.v[m.mask]) < 20. and np.nanmedian(m.width[m.mask]) < 20.

    def test_intermittency_raises_increment_kurtosis(self):
        kurt = {}
        for gamma in (0., 1.0):
            k = []
            for seed in range(6):
                g = TurbulenceField(1., 0.67, 20., gamma=gamma).draw((48, 48), 1., np.random.default_rng(seed))['g']
                d = (g[:, 2:] - g[:, :-2]).ravel()
                k.append(np.mean((d - d.mean()) ** 4) / np.var(d) ** 2)
            kurt[gamma] = np.mean(k)
        assert kurt[0.] == pytest.approx(3., abs=0.4) and kurt[1.0] > kurt[0.] + 1.

    def test_from_vmap(self):
        vm = _obs().simulate((10., 1.0, 1.0, 0.), np.random.default_rng(7))
        obs = Observation.from_vmap(vm)
        assert obs.kpp == 1.5 and obs.psf_sig == pytest.approx(1.2)
        np.testing.assert_array_equal(obs.mask, vm.mask)
        with pytest.raises(ValueError):
            Observation.from_vmap(VelocityMap(vm.v, vm.v_err, vm.mask))


class TestDemo:

    def test_velocity_field_truth(self):
        flux, vel, vt, gal = demo_velocity_field()
        m = flux > 3
        assert vt[m].std() == pytest.approx(DEMO_TRUTH['sig_turb'], rel=1e-10) and abs(vt[m].mean()) < 1e-10
        assert gal.sum() > 5 and flux.shape == (40, 50)

    def test_cube_end_to_end(self):
        cube, truth = demo_cube('OII', 0.785273, shape=(30, 36))
        lo, hi = truth['window']
        lm = LineModel('OII', 0.785273, dict(fit=(lo, hi)), air=True, lsf=lsf_muse)
        nm = NoiseModel.from_cube(cube)
        fit = GaussianFitter(lm).fit_cube(cube, nm, detect_snr=5., verbose=False)
        q = QualityMask(snr_min=10, verr_max=15, chi2_max=3).apply(fit)
        vm = VelocityMap.from_fit(fit, q, kpc_per_pix=1.49)
        assert vm.mask.sum() >= 150
        pull = (vm.v - truth['vel'])[vm.mask] / vm.v_err[vm.mask]
        assert abs(pull.mean()) < 0.3 and 0.8 < pull.std() < 1.25
        assert np.median(fit['sig'][q.mask]) == pytest.approx(truth['sig_int'], rel=0.05)
        assert np.median(fit['R'][q.mask]) == pytest.approx(1.35, rel=0.05)
        StructureFunction(vm, n_bins=12)                              # the chain runs through to S2
