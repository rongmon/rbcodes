"""Tests for rbcodes.turbulence.summaries."""
import numpy as np
import pytest

from rbcodes.ifu.vmap import VelocityMap
from rbcodes.turbulence.forward import Observation, TurbulenceField
from rbcodes.turbulence.summaries import Kurtosis, LineWidth, S2, Scattering, Summaries
from rbcodes.turbulence.vsf import StructureFunction


def _obs():
    yy, xx = np.mgrid[0:30, 0:36]
    flux = np.exp(-np.hypot(xx - 18, (yy - 15) * 1.3) / 9.)
    return Observation(flux > 0.15, flux, 1.5, 1.0, np.full(flux.shape, 3.))


def _data(seed=0, **kw):
    return _obs().simulate((15., 0.8, 1.2, 0.), np.random.default_rng(seed), **kw)


def test_scattering_scaling_and_labels():
    vm = _data()
    W = Scattering(vm.mask, J=4, L=4)
    assert len(W.labels) == 4 + 6
    u = np.where(vm.mask, vm.v, np.nan)
    a, b = W(u), W(3.7 * u)
    np.testing.assert_allclose(b[:4], a[:4] + np.log(3.7), atol=1e-10)      # log S1 shifts
    np.testing.assert_allclose(b[4:], a[4:], rtol=1e-10)                    # s2 = S2/S1 unchanged
    with pytest.raises(ValueError):
        Scattering(np.zeros((5, 5), bool))


def test_scattering_separates_intermittent_fields():
    obs = _obs()
    W = Scattering(obs.mask)
    G, I = [], []
    for seed in range(25):
        rng = np.random.default_rng(seed)
        G.append(W(obs.simulate((15., 0.67, 1.0, 0.), rng, noise=False).v)[4:])
        I.append(W(obs.simulate((15., 0.67, 1.0, 1.0), rng, noise=False, L_env=8.).v)[4:])
    G, I = np.array(G), np.array(I)
    C = 0.5 * (np.cov(G.T) + np.cov(I.T)) + 1e-9 * np.eye(G.shape[1])
    dm = I.mean(0) - G.mean(0)
    assert np.sqrt(dm @ np.linalg.solve(C, dm)) > 2.


def test_s2_block_matches_structure_function():
    data, mock = _data(1), _data(2)
    sf = StructureFunction(data, n_bins=10)
    summ = Summaries(sf, blocks=[S2()], removal=None)
    np.testing.assert_allclose(summ(data), sf.S[sf.ok])
    direct = StructureFunction(mock, n_bins=10)                # same mask, errors and bins
    np.testing.assert_allclose(summ(mock), direct.S[sf.ok])
    t = Summaries(sf, blocks=[S2()], removal=None, transform=True, s2_ref=50.)
    np.testing.assert_allclose(t(mock), np.arcsinh(direct.S[sf.ok] / 50.))


def test_bulk_removal_makes_summaries_plane_invariant():
    data = _data(3)
    sf = StructureFunction(data, n_bins=10)
    summ = Summaries(sf, blocks=['S2', 'WST', 'KURT'], removal='plane')
    tilted = data.copy()
    tilted.v = data.v + 4. * data.x - 2.5 * data.y + 30.
    np.testing.assert_allclose(summ(tilted), summ(data), rtol=1e-8, atol=1e-8)
    assert len(summ(data)) == len(summ.labels) == int(sf.ok.sum()) + 10 + 3


def test_kurtosis_gaussian_vs_heavy_tails():
    rng = np.random.default_rng(4)
    sf = StructureFunction(_data(), n_bins=10)
    k = Kurtosis(seps=(1, 2))
    summ = Summaries(sf, blocks=[k], removal=None)
    gauss = np.mean([summ(VelocityMap(np.where(sf.mask, rng.normal(0, 10, sf.mask.shape), np.nan), np.ones(sf.mask.shape),
                                      sf.mask)) for _ in range(10)], 0)
    heavy = np.mean([summ(VelocityMap(np.where(sf.mask, rng.laplace(0, 10, sf.mask.shape), np.nan), np.ones(sf.mask.shape),
                                      sf.mask)) for _ in range(10)], 0)
    assert np.all(np.abs(gauss) < 0.3) and np.all(heavy > 1.0)        # Laplace differences: excess kurtosis 3


def test_linewidth():
    data = _data()
    data.width = np.where(data.mask, np.linspace(20, 40, data.v.size).reshape(data.v.shape), np.nan)
    sf = StructureFunction(data, n_bins=8)
    w2 = data.width[sf.mask] ** 2
    lw = Summaries(sf, blocks=[LineWidth()], removal=None)
    np.testing.assert_allclose(lw(data), np.percentile(w2, [16, 50, 84]) / 1e3)
    np.testing.assert_allclose(Summaries(sf, blocks=['SIG'], removal=None, transform=True)(data),
                               np.log(np.percentile(w2, [16, 50, 84])))
    with pytest.raises(ValueError):
        lw(_data())


def test_3d_mocks_have_line_width_summary():
    obs = _obs()
    sf = StructureFunction(_data(), n_bins=8)
    summ = Summaries(sf, blocks=['S2', 'SIG'], transform=True)
    m = obs.simulate((25., 0.8, 1.2, 0.), np.random.default_rng(5), mode='3D', depth_kpc=10.)
    assert np.all(np.isfinite(summ(m)))
