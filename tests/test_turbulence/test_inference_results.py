"""Known-answer tests for rbcodes.turbulence.inference and rbcodes.turbulence.results."""
import numpy as np
import pytest
from scipy.stats import multivariate_normal

from rbcodes.turbulence.forward import Observation
from rbcodes.turbulence.inference import (Emulator, adaptive_metropolis, credible_interval, fisher,
                                          gelman_rubin, latin_hypercube, moped_loglike, moped_vectors,
                                          synthetic_loglike, validate_emulator)
from rbcodes.turbulence.results import Results, asym, num, pm
from rbcodes.turbulence.summaries import Summaries
from rbcodes.turbulence.vsf import StructureFunction


def test_synthetic_loglike_matches_gaussian():
    rng = np.random.default_rng(0)
    C = np.array([[2., 0.6], [0.6, 1.]])
    sims = rng.multivariate_normal([1., -1.], C, size=20000)
    a, b = np.array([1.5, -0.2]), np.array([0., -2.])
    got = synthetic_loglike(sims, a, hartlap=False) - synthetic_loglike(sims, b, hartlap=False)
    ref = multivariate_normal([1., -1.], C).logpdf(a) - multivariate_normal([1., -1.], C).logpdf(b)
    assert got == pytest.approx(ref, abs=0.02)
    # with few simulations: Hartlap factor (N - d - 2)/(N - 1) on the quadratic form
    small = sims[:20]
    S = np.cov(small.T)
    S = S + 1e-6 * np.eye(2) * np.mean(np.diag(S))
    r = a - small.mean(0)
    h = (20 - 2 - 2) / 19
    assert synthetic_loglike(small, a) == pytest.approx(-0.5 * h * r @ np.linalg.solve(S, r)
                                                        - 0.5 * np.linalg.slogdet(S)[1])


def test_latin_hypercube_strata():
    rng = np.random.default_rng(1)
    X = latin_hypercube(50, [(0, 1), (-3, 7)], rng)
    assert X.shape == (50, 2)
    for k, (lo, hi) in enumerate([(0, 1), (-3, 7)]):
        strata = np.floor((X[:, k] - lo) / (hi - lo) * 50).astype(int)
        assert sorted(strata) == list(range(50))
    np.testing.assert_array_equal(X, latin_hypercube(50, [(0, 1), (-3, 7)], np.random.default_rng(1)))


def _toy(theta, rng, nrep):
    """s = [t0 + t1^2, 2 t1, t0 t1] + noise whose sd grows with t0."""
    t0, t1 = theta
    mu = np.array([t0 + t1 ** 2, 2 * t1, t0 * t1])
    sd = np.array([0.1, 0.2, 0.15]) * (1 + 0.5 * (t0 + 1))
    return mu + rng.normal(size=(nrep, 3)) * sd, mu, sd


def test_emulator_mean_variance_jacobian():
    rng = np.random.default_rng(2)
    bounds = [(-1, 1), (-1, 1)]
    design = latin_hypercube(150, bounds, rng)
    sims = np.array([_toy(t, rng, 30)[0] for t in design])
    emu = Emulator(bounds).fit(design, sims)
    for th in ([0.2, -0.4], [-0.7, 0.6]):
        _, mu, sd = _toy(th, rng, 1)
        np.testing.assert_allclose(emu.mean(th)[0], mu, atol=0.03)
        np.testing.assert_allclose(np.sqrt(np.diag(emu.cov(th))), sd, rtol=0.15)
        # derivative of a polynomial fitted to noisy simulations: a few hundredths of noise
        np.testing.assert_allclose(emu.jacobian(th), [[1, 2 * th[1]], [0, 2], [th[1], th[0]]], atol=0.1)
    emu_var, rep = validate_emulator(design, sims, bounds, rng)
    assert rep['std_z'] < 2.0 and rep['err_ratio'] < 0.5 and emu_var.shape == (3,)


def test_adaptive_metropolis_and_gelman_rubin():
    C = np.array([[1., 0.8], [0.8, 1.]])
    Ci = np.linalg.inv(C)
    lp = lambda x: -0.5 * x @ Ci @ x
    chains = []
    for seed in range(4):
        c, _, acc = adaptive_metropolis(lp, np.array([2., -2.]) * (seed - 1.5), 8000, np.random.default_rng(seed))
        assert 0.1 < acc < 0.7
        chains.append(c[3000:])
    s = np.vstack(chains)
    np.testing.assert_allclose(s.mean(0), 0, atol=0.1)
    np.testing.assert_allclose(np.cov(s.T), C, atol=0.12)
    assert np.all(gelman_rubin(chains) < 1.05)
    stuck = [np.random.default_rng(9).normal(m, 0.1, (500, 2)) for m in (0., 3.)]
    assert np.all(gelman_rubin(stuck) > 1.5)


class _Linear:
    """Exact linear-Gaussian 'emulator': mu = A theta + b, fixed Sigma."""
    def __init__(self, A, b, S):
        self.A, self.b, self.S, self.d = A, b, S, len(b)
        self.emu_var = np.zeros(self.d)
        self.bounds = np.array([[-10, 10]] * A.shape[1], float)
    mean = lambda self, th: np.atleast_2d(self.A @ np.asarray(th) + self.b)
    cov = lambda self, th: self.S
    jacobian = Emulator.jacobian


def test_fisher_and_moped_on_linear_gaussian():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(6, 2))
    b = rng.normal(size=6)
    M = rng.normal(size=(6, 6))
    S = M @ M.T + 6 * np.eye(6)
    emu = _Linear(A, b, S)
    F, err = fisher(emu, np.zeros(2))
    np.testing.assert_allclose(F, A.T @ np.linalg.solve(S, A), rtol=1e-6)
    np.testing.assert_allclose(err, np.sqrt(np.diag(np.linalg.inv(F))), rtol=1e-6)
    s_obs = A @ np.array([0.3, -0.5]) + b + rng.multivariate_normal(np.zeros(6), S)
    Bm, mu_fid = moped_vectors(emu, np.zeros(2))
    t_obs = Bm.T @ (s_obs - mu_fid)
    full = lambda th: multivariate_normal(emu.mean(th)[0], S).logpdf(s_obs)
    comp = lambda th: moped_loglike(emu, Bm, mu_fid, t_obs, th)
    ths = [np.zeros(2), np.array([1., 0.]), np.array([-0.4, 0.9])]
    d_full = [full(t) - full(ths[0]) for t in ths]
    d_comp = [comp(t) - comp(ths[0]) for t in ths]
    np.testing.assert_allclose(d_comp, d_full, atol=1e-6)        # lossless for a linear model


def test_credible_interval():
    x = np.linspace(-5, 5, 2001)
    lo, md, hi = credible_interval(x, np.exp(-0.5 * x ** 2))
    assert (lo, md, hi) == pytest.approx((-0.994, 0., 0.994), abs=0.01)


def test_end_to_end_recovers_injected_parameters():
    yy, xx = np.mgrid[0:24, 0:28]
    flux = np.exp(-np.hypot(xx - 14, (yy - 12) * 1.2) / 7.)
    obs = Observation(flux > 0.12, flux, 1.5, 1.0, np.full(flux.shape, 2.))
    rng = np.random.default_rng(4)
    truth = np.array([np.log10(15.), 0.8])                       # (log10 sigma, alpha); L, gamma fixed
    sim = lambda th, r: obs.simulate((10 ** th[0], th[1], 1.5, 0.), r)
    data = sim(truth, rng)
    summ = Summaries(StructureFunction(data, n_bins=10), blocks=['S2'], transform=True)
    s_obs = summ(data)
    bounds = [(np.log10(5.), np.log10(40.)), (0.3, 1.5)]
    design = latin_hypercube(60, bounds, rng)
    sims = np.array([[summ(sim(t, rng)) for _ in range(6)] for t in design])
    emu_var, rep = validate_emulator(design, sims, bounds, rng)
    emu = Emulator(bounds).fit(design, sims).set_emulator_error(emu_var)
    lo, hi = np.array(bounds).T
    lp = lambda t: emu.loglike(t, s_obs) if np.all((t > lo) & (t < hi)) else -np.inf
    chain, _, _ = adaptive_metropolis(lp, design[np.argmax([lp(t) for t in design])], 6000, rng,
                                      cov0=np.diag(((hi - lo) / 30) ** 2))
    post = chain[2000:]
    q = np.percentile(post, [2.5, 50, 97.5], 0)
    assert q[0, 0] < truth[0] < q[2, 0]                          # sigma inside the 95% interval
    assert q[0, 1] < truth[1] < q[2, 1]                          # alpha inside the 95% interval
    assert np.std(post[:, 0]) < 0.25 * (hi[0] - lo[0])           # sigma is actually constrained


def test_results_json_tex(tmp_path):
    R = Results(tmp_path / 'paper' / 'results.json', tmp_path / 'paper' / 'results.tex')
    R.update({'alpha': pm(1.08, 0.46), 'rgeo': num(7.1)}, section='nb03')
    R.update({'post-alpha': asym(0.8, 0.6, 1.1)}, section='nb04')
    R.update({'alpha': pm(0.99, 0.43)}, section='nb03')         # rerun replaces the section
    allv = R.load()
    assert set(allv) == {'nb03', 'nb04'} and 'rgeo' not in allv['nb03']
    assert R['alpha'] == '\\ensuremath{0.99 \\pm 0.43}'
    assert R['post-alpha'] == '\\ensuremath{0.80^{+0.30}_{-0.20}}'
    tex = (tmp_path / 'paper' / 'results.tex').read_text()
    assert '\\expandafter\\def\\csname R@alpha\\endcsname{\\ensuremath{0.99 \\pm 0.43}}' in tex and 'rgeo' not in tex
