"""
Simulation-based inference with emulated synthetic likelihoods.

    from rbcodes.turbulence.inference import (latin_hypercube, Emulator, validate_emulator,
                                              adaptive_metropolis, gelman_rubin, fisher, moped_vectors,
                                              moped_loglike, synthetic_loglike)
    design = latin_hypercube(400, bounds, rng)                     # theta points
    sims   = ...                                                   # (n_design, n_rep, n_summaries)
    emu_var, report = validate_emulator(design, sims, bounds, rng) # held-out check, emulator error
    emu    = Emulator(bounds).fit(design, sims).set_emulator_error(emu_var)
    chain, lp, acc = adaptive_metropolis(lambda t: emu.loglike(t, s_obs) if inside(t) else -np.inf, x0, 20000, rng)

Ported from vsf_pipeline (``synthetic_loglike``, ``latin_hypercube``, ``Emulator``,
``adaptive_metropolis``, ``gelman_rubin``, ``moped_vectors``; and the emulator validation, Fisher
and MOPED steps of notebook 04). Random number generators are explicit arguments.
"""
from itertools import combinations_with_replacement

import numpy as np

__all__ = ['synthetic_loglike', 'latin_hypercube', 'Emulator', 'validate_emulator', 'adaptive_metropolis',
           'gelman_rubin', 'fisher', 'moped_vectors', 'moped_loglike', 'credible_interval']


def synthetic_loglike(sims, obs, hartlap=True, ridge=1e-6):
    """Gaussian synthetic log-likelihood (Wood 2010) from sims (n, d), with the Hartlap et al. (2007) factor."""
    N, d = sims.shape
    mu = sims.mean(0)
    S = np.cov(sims.T)
    S = S + ridge * np.eye(d) * np.mean(np.diag(S)) if d > 1 else np.atleast_2d(S)
    h = max((N - d - 2) / (N - 1), 0.05) if hartlap else 1.
    r = obs - mu
    return -0.5 * h * r @ np.linalg.solve(S, r) - 0.5 * np.linalg.slogdet(S)[1]


def latin_hypercube(n, bounds, rng):
    """Latin-hypercube design (McKay et al. 1979) of n points in the box ``bounds`` = [(lo, hi), ...]."""
    p = len(bounds)
    u = np.zeros((n, p))
    for k in range(p):
        u[:, k] = (rng.permutation(n) + rng.random(n)) / n
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    return lo + u * (hi - lo)


def _poly_features(X, degree):
    cols = [np.ones(len(X))]
    for d in range(1, degree + 1):
        for c in combinations_with_replacement(range(X.shape[1]), d):
            cols.append(np.prod(X[:, c], axis=1))
    return np.column_stack(cols)


class Emulator:
    """
    Smooth emulator of the summary distribution p(s | theta) ~ N(mu(theta), Sigma(theta)).

    mu(theta): ridge-regularised polynomial (``degree``) in theta scaled to [-1, 1];
    Sigma(theta) = D(theta) R D(theta): per-summary log-variance polynomial (degree 2) and a pooled
    correlation matrix R. (Heitmann et al. 2009 for emulation; Price et al. 2018 for synthetic
    likelihoods.) ``set_emulator_error`` adds a per-summary emulator variance to the likelihood.
    """

    def __init__(self, bounds, degree=3, ridge=1e-4):
        self.bounds = np.array(bounds, float)
        self.degree, self.ridge = degree, ridge

    def _x(self, theta):
        th = np.atleast_2d(theta)
        lo, hi = self.bounds[:, 0], self.bounds[:, 1]
        return 2 * (th - lo) / (hi - lo) - 1

    def fit(self, thetas, sims):
        """thetas (n, p); sims (n, nrep, d)."""
        X = _poly_features(self._x(thetas), self.degree)
        m = sims.mean(1)
        lv = np.log(sims.var(1, ddof=1) + 1e-12)
        A = X.T @ X + self.ridge * np.eye(X.shape[1])
        self.cm = np.linalg.solve(A, X.T @ m)
        X2 = _poly_features(self._x(thetas), 2)
        self.cv = np.linalg.solve(X2.T @ X2 + self.ridge * np.eye(X2.shape[1]), X2.T @ lv)
        z = (sims - m[:, None]) / np.sqrt(np.exp(lv))[:, None]
        z = z.reshape(-1, sims.shape[2])
        self.R = np.corrcoef(z.T)
        self.d = sims.shape[2]
        if self.d == 1:
            self.R = np.ones((1, 1))
        self.nrep = sims.shape[1]
        self.emu_var = np.zeros(self.d)
        return self

    def set_emulator_error(self, var):
        """Add an emulator-error variance (per summary) to the likelihood covariance."""
        self.emu_var = np.asarray(var, float)
        return self

    def mean(self, theta):
        return _poly_features(self._x(theta), self.degree) @ self.cm

    def cov(self, theta):
        sd = np.sqrt(np.exp(_poly_features(self._x(theta), 2) @ self.cv))[0]
        return self.R * np.outer(sd, sd)

    def loglike(self, theta, obs):
        mu = self.mean(theta)[0]
        S = self.cov(theta) + np.diag(self.emu_var) + 1e-9 * np.eye(self.d)
        r = obs - mu
        try:
            L = np.linalg.cholesky(S)
        except np.linalg.LinAlgError:
            return -np.inf
        z = np.linalg.solve(L, r)
        return -0.5 * z @ z - np.sum(np.log(np.diag(L)))

    def jacobian(self, theta, eps=1e-3):
        """d mu / d theta by central differences (d x p)."""
        theta = np.asarray(theta, float)
        J = np.zeros((self.d, len(theta)))
        for k in range(len(theta)):
            h = eps * (self.bounds[k, 1] - self.bounds[k, 0])
            tp, tm = theta.copy(), theta.copy()
            tp[k] += h
            tm[k] -= h
            J[:, k] = (self.mean(tp)[0] - self.mean(tm)[0]) / (2 * h)
        return J


def validate_emulator(design, sims, bounds, rng, frac_test=1 / 7, min_test=10, **emulator_kw):
    """
    Held-out check of an emulator and its error variance (notebook 04 procedure).

    A fraction of the design points is left out, the emulator is fitted to the rest, and the
    simulated means at the held-out points are compared with the emulated means:
    z = resid / (sd/sqrt(nrep)) (should have std ~1 if the emulator error is negligible) and
    zr = resid / sd (emulator error relative to the noise of one map).
    The emulator-error variance per summary is mean(resid^2) - mean(sd^2)/nrep (clipped at 0).

    Returns (emu_var, dict(z, zr, test, train, std_z, err_ratio, mu_test, sd_test)).
    """
    n, nrep = sims.shape[:2]
    perm = rng.permutation(n)
    nte = max(min_test, int(n * frac_test))
    te, tr = perm[:nte], perm[nte:]
    emu = Emulator(bounds, **emulator_kw).fit(design[tr], sims[tr])
    mu_te = emu.mean(design[te])
    sd_te = np.array([np.sqrt(np.diag(emu.cov(t))) for t in design[te]])
    resid = sims[te].mean(1) - mu_te
    z = resid / (sd_te / np.sqrt(nrep))
    zr = resid / sd_te
    emu_var = np.clip(np.mean(resid ** 2, 0) - np.mean(sd_te ** 2, 0) / nrep, 0, None)
    return emu_var, dict(z=z, zr=zr, test=te, train=tr, std_z=float(np.std(z)),
                         err_ratio=float(np.sqrt(np.mean(zr ** 2))), mu_test=mu_te, sd_test=sd_te)


def adaptive_metropolis(logpost, x0, n_steps, rng, cov0=None, adapt_start=500, scale=None):
    """
    Adaptive Metropolis (Haario, Saksman & Tamminen 2001): Gaussian proposal whose covariance is
    the running chain covariance x 2.38^2/p. Returns (chain, log-posterior values, acceptance).
    """
    p = len(x0)
    x = np.array(x0, float)
    lp = logpost(x)
    cov = np.diag((0.05 * np.abs(x) + 0.01) ** 2) if cov0 is None else cov0
    sc = (2.38 ** 2 / p) if scale is None else scale
    chain = np.zeros((n_steps, p))
    lps = np.zeros(n_steps)
    acc = 0
    for t in range(n_steps):
        if t >= adapt_start and t % 50 == 0:
            cov = np.cov(chain[t // 2:t].T) + 1e-8 * np.eye(p)
        prop = rng.multivariate_normal(x, sc * cov)
        lpp = logpost(prop)
        if np.log(rng.random()) < lpp - lp:
            x, lp = prop, lpp
            acc += 1
        chain[t], lps[t] = x, lp
    return chain, lps, acc / n_steps


def gelman_rubin(chains):
    """Potential scale reduction R-hat (Gelman & Rubin 1992) for a list of equal-length chains (n, p)."""
    c = np.array(chains)
    m, n, p = c.shape
    W = c.var(1, ddof=1).mean(0)
    B = n * c.mean(1).var(0, ddof=1)
    return np.sqrt(((n - 1) / n * W + B / n) / W)


def fisher(emu, theta, idx=None):
    """
    Fisher matrix F = J^T Sigma^-1 J at ``theta`` from the emulator (J = d mu / d theta; Sigma
    including the emulator error), for the summaries ``idx`` (default all).
    Returns (F, marginal errors sqrt(diag(F^-1))).
    """
    J = emu.jacobian(theta)
    S = emu.cov(theta) + np.diag(emu.emu_var) + 1e-9 * np.eye(emu.d)
    idx = np.arange(emu.d) if idx is None else np.asarray(idx)
    F = J[idx].T @ np.linalg.solve(S[np.ix_(idx, idx)], J[idx])
    return F, np.sqrt(np.diag(np.linalg.pinv(F)))


def moped_vectors(emu, theta_fid):
    """
    MOPED / score compression (Heavens et al. 2000; Alsing & Wandelt 2018):
    t = (d mu/d theta)^T Sigma^-1 (s - mu(theta_fid)) — one number per parameter, lossless at
    theta_fid for Gaussian s. Returns (B (d x p), mu(theta_fid)); t = B^T (s - mu_fid).
    """
    J = emu.jacobian(theta_fid)
    S = emu.cov(theta_fid) + 1e-9 * np.eye(emu.d)
    return np.linalg.solve(S, J), emu.mean(theta_fid)[0]


def moped_loglike(emu, B, mu_fid, t_obs, theta):
    """Gaussian log-likelihood of the compressed summaries t_obs = B^T (s_obs - mu_fid) at theta."""
    mt = B.T @ (emu.mean(theta)[0] - mu_fid)
    Ct = B.T @ (emu.cov(theta) + np.diag(emu.emu_var)) @ B + 1e-9 * np.eye(len(t_obs))
    r = t_obs - mt
    return -0.5 * r @ np.linalg.solve(Ct, r) - 0.5 * np.linalg.slogdet(Ct)[1]


def credible_interval(x, p, q=(0.16, 0.5, 0.84)):
    """Quantiles of a gridded 1-D posterior p(x) (notebook 04 ``cred``)."""
    c = np.cumsum(p) / np.sum(p)
    return np.interp(q, c, x)
