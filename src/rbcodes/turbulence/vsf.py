"""
Second-order velocity structure function S2(r) of a VelocityMap.

    from rbcodes.turbulence.vsf import StructureFunction
    sf = StructureFunction(res, n_bins=20, estimator='debiased')     # res: bulk-subtracted VelocityMap
    sf.rc, sf.S, sf.err                                              # bin centres (kpc), S2, jackknife errors
    fit = sf.fit_powerlaw(rmin=5, rmax=10)                           # alpha, alpha_err, A
    d = sf.directional(theta)                                        # parallel / perpendicular to an axis

Pairs are formed between the kept spaxels; separations come from the map's positions:
image plane (pixels x kpc_per_pix) or source plane (beta, after a lens model). With several
lensed images, pairs across images are left out unless ``cross_images=True``.

Estimators, per pair (i, j) with increment d = v_i - v_j and errors sigma_i, sigma_j, averaged
in logarithmic separation bins:

* ``raw``:      d^2
* ``debiased``: d^2 - sigma_i^2 - sigma_j^2
* ``ivw``:      debiased, weighted by w = 1/(sigma_i^2 + sigma_j^2)^2 (inverse variance of the
                pair estimate when noise dominates)
* ``cross``:    (vA_i - vA_j)(vB_i - vB_j) for two maps with independent noise (``other=``,
                e.g. from two halves of the exposures); no noise bias.

Errors: delete-one-block jackknife covariance, with square blocks of ``block`` pixels
(image plane) or of the same size in the source plane. The jackknife only sees variance
from structure smaller than the region: for fields with power on scales comparable to or
larger than the region it underestimates the scatter between independent regions (use
forward-modelled mocks for that).
"""
import numpy as np
from scipy.optimize import least_squares

__all__ = ['StructureFunction', 'ESTIMATORS']

ESTIMATORS = ('raw', 'debiased', 'ivw', 'cross')


class StructureFunction:
    """
    Parameters
    ----------
    vm : VelocityMap
    n_bins : int
        Logarithmic separation bins between ``rmin`` and the largest pair separation.
    rmin : float, optional
        Smallest separation (position units); default one pixel.
    min_pairs : int
        Bins with fewer pairs are not used.
    block : float, optional
        Jackknife block size in pixels; default max(3, round(sqrt(N / n_blocks_target))).
    n_blocks_target : int
    positions : {'image', 'source'}
    estimator : one of ESTIMATORS
    other : VelocityMap, for estimator='cross'
    cross_images : bool
        Keep pairs between different lensed images (default False).

    Attributes
    ----------
    rc, edges, ok, cnt : bin centres, edges, usable bins, pair counts
    S, C, err : S2, its jackknife covariance and errors (for ``estimator``)
    S_raw : raw S2
    major, minor, r_geo : extent of the kept region (4 sqrt of the eigenvalues of the position
        covariance) and r_geo = 0.3 x minor
    """

    def __init__(self, vm, n_bins=14, rmin=None, min_pairs=100, block=None, n_blocks_target=16,
                 positions='image', estimator='debiased', other=None, cross_images=False):
        if estimator not in ESTIMATORS:
            raise ValueError(f'estimator must be one of {ESTIMATORS}.')
        if (estimator == 'cross') != (other is not None):
            raise ValueError("estimator='cross' needs other= (a second map); other= only with 'cross'.")
        self.estimator, self.positions = estimator, positions
        mask = vm.mask & (other.mask if other is not None else True)
        if positions == 'image':
            unit = 1. if vm.kpc_per_pix is None else vm.kpc_per_pix
            X, Y = vm.x * unit, vm.y * unit
        elif positions == 'source':
            if not vm.has_source_positions:
                raise ValueError('No source-plane positions; apply a lens model first.')
            X, Y = vm.beta_x, vm.beta_y
            mask = mask & np.isfinite(X) & np.isfinite(Y)
        else:
            raise ValueError("positions must be 'image' or 'source'.")
        self.x, self.y = X[mask], Y[mask]
        self.u, self.e2 = vm.v[mask], vm.v_err[mask] ** 2
        self.u2 = other.v[mask] if other is not None else None
        img = vm.image[mask] if vm.image is not None else None
        n = len(self.u)
        if n < 3:
            raise ValueError(f'Only {n} spaxels in the map.')

        # pixel spacing in position units
        if positions == 'image':
            self.pixel = unit
        else:
            from scipy.spatial import cKDTree
            d, _ = cKDTree(np.c_[self.x, self.y]).query(np.c_[self.x, self.y], k=2)
            self.pixel = float(np.median(d[:, 1]))

        i, j = np.triu_indices(n, 1)
        if img is not None and not cross_images:
            same = img[i] == img[j]
            i, j = i[same], j[same]
        self.i, self.j = i, j
        dx, dy = self.x[j] - self.x[i], self.y[j] - self.y[i]
        self.r = np.hypot(dx, dy)
        self.ang = np.arctan2(dy, dx) % np.pi

        rmin = self.pixel if rmin is None else rmin
        edges = np.logspace(np.log10(rmin * 0.999), np.log10(self.r.max() * 1.001), n_bins + 1)
        b = np.digitize(self.r, edges) - 1
        b[(b < 0) | (b >= n_bins)] = -1
        cnt = np.bincount(b[b >= 0], minlength=n_bins)
        ok = cnt >= min_pairs
        self.edges, self.rc, self.ok, self.nb = edges, np.sqrt(edges[1:] * edges[:-1]), ok, n_bins
        self.b = np.where((b >= 0) & ok[np.clip(b, 0, n_bins - 1)], b, -1)
        self.cnt = np.bincount(self.b[self.b >= 0], minlength=n_bins)

        # jackknife blocks (square, `block` pixels on a side)
        if block is None:
            block = max(3, int(round(np.sqrt(n / n_blocks_target))))
        self.block = block
        size = block * self.pixel
        x0, y0 = (0., 0.) if positions == 'image' else (self.x.min(), self.y.min())
        bx = np.floor((self.x - x0) / size + 1e-9).astype(int)
        by = np.floor((self.y - y0) / size + 1e-9).astype(int)
        ub, inv = np.unique(by * 100000 + bx, return_inverse=True)
        members = np.bincount(inv)
        self.blocks = [np.where(inv == k)[0] for k in range(len(ub)) if members[k] >= 0.25 * block ** 2]

        cxy = np.cov(np.vstack([self.x, self.y]))
        ev = np.sort(np.linalg.eigvalsh(cxy))
        self.minor, self.major = 4 * np.sqrt(ev[0]), 4 * np.sqrt(ev[1])
        self.r_geo = 0.3 * self.minor

        self.S, self.C = self.jackknife()
        self.err = np.sqrt(np.diag(self.C))
        self.S_raw = self.s2('raw')

    def __repr__(self):
        return (f'StructureFunction({self.estimator}, {self.positions} plane, {len(self.u)} spaxels, '
                f'{len(self.i)} pairs, {int(self.ok.sum())} bins, {len(self.blocks)} jackknife blocks)')

    # --- estimators -----------------------------------------------------------

    def _pair_values(self, estimator):
        i, j = self.i, self.j
        if estimator == 'cross':
            return (self.u[i] - self.u[j]) * (self.u2[i] - self.u2[j]), None
        d2 = (self.u[i] - self.u[j]) ** 2
        if estimator == 'raw':
            return d2, None
        e2 = self.e2[i] + self.e2[j]
        return d2 - e2, (1 / e2 ** 2 if estimator == 'ivw' else None)

    def s2(self, estimator=None, sel=None):
        """S2 per bin; ``sel`` is a boolean over pairs (e.g. a direction or a jackknife subset)."""
        est = self.estimator if estimator is None else estimator
        if est == 'cross' and self.u2 is None:
            raise ValueError("The 'cross' estimator needs other=.")
        val, w = self._pair_values(est)
        b = self.b if sel is None else np.where(sel, self.b, -1)
        g = b >= 0
        w = np.ones_like(val) if w is None else w
        num = np.bincount(b[g], weights=(w * val)[g], minlength=self.nb)
        den = np.bincount(b[g], weights=w[g], minlength=self.nb)
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(den > 0, num / den, np.nan)

    def jackknife(self, estimator=None, sel=None):
        """S2 and its delete-one-block jackknife covariance."""
        full = self.s2(estimator, sel)
        keep = np.ones(len(self.u), bool)
        reps = []
        for blk in self.blocks:
            keep[:] = True
            keep[blk] = False
            s = keep[self.i] & keep[self.j]
            reps.append(self.s2(estimator, s if sel is None else (s & sel)))
        reps = np.array(reps)
        nb = len(reps)
        D = np.nan_to_num(reps - np.nanmean(reps, 0))
        return full, (nb - 1) / nb * D.T @ D

    def directional(self, theta, tol=np.pi / 6, estimator=None):
        """
        S2 from pairs parallel (|angle - theta| < tol) and perpendicular (> pi/2 - tol) to the
        axis at angle ``theta`` (radians, from +x towards +y). Returns
        dict(parallel=(S, C), perpendicular=(S, C)).
        """
        dang = np.abs(((self.ang - theta) + np.pi / 2) % np.pi - np.pi / 2)
        return dict(parallel=self.jackknife(estimator, dang < tol),
                    perpendicular=self.jackknife(estimator, dang > np.pi / 2 - tol))

    def increments(self, sep, tol=None):
        """Velocity increments v_i - v_j of pairs with |r - sep| < tol (default half a pixel)."""
        tol = 0.5 * self.pixel if tol is None else tol
        sel = np.abs(self.r - sep) < tol
        return (self.u[self.i] - self.u[self.j])[sel]

    # --- descriptive fits ---------------------------------------------------

    def fit_range(self, rmin_pix=2., rmax_frac=0.5):
        """Default separation range for slope fits: [rmin_pix pixels, rmax_frac x largest separation]."""
        return rmin_pix * self.pixel, rmax_frac * self.r.max()

    def fit_powerlaw(self, rmin, rmax, r0=10., estimator=None, noise_floor=False, psf_fwhm=None):
        """
        Descriptive fit S2 = A (r/r0)^alpha [+ N] over rmin <= r <= rmax, in linear space with
        the (shrunk) jackknife covariance.

        noise_floor : fit S2 = A (r/r0)^alpha + N (use with the raw estimator to compare N
            with 2 <sigma_v^2>).
        psf_fwhm : fit the power law smoothed by a Gaussian PSF of this FWHM (position units),
            for a field of uniform brightness: S2(r) = A I(r) / I_0(r0), with
            I(r) = int_0^inf k^-(alpha+1) exp(-k^2 s^2) [1 - J0(k r)] dk, s = FWHM/2.355,
            I_0 the same with s = 0 (so A is the unsmoothed amplitude at r0). 0 < alpha < 2.
            The smoothed model is averaged over the pair separations in each bin (8 sub-bins
            weighted by pair counts), since it curves strongly within a bin; the plain power
            law is evaluated at the bin centres.

        Returns dict(alpha, alpha_err, A, A_err, N, N_err, chi2, nbins).
        """
        est = self.estimator if estimator is None else estimator
        S, C = (self.S, self.C) if est == self.estimator else self.jackknife(est)
        rc = self.rc
        k = (rc >= rmin) & (rc <= rmax) & np.isfinite(S) & self.ok
        out = dict(alpha=np.nan, alpha_err=np.nan, A=np.nan, A_err=np.nan, N=np.nan, N_err=np.nan,
                   chi2=np.nan, nbins=int(k.sum()))
        if k.sum() < 3 + bool(noise_floor):
            return out
        Ck = C[np.ix_(k, k)]
        Ck = 0.5 * Ck + 0.5 * np.diag(np.diag(Ck))
        Ck = Ck + np.eye(k.sum()) * (1e-6 * np.max(np.diag(Ck)) + 1e-9)
        Li = np.linalg.inv(np.linalg.cholesky(Ck))
        r, y = rc[k], S[k]
        if psf_fwhm:
            s = psf_fwhm / 2.3548
            rs, ws = self._sub_bins(np.where(k)[0])
            shape = lambda a: (ws @ _smoothed_powerlaw(rs, a, s)) / _smoothed_powerlaw(np.array([r0]), a, 0.)[0]
            a_lo, a_hi = 0.05, 1.95
        else:
            shape = lambda a: (r / r0) ** a
            a_lo, a_hi = -1., 3.
        A0 = max(np.nanmedian(y), 1.)
        if noise_floor:
            f = lambda p: Li @ (p[0] * shape(p[1]) + p[2] - y)
            p0, lo, hi = [A0, 0.7, 0.], [0, a_lo, 0], [np.inf, a_hi, np.inf]
        else:
            f = lambda p: Li @ (p[0] * shape(p[1]) - y)
            p0, lo, hi = [A0, 0.7], [0, a_lo], [np.inf, a_hi]
        res = least_squares(f, p0, bounds=(lo, hi))
        cov = np.linalg.pinv(res.jac.T @ res.jac)
        e = np.sqrt(np.clip(np.diag(cov), 0, None))
        out.update(alpha=res.x[1], alpha_err=e[1], A=res.x[0], A_err=e[0], chi2=float(np.sum(res.fun ** 2)))
        if noise_floor:
            out.update(N=res.x[2], N_err=e[2])
        return out

    def _sub_bins(self, bins, nsub=8):
        """Sub-bin centres and, per bin, the pair-count weights (rows sum to 1)."""
        centres, W = [], np.zeros((len(bins), len(bins) * nsub))
        for n, b in enumerate(bins):
            e = np.geomspace(self.edges[b], self.edges[b + 1], nsub + 1)
            c, _ = np.histogram(self.r[self.b == b], e)
            centres.append(np.sqrt(e[1:] * e[:-1]))
            W[n, n * nsub:(n + 1) * nsub] = c / max(c.sum(), 1)
        return np.concatenate(centres), W

    def r_turn(self, rmax=None, frac=0.9, S=None):
        """First separation (<= rmax) where S2 reaches frac x its maximum; log-interpolated."""
        S = self.S if S is None else S
        rmax = self.r.max() if rmax is None else rmax
        k = np.isfinite(S) & (self.rc <= rmax) & self.ok
        r, y = self.rc[k], S[k]
        if len(y) < 2:
            return np.nan
        top = np.nanmax(y)
        i = np.where(y >= frac * top)[0][0]
        if i == 0:
            return r[0]
        f = (frac * top - y[i - 1]) / (y[i] - y[i - 1])
        return 10 ** (np.log10(r[i - 1]) + f * np.log10(r[i] / r[i - 1]))


def _smoothed_powerlaw(r, alpha, s, X=200.):
    """
    I(r) = int_0^inf k^-(alpha+1) exp(-k^2 s^2) [1 - J0(k r)] dk = r^alpha c(alpha, s/r), with
    c(alpha, t) = int_0^inf x^-(alpha+1) exp(-x^2 t^2) [1 - J0(x)] dx: integrated numerically to
    x = X, beyond which J0 is dropped (its contribution is < X^-(alpha+1.5)).
    """
    from scipy.integrate import quad
    from scipy.special import j0
    r = np.atleast_1d(np.asarray(r, dtype=float))
    out = np.empty_like(r)
    for n, rr in enumerate(r):
        t = s / rr
        g = (lambda x: np.exp(-(x * t) ** 2)) if t > 0 else (lambda x: 1.)
        head = quad(lambda x: x ** (-alpha - 1) * g(x) * (1 - j0(x)), 0, X, limit=500)[0]
        tail = (quad(lambda x: x ** (-alpha - 1) * g(x), X, np.inf)[0] if t > 0 else X ** (-alpha) / alpha)
        out[n] = rr ** alpha * (head + tail)
    return out
