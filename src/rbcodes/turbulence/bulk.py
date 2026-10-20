"""
Bulk-flow models for velocity maps.

    from rbcodes.turbulence.bulk import Plane, ArctanDisk, compare_bulk
    res = vm.subtract_bulk(Plane())                    # residual VelocityMap
    results, best = compare_bulk(vm)                   # all models, lowest BIC

Every model has ``fit(x, y, v, v_err)``, ``predict(x, y)``, ``describe()`` and ``npar``
(the protocol used by ``VelocityMap.subtract_bulk``). Positions are in the units of
``VelocityMap.points`` (kpc, or pixels without a pixel scale).

The polynomial models are least-squares fits centred on the mean position. The disk and
two-lobe models use a robust (soft L1, f_scale = 20 km/s) bounded fit. Length bounds and
starting values are given in pixels, as in the original pipeline, and multiplied by the
pixel spacing of the input positions (the smallest separation between them), so they do
not depend on the units; all of them are options.
"""
import numpy as np
from scipy.optimize import least_squares

__all__ = ['Constant', 'Plane', 'Quadratic', 'ArctanDisk', 'TwoLobe', 'BULK_MODELS', 'compare_bulk']


def _pixel(x, y):
    """Smallest non-zero separation between positions along x or y (one pixel on a grid)."""
    d = []
    for a in (x, y):
        u = np.unique(np.round(a, 9))
        if len(u) > 1:
            d.append(np.min(np.diff(u)))
    return min(d) if d else 1.


class _Poly:
    kind, npar = None, None

    def __init__(self, weighted=False):
        self.weighted = weighted

    def _design(self, x, y):
        raise NotImplementedError

    def fit(self, x, y, v, v_err=None):
        self.x0, self.y0 = float(np.mean(x)), float(np.mean(y))
        A = self._design(x - self.x0, y - self.y0)
        if self.weighted and v_err is not None:
            w = 1 / np.asarray(v_err)
            self.p, *_ = np.linalg.lstsq(A * w[:, None], v * w, rcond=None)
        else:
            self.p, *_ = np.linalg.lstsq(A, v, rcond=None)
        return self

    def predict(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        out = self._design((x - self.x0).ravel(), (y - self.y0).ravel()) @ self.p
        return out.reshape(x.shape)

    def describe(self):
        return dict(kind=self.kind, p=self.p.tolist(), centre=[self.x0, self.y0], weighted=self.weighted)


class Constant(_Poly):
    """v = c (no bulk flow; removes the mean)."""
    kind, npar = 'none', 1

    def _design(self, x, y):
        return np.c_[np.ones_like(x)]


class Plane(_Poly):
    """v = c0 + c1 x + c2 y (positions relative to the mean position)."""
    kind, npar = 'plane', 3

    def _design(self, x, y):
        return np.c_[np.ones_like(x), x, y]

    @property
    def gradient(self):
        """(|grad v| in km/s per position unit, position angle of the gradient in degrees, from +x to +y)."""
        return float(np.hypot(self.p[1], self.p[2])), float(np.degrees(np.arctan2(self.p[2], self.p[1])) % 360)


class Quadratic(_Poly):
    """v = c0 + c1 x + c2 y + c3 x^2 + c4 y^2 + c5 x y."""
    kind, npar = 'quad', 6

    def _design(self, x, y):
        return np.c_[np.ones_like(x), x, y, x * x, y * y, x * y]


def arctan_disk(p, x, y):
    """Thin rotating disk: v0 + (2/pi) V_s arctan(R/r_t) cos(phi), R = sqrt(x'^2 + (y'/q)^2)."""
    x0, y0, pa, q, vs, rt, v0 = p
    xr = (x - x0) * np.cos(pa) + (y - y0) * np.sin(pa)
    yr = -(x - x0) * np.sin(pa) + (y - y0) * np.cos(pa)
    R = np.hypot(xr, yr / q) + 1e-6
    return v0 + 2 / np.pi * vs * np.arctan(R / rt) * xr / R


def two_lobe(p, x, y):
    """Two kinematic lobes: v0 + dv tanh(s/w), s = distance along axis theta through (xc, yc)."""
    xc, yc, th, v0, dv, w = p
    return v0 + dv * np.tanh(((x - xc) * np.cos(th) + (y - yc) * np.sin(th)) / w)


class _Robust:
    kind, npar, names = None, None, ()

    def __init__(self, init=None, f_scale=20., v0_bound=400.):
        self.init, self.f_scale, self.v0_bound = init, f_scale, v0_bound

    def _plane_start(self, x, y, v):
        xm, ym = x.mean(), y.mean()
        A = np.c_[np.ones_like(x), x - xm, y - ym]
        c, *_ = np.linalg.lstsq(A, v, rcond=None)
        return xm, ym, np.arctan2(c[2], c[1]), np.ptp(A @ c) / 2 + 1., c[0]

    def _solve(self, func, p0, lo, hi, x, y, v):
        r = least_squares(lambda p: func(p, x, y) - v, np.clip(p0, lo, hi), bounds=(lo, hi),
                          loss='soft_l1', f_scale=self.f_scale)
        self.p, self.success = r.x, r.success
        return self

    def describe(self):
        return dict(kind=self.kind, p=self.p.tolist(), names=list(self.names))


class ArctanDisk(_Robust):
    """
    Thin rotating disk (``arctan_disk``); parameters x0, y0, pa, q, V_s, r_t, v0.

    rt_bounds_pix, rt0_pix: turnover-radius bounds and start in pixels (default 0.3-40, 3).
    q_bounds: axis-ratio bounds (default 0.15-1); q0 start (0.7).
    """
    kind, npar = 'arctan', 7
    names = ('x0', 'y0', 'pa', 'q', 'V_s', 'r_t', 'v0')

    def __init__(self, init=None, rt_bounds_pix=(0.3, 40.), rt0_pix=3., q_bounds=(0.15, 1.), q0=0.7,
                 f_scale=20., v0_bound=400.):
        super().__init__(init, f_scale, v0_bound)
        self.rt_bounds_pix, self.rt0_pix, self.q_bounds, self.q0 = rt_bounds_pix, rt0_pix, q_bounds, q0

    def fit(self, x, y, v, v_err=None):
        xm, ym, th0, span, c0 = self._plane_start(x, y, v)
        px = _pixel(x, y)
        p0 = self.init if self.init is not None else [xm, ym, th0, self.q0, span, self.rt0_pix * px, c0]
        lo = [x.min(), y.min(), -2 * np.pi, self.q_bounds[0], -5 * span - 50, self.rt_bounds_pix[0] * px, -self.v0_bound]
        hi = [x.max(), y.max(), 2 * np.pi, self.q_bounds[1], 5 * span + 50, self.rt_bounds_pix[1] * px, self.v0_bound]
        return self._solve(arctan_disk, p0, lo, hi, x, y, v)

    def predict(self, x, y):
        return arctan_disk(self.p, np.asarray(x, float), np.asarray(y, float))


class TwoLobe(_Robust):
    """
    Two kinematic lobes (``two_lobe``); parameters xc, yc, theta, v0, dv, w.

    w_bounds_pix, w0_pix: transition-width bounds and start in pixels (default 0.3-60, 4).
    dv_bound: |dv| bound in km/s (default 600).
    """
    kind, npar = 'twolobe', 6
    names = ('xc', 'yc', 'theta', 'v0', 'dv', 'w')

    def __init__(self, init=None, w_bounds_pix=(0.3, 60.), w0_pix=4., dv_bound=600., f_scale=20., v0_bound=400.):
        super().__init__(init, f_scale, v0_bound)
        self.w_bounds_pix, self.w0_pix, self.dv_bound = w_bounds_pix, w0_pix, dv_bound

    def fit(self, x, y, v, v_err=None):
        xm, ym, th0, span, c0 = self._plane_start(x, y, v)
        px = _pixel(x, y)
        p0 = self.init if self.init is not None else [xm, ym, th0, c0, span, self.w0_pix * px]
        lo = [x.min(), y.min(), -2 * np.pi, -self.v0_bound, -self.dv_bound, self.w_bounds_pix[0] * px]
        hi = [x.max(), y.max(), 2 * np.pi, self.v0_bound, self.dv_bound, self.w_bounds_pix[1] * px]
        return self._solve(two_lobe, p0, lo, hi, x, y, v)

    def predict(self, x, y):
        return two_lobe(self.p, np.asarray(x, float), np.asarray(y, float))


BULK_MODELS = dict(none=Constant, plane=Plane, quad=Quadratic, arctan=ArctanDisk, twolobe=TwoLobe)


def compare_bulk(vm, models=None, positions='image'):
    """
    Fit several bulk models to a VelocityMap and rank them by BIC.

    models : dict {name: model instance}; default one of each in BULK_MODELS.
    BIC = N ln(RSS/N) + k ln N over the kept spaxels (unweighted residuals).

    Returns
    -------
    results : {name: dict(model, resid (VelocityMap), rms, k, bic, dbic)}
    best : name with the lowest BIC
    """
    models = models or {k: c() for k, c in BULK_MODELS.items()}
    out = {}
    for name, mdl in models.items():
        res = vm.subtract_bulk(mdl, positions)
        r = res.v[vm.mask]
        n = len(r)
        rss = np.sum(r ** 2)
        out[name] = dict(model=mdl, resid=res, rms=np.sqrt(rss / n), k=mdl.npar,
                         bic=n * np.log(rss / n) + mdl.npar * np.log(n))
    best = min(out, key=lambda k: out[k]['bic'])
    for d in out.values():
        d['dbic'] = d['bic'] - out[best]['bic']
    return out, best
