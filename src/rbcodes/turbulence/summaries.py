"""
Summary statistics of velocity maps (data or mocks) for inference.

    from rbcodes.turbulence.summaries import S2, Scattering, Kurtosis, LineWidth, Summaries
    sf   = StructureFunction(vm_data, n_bins=20)                # defines pairs, bins and errors
    summ = Summaries(sf, blocks=[S2(), Scattering(sf.mask), Kurtosis(), LineWidth()],
                     removal='plane', transform=True)
    s_obs = summ(vm_data);  s_mock = summ(mock)                  # same length, same labels

Every map is measured on the data's spaxels, pairs and bins (``sf``), after fitting and removing
the same bulk-flow model (``removal``). Ported from vsf_pipeline (``Scattering``, ``Summaries``).

Blocks:

* ``S2``: debiased S2 in the bins ``kfit`` (default: all usable bins), with the data's velocity
  errors; ``transform``: asinh(S2 / s2_ref).
* ``Scattering``: reduced wavelet scattering transform (Mallat 2012; Allys et al. 2019): log S1(j),
  s2(j1, j2) = S2(j1, j2)/S1(j1); image-plane maps on a grid only.
* ``Kurtosis``: excess kurtosis of velocity increments at separations ``seps`` (pixels);
  ``transform``: asinh.
* ``LineWidth``: 16/50/84th percentiles of width^2 in the mask (units 10^3 km^2/s^2);
  ``transform``: log of width^2.
"""
import numpy as np
from scipy.ndimage import binary_erosion

from .bulk import BULK_MODELS

__all__ = ['S2', 'Scattering', 'Kurtosis', 'LineWidth', 'Summaries']


class Scattering:
    """
    Reduced wavelet scattering transform on a masked map.

    Morlet-like filters psi_{j,theta} (j = 1..J, scale 2^j px; L orientations), zero-padded to 2x the map.
    S1(j)     = < |u * psi_{j,theta}| >_{x in mask, theta}
    S2(j1,j2) = < ||u * psi_{j1,theta1}| * psi_{j2,theta2}| >_{x, theta1, theta2},  j2 > j1
    s2(j1,j2) = S2(j1,j2)/S1(j1)
    Output: [log S1(j) for j] + [s2(j1, j2) for j1 < j2]. The average uses the mask eroded by one pixel.
    """

    def __init__(self, mask, J=4, L=4, sigma0=0.8, xi0=3 * np.pi / 4):
        mask = np.asarray(mask, bool)
        ny, nx = mask.shape
        self.shape = (2 * ny, 2 * nx)
        self.J, self.L = J, L
        ky = 2 * np.pi * np.fft.fftfreq(2 * ny)[:, None]
        kx = 2 * np.pi * np.fft.fftfreq(2 * nx)[None, :]
        self.bank = np.zeros((J, L) + self.shape)
        for j in range(1, J + 1):
            xi = xi0 / 2 ** (j - 1)
            s = 1.0 / (sigma0 * 2 ** (j - 1))
            for t in range(L):
                th = np.pi * t / L
                cx, cy = xi * np.cos(th), xi * np.sin(th)
                g = np.exp(-((kx - cx) ** 2 + (ky - cy) ** 2) / (2 * s ** 2))
                g0 = np.exp(-(kx ** 2 + ky ** 2) / (2 * s ** 2))
                self.bank[j - 1, t] = g - np.exp(-(cx ** 2 + cy ** 2) / (2 * s ** 2)) * g0
        w = np.zeros(self.shape)
        w[:ny, :nx] = binary_erosion(mask, iterations=1)
        if w.sum() == 0:
            raise ValueError('The mask is empty after erosion by one pixel.')
        self.w = w / w.sum()
        self.ny, self.nx = ny, nx
        self.pairs = [(j1 + 1, j2 + 1) for j1 in range(J) for j2 in range(j1 + 1, J)]
        self.labels = [f'logS1 j={j}' for j in range(1, J + 1)] + [f's2 {p}' for p in self.pairs]

    def __call__(self, u):
        """``u``: 2-D map (NaN outside the mask)."""
        x = np.zeros(self.shape)
        x[:self.ny, :self.nx] = np.nan_to_num(u - np.nanmean(u))
        U1 = np.abs(np.fft.ifft2(np.fft.fft2(x)[None, None] * self.bank))
        S1 = (U1 * self.w).sum((-1, -2)).mean(1)
        U1h = np.fft.fft2(U1)
        s2 = []
        for j1 in range(self.J):
            for j2 in range(j1 + 1, self.J):
                U2 = np.abs(np.fft.ifft2(U1h[j1][:, None] * self.bank[j2][None]))
                s2.append((U2 * self.w).sum((-1, -2)).mean() / S1[j1])
        return np.concatenate([np.log(S1), s2])

    # block interface
    def setup(self, summ):
        return self.labels

    def compute(self, summ, res, u):
        return self(np.where(summ.sf.mask, res.v, np.nan))


class S2:
    """Debiased S2 in the bins ``kfit`` (bool over bins; default: the usable bins of ``sf``)."""

    def __init__(self, kfit=None, estimator='debiased'):
        self.kfit, self.estimator = kfit, estimator

    def setup(self, summ):
        sf = summ.sf
        self.k = sf.ok.copy() if self.kfit is None else np.asarray(self.kfit, bool)
        return [f'S2 r={r:.1f}' for r in sf.rc[self.k]]

    def compute(self, summ, res, u):
        x = summ.sf.s2(self.estimator, u=u)[self.k]
        return np.arcsinh(x / summ.s2_ref) if summ.transform else x


class Kurtosis:
    """Excess kurtosis of velocity increments at separations ``seps`` (pixels)."""

    def __init__(self, seps=(2, 4, 8)):
        self.seps = tuple(seps)

    def setup(self, summ):
        return [f'kurt s={s}px' for s in self.seps]

    def compute(self, summ, res, u):
        k = []
        for s in self.seps:
            d = summ.sf.increments(s * summ.sf.pixel, u=u)
            d = d - d.mean()
            k.append(np.mean(d ** 4) / np.mean(d ** 2) ** 2 - 3)
        k = np.array(k)
        return np.arcsinh(k) if summ.transform else k


class LineWidth:
    """16/50/84th percentiles of width^2 inside the mask (needs a width map, e.g. 3-D mocks)."""

    def setup(self, summ):
        return ['sig2 p16', 'sig2 p50', 'sig2 p84']

    def compute(self, summ, res, u):
        w = res.width[summ.sf.mask]
        if not np.isfinite(w).any():
            raise ValueError('LineWidth needs a width map.')
        x = np.nanpercentile(w ** 2, [16, 50, 84])
        return np.log(x) if summ.transform else x / 1e3


class Summaries:
    """
    Summary vector of a VelocityMap after removing a bulk-flow model.

    Parameters
    ----------
    sf : StructureFunction of the data (its spaxels, pairs, bins and errors are used for every map)
    blocks : list of block objects (S2, Scattering, Kurtosis, LineWidth), or names
        'S2', 'WST', 'KURT', 'SIG' (the vsf_pipeline names; 'WST' uses Scattering(sf.mask)).
    removal : bulk model name in ``BULK_MODELS`` (fitted to each map), a model instance, or None
    transform : variance-stabilising transforms (asinh for S2 and kurtosis, log for width^2), so that
        the summaries vary smoothly and nearly linearly with log-parameters (for emulators)
    s2_ref : scale of the S2 transform (km/s)^2
    """

    _NAMES = dict(S2=lambda sf: S2(), WST=lambda sf: Scattering(sf.mask), KURT=lambda sf: Kurtosis(),
                  SIG=lambda sf: LineWidth())

    def __init__(self, sf, blocks=('S2',), removal='plane', transform=False, s2_ref=100.):
        self.sf, self.removal, self.transform, self.s2_ref = sf, removal, transform, s2_ref
        self.blocks = [self._NAMES[b](sf) if isinstance(b, str) else b for b in blocks]
        self.labels = []
        for b in self.blocks:
            self.labels += b.setup(self)

    def residual(self, vm):
        if self.removal is None:
            return vm
        model = BULK_MODELS[self.removal]() if isinstance(self.removal, str) else self.removal
        return vm.subtract_bulk(model, self.sf.positions)

    def __call__(self, vm):
        res = self.residual(vm)
        u = self.sf.values(res)
        return np.concatenate([np.atleast_1d(b.compute(self, res, u)) for b in self.blocks])
