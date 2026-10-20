"""
Forward model: turbulent velocity fields seen through the instrument, as mock VelocityMaps.

    from rbcodes.turbulence.forward import TurbulenceField, Observation
    obs  = Observation.from_vmap(vm, psf_fwhm_kpc=5.2)            # mask, flux, errors of the data
    mock = obs(TurbulenceField(sigma=10., alpha=2/3, L_kpc=20.), rng)   # VelocityMap
    mock = obs.simulate((10., 2/3, np.log10(20.), 0.), rng)            # theta = (sigma, alpha, log10 L, gamma)

Field: Gaussian random field with P(k) ~ (k^2 + k0^2)^-(alpha+D)/2, k0 = 1/L, in D = 2 dimensions
(the map) or D = 3 (a slab of depth ``depth_kpc`` along the line of sight). Its structure function
scales as r^alpha for r << L. Toy intermittency: g -> g exp(gamma g2 - gamma^2/2), with g2 an
independent smooth field of scale ``L_env``.

Observation: flux-weighted PSF average, (F v (x) P) / (F (x) P); mean removed over the mask;
optional bulk-flow model added; Gaussian noise of the per-spaxel velocity errors.

Lensed systems (``beta=``): the field is made on a source-plane grid (``src_kpp`` kpc per pixel)
covering the source-plane positions beta of the spaxels (e.g. from
``rbcodes.lensing.lens_model.LensModel.raytrace``), sampled at each spaxel's beta (bilinear), and
then observed in the image plane as above (PSF, mask, noise). Spaxels of different images that map
to the same beta get the same field value. Spaxels without beta get zero weight.

* mode '2D': sigma is the rms of the *observed* (PSF-weighted) field in the mask.
* mode '3D': sigma is the rms of the 3-D field; the map is the PSF-weighted line-of-sight mean,
  and the width map is sigma_obs^2 = W(sigma_los^2 + v_c^2) - W(v_c)^2 + sigma_th^2 (W = the
  flux-weighted PSF average).

Ported from vsf_pipeline (``grf``, ``psf_weight``, ``MockMaker``).
"""
import numpy as np
from scipy.ndimage import gaussian_filter

from rbcodes.ifu.vmap import VelocityMap

__all__ = ['grf', 'psf_weight', 'TurbulenceField', 'Observation']

_KCACHE = {}


def _kgrid(shape, kpp):
    key = (shape, round(kpp, 6))
    if key not in _KCACHE:
        ks = np.meshgrid(*[np.fft.fftfreq(n) / kpp for n in shape], indexing='ij')
        _KCACHE[key] = sum(k ** 2 for k in ks)
    return _KCACHE[key]


def grf(alpha, L_kpc, shape, kpp, rng, pad=4):
    """
    Gaussian random field with P(k) ~ (k^2 + k0^2)^-(alpha+D)/2, k0 = 1/L, D = len(shape).

    Generated in a box ``pad`` times larger (per axis) and cut at a random offset, so large modes
    are present and the field is not periodic over the map. ``L_kpc`` None or inf: pure power law.
    """
    D = len(shape)
    big = tuple(n * pad for n in shape)
    k2 = _kgrid(big, kpp)
    k02 = 0. if (L_kpc is None or not np.isfinite(L_kpc)) else (1. / L_kpc) ** 2
    with np.errstate(divide='ignore'):
        P = (k2 + k02) ** (-(alpha + D) / 2.)
    P.flat[0] = 0.
    w = rng.normal(size=big) + 1j * rng.normal(size=big)
    f = np.fft.ifftn(np.sqrt(P) * w).real
    sl = tuple(slice(o, o + n) for o, n in zip([rng.integers(0, b - n + 1) for b, n in zip(big, shape)], shape))
    return f[sl]


def psf_weight(x, flux, psf_sig):
    """Flux-weighted PSF average (F x (x) P) / (F (x) P); ``psf_sig`` = Gaussian sigma in pixels."""
    num = gaussian_filter(np.nan_to_num(flux * x), psf_sig)
    den = gaussian_filter(np.nan_to_num(flux), psf_sig)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 0, num / den, 0.)


class TurbulenceField:
    """
    Parameters of a turbulent velocity field.

    sigma : rms (km/s; see the module docstring for what it refers to in each mode)
    alpha : structure-function slope
    L_kpc : outer (driving) scale; None/inf for a pure power law
    gamma : toy intermittency (0 = Gaussian)
    mode : '2D' or '3D'
    depth_kpc : line-of-sight depth of the 3-D slab
    L_env : scale of the intermittency envelope (kpc)
    pad : box padding factor for the 2-D field (the 3-D slab uses 2)
    """

    def __init__(self, sigma, alpha, L_kpc, gamma=0., mode='2D', depth_kpc=20., L_env=10., pad=4):
        if mode not in ('2D', '3D'):
            raise ValueError("mode must be '2D' or '3D'.")
        self.sigma, self.alpha, self.L_kpc, self.gamma = float(sigma), float(alpha), L_kpc, float(gamma)
        self.mode, self.depth_kpc, self.L_env, self.pad = mode, float(depth_kpc), float(L_env), pad

    @classmethod
    def from_theta(cls, theta, **kw):
        """theta = (sigma, alpha, log10 L, gamma)."""
        sig, alpha, logL, gamma = theta
        return cls(sig, alpha, 10 ** logL, gamma, **kw)

    def _modulate(self, g, rng, shape, kpp, pad):
        if self.gamma <= 0:
            return g
        g2 = grf(1.6, self.L_env, shape, kpp, rng, pad)
        g2 = (g2 - g2.mean()) / g2.std()
        return g * np.exp(self.gamma * g2 - self.gamma ** 2 / 2)

    def draw(self, shape, kpp, rng):
        """
        One realization on a map of ``shape`` with ``kpp`` kpc per pixel.

        '2D': dict(g) — the (unnormalised) field; '3D': dict(vc, s2los) — the line-of-sight
        mean and variance of the slab field scaled to rms ``sigma``.
        """
        if self.mode == '2D':
            g = grf(self.alpha, self.L_kpc, shape, kpp, rng, self.pad)
            return dict(g=self._modulate(g, rng, shape, kpp, self.pad))
        nz = max(8, int(round(self.depth_kpc / kpp)))
        sh = (nz,) + tuple(shape)
        g = grf(self.alpha, self.L_kpc, sh, kpp, rng, pad=2)
        g = self._modulate(g, rng, sh, kpp, 2)
        g = g * self.sigma / g.std()
        return dict(vc=g.mean(0), s2los=g.var(0))


class Observation:
    """
    How a field is observed.

    Parameters
    ----------
    mask : bool (ny, nx), spaxels with measurements
    flux : (ny, nx) line flux used for the flux-weighted PSF average
    kpp : kpc per pixel
    psf_sig : PSF Gaussian sigma in pixels
    verr : (ny, nx) velocity errors (km/s) for the noise
    bulk_model : (ny, nx), optional, added to the turbulent field
    sigerr : (ny, nx), optional, width errors (3-D mode)
    sig_th : thermal width added in quadrature in 3-D mode (km/s; default 2.3)
    line, z, wcs : metadata copied to the mock maps
    beta : (beta_x, beta_y), optional
        Source-plane positions (kpc, 2-D like the mask, NaN where unknown) of a lensed system.
        Spaxels without beta are removed from the mask.
    image : (ny, nx) int, optional, multiple-image ids (copied to the mocks)
    src_kpp : source-plane grid spacing (kpc); default: the median distance between neighbouring
        spaxels in the source plane
    """

    def __init__(self, mask, flux, kpp, psf_sig, verr, bulk_model=None, sigerr=None, sig_th=2.3,
                 line='mock', z=None, wcs=None, beta=None, image=None, src_kpp=None):
        self.mask, self.flux, self.kpp, self.psf_sig = np.asarray(mask, bool), np.nan_to_num(flux), kpp, psf_sig
        self.shape = self.mask.shape
        self.beta, self.image, self.src_kpp = None, image, src_kpp
        if beta is not None:
            bx, by = (np.asarray(b, dtype=float) for b in beta)
            if bx.shape != self.shape or by.shape != self.shape:
                raise ValueError('beta arrays must have the shape of the mask.')
            have = np.isfinite(bx) & np.isfinite(by)
            self.mask = self.mask & have
            self.flux = np.where(have, self.flux, 0.)
            self.beta = (bx, by)
            if src_kpp is None:
                from scipy.spatial import cKDTree
                # distinct positions only: spaxels of different images can share a beta
                pts = np.unique(np.round(np.c_[bx[self.mask], by[self.mask]], 9), axis=0)
                if len(pts) < 2:
                    raise ValueError('Fewer than two distinct source-plane positions in the mask.')
                d, _ = cKDTree(pts).query(pts, k=2)
                self.src_kpp = float(np.median(d[:, 1]))
            if not self.src_kpp > 0:
                raise ValueError(f'src_kpp must be positive, got {self.src_kpp}.')
        self.verr = np.nan_to_num(np.broadcast_to(verr, self.shape)).astype(float)
        self.sigerr = None if sigerr is None else np.nan_to_num(np.broadcast_to(sigerr, self.shape)).astype(float)
        self.bulk = np.zeros(self.shape) if bulk_model is None else np.nan_to_num(bulk_model)
        self.sig_th, self.line, self.z, self.wcs = sig_th, line, z, wcs

    @classmethod
    def from_vmap(cls, vm, psf_fwhm_kpc=None, flux=None, bulk_model=None, positions='image', **kw):
        """Mask, flux, velocity errors, scale and PSF of a data VelocityMap; ``positions='source'``
        also takes its source-plane positions and image ids (a lensed system)."""
        if vm.kpc_per_pix is None:
            raise ValueError('The map needs kpc_per_pix.')
        fwhm = vm.psf_fwhm_kpc if psf_fwhm_kpc is None else psf_fwhm_kpc
        if fwhm is None:
            raise ValueError('Give psf_fwhm_kpc (the map has none).')
        flux = vm.flux if flux is None else flux
        if positions == 'source':
            if not vm.has_source_positions:
                raise ValueError('The map has no source-plane positions; ray-trace it first.')
            kw.setdefault('beta', (vm.beta_x, vm.beta_y))
            kw.setdefault('image', vm.image)
        return cls(vm.mask, flux, vm.kpc_per_pix, fwhm / 2.3548 / vm.kpc_per_pix, vm.v_err,
                   bulk_model=bulk_model, sigerr=np.where(np.isfinite(vm.width_err), vm.width_err, 0.),
                   line=vm.line, z=vm.z, wcs=vm.wcs, **kw)

    def _draw(self, field, rng, kpp):
        """The field on the image grid: drawn there, or (lensed) drawn in the source plane and sampled at beta."""
        if self.beta is None:
            return field.draw(self.shape, self.kpp if kpp is None else kpp, rng)
        if kpp is not None and kpp != self.kpp:
            raise ValueError('kpp cannot be changed for a lensed observation (positions come from beta).')
        from scipy.ndimage import map_coordinates
        bx, by = self.beta
        have = np.isfinite(bx) & np.isfinite(by)
        h = self.src_kpp
        x0, y0 = np.nanmin(bx) - 2 * h, np.nanmin(by) - 2 * h
        nx = int(np.ceil((np.nanmax(bx) - x0) / h)) + 3
        ny = int(np.ceil((np.nanmax(by) - y0) / h)) + 3
        if ny * nx > 4_000_000:
            raise ValueError(f'Source-plane grid of {ny} x {nx} pixels at src_kpp={h:.3g} kpc is too large; '
                             'give a larger src_kpp.')
        d = field.draw((ny, nx), h, rng)
        coords = np.vstack([((by[have] - y0) / h), ((bx[have] - x0) / h)])
        out = {}
        for k, f in d.items():
            a = np.zeros(self.shape)
            a[have] = map_coordinates(f, coords, order=1, mode='nearest')
            out[k] = a
        return out

    def observe(self, field, rng, kpp=None):
        """Noise-free observed turbulent maps: (v, width or None)."""
        d = self._draw(field, rng, kpp)
        if field.mode == '2D':
            t = psf_weight(d['g'], self.flux, self.psf_sig)
            t = t - t[self.mask].mean()
            return t * field.sigma / max(t[self.mask].std(), 1e-12), None
        vc, s2los = d['vc'], d['s2los']
        v = psf_weight(vc, self.flux, self.psf_sig)
        m2 = psf_weight(s2los + vc ** 2, self.flux, self.psf_sig)
        sig_obs = np.sqrt(np.clip(m2 - v ** 2, 0, None) + self.sig_th ** 2)
        return v - v[self.mask].mean(), sig_obs

    def __call__(self, field, rng, noise=True, kpp=None):
        """Mock VelocityMap of ``field`` (values outside the mask are NaN)."""
        v, s = self.observe(field, rng, kpp)
        v = v + self.bulk
        if noise:
            v = v + rng.normal(size=self.shape) * self.verr
            if s is not None and self.sigerr is not None:
                s = np.abs(s + rng.normal(size=self.shape) * self.sigerr)
        v = np.where(self.mask, v, np.nan)
        s = None if s is None else np.where(self.mask, s, np.nan)
        prov = [f'mock {field.mode}: sigma={field.sigma:g} alpha={field.alpha:g} L={field.L_kpc} '
                f'gamma={field.gamma:g} noise={noise}']
        if self.beta is not None:
            prov[0] += f' (lensed, source grid {self.src_kpp:.3g} kpc)'
        vm = VelocityMap(v, self.verr, self.mask, width=s, kpc_per_pix=self.kpp if kpp is None else kpp,
                         psf_fwhm_kpc=self.psf_sig * 2.3548 * self.kpp, line=self.line, z=self.z,
                         wcs=self.wcs, provenance=prov,
                         meta=dict(mock=dict(sigma=field.sigma, alpha=field.alpha, L_kpc=field.L_kpc,
                                             gamma=field.gamma, mode=field.mode, lensed=self.beta is not None)))
        if self.beta is not None:
            vm.set_source_positions(self.beta[0], self.beta[1], self.image)
        return vm

    def simulate(self, theta, rng, noise=True, kpp=None, **field_kw):
        """Mock VelocityMap for theta = (sigma, alpha, log10 L, gamma)."""
        return self(TurbulenceField.from_theta(theta, **field_kw), rng, noise, kpp)
