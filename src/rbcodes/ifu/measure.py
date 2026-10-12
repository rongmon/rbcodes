"""
Velocity measurements in IFU cubes: Gaussian fits, moments, Monte Carlo error checks, quality cuts.

    from rbcodes.ifu.lines import LineModel, lsf_muse
    from rbcodes.ifu.noise import NoiseModel
    from rbcodes.ifu.measure import GaussianFitter, QualityMask, SpaxelDecisions

    lm  = LineModel('OII', z=0.785273, windows=dict(flux=(6640, 6680), fit=(6640, 6680)),
                    air=False, lsf=lsf_muse)
    nm  = NoiseModel.from_cube(cube, linefree=lm, sky_box=(0, 5, 0, 8))
    fit = GaussianFitter(lm, vmax=600, sig_bounds=(1, 500)).fit_cube(cube, nm, detect_snr=3)
    qm  = QualityMask(snr_min=5, verr_max=30, chi2_max=3).apply(fit)
    fit.save('OII_fit.fits')

Coordinates are 0-based (iy, ix) = (row, column). Rejected spaxels are never filled or
interpolated. Per-spaxel decisions (``SpaxelDecisions``) update the fit maps without refitting
the cube and are stored in a JSON file so they are re-applied on every run.
"""
import json
import os
import time

import numpy as np
from scipy.ndimage import convolve, gaussian_filter, label
from scipy.optimize import curve_fit

from .lines import C_KMS, fit_domain

__all__ = ['FIT_KEYS', 'REJECT', 'STATUS', 'FitResult', 'GaussianFitter', 'MomentMeasurer',
           'MCErrors', 'QualityMask', 'QualityResult', 'SpaxelDecisions', 'detection_map',
           'spaxel_report', 'plot_spaxel']

FIT_KEYS = ['v', 'v_err', 'sig', 'sig_err', 'F', 'F_err', 'R', 'R_err', 'chi2', 'C', 'status']
STATUS = {0: 'ok', 1: 'bad data', 2: 'fit failed'}
# rejection codes for the quality mask
REJECT = {0: 'kept', 1: 'not attempted / bad data', 2: 'fit failed', 3: 'low S/N', 4: 'large σ_v',
          5: 'high χ²_ν', 6: 'at velocity bound', 7: 'isolated pixel', 8: 'small island',
          9: 'forced include', 10: 'forced exclude'}


def _cube_arrays(cube, noise, wave=None):
    """(flux, sigma, wave, wcs) from an IFUCube or arrays, and a NoiseModel or array."""
    if hasattr(cube, 'flux'):
        flux, wave, wcs = cube.flux, cube.wave, cube.wcs
    else:
        flux, wcs = np.asarray(cube, dtype=float), None
        if wave is None:
            raise ValueError('With an array cube, give wave=.')
    sigma = noise.sigma if hasattr(noise, 'sigma') else np.asarray(noise, dtype=float)
    wave = np.asarray(wave, dtype=float)
    if flux.shape != sigma.shape or flux.shape[0] != len(wave):
        raise ValueError(f'Shapes do not match: flux {flux.shape}, noise {sigma.shape}, wave {wave.shape}.')
    return flux, sigma, wave, wcs


def detection_map(model, flux, sigma, wave, detect_snr=3.):
    """Integrated flux over the flux window(s), its error, and the spaxels with S/N > detect_snr."""
    fl = model.masks(wave)[0]
    dl = np.median(np.diff(wave))
    F = np.nansum(flux[fl], 0) * dl
    Fe = np.sqrt(np.nansum(sigma[fl] ** 2, 0)) * dl
    with np.errstate(invalid='ignore', divide='ignore'):
        detect = np.isfinite(F) & (F / Fe > detect_snr)
    return F, Fe, detect


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

class FitResult:
    """
    Maps from a fitter or moment measurer.

    ``fit['v']``, ``fit['v_err']``, ``fit['sig']``, ``fit['sig_err']`` (km/s), ``fit['F']``,
    ``fit['F_err']``, ``fit['R']``, ``fit['R_err']``, ``fit['C']``, ``fit['chi2']``,
    ``fit['status']`` (NaN = not attempted, else STATUS), plus per-line keys (``F_<line>`` ...)
    for tied models. ``fit.detect`` = spaxels attempted by default; ``fit.meta`` = settings.
    """

    def __init__(self, maps, meta=None, detect=None, wcs=None):
        self.maps, self.meta, self.detect, self.wcs = dict(maps), dict(meta or {}), detect, wcs

    def __getitem__(self, k):
        return self.maps[k]

    def __setitem__(self, k, v):
        self.maps[k] = v

    def __contains__(self, k):
        return k in self.maps

    def keys(self):
        return self.maps.keys()

    @property
    def shape(self):
        return self.maps['v'].shape

    def __repr__(self):
        return (f"FitResult({self.meta.get('line', '?')}, {self.meta.get('method', '')}, shape={self.shape}, "
                f"fitted={int(np.isfinite(self.maps['v']).sum())})")

    def params(self, model, iy, ix):
        """Parameter tuple of spaxel (iy, ix) in ``model.param_names`` order."""
        return tuple(self.maps[k][iy, ix] for k in model.param_names)

    def save(self, path, overwrite=True):
        """Multi-extension FITS: one image per map (with the 2D WCS), settings as JSON in the header."""
        from astropy.io import fits
        ph = fits.Header()
        ph['ORIGIN'] = 'rbcodes.ifu.measure'
        ph['META'] = json.dumps(self.meta, default=_json_default)
        wh = self.wcs.to_header() if self.wcs is not None else fits.Header()
        hdus = [fits.PrimaryHDU(header=ph)]
        maps = dict(self.maps)
        if self.detect is not None:
            maps['detect'] = self.detect.astype(float)
        for k, m in maps.items():
            h = wh.copy()
            h['EXTNAME'] = k
            hdus.append(fits.ImageHDU(np.asarray(m, dtype=float), header=h))
        fits.HDUList(hdus).writeto(path, overwrite=overwrite)

    @classmethod
    def load(cls, path):
        from astropy.io import fits
        from astropy.wcs import WCS
        with fits.open(path, memmap=False) as hdul:
            meta = json.loads(hdul[0].header.get('META', '{}'))
            names = {k.upper(): k for k in meta.get('keys', FIT_KEYS)}
            maps = {names.get(h.name, h.name): h.data.astype(float) for h in hdul[1:]}
            wcs = WCS(hdul[1].header, naxis=2) if 'CTYPE1' in hdul[1].header else None
        detect = maps.pop('DETECT', maps.pop('detect', None))
        return cls(maps, meta, None if detect is None else detect.astype(bool), wcs)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, tuple):
        return list(o)
    return str(o)


# ---------------------------------------------------------------------------
# Gaussian fitting
# ---------------------------------------------------------------------------

class GaussianFitter:
    """
    Per-spaxel least-squares fit of a ``LineModel``.

    Parameters
    ----------
    model : LineModel
    vmax : float
        Velocity bounds +-vmax (km/s).
    sig_bounds : (lo, hi)
        Intrinsic dispersion bounds (km/s).
    inflate : bool
        Multiply the errors by sqrt(max(1, chi2_nu)).
    cont_bound : float
        Continuum bounds, +-cont_bound x median noise around ``cont_p0`` (the default 5 with
        cont_p0=0 assumes continuum-subtracted cubes).
    cont_p0 : float
        Starting continuum level.

    ``curve_fit(sigma=noise, absolute_sigma=True)``: the covariance follows the noise model.
    """

    def __init__(self, model, vmax=600., sig_bounds=(1., 500.), inflate=True, cont_bound=5., cont_p0=0.):
        self.model, self.vmax, self.sig_bounds = model, float(vmax), tuple(sig_bounds)
        self.inflate, self.cont_bound, self.cont_p0 = inflate, float(cont_bound), float(cont_p0)

    def settings(self):
        return dict(vmax=self.vmax, sig_bounds=self.sig_bounds, inflate=self.inflate,
                    cont_bound=self.cont_bound, cont_p0=self.cont_p0)

    def fit_spectrum(self, lam, f, e, vmax=None, sig_bounds=None, ratio_bounds=None, p0_v=None, p0_sig=None):
        """
        Fit one spectrum (lam, f, e) already cut to the fit domain.

        Returns dict with ``status`` (0 ok, 1 bad data, 2 fit failed) and, if 0: ``p``, ``perr``,
        ``chi2`` (reduced), ``inflate`` (factor applied to the errors), ``v``, ``v_err``, ``sig``,
        ``sig_err``, ``F``, ``F_err``, ``R``, ``R_err``, ``C`` (first line), per-line values for
        tied models, ``lo``, ``hi`` (bounds), ``nchan``.
        """
        m = self.model
        vmax = self.vmax if vmax is None else vmax
        sig_bounds = self.sig_bounds if sig_bounds is None else sig_bounds
        lam, f, e = (np.asarray(a, dtype=float) for a in (lam, f, e))
        out = dict(status=1)
        if len(lam) < m.npar + 2 or not (np.all(np.isfinite(f)) and np.all(np.isfinite(e)) and np.all(e > 0)):
            return out
        # fit in units of the median noise (the bounds and starting values do not depend on flux units)
        ce = np.median(e)
        fs, es = f / ce, e / ce
        dl = np.median(np.diff(lam))
        gi = m._group_index(lam)
        g0 = gi == 0
        if p0_v is None:
            fp = np.clip(fs[g0], 0, None)
            g = m.groups[0]
            lc = np.sum(fp * lam[g0]) / fp.sum() if fp.sum() > 0 else g['obs'].mean()
            wts = m.weights(g, 1.3) if g['kind'] == 'free' else g['w']
            p0_v = C_KMS * (lc / np.average(g['obs'], weights=wts) - 1)
        p0_v = float(np.clip(p0_v, -vmax + 1, vmax - 1))
        s0 = float(np.clip(p0_sig or 60., sig_bounds[0] * 1.01, sig_bounds[1] * 0.99))
        c0 = self.cont_p0 / ce
        p0, lo, hi = [], [], []
        for i, g in enumerate(m.groups):
            sel = gi == i
            F0 = max(np.sum(fs[sel]) * dl, 3 * dl)
            p0 += [c0, F0]
            lo += [c0 - self.cont_bound, 0.]
            hi += [c0 + self.cont_bound, 50 * F0 + 1]
            if g['kind'] == 'free':
                rb = ratio_bounds if (ratio_bounds is not None and m.single) else g['ratio_bounds']
                p0.append(float(np.clip(1.2, rb[0] * 1.01, rb[1] * 0.99)))
                lo.append(rb[0])
                hi.append(rb[1])
        p0 += [p0_v, s0]
        lo += [-vmax, sig_bounds[0]]
        hi += [vmax, sig_bounds[1]]
        try:
            p, cov = curve_fit(m, lam, fs, p0=p0, sigma=es, absolute_sigma=True, bounds=(lo, hi), maxfev=5000)
        except Exception:
            out['status'] = 2
            return out
        pe = np.sqrt(np.clip(np.diag(cov), 0, None))
        chi2 = np.sum(((fs - m(lam, *p)) / es) ** 2) / max(len(lam) - m.npar, 1)
        k = np.sqrt(max(1., chi2)) if self.inflate else 1.
        # back to the input flux units: continua and fluxes scale with ce
        scale = np.array([ce if n.startswith(('C', 'F')) else 1. for n in m.param_names])
        p, pe = p * scale, pe * scale
        val = dict(zip(m.param_names, p))
        err = dict(zip(m.param_names, pe * k))
        first = m.groups[0]['name']
        key = (lambda b: b) if m.single else (lambda b: f'{b}_{first}')
        out.update(status=0, p=p, perr=pe * k, chi2=chi2, inflate=k, v=val['v'], v_err=err['v'],
                   sig=val['sig'], sig_err=err['sig'], F=val[key('F')], F_err=err[key('F')],
                   R=val.get(key('R'), np.nan), R_err=err.get(key('R'), np.nan), C=val[key('C')],
                   lo=np.array(lo) * scale, hi=np.array(hi) * scale, nchan=len(lam))
        if not m.single:
            for n in m.param_names[:-2]:
                out[n], out[n + '_err'] = val[n], err[n]
        return out

    def _map_keys(self):
        extra = [] if self.model.single else [k for n in self.model.param_names[:-2] for k in (n, n + '_err')]
        return FIT_KEYS + extra

    def fit_cube(self, cube, noise, detect=None, detect_snr=3., overrides=None, decisions=None, wave=None,
                 verbose=True):
        """
        Fit every spaxel in ``detect`` (default: integrated S/N over the flux window > detect_snr).

        overrides = {(iy, ix): dict(fit=(lo, hi), exclude=[(lo, hi)], vmax=..., sig_bounds=(..),
        ratio_bounds=(..), p0_v=..., p0_sig=...)} changes the settings for individual spaxels (and
        forces a fit even if not detected). ``decisions`` (SpaxelDecisions) adds the stored refits.
        Returns a FitResult.
        """
        flux, sigma, wave, wcs = _cube_arrays(cube, noise, wave)
        m = self.model
        F_int, F_int_err, det = detection_map(m, flux, sigma, wave, detect_snr)
        detect = det if detect is None else np.asarray(detect, bool)
        overrides = dict(overrides or {})
        if decisions is not None:
            overrides.update(decisions.overrides(m.name))
        ny, nx = flux.shape[1:]
        keys = self._map_keys()
        out = {k: np.full((ny, nx), np.nan) for k in keys}
        base = m.fit_mask(wave)
        todo = set(zip(*np.where(detect))) | set(overrides)
        t0, nfail = time.time(), 0
        for iy, ix in sorted(todo):
            o = overrides.get((iy, ix), {})
            dom = m.fit_mask(wave, o) if ('fit' in o or 'exclude' in o) else base
            r = self.fit_spectrum(wave[dom], flux[dom, iy, ix], sigma[dom, iy, ix], o.get('vmax'),
                                  o.get('sig_bounds'), o.get('ratio_bounds'), o.get('p0_v'), o.get('p0_sig'))
            out['status'][iy, ix] = r['status']
            if r['status'] == 2:
                nfail += 1
            if r['status'] != 0:
                continue
            for k in keys:
                if k != 'status':
                    out[k][iy, ix] = r[k]
        if verbose:
            print(f'  {m.name}: fitted {np.isfinite(out["v"]).sum()} spaxels in {time.time() - t0:.1f} s '
                  f'({nfail} failed, {len(overrides)} with per-spaxel overrides)')
        out['F_int'], out['F_int_err'] = F_int, F_int_err
        meta = dict(m.meta(), method='gaussian', detect_snr=detect_snr, keys=keys + ['F_int', 'F_int_err'],
                    overrides={f'{iy},{ix}': o for (iy, ix), o in overrides.items()}, **self.settings())
        return FitResult(out, meta, detect, wcs)

    def refit(self, cube, noise, iy, ix, wave=None, **override):
        """Fit one spaxel with the fitter settings changed by ``override`` (keys as in fit_cube).
        The result (with ``override`` attached) can be passed to ``SpaxelDecisions.accept``."""
        flux, sigma, wave, _ = _cube_arrays(cube, noise, wave)
        dom = self.model.fit_mask(wave, override)
        r = self.fit_spectrum(wave[dom], flux[dom, iy, ix], sigma[dom, iy, ix], override.get('vmax'),
                              override.get('sig_bounds'), override.get('ratio_bounds'), override.get('p0_v'),
                              override.get('p0_sig'))
        r['override'] = dict(override)
        return r


# ---------------------------------------------------------------------------
# Moments
# ---------------------------------------------------------------------------

class MomentMeasurer:
    """
    Velocity moments in a fixed velocity window, with analytic errors.

    The components of the model's lines (or ``components``, rest wavelengths in the model's
    air/vacuum convention) are put on a common velocity grid relative to each component's
    observed wavelength lambda_rest (1+z), and stacked with inverse-variance weights. Then, in
    ``vwin`` (km/s):

        M0 = sum f dv,  v = sum f u dv / M0,  w^2 = sum f (u - v)^2 dv / M0

    with errors from the noise: sigma_v^2 = sum[(u - v) e dv]^2 / M0^2 and
    sigma(w^2)^2 = sum[((u - v)^2 - w^2) e dv]^2 / M0^2. The width w includes the instrumental
    width; it is NaN where noise makes w^2 < 0 (choose ``vwin`` close to the line extent: the
    noise in w^2 grows quickly with the window). ``F`` = M0 x lambda_1 / c (flux x Angstrom per
    component). Only for components that do not blend with each other (for a blended doublet
    give one reference wavelength). The stack weights by noise only, so it gains S/N for lines
    of similar strength; adding a much weaker line lowers it.

    Parameters
    ----------
    model : LineModel
    vwin : (lo, hi) km/s
    components : sequence of float, optional
        Rest wavelengths to stack (default: every component of the model).
    dv : float, optional
        Velocity grid step (km/s); default: the channel width at the first component.
    continuum : {'none', 'linefree'}
        'linefree' subtracts the per-spaxel median over the model's line-free channels.
    """

    def __init__(self, model, vwin, components=None, dv=None, continuum='none'):
        if continuum not in ('none', 'linefree'):
            raise ValueError("continuum must be 'none' or 'linefree'.")
        self.model, self.vwin, self.dv, self.continuum = model, tuple(vwin), dv, continuum
        rest = np.concatenate([g['rest'] for g in model.groups]) if components is None else np.asarray(components, float)
        self.obs = np.atleast_1d(rest) * (1 + model.z)

    def stack(self, wave, f, e, exclude=None):
        """(u, f_stack, e_stack) on the velocity grid; excluded channels are dropped first."""
        keep = np.isfinite(f) & np.isfinite(e) & (e > 0)
        for lo, hi in (exclude or []):
            keep &= ~((wave > lo) & (wave < hi))
        if len(self.obs) == 1:
            u = C_KMS * (wave / self.obs[0] - 1)
            return u[keep], f[keep], e[keep]
        dv = self.dv or C_KMS * np.median(np.diff(wave)) / self.obs[0]
        lo, hi = self.vwin
        u = np.arange(lo, hi + dv / 2, dv)
        num, wsum = np.zeros_like(u), np.zeros_like(u)
        for lk in self.obs:
            uk = C_KMS * (wave[keep] / lk - 1)
            if uk.min() > lo or uk.max() < hi:
                continue                                  # window not covered for this component
            fk = np.interp(u, uk, f[keep])
            ek = np.interp(u, uk, e[keep])
            w = 1 / ek ** 2
            num += w * fk
            wsum += w
        ok = wsum > 0
        return u[ok], num[ok] / wsum[ok], 1 / np.sqrt(wsum[ok])

    def moments(self, u, f, e):
        sel = (u >= self.vwin[0]) & (u <= self.vwin[1])
        u, f, e = u[sel], f[sel], e[sel]
        out = dict(status=1)
        if len(u) < 3:
            return out
        dv = np.gradient(u)
        m0 = np.sum(f * dv)
        out['F_vel'], out['F_vel_err'] = m0, np.sqrt(np.sum((e * dv) ** 2))
        if not m0 > 0:
            out['status'] = 2
            return out
        v = np.sum(f * u * dv) / m0
        w2 = np.sum(f * (u - v) ** 2 * dv) / m0
        v_err = np.sqrt(np.sum(((u - v) * e * dv) ** 2)) / m0
        w2_err = np.sqrt(np.sum((((u - v) ** 2 - w2) * e * dv) ** 2)) / m0
        w = np.sqrt(w2) if w2 > 0 else np.nan
        out.update(status=0, v=v, v_err=v_err, sig=w, sig_err=w2_err / (2 * w) if w2 > 0 else np.nan)
        return out

    def measure_cube(self, cube, noise, detect=None, detect_snr=3., wave=None):
        flux, sigma, wave, wcs = _cube_arrays(cube, noise, wave)
        m = self.model
        F_int, F_int_err, det = detection_map(m, flux, sigma, wave, detect_snr)
        detect = det if detect is None else np.asarray(detect, bool)
        if self.continuum == 'linefree':
            lf = m.masks(wave)[2]
            flux = flux - np.nanmedian(flux[lf], axis=0)[None]
        exclude = [r for g in m.groups for r in (g['W'].get('exclude') or [])]
        ny, nx = flux.shape[1:]
        keys = ['v', 'v_err', 'sig', 'sig_err', 'F', 'F_err', 'status']
        out = {k: np.full((ny, nx), np.nan) for k in keys}
        to_A = self.obs[0] / C_KMS
        for iy, ix in zip(*np.where(detect)):
            r = self.moments(*self.stack(wave, flux[:, iy, ix], sigma[:, iy, ix], exclude))
            out['status'][iy, ix] = r['status']
            if 'F_vel' in r:
                out['F'][iy, ix], out['F_err'][iy, ix] = r['F_vel'] * to_A, r['F_vel_err'] * to_A
            if r['status'] == 0:
                for k in ('v', 'v_err', 'sig', 'sig_err'):
                    out[k][iy, ix] = r[k]
        out['F_int'], out['F_int_err'] = F_int, F_int_err
        meta = dict(m.meta(), method='moments', vwin=self.vwin, components_obs=self.obs.tolist(),
                    continuum=self.continuum, detect_snr=detect_snr, keys=keys + ['F_int', 'F_int_err'],
                    width='includes instrumental width')
        return FitResult(out, meta, detect, wcs)


# ---------------------------------------------------------------------------
# Monte Carlo error validation
# ---------------------------------------------------------------------------

class MCErrors:
    """
    Monte Carlo check of the reported velocity errors.

    Noise is injected into the best-fit model of each chosen spaxel, which is refitted; the
    scatter of the refitted velocities is compared with the reported sigma_v. Noise can be
    correlated: white noise smoothed with a Gaussian of ``corr_spatial`` spaxels and
    ``corr_spectral`` channels (sigma), renormalised to unit variance, times the noise model.

    Parameters
    ----------
    fitter : GaussianFitter
    n_real : int
    seed : int
        Random seed (explicit).
    corr_spatial, corr_spectral : float
        Gaussian smoothing sigma of the injected noise (0 = white).
    """

    def __init__(self, fitter, n_real, seed, corr_spatial=0., corr_spectral=0.):
        self.fitter, self.n_real, self.seed = fitter, int(n_real), seed
        self.corr_spatial, self.corr_spectral = float(corr_spatial), float(corr_spectral)

    def _noise_field(self, rng, shape):
        z = rng.standard_normal(shape)
        if self.corr_spatial > 0 or self.corr_spectral > 0:
            z = gaussian_filter(z, (self.corr_spectral, self.corr_spatial, self.corr_spatial), mode='wrap')
            z /= z.std()
        return z

    def run(self, fit, cube, noise, spaxels=None, n_spaxels=50, wave=None):
        """
        Parameters
        ----------
        fit : FitResult from ``fitter.fit_cube``
        spaxels : list of (iy, ix), optional
            Default: ``n_spaxels`` drawn at random (with ``seed``) from the fitted spaxels.

        Returns
        -------
        dict(spaxels, v_std (scatter of refitted v), v_err (reported, original fit),
        v_err_mc (median reported in the realizations), ratio = v_std / v_err, median_ratio)
        """
        flux, sigma, wave, _ = _cube_arrays(cube, noise, wave)
        m = self.fitter.model
        rng = np.random.default_rng(self.seed)
        if spaxels is None:
            ys, xs = np.where(np.isfinite(fit['v']))
            pick = rng.choice(len(ys), size=min(n_spaxels, len(ys)), replace=False)
            spaxels = [(int(ys[i]), int(xs[i])) for i in pick]
        overrides = {tuple(map(int, k.split(','))): o for k, o in fit.meta.get('overrides', {}).items()}
        base = m.fit_mask(wave)
        doms = {s: (m.fit_mask(wave, overrides[s]) if s in overrides else base) for s in spaxels}
        models = {s: m(wave[doms[s]], *fit.params(m, *s)) for s in spaxels}
        v = np.full((len(spaxels), self.n_real), np.nan)
        ve = v.copy()
        for j in range(self.n_real):
            z = self._noise_field(rng, flux.shape)
            for i, s in enumerate(spaxels):
                d = doms[s]
                e = sigma[d, s[0], s[1]]
                o = overrides.get(s, {})
                r = self.fitter.fit_spectrum(wave[d], models[s] + z[d, s[0], s[1]] * e, e, o.get('vmax'),
                                             o.get('sig_bounds'), o.get('ratio_bounds'), o.get('p0_v'),
                                             o.get('p0_sig'))
                if r['status'] == 0:
                    v[i, j], ve[i, j] = r['v'], r['v_err']
        v_std = np.nanstd(v, axis=1, ddof=1)
        v_err = np.array([fit['v_err'][s] for s in spaxels])
        ratio = v_std / v_err
        return dict(spaxels=spaxels, v_std=v_std, v_err=v_err, v_err_mc=np.nanmedian(ve, axis=1),
                    ratio=ratio, median_ratio=float(np.nanmedian(ratio)), v=v)


# ---------------------------------------------------------------------------
# Quality mask and per-spaxel decisions
# ---------------------------------------------------------------------------

class QualityResult:
    """``mask`` (bool, kept spaxels) and ``code`` (int, keys of REJECT) from ``QualityMask.apply``."""

    def __init__(self, mask, code, settings):
        self.mask, self.code, self.settings = mask, code, settings

    def counts(self):
        return {REJECT[k]: int((self.code == k).sum()) for k in REJECT if (self.code == k).any()}

    def __repr__(self):
        return f'QualityResult(kept={int(self.mask.sum())}, {self.counts()})'


class QualityMask:
    """
    Quality cuts, in order: S/N = F/F_err > snr_min; v_err < verr_max; chi2_nu < chi2_max;
    |v| < 0.98 vmax; then spaxels with fewer than ``min_neigh`` kept neighbours (of 8) and
    connected groups smaller than ``min_island`` are removed. ``force_include`` spaxels with a
    successful fit are kept regardless; ``force_exclude`` spaxels are always removed (codes 9/10).
    Cuts on quantities a result does not have (e.g. chi2 for moments) are skipped.
    No hole filling.
    """

    def __init__(self, snr_min=5., verr_max=30., chi2_max=3., vmax=None, min_neigh=3, min_island=15,
                 force_include=(), force_exclude=()):
        self.snr_min, self.verr_max, self.chi2_max, self.vmax = snr_min, verr_max, chi2_max, vmax
        self.min_neigh, self.min_island = min_neigh, min_island
        self.force_include, self.force_exclude = list(force_include), list(force_exclude)

    def apply(self, fit, decisions=None):
        """Returns a QualityResult. ``decisions`` (SpaxelDecisions) adds stored includes/excludes."""
        v = fit['v']
        ny, nx = v.shape
        detect = fit.detect if fit.detect is not None else np.isfinite(fit['status'])
        vmax = self.vmax if self.vmax is not None else fit.meta.get('vmax', np.inf)
        code = np.full((ny, nx), 1, int)
        ok = np.isfinite(v)
        code[(detect | ok) & ~ok] = 2
        with np.errstate(invalid='ignore', divide='ignore'):
            snr = fit['F'] / fit['F_err']
        c = ok.copy()
        cuts = [(snr > self.snr_min, 3), (fit['v_err'] < self.verr_max, 4)]
        if 'chi2' in fit and np.isfinite(fit['chi2']).any():
            cuts.append((fit['chi2'] < self.chi2_max, 5))
        cuts.append((np.abs(v) < 0.98 * vmax, 6))
        for cond, cd in cuts:
            bad = c & ~np.nan_to_num(cond, nan=False).astype(bool)
            code[bad] = cd
            c &= ~bad
        k8 = np.ones((3, 3))
        k8[1, 1] = 0
        for _ in range(2):   # iterate: removing a pixel can isolate its neighbour
            iso = c & (convolve(c.astype(int), k8, mode='constant') < self.min_neigh)
            code[iso] = 7
            c &= ~iso
        lab, nl = label(c)
        for i in range(1, nl + 1):
            if (lab == i).sum() < self.min_island:
                code[lab == i] = 8
                c[lab == i] = False
        code[c] = 0
        inc, exc = list(self.force_include), list(self.force_exclude)
        if decisions is not None:
            inc += decisions.includes(fit.meta.get('line'))
            exc += decisions.excludes(fit.meta.get('line'))
        for iy, ix in inc:
            if ok[iy, ix]:
                c[iy, ix] = True
                code[iy, ix] = 9
            else:
                print(f'  force_include ({iy},{ix}) ignored: spaxel has no successful fit '
                      '(refit it and accept the result)')
        for iy, ix in exc:
            c[iy, ix] = False
            code[iy, ix] = 10
        settings = dict(snr_min=self.snr_min, verr_max=self.verr_max, chi2_max=self.chi2_max, vmax=vmax,
                        min_neigh=self.min_neigh, min_island=self.min_island)
        return QualityResult(c, code, settings)


class SpaxelDecisions:
    """
    Per-spaxel decisions stored in a JSON file and re-applied on every run:

    * ``accept(fit, iy, ix)`` keeps a spaxel with its current fit ('include');
    * ``accept(fit, iy, ix, result)`` writes a refit (``GaussianFitter.refit``) into the fit maps
      and keeps it ('refit'; the override is stored so ``fit_cube`` reproduces it);
    * ``reject(line, iy, ix)`` always removes it ('exclude');
    * ``forget(line, iy, ix)`` removes the decision.

    ``line`` is the model name (or a FitResult). Changing a decision does not refit the cube.
    """

    def __init__(self, path=None):
        self.path = None if path is None else os.fspath(path)
        self.data = {}
        if self.path and os.path.exists(self.path):
            with open(self.path) as fh:
                self.data = json.load(fh)

    @staticmethod
    def _line(line):
        return line.meta.get('line') if hasattr(line, 'meta') else line

    def _record(self, line, iy, ix, entry):
        d = self.data.setdefault(self._line(line), {})
        if entry is None:
            d.pop(f'{iy},{ix}', None)
        else:
            d[f'{iy},{ix}'] = entry
        self.save()

    def save(self):
        if self.path:
            with open(self.path, 'w') as fh:
                json.dump(self.data, fh, indent=1, default=_json_default)

    def accept(self, fit, iy, ix, result=None):
        if result is not None:
            if result.get('status') != 0:
                raise ValueError('that refit failed; nothing changed')
            for k in fit.keys():
                if k in result and k != 'status':
                    fit[k][iy, ix] = result[k]
            fit['status'][iy, ix] = 0
            entry = dict(action='refit', override=result.get('override', {}))
        else:
            if not np.isfinite(fit['v'][iy, ix]):
                raise ValueError('no fit for this spaxel: refit it and pass the result')
            entry = dict(action='include')
        self._record(fit, iy, ix, entry)

    def reject(self, line, iy, ix):
        self._record(line, iy, ix, dict(action='exclude'))

    def forget(self, line, iy, ix):
        self._record(line, iy, ix, None)

    def _items(self, line, actions):
        return {tuple(map(int, k.split(','))): d for k, d in self.data.get(self._line(line), {}).items()
                if d['action'] in actions}

    def overrides(self, line):
        return {s: {k: (tuple(v) if isinstance(v, list) and k != 'exclude' else v) for k, v in d['override'].items()}
                for s, d in self._items(line, ('refit',)).items()}

    def includes(self, line):
        return list(self._items(line, ('include', 'refit')))

    def excludes(self, line):
        return list(self._items(line, ('exclude',)))


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def spaxel_report(fit, code, iy, ix, snr_min=5., verr_max=30., chi2_max=3.):
    """One-line diagnosis of a spaxel: parameters, which cuts it passes/fails, rejection reason."""
    g = lambda k: fit[k][iy, ix]
    name = fit.meta.get('line', '')
    if not np.isfinite(g('v')):
        st = g('status')
        return f'({iy},{ix}) {name}: no fit (status {st:.0f}) | reason: {REJECT[int(code[iy, ix])]}'
    snr = g('F') / g('F_err')
    chk = lambda ok: 'ok' if ok else 'FAIL'
    chi = f'  χ²ν={g("chi2"):.2f} [{chk(g("chi2") < chi2_max)}]' if 'chi2' in fit else ''
    rr = f'  R={g("R"):.2f}' if 'R' in fit and np.isfinite(g('R')) else ''
    return (f'({iy},{ix}) {name}: v={g("v"):.1f}±{g("v_err"):.1f}  σ={g("sig"):.1f}±{g("sig_err"):.1f} km/s  '
            f'F={g("F"):.3g}  S/N={snr:.1f} [{chk(snr > snr_min)}]  σ_v [{chk(g("v_err") < verr_max)}]'
            f'{chi}{rr} | reason: {REJECT[int(code[iy, ix])]}')


def plot_spaxel(axes, model, wave, flux, noise, iy, ix, result=None, alt=None, alt_label='refit', title='',
                fit_window=None, exclude=None):
    """
    Data, noise band, best fit and components (``axes[0]``) and normalised residuals (``axes[1]``).
    ``result`` and ``alt`` are ``fit_spectrum``/``refit`` dicts (or None). Windows default to the model's.
    """
    a, b = axes
    sigma = noise.sigma if hasattr(noise, 'sigma') else noise
    g0 = model.groups[0]['W']
    fit_window = fit_window or (min(g['W']['fit'][0] for g in model.groups),
                                max(g['W']['fit'][1] for g in model.groups))
    exclude = exclude if exclude is not None else [r for g in model.groups for r in (g['W'].get('exclude') or [])]
    dom = model.fit_mask(wave, dict(fit=fit_window, exclude=exclude)) if model.single else model.fit_mask(wave)
    pad = (wave > fit_window[0] - 15) & (wave < fit_window[1] + 15)
    f, e = flux[:, iy, ix], sigma[:, iy, ix]
    a.step(wave[pad], f[pad], where='mid', color='k', lw=0.8, label='data')
    a.fill_between(wave[pad], -e[pad], e[pad], step='mid', color='0.85', label='±1σ noise')
    for g in model.groups:
        a.axvspan(*g['W']['fit'], color='C0', alpha=0.06)
    for lo, hi in exclude:
        a.axvspan(lo, hi, color='r', alpha=0.15)
    lf = np.linspace(wave[pad][0], wave[pad][-1], 800)
    for res, col, lab in [(result, 'C3', 'fit'), (alt, 'C2', alt_label)]:
        if res is None or res.get('status') != 0:
            continue
        a.plot(lf, model(lf, *res['p']), color=col, lw=1.4,
               label=f'{lab}: χ²ν={res["chi2"]:.2f}, v={res["v"]:.0f}±{res["v_err"]:.0f}')
        for comp in model.components(lf, *res['p']):
            a.plot(lf, res['p'][0] + comp, ':', color=col, lw=1)
        b.step(wave[dom], (f[dom] - model(wave[dom], *res['p'])) / e[dom], where='mid', color=col, lw=0.9)
    for lo_ in model.obs:
        a.axvline(lo_, color='C1', ls=':', lw=0.8)
    b.axhline(0, color='k', lw=0.6)
    b.fill_between(wave[pad], -1, 1, color='0.9')
    b.set_ylabel('(data−fit)/σ')
    b.set_xlabel('λ [Å]')
    a.legend(fontsize=7)
    a.set_title(title or f'{model.name} ({iy},{ix}) {g0["fit"]}', fontsize=8)
