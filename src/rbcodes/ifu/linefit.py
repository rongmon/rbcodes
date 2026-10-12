"""
Emission-line fitting of IFU cubes: velocity, dispersion and flux maps.

    from rbcodes.ifu import load_fits
    from rbcodes.ifu.linefit import EmissionLine, fit_cube, quality_mask

    cube = load_fits('cube.fits')
    line = EmissionLine('OII', z=0.7853, air=True, lsf='muse')
    res = fit_cube(cube, line)                   # noise from the variance extension
    keep, code = quality_mask(res, snr_min=5, verr_max=30)
    res.save('OII_fit.fits')

Model (all components share one velocity v and one intrinsic dispersion s):

    f(lambda) = C + F * sum_k w_k G(lambda; lambda_k, sigma_k)
    lambda_k  = lambda_rest,k (1 + z) (1 + v/c)
    sigma_k   = sqrt[(lambda_k s / c)^2 + sigma_LSF(lambda_k)^2]

G is a unit-area Gaussian, so F is the total line flux (flux density x Angstrom). The weights
w_k are fixed (``kind='fixed'``), set by a free ratio R = F_2/F_1 (``kind='free'``, two
components), or 1 (``kind='single'``). Because the LSF is in the model, s is the intrinsic
dispersion. Velocities are relative to the redshift z.

Rest wavelengths of the presets come from the rbcodes line lists (vacuum); ``air=True``
converts them to air for cubes in air wavelengths (e.g. MUSE).
"""
import time

import numpy as np
from astropy import constants as const
from scipy.ndimage import convolve, label
from scipy.optimize import curve_fit

__all__ = ['C_KMS', 'vac_to_air', 'lsf_muse', 'lsf_resolution', 'lsf_constant',
           'LINE_PRESETS', 'EmissionLine', 'fit_domain', 'fit_spectrum', 'detection_snr',
           'fit_cube', 'LineFitResult', 'FIT_KEYS', 'STATUS', 'REJECT', 'quality_mask',
           'spaxel_report', 'plot_spaxel']

C_KMS = const.c.to('km/s').value
FWHM_PER_SIGMA = 2. * np.sqrt(2. * np.log(2.))


# ---------------------------------------------------------------------------
# Wavelengths and line spread functions
# ---------------------------------------------------------------------------

def vac_to_air(wave):
    """Vacuum -> air wavelength in Angstrom (IAU standard, Morton 2000, ApJS 130, 403)."""
    wave = np.asarray(wave, dtype=float)
    s2 = (1e4 / wave) ** 2
    n = 1. + 8.34254e-5 + 2.406147e-2 / (130. - s2) + 1.5998e-4 / (38.9 - s2)
    return wave / n


def lsf_muse(lam):
    """MUSE line spread function, Gaussian sigma in Angstrom (FWHM polynomial of Bacon et al. 2017)."""
    return (5.866e-8 * lam ** 2 - 9.187e-4 * lam + 6.040) / FWHM_PER_SIGMA


def lsf_resolution(R):
    """LSF for a constant resolving power R = lambda / FWHM (e.g. KCWI)."""
    def lsf(lam):
        return np.asarray(lam, dtype=float) / (R * FWHM_PER_SIGMA)
    lsf.description = f'R={R:g}'
    return lsf


def lsf_constant(sigma):
    """LSF with a constant Gaussian sigma in Angstrom."""
    def lsf(lam):
        return np.full_like(np.asarray(lam, dtype=float), sigma)
    lsf.description = f'sigma={sigma:g}A'
    return lsf


def _resolve_lsf(lsf):
    if lsf is None:
        return lsf_constant(0.), 'none'
    if isinstance(lsf, str):
        if lsf.lower() == 'muse':
            return lsf_muse, 'muse'
        raise ValueError(f"Unknown lsf '{lsf}'. Use 'muse', a number (sigma in A), or a function.")
    if np.isscalar(lsf):
        f = lsf_constant(float(lsf))
        return f, f.description
    if callable(lsf):
        return lsf, getattr(lsf, 'description', 'custom')
    raise ValueError("lsf must be None, 'muse', a number (sigma in A) or a function of wavelength.")


# ---------------------------------------------------------------------------
# Lines
# ---------------------------------------------------------------------------

# names are entries of the rbcodes 'Gal_Em' line list, or of 'atom' for MgII.
# For 'free' doublets R = F(second) / F(first).
LINE_PRESETS = {
    'OII': dict(names=('[OII] 3726', '[OII] 3729'), kind='free', ratio_bounds=(0.3, 1.6),
                label='[O II] 3726,3729'),
    'OIII': dict(names=('[OIII] 5007',), kind='single', label='[O III] 5007'),
    'OIII_doublet': dict(names=('[OIII] 5007', '[OIII] 4959'), kind='fixed', ratios=(1., 1 / 2.98),
                         label='[O III] 4959,5007'),          # Storey & Zeippen (2000)
    'Hbeta': dict(names=('Hbeta',), kind='single', label='Hbeta'),
    'Halpha': dict(names=('Halpha',), kind='single', label='Halpha'),
    'NII': dict(names=('[NII] 6584', '[NII] 6548'), kind='fixed', ratios=(1., 1 / 3.05),
                label='[N II] 6548,6584'),
    'MgII': dict(names=('MgII 2796', 'MgII 2803'), kind='free', ratio_bounds=(0.45, 1.05),
                 label='Mg II 2796,2803'),                   # 0.5 optically thin ... 1 thick
}


def _rest_vacuum(name):
    from rbcodes.IGM.rb_setline import read_line_list
    for listname in ('Gal_Em', 'atom'):
        for entry in read_line_list(listname):
            if str(entry['ion']).strip() == name:
                return float(entry['wrest'])
    raise KeyError(f"Line '{name}' not found in the rbcodes 'Gal_Em' or 'atom' line lists.")


class EmissionLine:
    """
    One or more Gaussian components tied in velocity and intrinsic dispersion.

    Parameters
    ----------
    name : str
        A key of LINE_PRESETS ('OII', 'OIII', 'OIII_doublet', 'Hbeta', 'Halpha', 'NII', 'MgII'),
        or any label when ``rest`` is given.
    z : float
        Systemic redshift; fitted velocities are relative to it.
    air : bool
        Convert the (vacuum) rest wavelengths to air. True for MUSE cubes.
    lsf : None, 'muse', float or callable
        Instrumental sigma (Angstrom) as a function of observed wavelength: 'muse', a constant
        sigma, ``lsf_resolution(R)``, or any function. None = no instrumental broadening.
    rest : sequence of float, optional
        Rest wavelengths (vacuum unless ``air`` is False and you give air values) for a custom line.
    kind : {'single', 'fixed', 'free'}, optional
        Default: 'single' for one wavelength, 'fixed' with ``ratios``, else 'free'.
    ratios : sequence of float, optional
        Fixed flux ratios relative to the first component (kind='fixed').
    ratio_bounds : (lo, hi), optional
        Bounds of R = F_2/F_1 (kind='free').

    Examples
    --------
    >>> EmissionLine('Halpha', z=0.1, lsf=lsf_resolution(1800))
    >>> EmissionLine('CIII', z=2.0, rest=(1906.68, 1908.73), ratio_bounds=(0.5, 2.0))
    """

    def __init__(self, name, z=0., air=False, lsf=None, rest=None, kind=None, ratios=None,
                 ratio_bounds=None, label=None):
        if rest is None:
            if name not in LINE_PRESETS:
                raise KeyError(f"Unknown line '{name}'. Presets: {sorted(LINE_PRESETS)}; "
                               "or give rest=(...).")
            p = LINE_PRESETS[name]
            rest_vac = np.array([_rest_vacuum(n) for n in p['names']])
            rest = vac_to_air(rest_vac) if air else rest_vac
            kind = kind or p['kind']
            ratios = ratios if ratios is not None else p.get('ratios')
            ratio_bounds = ratio_bounds or p.get('ratio_bounds')
            label = label or p['label']
        rest = np.atleast_1d(np.asarray(rest, dtype=float))
        if kind is None:
            kind = 'single' if len(rest) == 1 else ('fixed' if ratios is not None else 'free')
        if kind == 'single' and len(rest) != 1:
            raise ValueError("kind='single' needs one rest wavelength.")
        if kind == 'free':
            if len(rest) != 2:
                raise ValueError("kind='free' needs two rest wavelengths.")
            ratio_bounds = tuple(ratio_bounds or (0.1, 10.))
        if kind == 'fixed':
            if ratios is None or len(ratios) != len(rest):
                raise ValueError("kind='fixed' needs one ratio per rest wavelength.")
        if kind not in ('single', 'fixed', 'free'):
            raise ValueError(f"Unknown kind '{kind}'.")

        self.name, self.z, self.air, self.kind = name, float(z), bool(air), kind
        self.label = label or name
        self.rest = rest
        self.obs = rest * (1. + self.z)
        self.ratio_bounds = ratio_bounds
        self.lsf, self.lsf_name = _resolve_lsf(lsf)
        if kind == 'fixed':
            r = np.asarray(ratios, dtype=float)
            self._w = r / r.sum()
        elif kind == 'single':
            self._w = np.array([1.])
        self.param_names = ('C', 'F', 'R', 'v', 'sig') if kind == 'free' else ('C', 'F', 'v', 'sig')
        self.npar = len(self.param_names)

    def __repr__(self):
        return (f"EmissionLine('{self.name}', z={self.z}, rest={np.round(self.rest, 3).tolist()}, "
                f"kind='{self.kind}', air={self.air}, lsf='{self.lsf_name}')")

    def weights(self, R=None):
        """Fraction of the flux in each component."""
        if self.kind == 'free':
            return np.array([1. / (1. + R), R / (1. + R)])
        return self._w

    def _unpack(self, p):
        if self.kind == 'free':
            C, F, R, v, s = p
        else:
            (C, F, v, s), R = p, None
        return C, F, R, v, s

    def components(self, lam, *p):
        """List of the Gaussian components (without the constant C)."""
        C, F, R, v, s = self._unpack(p)
        lam = np.asarray(lam, dtype=float)
        out = []
        for lk, wk in zip(self.obs * (1. + v / C_KMS), self.weights(R)):
            sk = np.sqrt((lk * s / C_KMS) ** 2 + self.lsf(lk) ** 2)
            out.append(F * wk * np.exp(-0.5 * ((lam - lk) / sk) ** 2) / (np.sqrt(2. * np.pi) * sk))
        return out

    def __call__(self, lam, *p):
        return p[0] + np.sum(self.components(lam, *p), axis=0)

    def default_window(self, vwin=900.):
        """(lo, hi) in observed Angstrom: +-vwin km/s beyond the outermost component."""
        return (self.obs.min() * (1. - vwin / C_KMS), self.obs.max() * (1. + vwin / C_KMS))


# ---------------------------------------------------------------------------
# Fitting one spectrum
# ---------------------------------------------------------------------------

STATUS = {0: 'ok', 1: 'bad data', 2: 'fit failed'}


def fit_domain(wave, window, exclude=None):
    """Channels inside window=(lo, hi) and outside every exclude=(lo, hi) range (Angstrom)."""
    wave = np.asarray(wave)
    m = (wave > window[0]) & (wave < window[1])
    for lo, hi in (exclude or []):
        m &= ~((wave > lo) & (wave < hi))
    return m


def fit_spectrum(line, wave, flux, err, vmax=600., sig_bounds=(1., 500.), ratio_bounds=None,
                 p0_v=None, p0_sig=None, inflate=True):
    """
    Fit ``line`` to one spectrum (already cut to the fitting domain).

    Parameters
    ----------
    line : EmissionLine
    wave, flux, err : 1-D arrays (err = 1-sigma noise)
    vmax : float
        Velocity bounds +-vmax km/s.
    sig_bounds : (lo, hi)
        Intrinsic dispersion bounds in km/s.
    ratio_bounds : (lo, hi), optional
        Overrides line.ratio_bounds for 'free' doublets.
    p0_v, p0_sig : float, optional
        Starting velocity and dispersion (default: flux-weighted centroid, 60 km/s).
    inflate : bool
        Multiply the errors by sqrt(max(1, chi2_nu)).

    Returns
    -------
    dict with ``status`` (see STATUS) and, if status is 0: ``p`` (parameters in
    ``line.param_names`` order), ``perr``, ``v``, ``v_err``, ``sig``, ``sig_err``, ``F``, ``F_err``,
    ``R``, ``R_err`` (NaN unless kind='free'), ``C``, ``C_err``, ``chi2`` (reduced), ``inflate``
    (factor applied to the errors), ``nchan``.
    """
    wave, flux, err = (np.asarray(a, dtype=float) for a in (wave, flux, err))
    out = dict(status=1)
    if (len(wave) < line.npar + 2 or not np.all(np.isfinite(flux))
            or not np.all(np.isfinite(err)) or not np.all(err > 0)):
        return out

    # work in units of the median noise so the fit does not depend on the flux units
    scale = np.median(err)
    f, e = flux / scale, err / scale
    dl = np.median(np.diff(wave))
    c0 = np.percentile(f, 25)
    F0 = max(np.sum(f - c0) * dl, 3. * dl)
    if p0_v is None:
        fp = np.clip(f - c0, 0, None)
        lc = np.sum(fp * wave) / fp.sum() if fp.sum() > 0 else np.mean(line.obs)
        w = line.weights(1.) if line.kind == 'free' else line.weights()
        p0_v = C_KMS * (lc / np.average(line.obs, weights=w) - 1.)
    p0_v = float(np.clip(p0_v, -vmax * 0.99, vmax * 0.99))
    s0 = float(np.clip(p0_sig or 60., sig_bounds[0] * 1.01, sig_bounds[1] * 0.99))
    c_lo, c_hi = c0 - 5. - abs(c0), c0 + 5. + abs(c0)
    F_hi = 50. * F0 + 100. * dl
    if line.kind == 'free':
        rb = ratio_bounds or line.ratio_bounds
        p0 = [c0, F0, float(np.clip(1., rb[0] * 1.01, rb[1] * 0.99)), p0_v, s0]
        lo = [c_lo, 0., rb[0], -vmax, sig_bounds[0]]
        hi = [c_hi, F_hi, rb[1], vmax, sig_bounds[1]]
    else:
        p0 = [c0, F0, p0_v, s0]
        lo = [c_lo, 0., -vmax, sig_bounds[0]]
        hi = [c_hi, F_hi, vmax, sig_bounds[1]]
    try:
        p, cov = curve_fit(line, wave, f, p0=p0, sigma=e, absolute_sigma=True,
                           bounds=(lo, hi), maxfev=5000)
    except (RuntimeError, ValueError, np.linalg.LinAlgError):
        out['status'] = 2
        return out
    perr = np.sqrt(np.clip(np.diag(cov), 0, None))
    chi2 = np.sum(((f - line(wave, *p)) / e) ** 2) / max(len(wave) - line.npar, 1)
    k = np.sqrt(max(1., chi2)) if inflate else 1.

    # back to the input flux units (C and F scale, R/v/sig do not)
    p, perr = p.copy(), perr.copy()
    p[:2] *= scale
    perr[:2] *= scale
    vals = dict(zip(line.param_names, p))
    errs = dict(zip(line.param_names, perr * k))
    out.update(status=0, p=p, perr=perr * k, chi2=chi2, inflate=k, nchan=len(wave),
               v=vals['v'], v_err=errs['v'], sig=vals['sig'], sig_err=errs['sig'],
               F=vals['F'], F_err=errs['F'], C=vals['C'], C_err=errs['C'],
               R=vals.get('R', np.nan), R_err=errs.get('R', np.nan))
    return out


# ---------------------------------------------------------------------------
# Fitting a cube
# ---------------------------------------------------------------------------

FIT_KEYS = ('v', 'v_err', 'sig', 'sig_err', 'F', 'F_err', 'R', 'R_err', 'C', 'C_err', 'chi2',
            'status', 'det_snr')


def detection_snr(flux, noise, wave, window):
    """Integrated S/N map: sum of the flux over ``window`` divided by the summed-in-quadrature noise."""
    m = fit_domain(wave, window)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.nansum(flux[m], axis=0) / np.sqrt(np.nansum(noise[m] ** 2, axis=0))


def fit_cube(cube, line, noise=None, wave=None, wcs=None, window=None, exclude=None,
             flux_window=None, detect=None, detect_snr=3., vmax=600., sig_bounds=(1., 500.),
             ratio_bounds=None, inflate=True, overrides=None, verbose=True):
    """
    Fit ``line`` to every detected spaxel of a cube.

    Parameters
    ----------
    cube : IFUCube (from ``rbcodes.ifu.load_fits``) or ndarray (n_wave, ny, nx)
    line : EmissionLine
    noise : ndarray like the cube, optional
        1-sigma noise. Default: sqrt of the cube's variance (see ``rbcodes.ifu.noise``).
    wave, wcs : needed only when ``cube`` is an array (wcs is optional, stored with the maps).
    window : (lo, hi) Angstrom, optional
        Fitting range. Default: ``line.default_window(max(900, vmax + 300))``.
    exclude : list of (lo, hi), optional
        Ranges left out of the fit (sky-line residuals, bad channels).
    flux_window : (lo, hi), optional
        Range summed for the detection S/N map (default: ``window``).
    detect : bool (ny, nx), optional
        Spaxels to fit. Default: detection S/N > ``detect_snr``. The detection S/N counts
        continuum too, so for cubes that are not continuum-subtracted pass ``detect`` (or a
        ``flux_window`` and a high ``detect_snr``).
    vmax, sig_bounds, ratio_bounds, inflate : see ``fit_spectrum``.
    overrides : dict, optional
        {(iy, ix): dict(window=(lo, hi), exclude=[...], vmax=..., sig_bounds=(..),
        ratio_bounds=(..), p0_v=..., p0_sig=...)}: settings for single spaxels. These spaxels are
        fitted even if not detected.

    Returns
    -------
    LineFitResult
    """
    if hasattr(cube, 'flux'):
        flux, wave, wcs = cube.flux, cube.wave, cube.wcs if wcs is None else wcs
        if noise is None:
            if cube.var is None:
                raise ValueError('The cube has no variance; give noise= (see rbcodes.ifu.noise).')
            from .noise import noise_cube
            noise, _ = noise_cube(flux, var=cube.var)
    else:
        flux = np.asarray(cube, dtype=float)
        if wave is None or noise is None:
            raise ValueError('With an array cube, give wave= and noise=.')
    wave = np.asarray(wave, dtype=float)
    if flux.ndim != 3 or flux.shape != noise.shape or flux.shape[0] != len(wave):
        raise ValueError(f'Shapes do not match: flux {flux.shape}, noise {np.shape(noise)}, '
                         f'wave {wave.shape}.')

    if window is None:
        window = line.default_window(max(900., vmax + 300.))
    window = tuple(float(w) for w in window)
    flux_window = tuple(flux_window) if flux_window is not None else window
    overrides = overrides or {}
    ny, nx = flux.shape[1:]
    maps = {k: np.full((ny, nx), np.nan) for k in FIT_KEYS}
    maps['det_snr'] = detection_snr(flux, noise, wave, flux_window)
    if detect is None:
        detect = np.nan_to_num(maps['det_snr'], nan=-np.inf) > detect_snr

    base = fit_domain(wave, window, exclude)
    todo = sorted(set(zip(*np.nonzero(detect))) | set(overrides))
    t0, nfail = time.time(), 0
    for iy, ix in todo:
        o = overrides.get((iy, ix), {})
        dom = (fit_domain(wave, o.get('window', window), o.get('exclude', exclude))
               if ('window' in o or 'exclude' in o) else base)
        r = fit_spectrum(line, wave[dom], flux[dom, iy, ix], noise[dom, iy, ix],
                         vmax=o.get('vmax', vmax), sig_bounds=o.get('sig_bounds', sig_bounds),
                         ratio_bounds=o.get('ratio_bounds', ratio_bounds),
                         p0_v=o.get('p0_v'), p0_sig=o.get('p0_sig'), inflate=inflate)
        maps['status'][iy, ix] = r['status']
        nfail += r['status'] == 2
        if r['status'] == 0:
            for k in FIT_KEYS[:-2]:
                maps[k][iy, ix] = r[k]
    if verbose:
        print(f'[linefit] {line.name}: fitted {int(np.isfinite(maps["v"]).sum())} of {len(todo)} '
              f'spaxels in {time.time() - t0:.1f} s ({nfail} failed, {len(overrides)} overrides)')

    meta = dict(line=line.name, label=line.label, z=line.z, air=line.air, kind=line.kind,
                rest=tuple(line.rest), lsf=line.lsf_name, window=window, exclude=list(exclude or []),
                vmax=vmax, sig_bounds=tuple(sig_bounds), inflate=inflate, detect_snr=detect_snr)
    return LineFitResult(maps, meta, wcs)


class LineFitResult:
    """
    Maps from ``fit_cube``.

    ``res['v']``, ``res['v_err']``, ``res['sig']``, ``res['sig_err']`` (km/s), ``res['F']``,
    ``res['F_err']`` (line flux), ``res['R']``, ``res['R_err']`` (doublet ratio), ``res['C']``,
    ``res['C_err']`` (continuum), ``res['chi2']`` (reduced), ``res['status']`` (NaN = not
    fitted, else STATUS), ``res['det_snr']`` (detection S/N). ``res.snr`` = F/F_err.
    ``res.meta`` holds the line and fit settings, ``res.wcs`` the 2D WCS (or None).
    """

    def __init__(self, maps, meta=None, wcs=None):
        self.maps, self.meta, self.wcs = dict(maps), dict(meta or {}), wcs

    def __getitem__(self, key):
        return self.maps[key]

    def keys(self):
        return self.maps.keys()

    @property
    def shape(self):
        return self.maps['v'].shape

    @property
    def snr(self):
        with np.errstate(invalid='ignore', divide='ignore'):
            return self.maps['F'] / self.maps['F_err']

    def params(self, iy, ix, line):
        """Best-fit parameter tuple of spaxel (iy, ix), for ``line(wave, *params)``."""
        return tuple(self.maps[k][iy, ix] for k in line.param_names)

    def __repr__(self):
        n = int(np.isfinite(self.maps['v']).sum())
        return f"LineFitResult({self.meta.get('line', '?')}, shape={self.shape}, fitted={n})"

    # --- FITS I/O ---------------------------------------------------------

    _META_CARDS = dict(line='LINE', label='LABEL', z='Z', air='AIR', kind='KIND', lsf='LSF',
                       vmax='VMAX', inflate='INFLATE', detect_snr='DETSNR')

    def save(self, path, overwrite=False):
        """Write all maps to a multi-extension FITS file (one image per map, 2D WCS on each)."""
        from astropy.io import fits
        ph = fits.Header()
        ph['ORIGIN'] = 'rbcodes.ifu.linefit'
        for key, card in self._META_CARDS.items():
            if key in self.meta:
                ph[card] = self.meta[key]
        for i, r in enumerate(self.meta.get('rest', ())):
            ph[f'REST{i + 1}'] = (float(r), 'rest wavelength [A]')
        if 'window' in self.meta:
            ph['WIN_LO'], ph['WIN_HI'] = self.meta['window']
        if 'sig_bounds' in self.meta:
            ph['SIG_LO'], ph['SIG_HI'] = self.meta['sig_bounds']
        for i, (lo, hi) in enumerate(self.meta.get('exclude', [])):
            ph[f'EXCL{i + 1}'] = f'{lo:.3f},{hi:.3f}'
        wh = self.wcs.to_header() if self.wcs is not None else fits.Header()
        hdus = [fits.PrimaryHDU(header=ph)]
        for k, m in self.maps.items():
            h = wh.copy()
            h['EXTNAME'] = k.upper()
            hdus.append(fits.ImageHDU(np.asarray(m, dtype=float), header=h))
        fits.HDUList(hdus).writeto(path, overwrite=overwrite)

    @classmethod
    def load(cls, path):
        """Read a file written by ``save``."""
        from astropy.io import fits
        from astropy.wcs import WCS
        maps, wcs = {}, None
        with fits.open(path, memmap=False) as hdul:
            ph = hdul[0].header
            meta = {key: ph[card] for key, card in cls._META_CARDS.items() if card in ph}
            meta['rest'] = tuple(ph[f'REST{i}'] for i in range(1, 100) if f'REST{i}' in ph)
            if 'WIN_LO' in ph:
                meta['window'] = (ph['WIN_LO'], ph['WIN_HI'])
            if 'SIG_LO' in ph:
                meta['sig_bounds'] = (ph['SIG_LO'], ph['SIG_HI'])
            meta['exclude'] = [tuple(float(x) for x in ph[f'EXCL{i}'].split(','))
                               for i in range(1, 100) if f'EXCL{i}' in ph]
            names = {k.upper(): k for k in FIT_KEYS}
            for h in hdul[1:]:
                maps[names.get(h.name, h.name.lower())] = h.data.astype(float)
            if 'CTYPE1' in hdul[1].header:
                wcs = WCS(hdul[1].header, naxis=2)
        return cls(maps, meta, wcs)


# ---------------------------------------------------------------------------
# Quality cuts and diagnostics
# ---------------------------------------------------------------------------

REJECT = {0: 'kept', 1: 'not fitted', 2: 'bad data', 3: 'fit failed', 4: 'low S/N',
          5: 'large velocity error', 6: 'high chi2', 7: 'at velocity bound', 8: 'isolated pixel',
          9: 'small island', 10: 'forced include', 11: 'forced exclude'}


def quality_mask(res, snr_min=5., verr_max=30., chi2_max=3., vmax=None, min_neigh=3,
                 min_island=15, force_include=(), force_exclude=()):
    """
    Spaxels to keep, and the first reason each other spaxel was rejected.

    Cuts, in order: S/N = F/F_err > snr_min; v_err < verr_max; chi2_nu < chi2_max;
    |v| < 0.98 vmax (fits stuck at the bound; vmax defaults to the one used in the fit);
    then spaxels with fewer than ``min_neigh`` kept neighbours (of 8) and connected groups
    smaller than ``min_island`` are removed. ``force_include`` spaxels (iy, ix) with a
    successful fit are kept regardless; ``force_exclude`` spaxels are always removed.

    Returns
    -------
    keep : bool (ny, nx)
    code : int (ny, nx), keys of REJECT
    """
    v, status = res['v'], res['status']
    vmax = vmax if vmax is not None else res.meta.get('vmax', np.inf)
    code = np.ones(v.shape, dtype=int)
    code[status == 1] = 2
    code[status == 2] = 3
    ok = np.isfinite(v)
    c = ok.copy()
    with np.errstate(invalid='ignore', divide='ignore'):
        cuts = [(res.snr > snr_min, 4), (res['v_err'] < verr_max, 5), (res['chi2'] < chi2_max, 6),
                (np.abs(v) < 0.98 * vmax, 7)]
    for cond, cd in cuts:
        bad = c & ~cond
        code[bad] = cd
        c &= ~bad
    k8 = np.ones((3, 3))
    k8[1, 1] = 0
    for _ in range(2):             # removing a pixel can isolate its neighbour
        iso = c & (convolve(c.astype(int), k8, mode='constant') < min_neigh)
        code[iso] = 8
        c &= ~iso
    lab, nlab = label(c)
    sizes = np.bincount(lab.ravel())
    small = (sizes < min_island)[lab] & c
    code[small] = 9
    c &= ~small
    code[c] = 0
    for iy, ix in force_include:
        if ok[iy, ix]:
            c[iy, ix] = True
            code[iy, ix] = 10
        else:
            print(f'[linefit] force_include ({iy},{ix}) ignored: no successful fit '
                  '(use overrides in fit_cube to fit it)')
    for iy, ix in force_exclude:
        c[iy, ix] = False
        code[iy, ix] = 11
    return c, code


def spaxel_report(res, code, iy, ix, snr_min=5., verr_max=30., chi2_max=3.):
    """One-line summary of a spaxel: parameters, which cuts it passes, rejection reason."""
    g = lambda k: res[k][iy, ix]
    name = res.meta.get('line', '')
    reason = REJECT[int(code[iy, ix])]
    if not np.isfinite(g('v')):
        st = g('status')
        st = 'not fitted' if not np.isfinite(st) else STATUS[int(st)]
        return f'({iy},{ix}) {name}: no fit ({st}) | reason: {reason}'
    snr = res.snr[iy, ix]
    chk = lambda ok: 'ok' if ok else 'FAIL'
    ratio = f'  R={g("R"):.2f}' if np.isfinite(g('R')) else ''
    return (f'({iy},{ix}) {name}: v={g("v"):.1f}+-{g("v_err"):.1f}  sig={g("sig"):.1f}+-'
            f'{g("sig_err"):.1f} km/s  F={g("F"):.3g}  S/N={snr:.1f} [{chk(snr > snr_min)}]  '
            f'v_err [{chk(g("v_err") < verr_max)}]  chi2={g("chi2"):.2f} [{chk(g("chi2") < chi2_max)}]'
            f'{ratio} | reason: {reason}')


def plot_spaxel(axes, line, wave, flux, noise, iy, ix, res=None, window=None, exclude=None,
                alt=None, alt_label='refit', title=''):
    """
    Data, noise band, fit and components (``axes[0]``) and residuals / noise (``axes[1]``).

    ``res`` is a LineFitResult (the fit of this spaxel is drawn); ``alt`` an optional
    ``fit_spectrum`` result to overlay. ``window``/``exclude`` default to those of ``res``.
    """
    a, b = axes
    meta = res.meta if res is not None else {}
    window = window or meta.get('window') or line.default_window()
    exclude = exclude if exclude is not None else meta.get('exclude')
    dom = fit_domain(wave, window, exclude)
    pad = (wave > window[0] - 15) & (wave < window[1] + 15)
    f, e = flux[:, iy, ix], noise[:, iy, ix]
    a.step(wave[pad], f[pad], where='mid', color='k', lw=0.8, label='data')
    a.fill_between(wave[pad], -e[pad], e[pad], step='mid', color='0.85', label='noise')
    a.axvspan(*window, color='C0', alpha=0.06)
    for lo, hi in (exclude or []):
        a.axvspan(lo, hi, color='r', alpha=0.15)
    fits_to_draw = []
    if res is not None and np.isfinite(res['v'][iy, ix]):
        fits_to_draw.append((res.params(iy, ix, line), res['chi2'][iy, ix], res['v'][iy, ix],
                             res['v_err'][iy, ix], 'C3', 'fit'))
    if alt is not None and alt.get('status') == 0:
        fits_to_draw.append((alt['p'], alt['chi2'], alt['v'], alt['v_err'], 'C2', alt_label))
    lf = np.linspace(wave[pad][0], wave[pad][-1], 800)
    for p, chi2, v, ve, col, lab in fits_to_draw:
        a.plot(lf, line(lf, *p), color=col, lw=1.4, label=f'{lab}: chi2={chi2:.2f}, v={v:.0f}+-{ve:.0f}')
        if len(line.rest) > 1:
            for comp in line.components(lf, *p):
                a.plot(lf, p[0] + comp, ':', color=col, lw=1)
        b.step(wave[dom], (f[dom] - line(wave[dom], *p)) / e[dom], where='mid', color=col, lw=0.9)
    for lo in line.obs:
        a.axvline(lo, color='C1', ls=':', lw=0.8)
    b.axhline(0, color='k', lw=0.6)
    b.fill_between(wave[pad], -1, 1, color='0.9')
    b.set_ylabel('(data-fit)/noise')
    b.set_xlabel('Wavelength [A]')
    a.legend(fontsize=7)
    a.set_title(title or f'{line.label} ({iy},{ix})', fontsize=8)
