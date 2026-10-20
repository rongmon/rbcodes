"""
Emission-line catalogue and line models for IFU fitting.

    from rbcodes.ifu.lines import LineModel, lsf_muse
    lm = LineModel('OII', z=0.785273, windows=dict(flux=(6640, 6680), fit=(6640, 6680)),
                   air=False, lsf=lsf_muse)

Model (all lines share one velocity v and one intrinsic dispersion s):

    f(lambda) = C_g + F_g * sum_k w_k G(lambda; lambda_k, sigma_k)     in the fit window of line g
    lambda_k  = lambda_rest,k (1 + z) (1 + v/c)
    sigma_k   = sqrt[(lambda_k s / c)^2 + sigma_LSF(lambda_k)^2]

G is a unit-area Gaussian, so F_g is the line flux. Within one line the component weights w_k
are fixed (``ratio='fixed'``), set by a free ratio R = F_2/F_1 within bounds (``ratio='free'``),
or 1 (``ratio='single'``). Several lines (e.g. ``['OII', 'Hbeta']``) are tied in v and s and
have independent fluxes and continua; each line has its own fit window and may have its own LSF.

Everything instrument- or physics-specific is an option: air or vacuum rest wavelengths
(``air``), the LSF, the doublet ratio bounds, the windows. Add lines by passing a
``catalogue`` with entries like those of ``LINES``.
"""
import numpy as np

__all__ = ['C_KMS', 'LINES', 'LineModel', 'lsf_muse', 'lsf_resolution', 'lsf_constant',
           'fit_domain', 'window_masks']

C_KMS = 299792.458

# Rest wavelengths: air and vacuum (Angstrom). 'ratio':
#   'free'   -> two components, F2/F1 = R free within ratio_bounds (R = flux of the 2nd / 1st entry)
#   'fixed'  -> components with fixed flux ratios (relative to the first)
#   'single' -> one Gaussian
# These are defaults; override ratio_bounds / ratios in LineModel, or pass your own catalogue.
LINES = {
    'OII':    dict(label='[O II] 3726,3729', air=(3726.032, 3728.815), vac=(3727.092, 3729.875),
                   ratio='free', ratio_bounds=(0.3, 1.6), amu=16.00),
    'OIII':   dict(label='[O III] 5007', air=(5006.843,), vac=(5008.240,), ratio='single', amu=16.00),
    'OIII_doublet': dict(label='[O III] 4959,5007', air=(5006.843, 4958.911), vac=(5008.240, 4960.295),
                         ratio='fixed', ratios=(1.0, 1 / 2.98), amu=16.00),     # Storey & Zeippen (2000)
    'Hbeta':  dict(label='Hbeta', air=(4861.333,), vac=(4862.683,), ratio='single', amu=1.008),
    'Halpha': dict(label='Halpha', air=(6562.819,), vac=(6564.614,), ratio='single', amu=1.008),
    'NII':    dict(label='[N II] 6584', air=(6583.454, 6548.050), vac=(6585.270, 6549.860),
                   ratio='fixed', ratios=(1.0, 1 / 3.05), amu=14.01),
    'MgII':   dict(label='Mg II 2796,2803', air=(2795.528, 2802.705), vac=(2796.354, 2803.531),
                   ratio='free', ratio_bounds=(1.0, 2.0), amu=24.31),
    # below 2000 A air and vacuum wavelengths are the same by convention; R = F(1909)/F(1907)
    # (low-density limit F(1907)/F(1909) ~ 1.5, i.e. R ~ 0.65; higher densities give larger R)
    'CIII':   dict(label='[C III] 1907, C III] 1909', air=(1906.683, 1908.734), vac=(1906.683, 1908.734),
                   ratio='free', ratio_bounds=(0.5, 10.0), amu=12.01),
}


# ---------------------------------------------------------------------------
# Line spread functions (Gaussian sigma in Angstrom as a function of observed wavelength)
# ---------------------------------------------------------------------------

def lsf_muse(lam):
    """MUSE LSF Gaussian sigma in Angstrom; FWHM(lambda) polynomial from Bacon et al. (2017, A&A 608, A1)."""
    return (5.866e-8 * lam ** 2 - 9.187e-4 * lam + 6.040) / 2.3548


def lsf_resolution(R):
    """LSF for a constant resolving power R = lambda / FWHM (e.g. KCWI)."""
    def lsf(lam):
        return np.asarray(lam, dtype=float) / (R * 2.3548)
    lsf.description = f'R={R:g}'
    return lsf


def lsf_constant(sigma):
    """LSF with a constant Gaussian sigma (Angstrom)."""
    def lsf(lam):
        return np.zeros_like(np.asarray(lam, dtype=float)) + sigma
    lsf.description = f'sigma={sigma:g} A'
    return lsf


def _as_lsf(lsf):
    if lsf is None:
        return lsf_constant(0.)
    if callable(lsf):
        return lsf
    if np.isscalar(lsf):
        return lsf_constant(float(lsf))
    raise ValueError('lsf must be None, a number (sigma in A) or a function of observed wavelength.')


# ---------------------------------------------------------------------------
# Windows (observed-frame Angstrom)
# ---------------------------------------------------------------------------

def fit_domain(wave, fit_window, exclude=None):
    """Boolean channel mask: inside fit_window=(lo, hi) Angstrom and outside every exclude=(lo, hi) range
    (use exclude for sky-line residuals or bad channels)."""
    wave = np.asarray(wave)
    m = (wave > fit_window[0]) & (wave < fit_window[1])
    for lo, hi in (exclude or []):
        m &= ~((wave > lo) & (wave < hi))
    return m


def window_masks(wave, W, margin=5.0):
    """Channel masks from a window dict W = dict(flux=(lo, hi), fit=(lo, hi), exclude=[...], linefree=[...] or None).
    Returns (flux_mask, fit_mask, linefree_mask). linefree=None -> channels more than ``margin`` A outside ``fit``."""
    wave = np.asarray(wave)
    flux = W.get('flux', W['fit'])
    fl = (wave > flux[0]) & (wave < flux[1])
    fit = fit_domain(wave, W['fit'], W.get('exclude'))
    if W.get('linefree'):
        lf = np.zeros(len(wave), bool)
        for lo, hi in W['linefree']:
            lf |= (wave > lo) & (wave < hi)
    else:
        lf = (wave < W['fit'][0] - margin) | (wave > W['fit'][1] + margin)
    return fl, fit, lf


# ---------------------------------------------------------------------------
# Line model
# ---------------------------------------------------------------------------

def _per_line(value, name, single):
    """Option given once for a single line, or as {line name: value} for several lines."""
    if isinstance(value, dict) and not single:
        return value.get(name)
    if isinstance(value, dict) and single and name in value:
        return value[name]
    return value


class LineModel:
    """
    One line (single or multiplet) or several lines tied in velocity and intrinsic dispersion.

    Parameters
    ----------
    names : str or list of str
        Keys of ``catalogue`` (default ``LINES``). A list fits the lines jointly (shared v, s;
        independent fluxes and continua).
    z : float
        Systemic redshift; velocities are relative to it.
    windows : dict
        Observed-frame Angstrom windows. For one line: dict(fit=(lo, hi), flux=(lo, hi),
        exclude=[(lo, hi), ...], linefree=[(lo, hi), ...] or None). For several lines:
        {name: dict(...)}. ``fit`` is required; ``flux`` defaults to ``fit``.
    air : bool
        Use air (True) or vacuum (False) rest wavelengths; must match the cube.
    lsf : None, float, callable, or {name: ...}
        Instrumental Gaussian sigma (Angstrom) vs observed wavelength: e.g. ``lsf_muse``,
        ``lsf_resolution(R)``, a constant, or None (then ``sig`` includes the instrumental width).
    ratio_bounds : (lo, hi) or {name: (lo, hi)}, optional
        Bounds of R = F_2/F_1 for 'free' doublets (default: the catalogue value).
    ratios : sequence or {name: sequence}, optional
        Fixed flux ratios for 'fixed' multiplets (default: the catalogue value).
    catalogue : dict, optional
        Line definitions (default ``LINES``); e.g. ``{**LINES, 'NeIII': dict(label='[Ne III] 3869',
        air=(3868.760,), vac=(3869.860,), ratio='single')}``.
    continuum : {'constant', 'linear'}
        Continuum of each line in its fit window: C (default), or C + S (lambda - lambda_c) with
        lambda_c the centre of the fit window (for cubes that are not continuum-subtracted).
    """

    def __init__(self, names, z, windows, air=False, lsf=None, ratio_bounds=None, ratios=None,
                 catalogue=None, continuum='constant'):
        if continuum not in ('constant', 'linear'):
            raise ValueError("continuum must be 'constant' or 'linear'.")
        self.continuum = continuum
        cat = LINES if catalogue is None else catalogue
        single = isinstance(names, str)
        names = [names] if single else list(names)
        self.names, self.z, self.air, self.single = names, float(z), bool(air), single
        self.name = names[0] if single else '+'.join(names)
        self.groups = []
        for nm in names:
            if nm not in cat:
                raise KeyError(f"Line '{nm}' is not in the catalogue ({sorted(cat)}).")
            cfg = cat[nm]
            W = windows if single and 'fit' in windows else windows.get(nm) if isinstance(windows, dict) else None
            if W is None or 'fit' not in W:
                raise ValueError(f"Give the observed-frame fit window for '{nm}': windows=dict(fit=(lo, hi), ...)"
                                 + ('' if single else f" as windows={{'{nm}': dict(fit=...)}}"))
            rest = np.array(cfg['air'] if air else cfg['vac'], dtype=float)
            kind = cfg['ratio']
            g = dict(name=nm, label=cfg.get('label', nm), cfg=cfg, rest=rest, obs=rest * (1 + self.z),
                     kind=kind, W=dict(W), lsf=_as_lsf(_per_line(lsf, nm, single)))
            if kind == 'free':
                g['ratio_bounds'] = tuple(_per_line(ratio_bounds, nm, single) or cfg['ratio_bounds'])
            elif kind == 'fixed':
                r = np.array(_per_line(ratios, nm, single) or cfg['ratios'], dtype=float)
                g['w'] = r / r.sum()
            else:
                g['w'] = np.array([1.0])
            self.groups.append(g)
        wins = sorted((g['W']['fit'][0], g['W']['fit'][1], g['name']) for g in self.groups)
        for a, b in zip(wins, wins[1:]):
            if b[0] < a[1]:
                raise ValueError(f"Fit windows of '{a[2]}' and '{b[2]}' overlap; use one multiplet entry instead.")

        # parameter layout: per line (C[, S], F[, R]), then v, s
        pn = []
        for g in self.groups:
            sfx = '' if single else f'_{g["name"]}'
            g['lc'] = 0.5 * (g['W']['fit'][0] + g['W']['fit'][1])
            pn += ([f'C{sfx}'] + ([f'S{sfx}'] if continuum == 'linear' else []) + [f'F{sfx}']
                   + ([f'R{sfx}'] if g['kind'] == 'free' else []))
        self.param_names = tuple(pn + ['v', 'sig'])
        self.npar = len(self.param_names)
        self.obs = np.concatenate([g['obs'] for g in self.groups])
        self.kind = self.groups[0]['kind'] if single else 'tied'

    def __repr__(self):
        return (f"LineModel({self.names if not self.single else repr(self.name)}, z={self.z}, air={self.air}, "
                f"params={self.param_names})")

    # --- parameters ---------------------------------------------------------

    def _split(self, p):
        """[(C, S, F, R or None) per line], v, s  (S = 0 for a constant continuum)."""
        out, i = [], 0
        lin = self.continuum == 'linear'
        for g in self.groups:
            C = p[i]
            S = p[i + 1] if lin else 0.
            i += 2 if lin else 1
            F = p[i]
            R = p[i + 1] if g['kind'] == 'free' else None
            i += 2 if g['kind'] == 'free' else 1
            out.append((C, S, F, R))
        return out, p[i], p[i + 1]

    @staticmethod
    def weights(g, R=None):
        if g['kind'] == 'free':
            return np.array([1 / (1 + R), R / (1 + R)])
        return g['w']

    def _gauss(self, g, lam, F, R, v, s):
        out = []
        for lk, wk in zip(g['obs'] * (1 + v / C_KMS), self.weights(g, R)):
            sk = np.sqrt((lk * s / C_KMS) ** 2 + g['lsf'](lk) ** 2)
            out.append(F * wk * np.exp(-0.5 * ((lam - lk) / sk) ** 2) / (np.sqrt(2 * np.pi) * sk))
        return out

    def _group_index(self, lam):
        """Line each wavelength belongs to (by fit window; nearest window outside them)."""
        if len(self.groups) == 1:
            return np.zeros(len(lam), int)
        centres = np.array([np.mean(g['W']['fit']) for g in self.groups])
        idx = np.argmin(np.abs(lam[:, None] - centres[None]), axis=1)
        for i, g in enumerate(self.groups):
            lo, hi = g['W']['fit']
            idx[(lam >= lo) & (lam <= hi)] = i
        return idx

    def __call__(self, lam, *p):
        lam = np.asarray(lam, dtype=float)
        per, v, s = self._split(p)
        out = np.zeros_like(lam)
        gi = self._group_index(lam)
        for i, (g, (C, S, F, R)) in enumerate(zip(self.groups, per)):
            m = gi == i
            if m.any():
                out[m] = C + S * (lam[m] - g['lc']) + np.sum(self._gauss(g, lam[m], F, R, v, s), axis=0)
        return out

    def continuum_at(self, lam, *p):
        """The continuum part of the model at ``lam``."""
        lam = np.asarray(lam, dtype=float)
        per, v, s = self._split(p)
        gi = self._group_index(lam)
        out = np.zeros_like(lam)
        for i, (g, (C, S, F, R)) in enumerate(zip(self.groups, per)):
            m = gi == i
            out[m] = C + S * (lam[m] - g['lc'])
        return out

    def components(self, lam, *p):
        """Individual Gaussian components (without the continua), as arrays over ``lam``."""
        lam = np.asarray(lam, dtype=float)
        per, v, s = self._split(p)
        out = []
        for g, (C, S, F, R) in zip(self.groups, per):
            out += self._gauss(g, lam, F, R, v, s)
        return out

    # --- windows ------------------------------------------------------------

    def masks(self, wave, margin=5.0):
        """(flux_mask, fit_mask, linefree_mask) over all lines (see ``window_masks``)."""
        fl = np.zeros(len(wave), bool)
        ft, lf = fl.copy(), fl.copy()
        for g in self.groups:
            a, b, c = window_masks(wave, g['W'], margin)
            fl |= a
            ft |= b
            lf |= c
        return fl, ft, lf & ~ft

    def fit_mask(self, wave, override=None):
        """Fit-domain mask; ``override`` may replace 'fit' / 'exclude' (single line) for one spaxel."""
        m = np.zeros(len(wave), bool)
        for g in self.groups:
            W = dict(g['W'])
            if override and self.single:
                W.update({k: override[k] for k in ('fit', 'exclude') if k in override})
            m |= fit_domain(wave, W['fit'], W.get('exclude'))
        return m

    def meta(self):
        """Settings for file headers."""
        return dict(line=self.name, z=self.z, air=self.air, continuum=self.continuum,
                    rest=[float(r) for g in self.groups for r in g['rest']],
                    windows={g['name']: g['W'] for g in self.groups},
                    lsf={g['name']: getattr(g['lsf'], 'description', getattr(g['lsf'], '__name__', 'custom'))
                         for g in self.groups},
                    ratio_bounds={g['name']: g.get('ratio_bounds') for g in self.groups})
