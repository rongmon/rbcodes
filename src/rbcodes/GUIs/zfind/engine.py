"""
engine.py — Computation engine for rb_zfind.

Contains three search modes:
    line_search()      — chi2 scan using rbcodes linelists (emission or absorption)
    template_search()  — chi2 scan against 1D FITS template spectra
    pca_search()       — chi2 scan using PCA eigenvectors (redrock-style, no redrock needed)

IMPORTANT: This module must NEVER import PyQt5, matplotlib, or any GUI library.
All functions accept rb_spectrum objects and return ZFindResult or AbsorberResult.
"""

import numpy as np
import os
import warnings
from scipy.stats import median_abs_deviation
from astropy import units as u

from rbcodes.GUIs.zfind.io import (
    ZFindResult, ZSolution, AbsorberResult, AbsorberCandidate
)

# Path to bundled templates
TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates')


# ---------------------------------------------------------------------------
# Shared preprocessing
# ---------------------------------------------------------------------------

def _preprocess(spec, fit_continuum=True):
    """
    Normalize rb_spectrum to plain arrays ready for engine use.

    Parameters
    ----------
    spec : rb_spectrum
    fit_continuum : bool
        If True and no continuum in spec, run rb_iter_contfit.
        If False, continuum = zeros (fine for emission mode).

    Returns
    -------
    wave : ndarray  (Angstrom, vacuum)
    flux : ndarray
    ivar : ndarray  (inverse variance, 0 for bad pixels)
    continuum : ndarray
    warnings_list : list of str
    """
    warn_list = []

    # --- air → vacuum ---
    if spec.meta.get('airvac', 'vac') == 'air':
        spec = spec.copy()
        spec.air2vac()
        warn_list.append("Air wavelengths converted to vacuum.")

    wave = spec.wavelength.to(u.AA).value
    flux = spec.flux.value

    # --- IVAR ---
    if spec.sig_is_set:
        sig = spec.sig.value
        ivar = np.where(sig > 0, 1.0 / sig**2, 0.0)
    else:
        sigma = 1.4826 * median_abs_deviation(flux, nan_policy='omit')
        if sigma == 0:
            sigma = np.std(flux[np.isfinite(flux)]) or 1.0
        ivar = np.full_like(flux, 1.0 / sigma**2)
        warn_list.append(
            f"No error array — using MAD-STD IVAR (sigma={sigma:.4g})."
        )

    # set ivar=0 for non-finite pixels
    ivar[~np.isfinite(flux)] = 0.0
    flux = np.nan_to_num(flux, nan=0.0)

    # --- continuum ---
    if spec.co_is_set:
        continuum = spec.co.value
    elif fit_continuum:
        try:
            from rbcodes.IGM.rb_iter_contfit import fit_optimal_polynomial
            err_arr = spec.sig.value if spec.sig_is_set else None
            result = fit_optimal_polynomial(
                wave, flux, error=err_arr,
                use_weights=False, silent=True, plot=False
            )
            continuum = result['continuum']
        except Exception as e:
            continuum = np.zeros_like(flux)
            warn_list.append(f"Continuum fit failed ({e}); using zero continuum.")
    else:
        continuum = np.zeros_like(flux)

    return wave, flux, ivar, continuum, warn_list


def _mad_ivar(flux):
    """Uniform IVAR from MAD-STD. Used as a standalone fallback."""
    sigma = 1.4826 * median_abs_deviation(flux, nan_policy='omit')
    if sigma == 0:
        sigma = np.std(flux[np.isfinite(flux)]) or 1.0
    return np.full_like(flux, 1.0 / sigma**2)


def _find_top_minima(z_array, chi2_array, n=3, min_dz=0.02):
    """
    Find top N global minima in chi2_array with minimum separation min_dz.

    Iteratively picks the global minimum, then masks out ±min_dz around it,
    and repeats. More robust than argrelmin for broad or plateau-like minima.

    Returns list of indices into z_array.
    """
    if not np.any(np.isfinite(chi2_array)):
        return []

    working = chi2_array.copy()
    selected = []

    while len(selected) < n:
        if not np.any(np.isfinite(working)):
            break
        idx = np.nanargmin(working)
        if not np.isfinite(working[idx]):
            break
        selected.append(int(idx))
        # mask out the neighbourhood so the next iteration picks a different z
        mask = np.abs(z_array - z_array[idx]) <= min_dz
        working[mask] = np.nan

    return selected


def _z_err_from_curvature(z_array, chi2_array, idx):
    """Estimate z_err from chi2 curvature at index idx."""
    try:
        d2 = np.gradient(np.gradient(chi2_array, z_array), z_array)
        curv = abs(d2[idx])
        if curv > 0:
            return float(np.sqrt(1.0 / curv))
    except Exception:
        pass
    return float('nan')


# ---------------------------------------------------------------------------
# Mode 1: Line search
# ---------------------------------------------------------------------------

def line_search(spec, linelist_df, z_min=0.0, z_max=6.0, n_steps=10000,
                mode='emission', fwhm_ang=0.0, fit_continuum=True,
                window_pixels=5, smooth_pixels=3,
                data_norm='subtract',
                wave_min=None, wave_max=None):
    """
    Scan redshift by matching emission or absorption lines to a spectrum.

    Parameters
    ----------
    spec : rb_spectrum
    linelist_df : DataFrame with columns ['wave', 'name']  (rest wavelengths in Å)
    z_min, z_max : float
    n_steps : int
    mode : 'emission' or 'absorption'
    fwhm_ang : float  — LSF FWHM in Å (0 = no convolution)
    fit_continuum : bool  — run rb_iter_contfit if no continuum in spec
    window_pixels : int  — half-window around each line for chi2 evaluation
    smooth_pixels : int  — boxcar kernel width; 1 = off
    data_norm : str  — how to prepare the observed spectrum before chi2:
                       'subtract'  flux − continuum  (default)
                       'normalize' flux / continuum − 1
                       'raw'       raw flux (no continuum removal)
    wave_min, wave_max : float or None — observed-frame Å clipping;
                         pixels outside the range are masked (ivar→0)

    Returns
    -------
    ZFindResult  (mode='emission')
    AbsorberResult  (mode='absorption')
    """
    wave, flux, ivar, continuum, warn_list = _preprocess(spec, fit_continuum=fit_continuum)

    # Wavelength mask — zero out pixels outside the requested range
    if wave_min is not None:
        ivar[wave < wave_min] = 0.0
    if wave_max is not None:
        ivar[wave > wave_max] = 0.0

    # Pre-smooth flux (and continuum) to suppress pixel noise before the scan.
    # ivar is left unchanged — the statistic is a detection SNR, not a formal chi2.
    if smooth_pixels > 1:
        try:
            from astropy.convolution import convolve, Box1DKernel
            _kernel = Box1DKernel(smooth_pixels)
            flux      = convolve(flux,      _kernel)
            continuum = convolve(continuum, _kernel)
        except Exception:
            pass  # silently skip smoothing if astropy not available

    # Prepare observed spectrum according to data_norm
    cont_safe = np.where(np.abs(continuum) > 1e-10, continuum, 1.0)
    if data_norm == 'normalize':
        flux_work     = flux / cont_safe - 1.0   # fractional deviation
        baseline_work = np.zeros_like(flux_work)
    elif data_norm == 'raw':
        flux_work     = flux
        baseline_work = np.zeros_like(flux_work)
    else:   # 'subtract' (default)
        flux_work     = flux - continuum
        baseline_work = np.zeros_like(flux_work)

    rest_waves = linelist_df['wave'].values.astype(float)
    line_names = linelist_df['name'].values

    z_array = np.linspace(z_min, z_max, n_steps)
    chi2_array = np.full(n_steps, np.nan)

    wmin, wmax = wave.min(), wave.max()
    dwave = np.median(np.diff(wave))  # approximate pixel scale

    for i, z in enumerate(z_array):
        obs_waves = rest_waves * (1.0 + z)
        in_range = (obs_waves > wmin + window_pixels * dwave) & \
                   (obs_waves < wmax - window_pixels * dwave)

        if not np.any(in_range):
            continue

        obs_in = obs_waves[in_range]
        accum_chi2 = 0.0
        n_used = 0

        for obs_w in obs_in:
            # pixel indices for window around this line
            pix = np.argmin(np.abs(wave - obs_w))
            lo = max(0, pix - window_pixels)
            hi = min(len(wave), pix + window_pixels + 1)

            fw = flux_work[lo:hi]   # already baseline-subtracted / normalized
            iv = ivar[lo:hi]

            if np.sum(iv > 0) < 2:
                continue

            iv_sum = np.sum(iv)
            if iv_sum <= 0:
                continue

            if mode == 'emission':
                S = np.sum(fw * iv)   # positive when flux > baseline
                # Sign check: only count when feature is genuinely in emission.
                if S > 0:
                    accum_chi2 -= (S ** 2) / iv_sum
                    n_used += 1

            else:  # absorption
                S = np.sum(-fw * iv)  # positive when flux < baseline (absorbed)
                # Sign check: only count when feature is genuinely in absorption.
                if S > 0:
                        accum_chi2 -= (S ** 2) / iv_sum
                        n_used += 1

        if n_used > 0:
            # Normalise by sqrt(n_used) → SNR-like statistic.
            # Dividing by n_used (mean) would dilute real detections with
            # undetected lines; sqrt gives the right noise scaling for
            # independent line contributions.
            chi2_array[i] = accum_chi2 / np.sqrt(n_used)

    # --- extract solutions ---
    top_idx = _find_top_minima(z_array, chi2_array, n=10)

    method_label = f'LineSearch:{linelist_df.attrs.get("name", "custom")}'

    if mode == 'emission':
        solutions = []
        for idx in top_idx:
            # count lines in range at this z
            obs = rest_waves * (1.0 + z_array[idx])
            n_lines = int(np.sum((obs > wmin) & (obs < wmax)))
            solutions.append(ZSolution(
                z=float(z_array[idx]),
                z_err=_z_err_from_curvature(z_array, chi2_array, idx),
                chi2_dof=float(chi2_array[idx]),
                method=method_label,
                template_type='Unknown',
                n_features=n_lines
            ))

        return ZFindResult(
            z_array=z_array,
            chi2_curves=[{'label': method_label, 'chi2': chi2_array.copy()}],
            solutions=solutions,
            input_spec=spec,
            warnings=warn_list
        )

    else:  # absorption
        # convert chi2 to significance (invert and normalize)
        finite = chi2_array[np.isfinite(chi2_array)]
        if len(finite) == 0:
            candidates = []
        else:
            noise = np.std(finite)
            baseline = np.median(finite)
            significance_curve = (baseline - chi2_array) / (noise + 1e-10)

            top_idx_abs = _find_top_minima(z_array, chi2_array, n=20)
            candidates = []
            for idx in top_idx_abs:
                obs = rest_waves * (1.0 + z_array[idx])
                matched_mask = (obs > wmin) & (obs < wmax)
                matched_names = line_names[matched_mask].tolist()
                sig_val = float(significance_curve[idx]) if np.isfinite(significance_curve[idx]) else 0.0

                candidates.append(AbsorberCandidate(
                    z=float(z_array[idx]),
                    significance=sig_val,
                    n_lines=int(np.sum(matched_mask)),
                    is_doublet=False,       # TODO: doublet detection in future
                    linelist_name=linelist_df.attrs.get('name', 'custom'),
                    lines_matched=matched_names
                ))

            candidates.sort(key=lambda c: c.significance, reverse=True)

        return AbsorberResult(
            z_array=z_array,
            significance_curve=significance_curve if len(finite) > 0 else np.zeros_like(z_array),
            candidates=candidates,
            input_spec=spec,
            warnings=warn_list
        )


# ---------------------------------------------------------------------------
# Mode 2: Template search
# ---------------------------------------------------------------------------

# Mapping: user-facing name → FITS filename under templates/marz/
TEMPLATE_FILES = {
    'EarlyType':        'EarlyType.fits',
    'Intermediate':     'Intermediate.fits',
    'LateTypeEmission': 'LateTypeEmission.fits',
    'Composite':        'Composite.fits',
    'QSO':              'QSO.fits',
    'HighZSFG':         'HighZSFG.fits',
}
TEMPLATE_NAMES = list(TEMPLATE_FILES.keys())

# Sensible z ranges per template type (z_min, z_max)
_TEMPLATE_Z_RANGE = {
    'EarlyType':        (0.0, 1.5),
    'Intermediate':     (0.0, 1.5),
    'LateTypeEmission': (0.0, 2.0),
    'Composite':        (0.0, 1.5),
    'QSO':              (0.0, 5.5),
    'HighZSFG':         (1.5, 6.0),
}


def _load_template(template_name: str):
    """
    Load a bundled MARZ template FITS file.

    Returns (wave_Å, flux) as plain ndarrays.  flux is median-normalised so
    amplitude differences between templates don't dominate the chi2.
    """
    fname = TEMPLATE_FILES.get(template_name)
    if fname is None:
        raise ValueError(
            f"Unknown template '{template_name}'. "
            f"Available: {TEMPLATE_NAMES}"
        )
    path = os.path.join(TEMPLATE_DIR, 'marz', fname)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Template file not found: {path}\n"
            f"Run: python {TEMPLATE_DIR}/download_marz_templates.py"
        )
    from astropy.io import fits as pyfits
    with pyfits.open(path) as hdul:
        wave = hdul['WAVE'].data.astype(float)
        flux = hdul['FLUX'].data.astype(float)

    # Compute a smooth pseudo-continuum of the template (wide boxcar).
    # Used when continuum_mode='normalize': template_norm = flux / pseudo_cont.
    smooth_width = max(51, len(flux) // 15)
    from astropy.convolution import convolve, Box1DKernel
    pseudo_cont = convolve(flux, Box1DKernel(smooth_width))
    # Avoid divide-by-zero; fall back to median where pseudo_cont is tiny
    median_abs = np.median(np.abs(flux[flux != 0])) if np.any(flux != 0) else 1.0
    pseudo_cont = np.where(np.abs(pseudo_cont) > 0.01 * median_abs,
                           pseudo_cont, median_abs)

    return wave, flux, pseudo_cont


def template_search(spec, template_name='LateTypeEmission',
                    z_min=None, z_max=None, n_steps=5000,
                    fit_continuum=True,
                    data_norm='normalize', model_norm='normalize',
                    smooth_pixels=1,
                    wave_min=None, wave_max=None):
    """
    Chi-square scan of *spec* against a bundled MARZ template.

    At each trial redshift z the template is interpolated onto the spectrum
    wavelength grid, scaled to the optimal amplitude, and the residual chi2
    is computed:

        A_opt   = Σ(F · T · ivar) / Σ(T² · ivar)
        chi2(z) = Σ((F - A·T)² · ivar) / N_overlap

    continuum_mode controls how F and T are prepared before comparison:

        'normalize'  (default) — F = flux / continuum,  T = t_flux / t_pseudo_cont
                                  Both are dimensionless ~1 + fractional features.
                                  Best for MARZ templates which carry a continuum.
        'subtract'             — F = flux - continuum,  T = t_flux - t_pseudo_cont
                                  Both in flux units relative to their continua.
        'raw'                  — F = flux,  T = t_flux (no continuum processing).
                                  Only useful when flux calibration matches template.

    Parameters
    ----------
    spec             : rb_spectrum
    template_name    : str   — one of TEMPLATE_NAMES
    z_min, z_max     : float or None  — defaults from _TEMPLATE_Z_RANGE
    n_steps          : int
    fit_continuum : bool  — fit a polynomial continuum if spec has none
    data_norm     : str  — how to prepare the observed spectrum:
                           'normalize'  flux / continuum  (default)
                           'subtract'   flux − continuum
                           'raw'        flux as-is
    model_norm    : str  — how to prepare the template:
                           'normalize'  t_flux / t_pseudo_cont  (default)
                           'subtract'   t_flux − t_pseudo_cont
                           'raw'        t_flux as-is
    smooth_pixels : int  — boxcar width applied to flux_work before scan; 1=off
    wave_min, wave_max : float or None — observed-frame Å mask (ivar→0 outside)

    Returns
    -------
    ZFindResult
    """
    from rbcodes.GUIs.zfind.io import ZFindResult, ZSolution

    wave, flux, ivar, continuum, warn_list = _preprocess(
        spec, fit_continuum=fit_continuum
    )

    # Wavelength mask
    if wave_min is not None:
        ivar[wave < wave_min] = 0.0
    if wave_max is not None:
        ivar[wave > wave_max] = 0.0

    t_wave, t_flux, t_pseudo_cont = _load_template(template_name)

    # --- prepare observed spectrum ---
    cont_safe = np.where(np.abs(continuum) > 1e-10, continuum, 1.0)
    if data_norm == 'normalize':
        flux_work = flux / cont_safe
    elif data_norm == 'subtract':
        flux_work = flux - continuum
    else:   # 'raw'
        flux_work = flux

    # --- prepare template ---
    if model_norm == 'normalize':
        t_flux_work = t_flux / t_pseudo_cont
    elif model_norm == 'subtract':
        t_flux_work = t_flux - t_pseudo_cont
    else:
        t_flux_work = t_flux

    # Optional pre-smoothing of the observed flux_work
    if smooth_pixels > 1:
        try:
            from astropy.convolution import convolve, Box1DKernel
            flux_work = convolve(flux_work, Box1DKernel(smooth_pixels))
        except Exception:
            pass

    # z range defaults
    z_lo, z_hi = _TEMPLATE_Z_RANGE.get(template_name, (0.0, 6.0))
    if z_min is None:
        z_min = z_lo
    if z_max is None:
        z_max = z_hi

    z_array   = np.linspace(z_min, z_max, n_steps)
    chi2_array = np.full(n_steps, np.nan)

    wmin, wmax = wave.min(), wave.max()

    for i, z in enumerate(z_array):
        # Redshift template into observed frame
        t_wave_obs = t_wave * (1.0 + z)

        # Check overlap: template must cover at least 20 % of spectrum range
        overlap_lo = max(wmin, t_wave_obs.min())
        overlap_hi = min(wmax, t_wave_obs.max())
        if (overlap_hi - overlap_lo) < 0.2 * (wmax - wmin):
            continue

        # Interpolate template onto spectrum wavelength grid
        T = np.interp(wave, t_wave_obs, t_flux_work, left=0.0, right=0.0)

        # Valid pixels: ivar > 0 and template non-zero
        good = (ivar > 0) & (T != 0.0)
        if good.sum() < 50:
            continue

        F  = flux_work[good]
        Tg = T[good]
        IV = ivar[good]

        # Optimal amplitude (least-squares)
        denom = np.sum(Tg ** 2 * IV)
        if denom == 0:
            continue
        A = np.sum(F * Tg * IV) / denom

        # Reduced chi2
        resid = F - A * Tg
        chi2_array[i] = np.sum(resid ** 2 * IV) / good.sum()

    # --- extract solutions ---
    top_idx = _find_top_minima(z_array, chi2_array, n=10, min_dz=0.02)

    method_label = f'Template:{template_name}'
    solutions = []
    for idx in top_idx:
        solutions.append(ZSolution(
            z=float(z_array[idx]),
            z_err=_z_err_from_curvature(z_array, chi2_array, idx),
            chi2_dof=float(chi2_array[idx]),
            method=method_label,
            template_type=template_name,
            n_features=0,
        ))

    return ZFindResult(
        z_array=z_array,
        chi2_curves=[{'label': method_label, 'chi2': chi2_array.copy()}],
        solutions=solutions,
        input_spec=spec,
        warnings=warn_list,
    )


def multi_template_search(spec, templates=None, z_min=None, z_max=None,
                           n_steps=5000, fit_continuum=True,
                           data_norm='normalize', model_norm='normalize',
                           smooth_pixels=1, wave_min=None, wave_max=None):
    """
    Run template_search for every template in *templates* (default: all),
    combine chi2 curves onto a common z grid, return the best-fitting
    solution set and all per-template curves for plotting.

    Returns ZFindResult where:
      chi2_curves  — one entry per template (for plotting overlay)
      solutions    — drawn from whichever template gave the lowest chi2/dof
    """
    from rbcodes.GUIs.zfind.io import ZFindResult

    if templates is None:
        templates = TEMPLATE_NAMES

    all_results = {}
    warn_all = []
    for tname in templates:
        try:
            r = template_search(spec, template_name=tname,
                                z_min=z_min, z_max=z_max,
                                n_steps=n_steps,
                                fit_continuum=fit_continuum,
                                data_norm=data_norm, model_norm=model_norm,
                                smooth_pixels=smooth_pixels,
                                wave_min=wave_min, wave_max=wave_max)
            all_results[tname] = r
            warn_all.extend(r.warnings)
        except Exception as e:
            warn_all.append(f'{tname} failed: {e}')

    if not all_results:
        raise RuntimeError('All templates failed.')

    # Common z grid from first result
    first = next(iter(all_results.values()))
    z_array = first.z_array

    # Build combined chi2_curves list
    chi2_curves = []
    for tname, r in all_results.items():
        chi2_curves.append({'label': tname, 'chi2': r.chi2_curves[0]['chi2']})

    # Best solution = lowest chi2/dof across all templates
    best_tname = min(
        all_results,
        key=lambda t: (all_results[t].solutions[0].chi2_dof
                       if all_results[t].solutions else np.inf)
    )
    solutions = all_results[best_tname].solutions

    return ZFindResult(
        z_array=z_array,
        chi2_curves=chi2_curves,
        solutions=solutions,
        input_spec=spec,
        warnings=list(dict.fromkeys(warn_all)),   # deduplicate
    )


# ---------------------------------------------------------------------------
# Mode 3: PCA search  (DESI redrock-style, no redrock required)
# ---------------------------------------------------------------------------

PCA_DIR = os.path.join(TEMPLATE_DIR, 'pca')

PCA_FILES = {
    'galaxy':  'rrtemplate-GALAXY-None-v2.6.fits',
    'qso_loz': 'rrtemplate-QSO-LOZ-v1.1.fits',
    'qso_hiz': 'rrtemplate-QSO-HIZ-v1.1.fits',
}
PCA_NAMES = list(PCA_FILES.keys())

_PCA_Z_RANGE = {
    'galaxy':  (0.0, 1.6),   # rest 1228–11000 Å → DESI/SDSS optical range
    'qso_loz': (0.0, 2.5),   # rest 1329–9634 Å, low-z QSO
    'qso_hiz': (1.0, 6.0),   # rest 444–4499 Å, high-z QSO (Lyα etc.)
}


def _load_pca(template_set: str):
    """
    Load DESI redrock PCA eigenvectors.

    Handles both wavelength grid conventions found in DESI templates:
      - log10-linear: CTYPE1 contains 'log10' or CRVAL1 < 10
      - linear (Å):   otherwise (e.g. GALAXY template)

    Returns
    -------
    wave       : ndarray (Å, rest frame)
    eigenvecs  : ndarray (n_components, n_pixels)
    """
    fname = PCA_FILES.get(template_set)
    if fname is None:
        raise ValueError(
            f"Unknown PCA template set '{template_set}'. "
            f"Available: {PCA_NAMES}"
        )
    path = os.path.join(PCA_DIR, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"PCA template not found: {path}\n"
            f"Run: python {TEMPLATE_DIR}/download_desi_pca.py"
        )
    from astropy.io import fits as pyfits
    with pyfits.open(path) as hdul:
        # Find the eigenvector HDU (named BASIS_VECTORS or first with 2-D data)
        data, hdr = None, None
        for h in hdul:
            if h.data is not None and np.ndim(h.data) == 2:
                data = h.data.astype(float)
                hdr  = h.header
                break
        if data is None:
            raise ValueError(f"Cannot find 2-D eigenvector data in {path}")

        crval1 = float(hdr['CRVAL1'])
        cdelt1 = float(hdr['CDELT1'])
        naxis1 = int(hdr['NAXIS1'])
        ctype1 = str(hdr.get('CTYPE1', '')).lower()

        pix = np.arange(naxis1, dtype=float)
        if 'log' in ctype1 or crval1 < 10.0:
            # log10-linear grid (QSO templates)
            wave = 10.0 ** (crval1 + cdelt1 * pix)
        else:
            # linear Å grid (GALAXY template: CRVAL1=1228, CDELT1=0.1)
            wave = crval1 + cdelt1 * pix

    # data shape (NAXIS2, NAXIS1) = (n_components, n_pixels) in Python/astropy
    if data.ndim == 1:
        data = data[np.newaxis, :]

    # Downsample over-resolved templates to at most ~10 000 pixels.
    # The galaxy template has 97 720 px at 0.1 Å/px — far more than needed
    # for a redshift search.  Downsampling 10× gives 1 Å/px, still adequate.
    max_pix = 10_000
    n_pix = data.shape[1]
    if n_pix > max_pix:
        step = n_pix // max_pix
        data = data[:, ::step]
        wave = wave[::step]

    return wave, data   # (n_comp, n_pix)


def pca_search(spec, template_set='galaxy', z_min=None, z_max=None, n_steps=5000,
               fit_continuum=True,
               data_norm='normalize', model_norm='normalize',
               smooth_pixels=1, wave_min=None, wave_max=None):
    """
    Chi-square scan using DESI redrock PCA eigenvectors (no redrock required).

    At each trial redshift the eigenvectors are shifted to the observed frame,
    interpolated onto the spectrum wavelength grid, and the optimal linear
    combination is solved by weighted least squares:

        A c = b   where  A = M^T diag(ivar) M,  b = M^T diag(ivar) F
        model = M c
        chi2(z) = Σ((F - model)² · ivar) / N_good

    This is equivalent to what redrock does but without the GPU backend.

    Parameters
    ----------
    spec          : rb_spectrum
    template_set  : str  — 'galaxy', 'qso_loz', or 'qso_hiz'
    z_min, z_max  : float or None  — defaults from _PCA_Z_RANGE
    n_steps       : int
    fit_continuum : bool
    data_norm     : str  — 'normalize' | 'subtract' | 'raw'
                          how to preprocess the observed spectrum
    model_norm    : str  — 'normalize' | 'subtract' | 'raw'
                          how to preprocess the PCA eigenvectors
                          ('normalize' = L2-normalise; 'subtract' = mean-center;
                           'raw' = use as-is)
    smooth_pixels : int  — boxcar pre-smoothing of flux_work; 1=off
    wave_min, wave_max : float or None — observed-frame Å mask (ivar→0 outside)

    Returns
    -------
    ZFindResult
    """
    wave, flux, ivar, continuum, warn_list = _preprocess(
        spec, fit_continuum=fit_continuum
    )

    # Wavelength mask
    if wave_min is not None:
        ivar[wave < wave_min] = 0.0
    if wave_max is not None:
        ivar[wave > wave_max] = 0.0

    # Prepare observed spectrum
    cont_safe = np.where(np.abs(continuum) > 1e-10, continuum, 1.0)
    if data_norm == 'normalize':
        flux_work = flux / cont_safe
    elif data_norm == 'subtract':
        flux_work = flux - continuum
    else:   # 'raw'
        flux_work = flux

    if smooth_pixels > 1:
        try:
            from astropy.convolution import convolve, Box1DKernel
            flux_work = convolve(flux_work, Box1DKernel(smooth_pixels))
        except Exception:
            pass

    t_wave, eigenvecs = _load_pca(template_set)
    n_comp = eigenvecs.shape[0]

    # Prepare eigenvectors according to model_norm
    if model_norm == 'normalize':
        # L2-normalise each eigenvector (scale-invariant matching)
        norms = np.sqrt(np.sum(eigenvecs ** 2, axis=1, keepdims=True))
        eigenvecs = eigenvecs / np.where(norms > 0, norms, 1.0)
    elif model_norm == 'subtract':
        # Mean-center each eigenvector (remove DC offset)
        eigenvecs = eigenvecs - eigenvecs.mean(axis=1, keepdims=True)
    # 'raw': use eigenvectors as loaded

    z_lo, z_hi = _PCA_Z_RANGE.get(template_set, (0.0, 6.0))
    if z_min is None:
        z_min = z_lo
    if z_max is None:
        z_max = z_hi

    z_array    = np.linspace(z_min, z_max, n_steps)
    chi2_array = np.full(n_steps, np.nan)
    wmin, wmax = wave.min(), wave.max()

    for i, z in enumerate(z_array):
        t_wave_obs = t_wave * (1.0 + z)

        # Require at least 20 % wavelength overlap
        overlap_lo = max(wmin, t_wave_obs.min())
        overlap_hi = min(wmax, t_wave_obs.max())
        if (overlap_hi - overlap_lo) < 0.2 * (wmax - wmin):
            continue

        # Interpolate each eigenvector onto the spectrum wavelength grid
        # M : (n_obs_pix, n_comp)
        M = np.zeros((len(wave), n_comp), dtype=float)
        for j in range(n_comp):
            M[:, j] = np.interp(wave, t_wave_obs, eigenvecs[j],
                                 left=0.0, right=0.0)

        good = ivar > 0
        if good.sum() < 50:
            continue

        Mg  = M[good]           # (n_good, n_comp)
        Fg  = flux_work[good]   # (n_good,)
        IVg = ivar[good]        # (n_good,)

        # Weighted normal equations
        WMg = Mg * IVg[:, np.newaxis]   # (n_good, n_comp)
        A   = WMg.T @ Mg                 # (n_comp, n_comp)
        b   = WMg.T @ Fg                 # (n_comp,)

        try:
            # Use lstsq for numerical robustness (handles near-singular A)
            c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        except np.linalg.LinAlgError:
            continue

        resid          = Fg - Mg @ c
        chi2_array[i]  = float(np.sum(resid ** 2 * IVg) / good.sum())

    # Extract solutions
    top_idx = _find_top_minima(z_array, chi2_array, n=10, min_dz=0.02)
    method_label = f'PCA:{template_set}'
    solutions = []
    for idx in top_idx:
        solutions.append(ZSolution(
            z=float(z_array[idx]),
            z_err=_z_err_from_curvature(z_array, chi2_array, idx),
            chi2_dof=float(chi2_array[idx]),
            method=method_label,
            template_type=template_set,
            n_features=n_comp,
        ))

    return ZFindResult(
        z_array=z_array,
        chi2_curves=[{'label': method_label, 'chi2': chi2_array.copy()}],
        solutions=solutions,
        input_spec=spec,
        warnings=warn_list,
    )


def multi_pca_search(spec, template_sets=None, z_min=None, z_max=None,
                     n_steps=5000, fit_continuum=True,
                     data_norm='normalize', model_norm='normalize',
                     smooth_pixels=1, wave_min=None, wave_max=None):
    """
    Run pca_search for every set in *template_sets* (default: all),
    return combined ZFindResult with one chi2 curve per set.

    Each template is run on its own default z range unless z_min/z_max are
    provided explicitly. Chi2 curves are interpolated onto a common grid
    spanning the union of all individual z ranges before returning.
    """
    if template_sets is None:
        template_sets = PCA_NAMES

    all_results = {}
    warn_all    = []
    for tset in template_sets:
        try:
            r = pca_search(spec, template_set=tset,
                           z_min=z_min, z_max=z_max,
                           n_steps=n_steps,
                           fit_continuum=fit_continuum,
                           data_norm=data_norm, model_norm=model_norm,
                           smooth_pixels=smooth_pixels,
                           wave_min=wave_min, wave_max=wave_max)
            all_results[tset] = r
            warn_all.extend(r.warnings)
        except Exception as e:
            warn_all.append(f'PCA:{tset} failed: {e}')

    if not all_results:
        raise RuntimeError('All PCA template sets failed.')

    # Build a common z grid spanning the union of all individual ranges
    z_lo = min(r.z_array[0]  for r in all_results.values())
    z_hi = max(r.z_array[-1] for r in all_results.values())
    z_array = np.linspace(z_lo, z_hi, n_steps)

    chi2_curves = []
    for tset, r in all_results.items():
        # Interpolate each chi2 curve onto the common grid (NaN outside range)
        interp_chi2 = np.interp(z_array, r.z_array, r.chi2_curves[0]['chi2'],
                                left=np.nan, right=np.nan)
        chi2_curves.append({'label': f'PCA:{tset}', 'chi2': interp_chi2})

    best_tset = min(
        all_results,
        key=lambda t: (all_results[t].solutions[0].chi2_dof
                       if all_results[t].solutions else np.inf)
    )
    solutions = all_results[best_tset].solutions

    return ZFindResult(
        z_array=z_array,
        chi2_curves=chi2_curves,
        solutions=solutions,
        input_spec=spec,
        warnings=list(dict.fromkeys(warn_all)),
    )
