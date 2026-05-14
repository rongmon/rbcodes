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
    Find top N local minima in chi2_array with minimum separation min_dz.

    Returns list of indices into z_array.
    """
    from scipy.signal import argrelmin
    finite_mask = np.isfinite(chi2_array)
    if not np.any(finite_mask):
        return []

    # find all local minima
    min_idx = argrelmin(chi2_array, order=5)[0]
    min_idx = min_idx[finite_mask[min_idx]]

    if len(min_idx) == 0:
        # fallback: global minimum
        min_idx = np.array([np.nanargmin(chi2_array)])

    # sort by chi2 value
    min_idx = min_idx[np.argsort(chi2_array[min_idx])]

    # enforce minimum z separation
    selected = []
    for idx in min_idx:
        z = z_array[idx]
        if all(abs(z - z_array[s]) > min_dz for s in selected):
            selected.append(idx)
        if len(selected) >= n:
            break

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
                window_pixels=5):
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

    Returns
    -------
    ZFindResult  (mode='emission')
    AbsorberResult  (mode='absorption')
    """
    wave, flux, ivar, continuum, warn_list = _preprocess(spec, fit_continuum=fit_continuum)

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

            fw = flux[lo:hi]
            cw = continuum[lo:hi]
            iv = ivar[lo:hi]

            if np.sum(iv > 0) < 2:
                continue

            if mode == 'emission':
                residual = fw - cw
                # fit Gaussian amplitude analytically: A = sum(residual * gauss * ivar) / sum(gauss^2 * ivar)
                # simplified: just use sum of (flux - continuum) weighted by ivar
                signal = np.sum(residual * iv)
                norm = np.sum(iv)
                if norm > 0:
                    A = signal / norm
                    accum_chi2 += np.sum((residual - A) ** 2 * iv) / max(1, np.sum(iv > 0))
                    n_used += 1

            else:  # absorption
                decrement = cw - fw
                signal = np.sum(decrement * iv)
                norm = np.sum(iv)
                if norm > 0:
                    A = signal / norm
                    accum_chi2 += np.sum((decrement - A) ** 2 * iv) / max(1, np.sum(iv > 0))
                    n_used += 1

        if n_used > 0:
            chi2_array[i] = accum_chi2 / n_used

    # --- extract solutions ---
    top_idx = _find_top_minima(z_array, chi2_array, n=3)

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

            top_idx_abs = _find_top_minima(z_array, chi2_array, n=10)
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
# Mode 2: Template search  (stub — implement in Phase 7)
# ---------------------------------------------------------------------------

def template_search(spec, template_name, z_min=0.0, z_max=6.0, n_steps=10000,
                    fwhm_ang=0.0, fit_continuum=True, template_set='marz'):
    """
    Chi-square scan against a bundled 1D FITS template spectrum.

    Parameters
    ----------
    spec : rb_spectrum
    template_name : str  — e.g. 'ELG', 'Passive', 'QSO', 'LRG', 'LBG'
    template_set : str   — 'marz' or 'sdss'

    Returns
    -------
    ZFindResult
    """
    raise NotImplementedError(
        "template_search() not yet implemented. See Phase 7 in ZFIND_TODO.md."
    )


# ---------------------------------------------------------------------------
# Mode 3: PCA search  (stub — implement in Phase 8)
# ---------------------------------------------------------------------------

def pca_search(spec, template_set='galaxy', z_min=0.0, z_max=6.0, n_steps=10000,
               fwhm_ang=0.0, fit_continuum=True):
    """
    Chi-square scan using PCA eigenvectors (redrock-style, no redrock required).
    Eigenvectors stored as FITS in templates/pca/.

    Parameters
    ----------
    spec : rb_spectrum
    template_set : str  — 'galaxy' or 'qso'

    Returns
    -------
    ZFindResult
    """
    raise NotImplementedError(
        "pca_search() not yet implemented. See Phase 8 in ZFIND_TODO.md."
    )
