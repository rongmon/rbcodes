"""
LineFitter.py — Quick line-fitting engine for rb_multispec viewer.

Two fitting modes
-----------------
fit_gaussian(wave, flux, x1, y1, x2, y2)  ->  dict
fit_com(wave, flux, x1, y1, x2, y2)       ->  dict

The two (x, y) anchor points define both the fit window [x1, x2] and
a linear continuum:

    cont(λ) = y1 + (y2 - y1) / (x2 - x1) * (λ - x1)

This means the user controls the continuum simply by clicking at the
appropriate flux level on each side of the line — a tilted continuum
is fully supported, as is fitting a half-profile.

Absorption vs emission is detected automatically from the sign of the
integrated residual: negative sum → absorption, positive → emission.

No error arrays are used; both fits are unweighted least-squares / moments.
"""

import numpy as np
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _linear_continuum(wave, x1, y1, x2, y2):
    """Evaluate the linear continuum at each wavelength in wave."""
    return y1 + (y2 - y1) / (x2 - x1) * (wave - x1)


def _ensure_left_right(x1, y1, x2, y2):
    """Return (x1,y1,x2,y2) with x1 < x2."""
    if x1 > x2:
        return x2, y2, x1, y1
    return x1, y1, x2, y2


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fit_gaussian(wave, flux, x1, y1, x2, y2):
    """
    Fit a single Gaussian to (flux - continuum) between x1 and x2.

    Parameters
    ----------
    wave, flux : array-like
        Full spectrum arrays (only the window [x1, x2] is used).
    x1, y1 : float
        Left anchor — wavelength and flux at the left edge of the fit window.
        This point also sets the left continuum level.
    x2, y2 : float
        Right anchor — same for the right edge.

    Returns
    -------
    dict
        centroid   : float  — fitted line centre (Å)
        fwhm_ang   : float  — FWHM in Å
        fwhm_kms   : float  — FWHM in km/s
        amplitude  : float  — peak height above/below continuum
                              (negative = absorption, positive = emission)
        direction  : int    — -1 absorption, +1 emission
        asymmetric : bool   — True if centroid is within 20 % of either edge
                              (suggests truncated profile; FWHM may be unreliable)
        fit_wave   : ndarray — wavelength grid for the model overlay
        fit_flux   : ndarray — model flux (continuum + Gaussian)
        cont_wave  : ndarray — [x1, x2] for plotting the continuum segment
        cont_flux  : ndarray — [y1, y2]

    Raises
    ------
    ValueError   : fewer than 5 finite pixels inside the window
    RuntimeError : scipy curve_fit did not converge
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    x1, y1, x2, y2 = _ensure_left_right(x1, y1, x2, y2)

    mask = (wave >= x1) & (wave <= x2) & np.isfinite(flux)
    n_pix = mask.sum()
    if n_pix < 5:
        raise ValueError(
            f"Only {n_pix} finite pixels in window [{x1:.1f}, {x2:.1f}] Å — need ≥ 5"
        )

    w = wave[mask]
    f = flux[mask]
    cont = _linear_continuum(w, x1, y1, x2, y2)
    residual = f - cont

    direction = -1 if residual.sum() < 0 else 1
    data = direction * residual  # always positive for the Gaussian fitter

    # Initial guesses
    amp_guess = float(np.max(data))
    if amp_guess <= 0:
        amp_guess = float(np.abs(residual).max())
    cen_guess = float(w[np.argmax(data)])
    sig_guess = (x2 - x1) / 4.0

    def _gauss(x, amp, cen, sig):
        return amp * np.exp(-0.5 * ((x - cen) / sig) ** 2)

    try:
        popt, _ = curve_fit(
            _gauss, w, data,
            p0=[amp_guess, cen_guess, sig_guess],
            bounds=([0.0, x1, 1e-3], [np.inf, x2, (x2 - x1)]),
            maxfev=3000,
        )
    except Exception as exc:
        raise RuntimeError(f"Gaussian curve_fit failed: {exc}")

    amp, cen, sig = popt
    fwhm_ang = 2.3548 * abs(sig)
    c_kms = 2.998e5
    fwhm_kms = fwhm_ang / cen * c_kms if cen > 0 else 0.0

    window = x2 - x1
    asymmetric = ((cen - x1) < 0.20 * window) or ((x2 - cen) < 0.20 * window)

    # Model curve for overlay
    fit_wave = np.linspace(x1, x2, 300)
    fit_cont = _linear_continuum(fit_wave, x1, y1, x2, y2)
    fit_flux = fit_cont + direction * _gauss(fit_wave, amp, cen, sig)

    return {
        'centroid':   float(cen),
        'fwhm_ang':   float(fwhm_ang),
        'fwhm_kms':   float(fwhm_kms),
        'amplitude':  float(direction * amp),
        'direction':  direction,
        'asymmetric': asymmetric,
        'fit_wave':   fit_wave,
        'fit_flux':   fit_flux,
        'cont_wave':  np.array([x1, x2]),
        'cont_flux':  np.array([y1, y2]),
    }


def fit_com(wave, flux, x1, y1, x2, y2):
    """
    Centre-of-mass centroiding between x1 and x2.

    The continuum is the same linear model as in fit_gaussian.  Weights are
    the (continuum-subtracted) flux values, clipped to zero for pixels where
    the residual has the wrong sign (noise).

    Parameters
    ----------
    wave, flux : array-like  (full spectrum)
    x1, y1, x2, y2 : float  (same meaning as fit_gaussian)

    Returns
    -------
    dict
        centroid       : float — CoM wavelength (Å)
        sigma_ang      : float — RMS width in Å  (the actual measured quantity)
        sigma_kms      : float — RMS width in km/s
        fwhm_equiv_ang : float — 2.355 × sigma (Å);  assumes Gaussian shape
        fwhm_kms       : float — 2.355 × sigma in km/s
        amplitude      : float — peak of |residual| (sign follows direction)
        direction      : int   — -1 absorption, +1 emission
        asymmetric     : bool  — True if centroid within 20 % of either edge

    Raises
    ------
    ValueError : fewer than 3 finite pixels, or all weights ≤ 0
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    x1, y1, x2, y2 = _ensure_left_right(x1, y1, x2, y2)

    mask = (wave >= x1) & (wave <= x2) & np.isfinite(flux)
    n_pix = mask.sum()
    if n_pix < 3:
        raise ValueError(
            f"Only {n_pix} finite pixels in window [{x1:.1f}, {x2:.1f}] Å — need ≥ 3"
        )

    w = wave[mask]
    f = flux[mask]
    cont = _linear_continuum(w, x1, y1, x2, y2)
    residual = f - cont

    direction = -1 if residual.sum() < 0 else 1
    weights = direction * residual
    weights = np.maximum(weights, 0.0)  # clip noise pixels below zero

    total_weight = weights.sum()
    if total_weight == 0.0:
        raise ValueError("Sum of weights is zero — no signal above continuum in window")

    centroid = float(np.sum(weights * w) / total_weight)
    sigma = float(np.sqrt(np.sum(weights * (w - centroid) ** 2) / total_weight))
    fwhm_equiv = 2.3548 * sigma

    c_kms = 2.998e5
    sigma_kms = sigma / centroid * c_kms if centroid > 0 else 0.0
    fwhm_kms = fwhm_equiv / centroid * c_kms if centroid > 0 else 0.0

    amplitude = float(direction * np.max(weights))

    window = x2 - x1
    asymmetric = ((centroid - x1) < 0.20 * window) or ((x2 - centroid) < 0.20 * window)

    return {
        'centroid':       centroid,
        'sigma_ang':      sigma,
        'sigma_kms':      sigma_kms,
        'fwhm_equiv_ang': fwhm_equiv,
        'fwhm_kms':       fwhm_kms,
        'amplitude':      amplitude,
        'direction':      direction,
        'asymmetric':     asymmetric,
    }
