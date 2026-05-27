"""
io.py — Input normalization and output dataclasses for rb_zfind.

All engine functions receive rb_spectrum objects.
All engine functions return ZFindResult or AbsorberResult objects.

This module has NO GUI dependencies.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional

from rbcodes.utils.rb_spectrum import rb_spectrum


# ---------------------------------------------------------------------------
# Input normalization
# ---------------------------------------------------------------------------

def to_rb_spectrum(wave_or_spec, flux=None, error=None, continuum=None):
    """
    Normalize any supported input to an rb_spectrum object.

    Accepted forms
    --------------
    to_rb_spectrum('spectrum.fits')
    to_rb_spectrum(rb_spectrum_obj)
    to_rb_spectrum(wave_array, flux_array)
    to_rb_spectrum(wave_array, flux_array, error_array)
    to_rb_spectrum(wave_array, flux_array, error_array, continuum_array)
    to_rb_spectrum((wave, flux))
    to_rb_spectrum((wave, flux, error))
    to_rb_spectrum((wave, flux, error, continuum))

    Returns
    -------
    rb_spectrum
    """
    # filepath string
    if isinstance(wave_or_spec, str):
        return rb_spectrum.from_file(wave_or_spec)

    # already an rb_spectrum
    if isinstance(wave_or_spec, rb_spectrum):
        return wave_or_spec

    # tuple input
    if isinstance(wave_or_spec, (tuple, list)):
        return rb_spectrum.from_tuple(tuple(wave_or_spec))

    # numpy array — needs at least flux
    if isinstance(wave_or_spec, np.ndarray):
        if flux is None:
            raise ValueError("flux array required when wave is passed as ndarray")
        parts = (wave_or_spec, flux)
        if error is not None:
            parts += (error,)
        if continuum is not None:
            parts += (continuum,)
        return rb_spectrum.from_tuple(parts)

    raise TypeError(
        f"Unsupported input type: {type(wave_or_spec)}. "
        "Expected str (filepath), rb_spectrum, tuple, or numpy array."
    )


# ---------------------------------------------------------------------------
# Output dataclasses
# ---------------------------------------------------------------------------

@dataclass
class ZSolution:
    """Single redshift solution from any engine mode."""
    z            : float
    z_err        : float        # from chi2 curvature: sqrt(1/|d²chi2/dz²|)
    chi2_dof     : float
    method       : str          # 'LineSearch', 'Template:Passive', 'PCA:Galaxy', etc.
    template_type: str          # 'Galaxy', 'QSO', 'Star', 'Unknown'
    n_features   : int          # lines or template pixels contributing to fit

    def __repr__(self):
        return (f"ZSolution(z={self.z:.5f}, z_err={self.z_err:.5f}, "
                f"chi2/dof={self.chi2_dof:.3f}, method={self.method})")


@dataclass
class ZFindResult:
    """
    Output from emission-mode search (galaxy/QSO redshift).

    chi2_curves is a list of dicts so multiple engine modes can be overlaid:
        [{'label': 'LineSearch:ISM', 'chi2': ndarray}, ...]
    """
    z_array     : np.ndarray
    chi2_curves : List[dict]            # one entry per engine mode that ran
    solutions   : List[ZSolution]       # top N solutions, sorted by chi2_dof asc
    input_spec  : rb_spectrum           # stored for re-plotting without recomputing
    warnings    : List[str] = field(default_factory=list)

    def best(self):
        """Return the best solution (lowest chi2_dof)."""
        if not self.solutions:
            return None
        return self.solutions[0]

    def __repr__(self):
        n = len(self.solutions)
        best = self.best()
        return (f"ZFindResult({n} solutions, "
                f"best z={best.z:.5f} [{best.method}])" if best else "ZFindResult(empty)")


@dataclass
class AbsorberCandidate:
    """Single absorber candidate from absorption-mode search."""
    z            : float
    significance : float        # sigma above noise floor
    n_lines      : int
    is_doublet   : bool         # doublet match = far more reliable than single line
    linelist_name: str          # e.g. 'CIV', 'MgII', 'SiIV'
    lines_matched: List[str]    # e.g. ['CIV 1548', 'CIV 1550']

    def __repr__(self):
        dtype = "doublet" if self.is_doublet else "single"
        return (f"AbsorberCandidate(z={self.z:.5f}, {self.significance:.1f}σ, "
                f"{self.linelist_name} {dtype})")


@dataclass
class AbsorberResult:
    """
    Output from absorption-mode search (QSO sightline absorbers).

    Returns a list of candidates, not a single z — multiple absorbers
    at different redshifts are the expected output.
    """
    z_array           : np.ndarray
    significance_curve: np.ndarray      # significance vs z (not chi2)
    candidates        : List[AbsorberCandidate]  # sorted by significance desc
    input_spec        : rb_spectrum
    warnings          : List[str] = field(default_factory=list)

    def __repr__(self):
        return f"AbsorberResult({len(self.candidates)} candidates)"
