"""
Phase 2 tests: line_search() emission mode.

Tests cover:
  - Synthetic spectrum at known z → best solution within tolerance
  - Result type and structure (ZFindResult, chi2_curves, solutions)
  - MAD-STD IVAR fallback when no error array is provided
  - Solutions ordered by chi2 (most negative first)
  - z_err is finite and small
  - n_features count
  - All lines out of range → empty/nan chi2 handled gracefully
  - Single-line search still finds true z
  - High-z test (z=1.5)
  - Air-to-vacuum conversion warning
  - n_steps parameter respected
"""

import numpy as np
import pandas as pd
import pytest
import astropy.units as u

from rbcodes.utils.rb_spectrum import rb_spectrum
from rbcodes.GUIs.zfind.engine import line_search
from rbcodes.GUIs.zfind.io import ZFindResult, ZSolution


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Common galaxy emission lines (rest wavelengths in Å, vacuum)
GAL_EM = {
    'Halpha':     6562.80,
    'Hbeta':      4861.33,
    '[OIII]5007': 5006.84,
    '[OIII]4959': 4958.91,
    '[OII]3727':  3727.09,
}


def _make_linelist(line_dict, name='TestEM'):
    """Build a DataFrame linelist from a {name: wave_rest} dict."""
    df = pd.DataFrame({
        'wave': list(line_dict.values()),
        'name': list(line_dict.keys()),
    })
    df.attrs['name'] = name
    return df


def _make_em_spectrum(true_z, lines, n_pix=3000, snr=100.0, rng_seed=42,
                      provide_error=True, provide_continuum=True,
                      airvac='vac'):
    """
    Synthetic spectrum with narrow emission lines at true_z.

    Returns rb_spectrum with wavelength in observed frame.
    """
    wave_min = 3000.0 * (1.0 + true_z)
    wave_max = 9000.0 * (1.0 + true_z)
    wave_obs = np.linspace(wave_min, wave_max, n_pix)

    sigma_noise = 1.0 / snr
    rng = np.random.default_rng(rng_seed)
    flux = 1.0 + rng.normal(0, sigma_noise, n_pix)

    for lam_rest in lines.values():
        lam_obs = lam_rest * (1.0 + true_z)
        idx = np.argmin(np.abs(wave_obs - lam_obs))
        flux[idx - 1:idx + 2] += 5.0  # 3-pixel injection

    error = np.full_like(flux, sigma_noise) if provide_error else None
    continuum = np.ones_like(flux) if provide_continuum else None

    kw = dict(meta={'airvac': airvac})
    return rb_spectrum(
        wave_obs * u.AA,
        flux * u.dimensionless_unscaled,
        error=error * u.dimensionless_unscaled if error is not None else None,
        continuum=continuum * u.dimensionless_unscaled if continuum is not None else None,
        **kw,
    )


# ---------------------------------------------------------------------------
# Tests: result type and structure
# ---------------------------------------------------------------------------

class TestLineSearchEmissionStructure:

    def setup_method(self):
        self.true_z = 0.5
        self.lines = GAL_EM
        self.df = _make_linelist(self.lines)
        self.spec = _make_em_spectrum(self.true_z, self.lines)
        self.result = line_search(
            self.spec, self.df,
            z_min=0.3, z_max=0.7, n_steps=2000,
            mode='emission', fit_continuum=False,
        )

    def test_returns_zfind_result(self):
        assert isinstance(self.result, ZFindResult)

    def test_z_array_shape(self):
        assert len(self.result.z_array) == 2000

    def test_z_array_bounds(self):
        assert self.result.z_array[0] == pytest.approx(0.3, abs=1e-9)
        assert self.result.z_array[-1] == pytest.approx(0.7, abs=1e-9)

    def test_chi2_curves_present(self):
        assert len(self.result.chi2_curves) == 1
        curve = self.result.chi2_curves[0]
        assert 'label' in curve
        assert 'chi2' in curve
        assert len(curve['chi2']) == 2000

    def test_chi2_label_contains_linelist_name(self):
        assert 'TestEM' in self.result.chi2_curves[0]['label']

    def test_solutions_list(self):
        assert isinstance(self.result.solutions, list)
        assert len(self.result.solutions) >= 1

    def test_solutions_are_zsolution(self):
        for sol in self.result.solutions:
            assert isinstance(sol, ZSolution)

    def test_input_spec_preserved(self):
        assert self.result.input_spec is self.spec

    def test_no_warnings_when_error_and_continuum_given(self):
        assert self.result.warnings == []


# ---------------------------------------------------------------------------
# Tests: accuracy
# ---------------------------------------------------------------------------

class TestLineSearchEmissionAccuracy:

    def _run(self, true_z, lines, z_min, z_max, **kw):
        df = _make_linelist(lines)
        spec = _make_em_spectrum(true_z, lines)
        return line_search(spec, df, z_min=z_min, z_max=z_max,
                           n_steps=2000, mode='emission',
                           fit_continuum=False, **kw), true_z

    def test_three_lines_z05(self):
        # Tolerance 0.002 — boxcar pre-smoothing + z-grid quantization over
        # [0.3, 0.7] at 2000 steps can place the minimum ~0.001-0.002 from truth.
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta', '[OIII]5007')}
        result, true_z = self._run(0.5, lines, z_min=0.3, z_max=0.7)
        assert abs(result.best().z - true_z) < 0.002

    def test_five_lines_z02(self):
        # Tolerance is 0.002 — grid quantization at 2000 steps over [0,0.4] can
        # put the best-fit z up to ~0.5 pixel_scale / rest_wave ≈ 0.001 from truth.
        result, true_z = self._run(0.2, GAL_EM, z_min=0.0, z_max=0.4)
        assert abs(result.best().z - true_z) < 0.002

    def test_high_redshift_z15(self):
        # Use optical lines — at z=1.5 they land at 16408, 12153, 12528 Å,
        # well within the 3000*(1+1.5)–9000*(1+1.5) = 7500–22500 Å range.
        optical_lines = {'Halpha': 6562.80, 'Hbeta': 4861.33, '[OIII]5007': 5006.84}
        df = _make_linelist(optical_lines, name='OptHiZ')
        spec = _make_em_spectrum(1.5, optical_lines, n_pix=4000)
        result = line_search(spec, df, z_min=1.3, z_max=1.7, n_steps=2000,
                             mode='emission', fit_continuum=False)
        assert result.best() is not None
        assert abs(result.best().z - 1.5) < 0.005

    def test_single_line_still_finds_z(self):
        lines = {'Halpha': 6562.8}
        df = _make_linelist(lines, name='HalphaOnly')
        spec = _make_em_spectrum(0.4, lines)
        result = line_search(spec, df, z_min=0.2, z_max=0.6, n_steps=2000,
                             mode='emission', fit_continuum=False)
        assert abs(result.best().z - 0.4) < 0.002

    def test_best_z_is_most_negative_chi2(self):
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta', '[OIII]5007')}
        result, _ = self._run(0.5, lines, z_min=0.3, z_max=0.7)
        best = result.best()
        # best solution must have chi2_dof ≤ all other solutions
        for sol in result.solutions:
            assert best.chi2_dof <= sol.chi2_dof


# ---------------------------------------------------------------------------
# Tests: z_err
# ---------------------------------------------------------------------------

class TestLineSearchEmissionZErr:

    def test_z_err_is_finite(self):
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta', '[OIII]5007')}
        df = _make_linelist(lines)
        spec = _make_em_spectrum(0.5, lines)
        result = line_search(spec, df, z_min=0.3, z_max=0.7,
                             n_steps=2000, mode='emission', fit_continuum=False)
        assert np.isfinite(result.best().z_err)

    def test_z_err_small_for_strong_lines(self):
        # high SNR, many lines → z_err should be very small
        lines = GAL_EM
        df = _make_linelist(lines)
        spec = _make_em_spectrum(0.5, lines, snr=500.0)
        result = line_search(spec, df, z_min=0.3, z_max=0.7,
                             n_steps=4000, mode='emission', fit_continuum=False)
        assert result.best().z_err < 0.001


# ---------------------------------------------------------------------------
# Tests: n_features
# ---------------------------------------------------------------------------

class TestLineSearchNFeatures:

    def test_n_features_at_least_one(self):
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta')}
        df = _make_linelist(lines)
        spec = _make_em_spectrum(0.5, lines)
        result = line_search(spec, df, z_min=0.3, z_max=0.7,
                             n_steps=2000, mode='emission', fit_continuum=False)
        assert result.best().n_features >= 1

    def test_n_features_bounded_by_linelist_size(self):
        df = _make_linelist(GAL_EM)
        spec = _make_em_spectrum(0.3, GAL_EM)
        result = line_search(spec, df, z_min=0.1, z_max=0.5,
                             n_steps=2000, mode='emission', fit_continuum=False)
        assert result.best().n_features <= len(GAL_EM)


# ---------------------------------------------------------------------------
# Tests: warnings and fallbacks
# ---------------------------------------------------------------------------

class TestLineSearchEmissionWarnings:

    def test_mad_std_warning_when_no_error(self):
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta', '[OIII]5007')}
        df = _make_linelist(lines)
        spec = _make_em_spectrum(0.5, lines, provide_error=False)
        result = line_search(spec, df, z_min=0.3, z_max=0.7,
                             n_steps=1000, mode='emission', fit_continuum=False)
        assert any('MAD' in w for w in result.warnings)

    def test_mad_std_fallback_still_finds_z(self):
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta', '[OIII]5007')}
        df = _make_linelist(lines)
        spec = _make_em_spectrum(0.5, lines, provide_error=False)
        result = line_search(spec, df, z_min=0.3, z_max=0.7,
                             n_steps=2000, mode='emission', fit_continuum=False)
        assert abs(result.best().z - 0.5) < 0.002

    def test_air_wavelength_warning(self):
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta')}
        df = _make_linelist(lines)
        # Construct an air-wavelength spectrum (wavelengths slightly shorter)
        spec = _make_em_spectrum(0.5, lines, airvac='air')
        result = line_search(spec, df, z_min=0.3, z_max=0.7,
                             n_steps=1000, mode='emission', fit_continuum=False)
        assert any('vacuum' in w.lower() or 'air' in w.lower()
                   for w in result.warnings)


# ---------------------------------------------------------------------------
# Tests: edge cases
# ---------------------------------------------------------------------------

class TestLineSearchEmissionEdgeCases:

    def test_all_lines_out_of_range_gives_all_nan_chi2(self):
        # Spectrum covers 3000–5000 Å; search at z=0.0 → Halpha (6563 Å) is out of range
        wave = np.linspace(3000, 5000, 500)
        flux = np.ones(500)
        spec = rb_spectrum(wave * u.AA, flux * u.dimensionless_unscaled,
                           continuum=flux * u.dimensionless_unscaled)
        df = pd.DataFrame({'wave': [6562.8], 'name': ['Halpha']})
        df.attrs['name'] = 'Halpha'
        result = line_search(spec, df, z_min=0.0, z_max=0.1,
                             n_steps=100, mode='emission', fit_continuum=False)
        assert isinstance(result, ZFindResult)
        # All chi2 should be NaN (no line in range)
        assert np.all(np.isnan(result.chi2_curves[0]['chi2']))

    def test_all_lines_out_of_range_best_is_none(self):
        wave = np.linspace(3000, 5000, 500)
        flux = np.ones(500)
        spec = rb_spectrum(wave * u.AA, flux * u.dimensionless_unscaled,
                           continuum=flux * u.dimensionless_unscaled)
        df = pd.DataFrame({'wave': [6562.8], 'name': ['Halpha']})
        df.attrs['name'] = 'Halpha'
        result = line_search(spec, df, z_min=0.0, z_max=0.1,
                             n_steps=100, mode='emission', fit_continuum=False)
        assert result.best() is None

    def test_n_steps_respected(self):
        lines = {k: GAL_EM[k] for k in ('Halpha',)}
        df = _make_linelist(lines)
        spec = _make_em_spectrum(0.5, lines)
        for n in (500, 1000, 3000):
            result = line_search(spec, df, z_min=0.3, z_max=0.7,
                                 n_steps=n, mode='emission', fit_continuum=False)
            assert len(result.z_array) == n
            assert len(result.chi2_curves[0]['chi2']) == n

    def test_method_label_in_solution(self):
        lines = {k: GAL_EM[k] for k in ('Halpha', 'Hbeta')}
        df = _make_linelist(lines, name='MyList')
        spec = _make_em_spectrum(0.5, lines)
        result = line_search(spec, df, z_min=0.3, z_max=0.7,
                             n_steps=1000, mode='emission', fit_continuum=False)
        assert 'MyList' in result.best().method
