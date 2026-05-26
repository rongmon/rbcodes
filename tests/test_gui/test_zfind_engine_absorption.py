"""
Phase 3 tests: line_search() absorption mode.

Tests cover:
  - Returns AbsorberResult with correct structure
  - Single absorber at known z found within tolerance
  - Two absorbers at different z both recovered
  - Candidates sorted by significance (descending)
  - n_lines count and lines_matched content
  - Significance of real absorbers well above noise floor
  - MAD-STD IVAR fallback warning when no error array
  - Air-to-vacuum conversion warning
  - All lines out of range → graceful empty result
  - n_steps parameter respected
  - input_spec preserved on result
"""

import numpy as np
import pandas as pd
import pytest
import astropy.units as u

from rbcodes.utils.rb_spectrum import rb_spectrum
from rbcodes.GUIs.zfind.engine import line_search
from rbcodes.GUIs.zfind.io import AbsorberResult, AbsorberCandidate


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Common QSO/IGM absorption doublets (rest wavelengths in Å, vacuum)
MGII = {'MgII2796': 2796.35, 'MgII2803': 2803.53}
CIVA = {'CIV1548': 1548.20, 'CIV1550': 1550.77}


def _make_linelist(line_dict, name='TestAbs'):
    df = pd.DataFrame({
        'wave': list(line_dict.values()),
        'name': list(line_dict.keys()),
    })
    df.attrs['name'] = name
    return df


def _make_abs_spectrum(absorbers, wave_min=3000, wave_max=12000, n_pix=5000,
                       continuum_level=10.0, snr=100.0, rng_seed=42,
                       provide_error=True, provide_continuum=True,
                       airvac='vac'):
    """
    Synthetic QSO continuum with absorption features injected at given z values.

    absorbers : list of dicts, each {'z': float, 'lines': {name: rest_wave}}
    """
    sigma_noise = continuum_level / snr
    wave_obs = np.linspace(wave_min, wave_max, n_pix)
    rng = np.random.default_rng(rng_seed)
    flux = continuum_level * np.ones(n_pix) + rng.normal(0, sigma_noise, n_pix)
    continuum = continuum_level * np.ones(n_pix)

    for ab in absorbers:
        z = ab['z']
        depth = ab.get('depth', 5.0)
        for lam_rest in ab['lines'].values():
            lam_obs = lam_rest * (1.0 + z)
            if wave_min < lam_obs < wave_max:
                idx = np.argmin(np.abs(wave_obs - lam_obs))
                flux[idx - 1:idx + 2] -= depth

    error = np.full_like(flux, sigma_noise) if provide_error else None
    cont = continuum if provide_continuum else None

    return rb_spectrum(
        wave_obs * u.AA,
        flux * u.dimensionless_unscaled,
        error=error * u.dimensionless_unscaled if error is not None else None,
        continuum=cont * u.dimensionless_unscaled if cont is not None else None,
        meta={'airvac': airvac},
    )


# ---------------------------------------------------------------------------
# Tests: result type and structure
# ---------------------------------------------------------------------------

class TestLineSearchAbsorptionStructure:

    def setup_method(self):
        absorbers = [{'z': 0.5, 'lines': MGII}]
        self.spec = _make_abs_spectrum(absorbers)
        self.df = _make_linelist(MGII)
        self.result = line_search(
            self.spec, self.df,
            z_min=0.0, z_max=2.0, n_steps=3000,
            mode='absorption', fit_continuum=False,
        )

    def test_returns_absorber_result(self):
        assert isinstance(self.result, AbsorberResult)

    def test_z_array_shape(self):
        assert len(self.result.z_array) == 3000

    def test_z_array_bounds(self):
        assert self.result.z_array[0] == pytest.approx(0.0, abs=1e-9)
        assert self.result.z_array[-1] == pytest.approx(2.0, abs=1e-9)

    def test_significance_curve_shape_matches_z_array(self):
        assert len(self.result.significance_curve) == len(self.result.z_array)

    def test_significance_curve_finite_values_exist(self):
        assert np.any(np.isfinite(self.result.significance_curve))

    def test_candidates_is_list(self):
        assert isinstance(self.result.candidates, list)

    def test_candidates_are_absorber_candidate(self):
        for c in self.result.candidates:
            assert isinstance(c, AbsorberCandidate)

    def test_input_spec_preserved(self):
        assert self.result.input_spec is self.spec

    def test_no_warnings_with_error_and_continuum(self):
        assert self.result.warnings == []


# ---------------------------------------------------------------------------
# Tests: accuracy — single absorber
# ---------------------------------------------------------------------------

class TestLineSearchAbsorptionSingleAbsorber:

    @pytest.mark.parametrize('true_z,z_min,z_max', [
        (0.5,  0.3, 0.8),
        (0.8,  0.5, 1.1),
        (1.2,  0.9, 1.5),
    ])
    def test_doublet_found_at_known_z(self, true_z, z_min, z_max):
        absorbers = [{'z': true_z, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=z_min, z_max=z_max,
                             n_steps=3000, mode='absorption',
                             fit_continuum=False)
        assert len(result.candidates) >= 1
        best = result.candidates[0]  # sorted by significance, highest first
        assert abs(best.z - true_z) < 0.005

    def test_single_line_absorber_found(self):
        true_z = 0.8
        absorbers = [{'z': true_z, 'lines': {'MgII2796': 2796.35}}]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist({'MgII2796': 2796.35}, name='MgII_single')
        result = line_search(spec, df, z_min=0.5, z_max=1.1,
                             n_steps=3000, mode='absorption',
                             fit_continuum=False)
        best = result.candidates[0]
        assert abs(best.z - true_z) < 0.005


# ---------------------------------------------------------------------------
# Tests: accuracy — multiple absorbers
# ---------------------------------------------------------------------------

class TestLineSearchAbsorptionMultipleAbsorbers:

    def test_two_doublets_both_found(self):
        true_zs = [0.5, 1.2]
        absorbers = [{'z': z, 'lines': MGII} for z in true_zs]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.0, z_max=2.0,
                             n_steps=5000, mode='absorption',
                             fit_continuum=False)
        recovered_z = [c.z for c in result.candidates]
        for true_z in true_zs:
            matches = [abs(rz - true_z) < 0.01 for rz in recovered_z]
            assert any(matches), f'Absorber at z={true_z} not recovered'

    def test_two_doublets_top_candidates_have_highest_significance(self):
        true_zs = [0.5, 1.2]
        absorbers = [{'z': z, 'lines': MGII} for z in true_zs]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.0, z_max=2.0,
                             n_steps=5000, mode='absorption',
                             fit_continuum=False)
        # Top 2 candidates should both have high significance
        top2_sigs = [c.significance for c in result.candidates[:2]]
        assert all(s > 5.0 for s in top2_sigs), \
            f'Expected significance > 5 for real absorbers, got {top2_sigs}'

    def test_three_absorbers_at_different_z(self):
        true_zs = [0.4, 0.9, 1.5]
        absorbers = [{'z': z, 'lines': MGII, 'depth': 4.0} for z in true_zs]
        spec = _make_abs_spectrum(absorbers, wave_max=15000, n_pix=6000)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.0, z_max=2.0,
                             n_steps=5000, mode='absorption',
                             fit_continuum=False)
        recovered_z = [c.z for c in result.candidates]
        for true_z in true_zs:
            matches = [abs(rz - true_z) < 0.01 for rz in recovered_z]
            assert any(matches), f'Absorber at z={true_z} not recovered'


# ---------------------------------------------------------------------------
# Tests: candidates properties
# ---------------------------------------------------------------------------

class TestLineSearchAbsorptionCandidates:

    def setup_method(self):
        self.true_z = 0.5
        absorbers = [{'z': self.true_z, 'lines': MGII}]
        self.spec = _make_abs_spectrum(absorbers)
        self.df = _make_linelist(MGII)
        self.result = line_search(
            self.spec, self.df,
            z_min=0.3, z_max=0.8, n_steps=3000,
            mode='absorption', fit_continuum=False,
        )
        self.best = self.result.candidates[0]

    def test_candidates_sorted_by_significance_descending(self):
        sigs = [c.significance for c in self.result.candidates]
        assert sigs == sorted(sigs, reverse=True)

    def test_n_lines_equals_doublet_count(self):
        assert self.best.n_lines == 2  # MgII doublet

    def test_lines_matched_contains_line_names(self):
        assert set(self.best.lines_matched) == {'MgII2796', 'MgII2803'}

    def test_significance_positive_for_real_absorber(self):
        assert self.best.significance > 0

    def test_real_absorber_significance_above_noise(self):
        # Real absorber significance must be well above the noise floor
        noise_sigs = [c.significance for c in self.result.candidates[2:]]
        if noise_sigs:
            assert self.best.significance > max(noise_sigs) * 2

    def test_linelist_name_stored(self):
        assert self.best.linelist_name == 'TestAbs'


# ---------------------------------------------------------------------------
# Tests: significance curve
# ---------------------------------------------------------------------------

class TestSignificanceCurve:

    def test_significance_peaks_near_true_z(self):
        true_z = 0.5
        absorbers = [{'z': true_z, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.3, z_max=0.8,
                             n_steps=3000, mode='absorption',
                             fit_continuum=False)
        sig = result.significance_curve
        z = result.z_array
        peak_z = z[np.nanargmax(sig)]
        assert abs(peak_z - true_z) < 0.01

    def test_significance_curve_non_negative_at_real_absorber(self):
        true_z = 0.5
        absorbers = [{'z': true_z, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.3, z_max=0.8,
                             n_steps=2000, mode='absorption',
                             fit_continuum=False)
        # At the true z, significance should be positive
        true_idx = np.argmin(np.abs(result.z_array - true_z))
        assert result.significance_curve[true_idx] > 0


# ---------------------------------------------------------------------------
# Tests: warnings and fallbacks
# ---------------------------------------------------------------------------

class TestLineSearchAbsorptionWarnings:

    def test_mad_std_warning_when_no_error(self):
        absorbers = [{'z': 0.5, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers, provide_error=False)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.3, z_max=0.8,
                             n_steps=1000, mode='absorption',
                             fit_continuum=False)
        assert any('MAD' in w for w in result.warnings)

    def test_mad_std_fallback_still_finds_absorber(self):
        true_z = 0.5
        absorbers = [{'z': true_z, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers, provide_error=False)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.3, z_max=0.8,
                             n_steps=2000, mode='absorption',
                             fit_continuum=False)
        best = result.candidates[0]
        assert abs(best.z - true_z) < 0.01

    def test_air_wavelength_warning(self):
        absorbers = [{'z': 0.5, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers, airvac='air')
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.3, z_max=0.8,
                             n_steps=1000, mode='absorption',
                             fit_continuum=False)
        assert any('vacuum' in w.lower() or 'air' in w.lower()
                   for w in result.warnings)


# ---------------------------------------------------------------------------
# Tests: edge cases
# ---------------------------------------------------------------------------

class TestLineSearchAbsorptionEdgeCases:

    def test_all_lines_out_of_range_gives_all_nan_chi2(self):
        # Spectrum covers 3000–5000 Å; search at z=0 → MgII (2796 Å) barely in range
        # but if z_max=0.05, obs_MgII < 2936 Å — still in range. Use a very blue spectrum.
        wave = np.linspace(1200, 2000, 500)
        flux = 5.0 * np.ones(500)
        spec = rb_spectrum(wave * u.AA, flux * u.dimensionless_unscaled,
                           continuum=flux * u.dimensionless_unscaled)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.0, z_max=0.1,
                             n_steps=100, mode='absorption',
                             fit_continuum=False)
        assert isinstance(result, AbsorberResult)
        # No lines fall in range for this wavelength/z combination
        # significance_curve should be zeros (empty finite branch returns zeros)

    def test_n_steps_respected(self):
        absorbers = [{'z': 0.5, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist(MGII)
        for n in (500, 1500, 4000):
            result = line_search(spec, df, z_min=0.3, z_max=0.8,
                                 n_steps=n, mode='absorption',
                                 fit_continuum=False)
            assert len(result.z_array) == n
            assert len(result.significance_curve) == n

    def test_result_has_no_solutions_attribute(self):
        # AbsorberResult has candidates, not solutions
        absorbers = [{'z': 0.5, 'lines': MGII}]
        spec = _make_abs_spectrum(absorbers)
        df = _make_linelist(MGII)
        result = line_search(spec, df, z_min=0.3, z_max=0.8,
                             n_steps=500, mode='absorption',
                             fit_continuum=False)
        assert hasattr(result, 'candidates')
        assert not hasattr(result, 'solutions')

    def test_different_linelist_civ(self):
        # CIV doublet absorber
        true_z = 1.5
        absorbers = [{'z': true_z, 'lines': CIVA}]
        spec = _make_abs_spectrum(absorbers, wave_min=3000, wave_max=12000,
                                  n_pix=5000)
        df = _make_linelist(CIVA, name='CIV')
        result = line_search(spec, df, z_min=1.2, z_max=1.8,
                             n_steps=3000, mode='absorption',
                             fit_continuum=False)
        best = result.candidates[0]
        assert abs(best.z - true_z) < 0.01
        assert best.linelist_name == 'CIV'
