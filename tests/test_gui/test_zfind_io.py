"""
Phase 1 tests for rbcodes/GUIs/zfind/io.py

Tests:
  - to_rb_spectrum(): all input forms
  - ZFindResult / ZSolution dataclasses
  - AbsorberResult / AbsorberCandidate dataclasses
"""

import unittest
import numpy as np
import os
import math

from rbcodes.GUIs.zfind.io import (
    to_rb_spectrum,
    ZFindResult, ZSolution,
    AbsorberResult, AbsorberCandidate,
)
from rbcodes.utils.rb_spectrum import rb_spectrum

# path to a known example FITS file
EXAMPLE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', '..', 'example-data'
)
TEST_FITS = os.path.join(EXAMPLE_DIR, 'test.fits')     # FLUX/ERROR/WAVELENGTH/CONTINUUM
SDSS_FITS = os.path.join(EXAMPLE_DIR, 'sdss1.fits')   # COADD binary table


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_arrays(n=200):
    """Simple synthetic wave/flux/error arrays."""
    wave = np.linspace(4000.0, 9000.0, n)
    flux = np.random.default_rng(42).normal(1.0, 0.1, n).astype(np.float32)
    error = np.full(n, 0.1, dtype=np.float32)
    return wave, flux, error


def _make_rb_spectrum():
    wave, flux, error = _make_arrays()
    return rb_spectrum.from_tuple((wave, flux, error))


# ---------------------------------------------------------------------------
# to_rb_spectrum tests
# ---------------------------------------------------------------------------

class TestToRbSpectrum(unittest.TestCase):

    # --- filepath inputs ---

    def test_from_filepath_test_fits(self):
        """Load from a known FITS filepath."""
        if not os.path.exists(TEST_FITS):
            self.skipTest(f"test.fits not found at {TEST_FITS}")
        spec = to_rb_spectrum(TEST_FITS)
        self.assertIsInstance(spec, rb_spectrum)
        self.assertFalse(spec._read_failed)
        self.assertGreater(len(spec.flux), 0)
        self.assertTrue(spec.sig_is_set)
        self.assertTrue(spec.co_is_set)   # test.fits has CONTINUUM

    def test_from_filepath_sdss_fits(self):
        """Load from an SDSS-format FITS filepath."""
        if not os.path.exists(SDSS_FITS):
            self.skipTest(f"sdss1.fits not found at {SDSS_FITS}")
        spec = to_rb_spectrum(SDSS_FITS)
        self.assertIsInstance(spec, rb_spectrum)
        self.assertFalse(spec._read_failed)
        self.assertGreater(len(spec.flux), 0)

    def test_from_bad_filepath(self):
        """Non-existent filepath returns a failed rb_spectrum, not exception."""
        spec = to_rb_spectrum('/nonexistent/path/fake.fits')
        self.assertIsInstance(spec, rb_spectrum)
        self.assertTrue(spec._read_failed)

    # --- rb_spectrum passthrough ---

    def test_from_rb_spectrum(self):
        """Passing an rb_spectrum returns it unchanged."""
        orig = _make_rb_spectrum()
        result = to_rb_spectrum(orig)
        self.assertIs(result, orig)

    # --- tuple inputs ---

    def test_from_tuple_2(self):
        """(wave, flux) tuple — no error."""
        wave, flux, _ = _make_arrays()
        spec = to_rb_spectrum((wave, flux))
        self.assertIsInstance(spec, rb_spectrum)
        self.assertFalse(spec._read_failed)
        self.assertFalse(spec.sig_is_set)
        np.testing.assert_array_almost_equal(spec.flux.value, flux)

    def test_from_tuple_3(self):
        """(wave, flux, error) tuple."""
        wave, flux, error = _make_arrays()
        spec = to_rb_spectrum((wave, flux, error))
        self.assertIsInstance(spec, rb_spectrum)
        self.assertTrue(spec.sig_is_set)
        np.testing.assert_array_almost_equal(spec.sig.value, error)

    def test_from_tuple_4(self):
        """(wave, flux, error, continuum) tuple."""
        wave, flux, error = _make_arrays()
        cont = np.ones_like(flux)
        spec = to_rb_spectrum((wave, flux, error, cont))
        self.assertIsInstance(spec, rb_spectrum)
        self.assertTrue(spec.co_is_set)
        np.testing.assert_array_almost_equal(spec.co.value, cont)

    def test_from_list_tuple(self):
        """List form is also accepted."""
        wave, flux, error = _make_arrays()
        spec = to_rb_spectrum([wave, flux, error])
        self.assertIsInstance(spec, rb_spectrum)
        self.assertFalse(spec._read_failed)

    # --- numpy array inputs ---

    def test_from_arrays_wave_flux(self):
        """wave ndarray + flux keyword."""
        wave, flux, _ = _make_arrays()
        spec = to_rb_spectrum(wave, flux=flux)
        self.assertIsInstance(spec, rb_spectrum)
        self.assertFalse(spec.sig_is_set)

    def test_from_arrays_wave_flux_error(self):
        """wave ndarray + flux + error keywords."""
        wave, flux, error = _make_arrays()
        spec = to_rb_spectrum(wave, flux=flux, error=error)
        self.assertIsInstance(spec, rb_spectrum)
        self.assertTrue(spec.sig_is_set)

    def test_from_arrays_wave_flux_error_continuum(self):
        """wave ndarray + flux + error + continuum keywords."""
        wave, flux, error = _make_arrays()
        cont = np.ones_like(flux)
        spec = to_rb_spectrum(wave, flux=flux, error=error, continuum=cont)
        self.assertIsInstance(spec, rb_spectrum)
        self.assertTrue(spec.co_is_set)

    def test_from_arrays_missing_flux_raises(self):
        """ndarray wave with no flux should raise ValueError."""
        wave, _, _ = _make_arrays()
        with self.assertRaises(ValueError):
            to_rb_spectrum(wave)

    # --- type errors ---

    def test_unsupported_type_raises(self):
        """Completely unsupported type should raise TypeError."""
        with self.assertRaises(TypeError):
            to_rb_spectrum(12345)

    def test_none_raises(self):
        """None input should raise TypeError."""
        with self.assertRaises(TypeError):
            to_rb_spectrum(None)

    # --- wavelength range preserved ---

    def test_wavelength_range_preserved(self):
        """Wave array values survive round-trip through to_rb_spectrum."""
        wave, flux, error = _make_arrays()
        spec = to_rb_spectrum(wave, flux=flux, error=error)
        np.testing.assert_array_almost_equal(
            spec.wavelength.value, np.sort(wave)
        )

    # --- air/vac metadata default ---

    def test_default_airvac_is_vac(self):
        """rb_spectrum created from arrays should default to vacuum."""
        wave, flux, _ = _make_arrays()
        spec = to_rb_spectrum((wave, flux))
        self.assertEqual(spec.meta.get('airvac', 'vac'), 'vac')


# ---------------------------------------------------------------------------
# ZSolution tests
# ---------------------------------------------------------------------------

class TestZSolution(unittest.TestCase):

    def _make(self, z=1.2345, z_err=0.0003, chi2=1.12, method='LineSearch',
              ttype='Galaxy', n=5):
        return ZSolution(z=z, z_err=z_err, chi2_dof=chi2,
                         method=method, template_type=ttype, n_features=n)

    def test_creation(self):
        s = self._make()
        self.assertAlmostEqual(s.z, 1.2345)
        self.assertAlmostEqual(s.z_err, 0.0003)
        self.assertAlmostEqual(s.chi2_dof, 1.12)
        self.assertEqual(s.method, 'LineSearch')
        self.assertEqual(s.template_type, 'Galaxy')
        self.assertEqual(s.n_features, 5)

    def test_repr(self):
        s = self._make()
        r = repr(s)
        self.assertIn('1.2345', r)
        self.assertIn('LineSearch', r)

    def test_nan_z_err(self):
        """z_err can be NaN (curvature estimation may fail)."""
        s = self._make(z_err=float('nan'))
        self.assertTrue(math.isnan(s.z_err))


# ---------------------------------------------------------------------------
# ZFindResult tests
# ---------------------------------------------------------------------------

class TestZFindResult(unittest.TestCase):

    def _make(self):
        z_array = np.linspace(0, 3, 1000)
        chi2 = np.random.default_rng(0).uniform(1, 5, 1000)
        chi2[400] = 0.8   # plant a minimum at z≈1.2
        chi2[700] = 0.9   # and another at z≈2.1
        solutions = [
            ZSolution(z=1.2, z_err=0.001, chi2_dof=0.8,
                      method='LineSearch', template_type='Galaxy', n_features=4),
            ZSolution(z=2.1, z_err=0.002, chi2_dof=0.9,
                      method='LineSearch', template_type='QSO', n_features=3),
        ]
        spec = _make_rb_spectrum()
        return ZFindResult(
            z_array=z_array,
            chi2_curves=[{'label': 'LineSearch:ISM', 'chi2': chi2}],
            solutions=solutions,
            input_spec=spec,
            warnings=['No error array — using MAD-STD IVAR.']
        )

    def test_creation(self):
        r = self._make()
        self.assertEqual(len(r.solutions), 2)
        self.assertEqual(len(r.chi2_curves), 1)
        self.assertEqual(len(r.z_array), 1000)

    def test_best(self):
        r = self._make()
        best = r.best()
        self.assertIsNotNone(best)
        self.assertAlmostEqual(best.z, 1.2)
        self.assertAlmostEqual(best.chi2_dof, 0.8)

    def test_best_empty(self):
        """best() returns None when no solutions."""
        spec = _make_rb_spectrum()
        r = ZFindResult(
            z_array=np.linspace(0, 3, 100),
            chi2_curves=[],
            solutions=[],
            input_spec=spec
        )
        self.assertIsNone(r.best())

    def test_repr(self):
        r = self._make()
        rep = repr(r)
        self.assertIn('2 solutions', rep)
        self.assertIn('1.2', rep)

    def test_warnings_stored(self):
        r = self._make()
        self.assertEqual(len(r.warnings), 1)
        self.assertIn('MAD-STD', r.warnings[0])

    def test_input_spec_preserved(self):
        """The input rb_spectrum is stored and accessible."""
        r = self._make()
        self.assertIsInstance(r.input_spec, rb_spectrum)
        self.assertFalse(r.input_spec._read_failed)

    def test_multiple_chi2_curves(self):
        """Multiple engine modes produce multiple curves."""
        z_array = np.linspace(0, 3, 100)
        spec = _make_rb_spectrum()
        r = ZFindResult(
            z_array=z_array,
            chi2_curves=[
                {'label': 'LineSearch:ISM', 'chi2': np.ones(100)},
                {'label': 'Template:Passive', 'chi2': np.ones(100) * 1.2},
            ],
            solutions=[],
            input_spec=spec
        )
        self.assertEqual(len(r.chi2_curves), 2)
        self.assertEqual(r.chi2_curves[0]['label'], 'LineSearch:ISM')
        self.assertEqual(r.chi2_curves[1]['label'], 'Template:Passive')


# ---------------------------------------------------------------------------
# AbsorberCandidate tests
# ---------------------------------------------------------------------------

class TestAbsorberCandidate(unittest.TestCase):

    def _make(self, z=0.512, sig=5.2, doublet=True):
        return AbsorberCandidate(
            z=z, significance=sig, n_lines=2,
            is_doublet=doublet, linelist_name='MgII',
            lines_matched=['MgII 2796', 'MgII 2803']
        )

    def test_creation(self):
        c = self._make()
        self.assertAlmostEqual(c.z, 0.512)
        self.assertAlmostEqual(c.significance, 5.2)
        self.assertEqual(c.n_lines, 2)
        self.assertTrue(c.is_doublet)
        self.assertEqual(c.linelist_name, 'MgII')
        self.assertEqual(len(c.lines_matched), 2)

    def test_single_line(self):
        c = self._make(doublet=False)
        self.assertFalse(c.is_doublet)

    def test_repr(self):
        c = self._make()
        r = repr(c)
        self.assertIn('0.512', r)
        self.assertIn('MgII', r)
        self.assertIn('doublet', r)


# ---------------------------------------------------------------------------
# AbsorberResult tests
# ---------------------------------------------------------------------------

class TestAbsorberResult(unittest.TestCase):

    def _make(self):
        z_array = np.linspace(0, 3, 1000)
        sig_curve = np.random.default_rng(1).normal(0, 1, 1000)
        sig_curve[300] = 5.2   # absorber at z≈0.9
        sig_curve[600] = 4.1   # absorber at z≈1.8
        candidates = [
            AbsorberCandidate(z=0.9, significance=5.2, n_lines=2,
                              is_doublet=True, linelist_name='MgII',
                              lines_matched=['MgII 2796', 'MgII 2803']),
            AbsorberCandidate(z=1.8, significance=4.1, n_lines=2,
                              is_doublet=True, linelist_name='CIV',
                              lines_matched=['CIV 1548', 'CIV 1550']),
        ]
        spec = _make_rb_spectrum()
        return AbsorberResult(
            z_array=z_array,
            significance_curve=sig_curve,
            candidates=candidates,
            input_spec=spec,
            warnings=[]
        )

    def test_creation(self):
        r = self._make()
        self.assertEqual(len(r.candidates), 2)
        self.assertEqual(len(r.z_array), 1000)

    def test_candidates_accessible(self):
        r = self._make()
        self.assertAlmostEqual(r.candidates[0].z, 0.9)
        self.assertEqual(r.candidates[0].linelist_name, 'MgII')
        self.assertAlmostEqual(r.candidates[1].z, 1.8)
        self.assertEqual(r.candidates[1].linelist_name, 'CIV')

    def test_repr(self):
        r = self._make()
        self.assertIn('2 candidates', repr(r))

    def test_empty_candidates(self):
        spec = _make_rb_spectrum()
        r = AbsorberResult(
            z_array=np.linspace(0, 3, 100),
            significance_curve=np.zeros(100),
            candidates=[],
            input_spec=spec
        )
        self.assertEqual(len(r.candidates), 0)


# ---------------------------------------------------------------------------
# Adapter smoke tests (quick check adapters read from ZFindResult correctly)
# ---------------------------------------------------------------------------

class TestAdapters(unittest.TestCase):

    def test_zfind_to_zgui_z(self):
        from rbcodes.GUIs.zfind.adapters import zfind_to_zgui_z
        solutions = [ZSolution(z=1.234, z_err=0.001, chi2_dof=1.1,
                               method='LineSearch', template_type='Galaxy', n_features=3)]
        spec = _make_rb_spectrum()
        result = ZFindResult(z_array=np.linspace(0,3,100),
                             chi2_curves=[], solutions=solutions, input_spec=spec)
        out = zfind_to_zgui_z(result)
        self.assertEqual(len(out), 2)
        self.assertAlmostEqual(out[0], 1.234)
        self.assertAlmostEqual(out[1], 0.001)

    def test_zfind_to_zgui_z_empty(self):
        from rbcodes.GUIs.zfind.adapters import zfind_to_zgui_z
        spec = _make_rb_spectrum()
        result = ZFindResult(z_array=np.linspace(0,3,100),
                             chi2_curves=[], solutions=[], input_spec=spec)
        out = zfind_to_zgui_z(result)
        self.assertTrue(math.isnan(out[0]))
        self.assertTrue(math.isnan(out[1]))

    def test_absorbers_to_multispec(self):
        from rbcodes.GUIs.zfind.adapters import absorbers_to_multispec
        candidates = [
            AbsorberCandidate(z=0.5, significance=5.0, n_lines=2, is_doublet=True,
                              linelist_name='MgII', lines_matched=['MgII 2796', 'MgII 2803']),
            AbsorberCandidate(z=1.2, significance=3.0, n_lines=2, is_doublet=True,
                              linelist_name='CIV', lines_matched=['CIV 1548', 'CIV 1550']),
        ]
        spec = _make_rb_spectrum()
        result = AbsorberResult(z_array=np.linspace(0,3,100),
                                significance_curve=np.zeros(100),
                                candidates=candidates, input_spec=spec)
        out = absorbers_to_multispec(result, accepted_indices=[0])
        self.assertEqual(len(out), 1)
        self.assertAlmostEqual(out[0]['zabs'], 0.5)
        self.assertEqual(out[0]['name'], 'MgII')
        self.assertIn('z=0.5', out[0]['label'])

    def test_zfind_to_multispec_z(self):
        from rbcodes.GUIs.zfind.adapters import zfind_to_multispec_z
        solutions = [ZSolution(z=2.567, z_err=0.002, chi2_dof=1.0,
                               method='PCA:Galaxy', template_type='Galaxy', n_features=10)]
        spec = _make_rb_spectrum()
        result = ZFindResult(z_array=np.linspace(0,3,100),
                             chi2_curves=[], solutions=solutions, input_spec=spec)
        self.assertAlmostEqual(zfind_to_multispec_z(result), 2.567)


if __name__ == '__main__':
    unittest.main()
