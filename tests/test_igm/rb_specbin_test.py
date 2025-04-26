"""
Test module for rb_specbin function in rbcodes.IGM.rb_specbin
"""

import unittest
import numpy as np
from numpy.testing import assert_allclose

# Import the function to test
from rbcodes.IGM.rb_specbin import rb_specbin

class TestRbSpecbin(unittest.TestCase):
    """Test cases for rb_specbin function"""

    def setUp(self):
        """Set up test data"""
        # Create a simple flux array for testing
        self.flux = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
        
        # Create corresponding wavelength and variance arrays
        self.wavelength = np.arange(1000, 1000 + len(self.flux))
        self.variance = np.ones_like(self.flux) * 0.1

    def test_basic_binning(self):
        """Test basic binning functionality without optional parameters."""
        nbin = 3
        result = rb_specbin(self.flux, nbin)
        
        # Expected result: mean of each bin
        expected_flux = np.array([2.0, 5.0, 8.0])
        
        # Check keys in output dictionary
        self.assertIn('flux', result)
        self.assertNotIn('wave', result)
        self.assertNotIn('error', result)
        
        # Check values
        assert_allclose(result['flux'], expected_flux)

    def test_binning_with_variance(self):
        """Test binning with variance parameter."""
        nbin = 3
        result = rb_specbin(self.flux, nbin, var=self.variance)
        
        # Expected results
        expected_flux = np.array([2.0, 5.0, 8.0])
        expected_error = np.sqrt(0.1 / nbin) * np.ones(3)
        
        # Check keys
        self.assertIn('flux', result)
        self.assertIn('error', result)
        self.assertNotIn('wave', result)
        
        # Check values
        assert_allclose(result['flux'], expected_flux)
        assert_allclose(result['error'], expected_error)

    def test_binning_with_wavelength(self):
        """Test binning with wavelength parameter."""
        nbin = 3
        result = rb_specbin(self.flux, nbin, wave=self.wavelength)
        
        # Expected results
        expected_flux = np.array([2.0, 5.0, 8.0])
        expected_wave = np.array([1001.0, 1004.0, 1007.0])
        
        # Check keys
        self.assertIn('flux', result)
        self.assertIn('wave', result)
        self.assertNotIn('error', result)
        
        # Check values
        assert_allclose(result['flux'], expected_flux)
        assert_allclose(result['wave'], expected_wave)

    def test_binning_with_all_parameters(self):
        """Test binning with both variance and wavelength parameters."""
        nbin = 3
        result = rb_specbin(self.flux, nbin, var=self.variance, wave=self.wavelength)
        
        # Expected results
        expected_flux = np.array([2.0, 5.0, 8.0])
        expected_error = np.sqrt(0.1 / nbin) * np.ones(3)
        expected_wave = np.array([1001.0, 1004.0, 1007.0])
        
        # Check keys
        self.assertIn('flux', result)
        self.assertIn('wave', result)
        self.assertIn('error', result)
        
        # Check values
        assert_allclose(result['flux'], expected_flux)
        assert_allclose(result['error'], expected_error)
        assert_allclose(result['wave'], expected_wave)

    def test_non_integer_bin_size(self):
        """Test handling of non-integer bin size."""
        with self.assertRaises(ValueError):
            rb_specbin(self.flux, 2.5)

    def test_negative_bin_size(self):
        """Test handling of negative bin size."""
        with self.assertRaises(ValueError):
            rb_specbin(self.flux, -3)

    def test_bin_size_too_large(self):
        """Test handling of bin size larger than array."""
        with self.assertRaises(ValueError):
            rb_specbin(self.flux, len(self.flux) + 1)

    def test_uneven_division(self):
        """Test binning when array length is not divisible by bin size."""
        # Array of length 10, bin size 3 => 4 bins (3+3+3+1)
        flux = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        result = rb_specbin(flux, 3)
        
        expected_flux = np.array([2.0, 5.0, 8.0, 10.0])
        assert_allclose(result['flux'], expected_flux)

    def test_non_numpy_input(self):
        """Test with Python list inputs instead of numpy arrays."""
        flux_list = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        wavelength_list = [1000.0, 1001.0, 1002.0, 1003.0, 1004.0, 1005.0]
        variance_list = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
        
        result = rb_specbin(flux_list, 2, var=variance_list, wave=wavelength_list)
        
        expected_flux = np.array([1.5, 3.5, 5.5])
        expected_error = np.sqrt(0.1 / 2) * np.ones(3)
        expected_wave = np.array([1000.5, 1002.5, 1004.5])
        
        assert_allclose(result['flux'], expected_flux)
        assert_allclose(result['error'], expected_error)
        assert_allclose(result['wave'], expected_wave)

    def test_edge_case_bin_size_one(self):
        """Test with bin size of 1 (should return original arrays)."""
        result = rb_specbin(self.flux, 1, var=self.variance, wave=self.wavelength)
        
        assert_allclose(result['flux'], self.flux)
        assert_allclose(result['error'], np.sqrt(self.variance))
        assert_allclose(result['wave'], self.wavelength)

    def test_mismatch_length_arrays(self):
        """Test error handling when arrays have different lengths."""
        short_wave = np.arange(1000, 1000 + len(self.flux) - 1)
        
        with self.assertRaises(ValueError):
            rb_specbin(self.flux, 3, wave=short_wave)
        
        short_var = np.ones(len(self.flux) - 1) * 0.1
        
        with self.assertRaises(ValueError):
            rb_specbin(self.flux, 3, var=short_var)


if __name__ == '__main__':
    unittest.main()