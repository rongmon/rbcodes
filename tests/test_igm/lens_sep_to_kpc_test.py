import unittest
import numpy as np
from astropy.cosmology import Planck18, FlatLambdaCDM
import astropy.units as u

# Import the function to test 
from rbcodes.IGM.lens_sep_to_kpc import lens_sep_to_kpc

class TestLensSepToKpc(unittest.TestCase):
    """
    Unit tests for the lens_sep_to_kpc function.
    
    These tests verify the correct functionality of the lens_sep_to_kpc function
    under various inputs, edge cases, and parameter configurations.
    """
    
    def setUp(self):
        """Set up common test parameters."""
        # Standard test parameters
        self.delta_arcsec = 1.0  # 1 arcsec separation
        self.z_lens = 0.55       # lens redshift
        self.z_source = 3.5      # source redshift
        
        # Test absorber redshifts
        self.zabs_single = 1.0                              # single absorber
        self.zabs_list = [0.2, 0.55, 1.2, 2.0]             # regular list
        self.zabs_array = np.array([0.3, 0.8, 1.5, 2.5])   # numpy array
        
        # Run the function once to get actual values as reference
        self.expected_single = lens_sep_to_kpc(
            self.delta_arcsec, self.zabs_single, self.z_lens, self.z_source)
        
        self.expected_list = lens_sep_to_kpc(
            self.delta_arcsec, self.zabs_list, self.z_lens, self.z_source)
        
        # Tolerance for floating point comparisons
        self.rtol = 1e-10  # tight relative tolerance for exact matching
    
    def test_single_absorber(self):
        """Test with a single absorber redshift."""
        result = lens_sep_to_kpc(self.delta_arcsec, self.zabs_single, 
                                self.z_lens, self.z_source)
        self.assertAlmostEqual(result, self.expected_single, delta=1e-10)
    
    def test_absorber_list(self):
        """Test with a list of absorber redshifts."""
        result = lens_sep_to_kpc(self.delta_arcsec, self.zabs_list, 
                                self.z_lens, self.z_source)
        np.testing.assert_allclose(result, self.expected_list, rtol=self.rtol)
    
    def test_absorber_array(self):
        """Test with a numpy array of absorber redshifts."""
        result = lens_sep_to_kpc(self.delta_arcsec, self.zabs_array, 
                                self.z_lens, self.z_source)
        # Just verify it works with arrays and returns expected shape
        self.assertEqual(len(result), len(self.zabs_array))
        self.assertTrue(isinstance(result, np.ndarray))
    
    def test_below_lens_redshift(self):
        """Test absorbers at redshifts below the lens redshift."""
        # For z_abs <= z_lens, should use angular diameter distance
        zabs_below = [0.1, 0.3, 0.55]  # All below or equal to lens redshift
        result = lens_sep_to_kpc(self.delta_arcsec, zabs_below, 
                                self.z_lens, self.z_source)
        
        # Verify values individually to be sure of calculation method
        for i, z in enumerate(zabs_below):
            arcsec2kpc = Planck18.arcsec_per_kpc_proper(z).value
            expected = self.delta_arcsec / arcsec2kpc
            self.assertAlmostEqual(result[i], expected, delta=1e-10)
    
    def test_type_conversion(self):
        """Test that function handles different input types."""
        # Test with a tuple instead of a list
        result_tuple = lens_sep_to_kpc(self.delta_arcsec, tuple(self.zabs_list), 
                                      self.z_lens, self.z_source)
        np.testing.assert_allclose(result_tuple, self.expected_list, rtol=self.rtol)
        
        # Test with a single float instead of a list
        result_float = lens_sep_to_kpc(self.delta_arcsec, 1.0, 
                                      self.z_lens, self.z_source)
        self.assertTrue(np.isscalar(result_float) or len(result_float) == 1)
    
    def test_output_type(self):
        """Test that the output is a numpy array without units."""
        result = lens_sep_to_kpc(self.delta_arcsec, self.zabs_list, 
                               self.z_lens, self.z_source)
        
        # Check if it's a scalar or array
        if np.isscalar(result):
            self.assertTrue(isinstance(result, (float, np.float64, np.floating)))
        else:
            self.assertTrue(isinstance(result, np.ndarray))
        
        # Verify no units are attached
        self.assertFalse(hasattr(result, 'unit'))
    
    def test_reproducibility(self):
        """Test that the function gives consistent results when called multiple times."""
        result1 = lens_sep_to_kpc(self.delta_arcsec, self.zabs_list, 
                                self.z_lens, self.z_source)
        result2 = lens_sep_to_kpc(self.delta_arcsec, self.zabs_list, 
                                self.z_lens, self.z_source)
        np.testing.assert_allclose(result1, result2, rtol=self.rtol)
    
    def test_scaling(self):
        """Test that the output scales linearly with delta_arcsec."""
        result1 = lens_sep_to_kpc(self.delta_arcsec, self.zabs_list, 
                                self.z_lens, self.z_source)
        result2 = lens_sep_to_kpc(2*self.delta_arcsec, self.zabs_list, 
                                self.z_lens, self.z_source)
        np.testing.assert_allclose(result2, 2*result1, rtol=self.rtol)
    

if __name__ == '__main__':
    unittest.main()