"""
Tests for the compute_EW module comparing new and old implementations.

This module tests the compute_EW functions from the rbcodes.IGM package,
ensuring that the new implementation is compatible with the old implementation
while also verifying the improved robustness and functionality.

# From the project root directory
python -m unittest tests.test_compute_EW

# Or for more verbose output
python -m unittest -v tests.test_compute_EW

"""

import unittest
import numpy as np
from rbcodes.IGM.compute_EW import compute_EW
from rbcodes.IGM.compute_EW_old import compute_EW as compute_EW_old

class TestComputeEW(unittest.TestCase):
    """Test case for compute_EW module."""
    
    def setUp(self):
        """Set up test fixtures for each test."""
        # Create test data
        self.wavelength = np.linspace(1200, 1400, 1000)
        
        # Standard spectrum with absorption feature
        self.flux_standard = np.ones_like(self.wavelength)
        center_idx = 500
        width = 50
        self.flux_standard[center_idx-width:center_idx+width] = 0.5
        self.error_standard = np.ones_like(self.wavelength) * 0.05
        
        # Common parameters
        self.wrest = 1300.0
        self.lmts = [-100, 100]  # Velocity limits in km/s
        self.zabs = 0.0
        self.f0 = 0.1  # Oscillator strength
        
    def test_standard_case(self):
        """Test that both functions produce similar results for standard input."""
        # Run both implementations
        result_new = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        result_old = compute_EW_old(
            self.wavelength, self.flux_standard, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Compare equivalent width values with appropriate tolerance
        self.assertAlmostEqual(
            result_new['ew_tot'], result_old['ew_tot'], 
            delta=abs(result_old['ew_tot'] * 0.05),  # 5% tolerance
            msg="Equivalent width values differ significantly"
        )
        
        # Compare equivalent width errors
        self.assertAlmostEqual(
            result_new['err_ew_tot'], result_old['err_ew_tot'], 
            delta=abs(result_old['err_ew_tot'] * 0.1),  # 10% tolerance for errors
            msg="Equivalent width error values differ significantly"
        )
        
        # Check column density if available (in log space)
        if 'col' in result_new and 'col' in result_old:
            # Handle potential zero or negative values
            log_col_new = np.log10(max(result_new['col'], 1e-30))
            log_col_old = np.log10(max(result_old['col'], 1e-30))
            
            self.assertAlmostEqual(
                log_col_new, log_col_old, 
                delta=0.1,  # 0.1 dex tolerance (common in astronomy)
                msg=f"Column density values differ significantly: {log_col_new} vs {log_col_old}"
            )
            
            # Compare column density errors
            # Using relative error comparison for column density errors
            if result_old['col'] > 0 and result_old['colerr'] > 0:
                rel_err_new = result_new['colerr'] / result_new['col']
                rel_err_old = result_old['colerr'] / result_old['col']
                
                self.assertAlmostEqual(
                    rel_err_new, rel_err_old,
                    delta=0.2,  # 20% tolerance for relative errors
                    msg=f"Column density relative errors differ significantly: {rel_err_new} vs {rel_err_old}"
                )
        
        # Compare velocity dispersion and centroid
        if 'vel_disp' in result_new and 'vel_disp' in result_old:
            self.assertAlmostEqual(
                result_new['vel_disp'], result_old['vel_disp'],
                delta=abs(result_old['vel_disp'] * 0.15),  # 15% tolerance
                msg="Velocity dispersion values differ significantly"
            )
        
        if 'med_vel' in result_new and 'med_vel' in result_old:
            self.assertAlmostEqual(
                result_new['med_vel'], result_old['med_vel'],
                delta=abs(10.0),  # 10 km/s tolerance
                msg="Velocity centroid values differ significantly"
            )    
    def test_with_nan_values(self):
        """Test new implementation's handling of NaN values compared to old implementation."""
        # Create a spectrum with some NaN values
        flux_with_nans = self.flux_standard.copy()
        flux_with_nans[100:120] = np.nan
        
        # New implementation should handle NaNs
        result_new = compute_EW(
            self.wavelength, flux_with_nans, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Verify result is not NaN
        self.assertFalse(
            np.isnan(result_new['ew_tot']), 
            "New implementation returned NaN for equivalent width with NaN in flux"
        )
        
        # Try with old implementation for comparison
        try:
            result_old = compute_EW_old(
                self.wavelength, flux_with_nans, self.wrest, self.lmts, 
                self.error_standard, f0=self.f0, zabs=self.zabs
            )
            # If the old implementation handles it, check consistency
            if not np.isnan(result_old['ew_tot']):
                self.assertAlmostEqual(
                    result_new['ew_tot'], result_old['ew_tot'], 
                    delta=abs(result_old['ew_tot'] * 0.1),  # 10% tolerance
                    msg="Results differ for NaN handling"
                )
        except Exception as e:
            # It's okay if the old implementation fails - just log it
            print(f"Old implementation failed with NaN values: {str(e)}")
    
    def test_saturated_flux(self):
        """Test handling of saturated flux values."""
        # Create a spectrum with saturated values (zero flux)
        flux_saturated = self.flux_standard.copy()
        flux_saturated[480:520] = 0.0  # Saturated in the middle of the absorption
        
        # Both implementations should handle this
        result_new = compute_EW(
            self.wavelength, flux_saturated, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs, sat_limit=0.1
        )
        
        # Verify the new implementation produces finite values
        self.assertTrue(
            np.isfinite(result_new['ew_tot']),
            "New implementation failed to handle saturated flux"
        )
        
        # Try with old implementation
        try:
            result_old = compute_EW_old(
                self.wavelength, flux_saturated, self.wrest, self.lmts, 
                self.error_standard, f0=self.f0, zabs=self.zabs, sat_limit=0.1
            )
            
            # Compare if old implementation also handled it
            if np.isfinite(result_old['ew_tot']):
                # Results should be reasonably similar but may differ due to improved handling
                self.assertAlmostEqual(
                    result_new['ew_tot'], result_old['ew_tot'], 
                    delta=abs(result_old['ew_tot'] * 0.15),  # 15% tolerance for saturated case
                    msg="Results differ significantly for saturated flux"
                )
        except Exception as e:
            # It's okay if the old implementation fails - just log it
            print(f"Old implementation failed with saturated values: {str(e)}")
    
    def test_normalization_option(self):
        """Test the new normalization option in the updated implementation."""
        # Create a spectrum that needs normalization
        flux_unnormalized = self.flux_standard.copy() * 2.0  # Scale by factor of 2
        
        # New implementation with median normalization
        result_new_normalized = compute_EW(
            self.wavelength, flux_unnormalized, self.wrest, self.lmts, 
            self.error_standard * 2.0, f0=self.f0, zabs=self.zabs,
            normalization='median'
        )
        
        # New implementation without normalization
        result_new_unnormalized = compute_EW(
            self.wavelength, flux_unnormalized, self.wrest, self.lmts, 
            self.error_standard * 2.0, f0=self.f0, zabs=self.zabs,
            normalization='none'
        )
        
        # The normalized result should be close to the original standard case
        result_standard = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Normalized result should match standard result
        self.assertAlmostEqual(
            result_new_normalized['ew_tot'], result_standard['ew_tot'], 
            delta=abs(result_standard['ew_tot'] * 0.05),  # 5% tolerance
            msg="Normalized result doesn't match expected value"
        )
        
        # Unnormalized result should differ (approximately by the scale factor)
        self.assertNotAlmostEqual(
            result_new_unnormalized['ew_tot'], result_standard['ew_tot'], 
            delta=abs(result_standard['ew_tot'] * 0.5),  # Should differ by >50%
            msg="Unnormalized result should differ substantially from standard"
        )
    
    def test_empty_velocity_range(self):
        """Test behavior when no pixels fall within the velocity range."""
        # Set velocity range outside the spectrum
        out_of_range_lmts = [1e6, 1.5e6]
        
        # New implementation should handle this gracefully
        result_new = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, out_of_range_lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Print actual result for debugging
        print(f"Result for empty range: {result_new['ew_tot']}, isNaN: {np.isnan(result_new['ew_tot'])}")
        
        # Verify that the result is actually NaN
        self.assertTrue(
            np.isnan(result_new['ew_tot']),
            f"Should return NaN for equivalent width when no pixels in velocity range. Got: {result_new['ew_tot']}"
        )        
        # Try with old implementation
        try:
            result_old = compute_EW_old(
                self.wavelength, self.flux_standard, self.wrest, out_of_range_lmts, 
                self.error_standard, f0=self.f0, zabs=self.zabs
            )
            # Record behavior - may be NaN or error
            print(f"Old implementation with empty range returned: {result_old['ew_tot']}")
        except Exception as e:
            # It's okay if the old implementation fails - just log it
            print(f"Old implementation failed with empty range: {str(e)}")


if __name__ == '__main__':
    unittest.main()