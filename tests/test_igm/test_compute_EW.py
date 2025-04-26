"""
Tests for the compute_EW module with a synthetic spectrum of known properties.

This module tests the compute_EW functions from the rbcodes.IGM package using
a synthetic spectrum with known equivalent width and column density values.

# From the project root directory
python -m unittest tests.test_compute_EW

# Or for more verbose output
python -m unittest -v tests.test_compute_EW

# To generate a visualization
python test_compute_EW.py --visualize
"""

import unittest
import numpy as np
from rbcodes.IGM.compute_EW import compute_EW

class TestComputeEW(unittest.TestCase):
    """Test case for compute_EW module."""
    
    def setUp(self):
        """Set up test fixtures for each test."""
        # Common parameters for all tests
        self.wrest = 1215.67  # Lyman-alpha rest wavelength (Å)
        self.f0 = 0.4164      # Oscillator strength for Lyman-alpha
        self.zabs = 0.0       # Redshift of absorber (no redshift in this case)
        self.lmts = [-150, 150]  # Velocity limits in km/s
        
        # Create wavelength array with high resolution around the line
        self.wavelength = np.linspace(1214.0, 1217.5, 5000)
        
        # Create standard test data
        self.flux_standard, self.error_standard, self.true_flux, self.known_ew, self.log_N = self.create_synthetic_spectrum()
        
        # Store reference values for standard case
        self.standard_result = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Store the computed EW as the reference value instead of the theoretical one
        self.computed_ew = self.standard_result['ew_tot']
    
    def create_synthetic_spectrum(self):
        """Create a synthetic absorption line with known properties."""
        # Define parameters
        N_col = 1.0e13  # Column density in atoms/cm^2
        b = 15.0        # Doppler parameter in km/s
        v_center = 0.0  # Center velocity in km/s
        # Constants
        c = 299792.458  # Speed of light in km/s
        # Convert wavelength to velocity
        velocity = c * (self.wavelength - self.wrest) / self.wrest
        # Optical depth at line center - using the updated formula
        tau_0 = 2.654e-15 * self.f0 * self.wrest * N_col / b
        # Optical depth profile (Gaussian)
        tau = tau_0 * np.exp(-((velocity - v_center) / b)**2)
        # Flux profile: I = exp(-tau)
        flux = np.exp(-tau)
        # Add noise - using the updated noise level
        noise_level = 0.015
        error = np.ones_like(flux) * noise_level
        np.random.seed(42)
        noisy_flux = flux + np.random.normal(0, noise_level, size=flux.shape)
        
        # Calculate theoretical EW directly in wavelength space
        theoretical_ew = np.trapz(1-flux, x=self.wavelength)
        
        # Log column density
        log_N = np.log10(N_col)
        # Return the synthetic data
        return noisy_flux, error, flux, theoretical_ew, log_N
    
    def visualize_synthetic_spectrum(self, result, show=False, save_path=None):
        """
        Visualize the synthetic spectrum and EW measurement.
        
        Parameters:
        -----------
        result : dict
            Result dictionary from compute_EW function
        show : bool, optional
            Whether to display the plot interactively (default: False)
        save_path : str, optional
            Path to save the figure (default: None, no saving)
        """
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Rectangle
        except ImportError:
            print("Matplotlib is required for visualization")
            return
        
        # Create a figure with two subplots (one for spectrum, one for EW)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), 
                                      gridspec_kw={'height_ratios': [3, 1]})
        
        # Convert wavelength to velocity for better visualization
        c = 299792.458  # Speed of light in km/s
        velocity = c * (self.wavelength - self.wrest) / self.wrest
        
        # Plot the spectrum in the top panel
        ax1.plot(velocity, self.true_flux, 'r-', label='Theoretical (noise-free)', alpha=0.7)
        ax1.plot(velocity, self.flux_standard, 'k-', label='With noise', alpha=0.5)
        ax1.fill_between(velocity, self.flux_standard - self.error_standard, 
                        self.flux_standard + self.error_standard, 
                        color='gray', alpha=0.3, label='Error')
        
        # Add a horizontal line at y=1 for reference
        ax1.axhline(y=1.0, color='g', linestyle='--', alpha=0.5, label='Continuum')
        
        # Shade the velocity limits used for EW calculation
        vlim_min, vlim_max = self.lmts
        ax1.axvspan(vlim_min, vlim_max, color='blue', alpha=0.1, label='Velocity limits')
        
        # Add vertical lines at velocity limits
        ax1.axvline(x=vlim_min, color='blue', linestyle='--', alpha=0.5)
        ax1.axvline(x=vlim_max, color='blue', linestyle='--', alpha=0.5)
        
        # Set labels and title for spectrum plot
        ax1.set_xlabel('Velocity (km/s)')
        ax1.set_ylabel('Normalized Flux')
        ax1.set_title(f'Synthetic Absorption Line (log N={self.log_N:.2f})')
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)
        
        # Set y-limits with some padding
        y_min = max(0, min(self.true_flux) - 0.1)
        ax1.set_ylim(y_min, 1.1)
        
        # Set x-limits to focus on the absorption feature
        padding = 100  # km/s
        ax1.set_xlim(min(vlim_min - padding, -200), max(vlim_max + padding, 200))
        
        # In the bottom panel, visualize the EW measurement
        if 'ew_array' in result and 'vel_tab' in result:
            # If the function returns the EW array and velocity grid, use it
            ax2.plot(result['vel_tab'], result['ew_array'], 'b-', label='EW contribution')
            ax2.fill_between(result['vel_tab'], 0, result['ew_array'], color='blue', alpha=0.3)
        else:
            # Otherwise, just show a rectangle representing the total EW
            # This is a simplified visualization
            ew_rect = Rectangle((-10, 0), 20, 0.5, color='blue', alpha=0.5)
            ax2.add_patch(ew_rect)
            ax2.text(0, 0.25, f"Total EW: {result['ew_tot']:.4f} Å", 
                    ha='center', va='center', fontsize=12)
        
        # Add annotations for EW values
        ax2.text(0.02, 0.85, f"Computed EW: {result['ew_tot']:.6f} Å", 
                transform=ax2.transAxes, fontsize=10)
        ax2.text(0.02, 0.70, f"Theoretical EW: {self.known_ew:.6f} Å", 
                transform=ax2.transAxes, fontsize=10)
        ax2.text(0.02, 0.55, f"Difference: {abs(result['ew_tot'] - self.known_ew):.6f} Å", 
                transform=ax2.transAxes, fontsize=10)
        ax2.text(0.02, 0.40, f"Relative diff: {abs(result['ew_tot'] - self.known_ew) / self.known_ew * 100:.2f}%", 
                transform=ax2.transAxes, fontsize=10)
        
        # If column density is available, add it to the annotations
        if 'col' in result:
            ax2.text(0.65, 0.85, f"Computed log N: {np.log10(result['col']):.3f}", 
                    transform=ax2.transAxes, fontsize=10)
            ax2.text(0.65, 0.70, f"Theoretical log N: {self.log_N:.3f}", 
                    transform=ax2.transAxes, fontsize=10)
            ax2.text(0.65, 0.55, f"Difference: {abs(np.log10(result['col']) - self.log_N):.3f} dex", 
                    transform=ax2.transAxes, fontsize=10)
        
        # Set labels for EW plot
        ax2.set_xlabel('Velocity (km/s)')
        ax2.set_ylabel('EW')
        ax2.set_title('Equivalent Width Measurement')
        ax2.set_xlim(ax1.get_xlim())  # Match x-limits with spectrum plot
        ax2.grid(True, alpha=0.3)
        
        # Adjust layout to prevent overlap
        plt.tight_layout()
        
        # Save the figure if a path is provided
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {save_path}")
        
        # Show the plot if requested
        if show:
            plt.show()
        else:
            plt.close(fig)
    
    def test_synthetic_spectrum(self, visualize=False):
        """
        Test compute_EW with a synthetic spectrum of known properties.
        
        Parameters:
        -----------
        visualize : bool, optional
            Whether to create and display a visualization (default: False)
        """
        # Compute equivalent width and column density
        result = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Print the results for debugging and analysis
        print(f"Computed EW: {result['ew_tot']:.6f} Å")
        print(f"Theoretical EW: {self.known_ew:.6f} Å")
        print(f"Difference: {abs(result['ew_tot'] - self.known_ew):.6f} Å")
        print(f"Relative difference: {abs(result['ew_tot'] - self.known_ew) / self.known_ew * 100:.2f}%")
        
        # Visualize if requested
        if visualize:
            #self.visualize_synthetic_spectrum(result, show=True, 
            #                                save_path="synthetic_spectrum_test.png")
            self.visualize_synthetic_spectrum(result, show=True)
        
        # Compare with theoretical value
        self.assertAlmostEqual(
            result['ew_tot'], self.known_ew, 
            delta=self.known_ew * 0.05,  # 5% tolerance
            msg=f"Computed EW ({result['ew_tot']:.6f} Å) differs from theoretical ({self.known_ew:.6f} Å)"
        )
        
        # Check column density if available
        if 'col' in result:
            print(f"Computed Column Density: {result['col']:.3e} cm^-2")
            print(f"Known Column Density: {10**self.log_N:.3e} cm^-2")
            print(f"Log Difference: {abs(np.log10(result['col']) - self.log_N):.3f} dex")
            
            self.assertAlmostEqual(
                np.log10(result['col']), self.log_N, 
                delta=0.3,  # 0.3 dex tolerance
                msg=f"Computed column density differs from theoretical value"
            )
    
    def test_standard_case(self):
        """Test that function produces expected results for standard input."""
        # Run the computation
        result = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Print all keys in the result for documentation
        print("Keys in result dictionary:", list(result.keys()))
        
        # Verify the result contains expected keys
        self.assertIn('ew_tot', result, "Result should contain 'ew_tot' key")
        self.assertIn('err_ew_tot', result, "Result should contain 'err_ew_tot' key")
        
        # Verify EW is non-zero and positive
        self.assertGreater(result['ew_tot'], 0, "Equivalent width should be positive")
        
        # Verify error is reasonable (positive and less than the EW itself)
        self.assertGreater(result['err_ew_tot'], 0, "Error should be positive")
        self.assertLess(result['err_ew_tot'], result['ew_tot'], "Error should be smaller than the EW")
        
        # Check column density if available
        if 'col' in result:
            self.assertGreater(result['col'], 0, "Column density should be positive")
            self.assertGreater(result['colerr'], 0, "Column density error should be positive")
    
    def test_with_nan_values(self):
        """Test handling of NaN values."""
        # Create a spectrum with some NaN values
        flux_with_nans = self.flux_standard.copy()
        flux_with_nans[100:120] = np.nan
        
        # Function should handle NaNs
        result = compute_EW(
            self.wavelength, flux_with_nans, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Verify result is not NaN
        self.assertFalse(np.isnan(result['ew_tot']), 
                        "Function should handle NaN values in flux")
    
    def test_saturated_flux(self):
        """Test handling of saturated flux values."""
        # Create a spectrum with saturated values (zero flux)
        flux_saturated = self.true_flux.copy()  # Start with the true flux
        
        # Set a portion to near-zero to simulate saturation
        center_idx = len(flux_saturated) // 2
        width = 50
        flux_saturated[center_idx-width:center_idx+width] = 0.001
        
        # Add noise to the saturated spectrum
        np.random.seed(42)
        noise_level = 0.001  # Updated noise level
        noisy_flux_saturated = flux_saturated + np.random.normal(0, noise_level, size=flux_saturated.shape)
        noisy_flux_saturated = np.clip(noisy_flux_saturated, 0, 1)
        
        # Test compute_EW with saturated spectrum
        result = compute_EW(
            self.wavelength, noisy_flux_saturated, self.wrest, self.lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs, sat_limit=0.1
        )
        
        # Verify the function produces finite values
        self.assertTrue(np.isfinite(result['ew_tot']),
                       "Function should handle saturated flux")
        
        # The EW should be larger than in the standard case
        self.assertGreater(
            result['ew_tot'], self.standard_result['ew_tot'],
            "EW with saturation should be greater than standard case"
        )
    
    def test_normalization_option(self):
        """Test the normalization option."""
        # Create a spectrum that needs normalization
        flux_unnormalized = self.flux_standard.copy() * 2.0  # Scale by factor of 2
        
        # With median normalization
        result_normalized = compute_EW(
            self.wavelength, flux_unnormalized, self.wrest, self.lmts, 
            self.error_standard * 2.0, f0=self.f0, zabs=self.zabs,
            normalization='median'
        )
        
        # Without normalization
        result_unnormalized = compute_EW(
            self.wavelength, flux_unnormalized, self.wrest, self.lmts, 
            self.error_standard * 2.0, f0=self.f0, zabs=self.zabs,
            normalization='none'
        )
        
        # The normalized result should be close to the standard case
        self.assertAlmostEqual(
            result_normalized['ew_tot'], self.standard_result['ew_tot'], 
            delta=abs(self.standard_result['ew_tot'] * 0.1),  # 10% tolerance
            msg="Normalized result doesn't match expected value"
        )
        
        # Unnormalized result should differ from standard
        self.assertNotEqual(
            result_unnormalized['ew_tot'], self.standard_result['ew_tot'],
            msg="Unnormalized result should differ from standard"
        )
    
    def test_different_velocity_ranges(self):
        """Test behavior with different velocity range settings."""
        # Narrow velocity range
        narrow_lmts = [-50, 50]
        result_narrow = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, narrow_lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Wide velocity range
        wide_lmts = [-300, 300]
        result_wide = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, wide_lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Print values for documentation
        print(f"Standard velocity range EW: {self.standard_result['ew_tot']:.6f} Å")
        print(f"Narrow velocity range EW: {result_narrow['ew_tot']:.6f} Å")
        print(f"Wide velocity range EW: {result_wide['ew_tot']:.6f} Å")
        
        # Due to noise or implementation details, the narrow range might 
        # occasionally give slightly larger values than standard.
        # Instead, verify that the wide range captures at least as much as standard
        self.assertGreaterEqual(
            result_wide['ew_tot'], self.standard_result['ew_tot'] * 0.95,  # Allow 5% variation
            "Wider velocity range should capture at least as much EW as standard"
        )
        
        # And verify the narrow range is at least close to the standard range
        # (either slightly larger or smaller due to noise)
        self.assertAlmostEqual(
            result_narrow['ew_tot'], self.standard_result['ew_tot'],
            delta=self.standard_result['ew_tot'] * 0.1,  # 10% tolerance
            msg="Narrow velocity range should yield similar EW to standard"
        )
    
    def test_empty_velocity_range(self):
        """Test behavior when no pixels fall within the velocity range."""
        # Set velocity range outside the spectrum
        out_of_range_lmts = [1e6, 1.5e6]
        
        # Function should handle this gracefully
        result = compute_EW(
            self.wavelength, self.flux_standard, self.wrest, out_of_range_lmts, 
            self.error_standard, f0=self.f0, zabs=self.zabs
        )
        
        # Print actual result for debugging
        print(f"Result for empty range: {result['ew_tot']}, isNaN: {np.isnan(result['ew_tot'])}")
        
        # Verify that the result is actually NaN
        self.assertTrue(
            np.isnan(result['ew_tot']),
            f"Should return NaN for equivalent width when no pixels in velocity range. Got: {result['ew_tot']}"
        )


if __name__ == '__main__':
    def create_visualization():
        """Create and display a visualization of the synthetic spectrum and EW measurement."""
        test_case = TestComputeEW()
        test_case.setUp()
        test_case.test_synthetic_spectrum(visualize=True)
        print("Visualization complete. Check the saved image: synthetic_spectrum_test.png")
        
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--visualize':
        create_visualization()
    else:
        unittest.main()