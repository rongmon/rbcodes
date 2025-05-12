"""
Spectral Analysis Pipeline Unit Tests

This test suite is designed to comprehensively validate the functionality 
of the spectral analysis pipeline, focusing on key astronomical spectroscopic 
data processing steps.

Test Coverage:
1. Spectrum Initialization
   - Verifies correct loading of spectral data
   - Checks basic attributes and data integrity

2. Spectrum Shifting
   - Tests conversion of observed wavelengths to rest frame
   - Validates redshift application

3. Spectrum Slicing
   - Examines ability to extract specific wavelength/velocity ranges
   - Ensures proper selection of spectral segments

4. Continuum Fitting
   - Validates different continuum fitting methods
   - Checks normalization quality
   - Handles real-world spectral variations

5. Equivalent Width Computation
   - Verifies column density and equivalent width calculations
   - Ensures physical constraints are met

6. Spectrum Saving
   - Tests data persistence 
   - Validates JSON and Pickle export capabilities

Key Parameters:
- Redshift: 0.348
- Transition Wavelength: 1037.5 Å
- Velocity Range: [-600, 600] km/s
- Masking Regions: Specific velocity windows to exclude from fitting

Note: These tests are designed to be robust and accommodate 
the natural variations found in astronomical spectral data.
python -m unittest -v rb_spec_unit_tests.py
"""

import os
import unittest
import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_filename

# Correct import paths
from rbcodes.GUIs import rb_spec as rs
from rbcodes.IGM import rb_setline as rsl
from rbcodes.IGM import compute_EW as cew

class TestSpectralAnalysisPipelineWithTestFits(unittest.TestCase):
    """
    Comprehensive unit tests for spectral analysis pipeline.
    
    This test suite exercises various methods of spectral data processing,
    ensuring robustness and accuracy of the analysis toolkit.
    """

    @classmethod
    def setUpClass(cls):
        """
        Prepare test environment:
        - Locate test FITS file
        - Set up spectroscopic analysis parameters
        - Load test spectrum
        """
        # Use resource_filename to locate the test.fits file
        cls.filename = resource_filename('rbcodes', 'example-data/test.fits')
        
        # Specific parameters for spectral analysis
        cls.zabs = 0.348  # Absorber redshift
        cls.transition = 1037.5  # Transition wavelength in Angstroms
        
        # Velocity limits for spectral analysis
        cls.xlim = [-600, 600]  # km/s
        
        # Masking regions to exclude from continuum fitting
        # These are velocity ranges that may contain absorption features
        # or other artifacts not suitable for continuum estimation
        cls.mask = [-560,-500,-330,-230,-170,120.,140.,290.]
        
        # Verify file exists
        if not os.path.exists(cls.filename):
            raise FileNotFoundError(f"Test FITS file not found: {cls.filename}")
        
        # Load the spectrum
        cls.spec = rs.rb_spec.from_file(cls.filename, filetype='xfits')
        print("\n--- Spectral Analysis Test Suite ---")
        print(f"Test File: {cls.filename}")
        print(f"Absorber Redshift: {cls.zabs}")
        print(f"Transition Wavelength: {cls.transition} Å")
        print(f"Velocity Range: {cls.xlim[0]} to {cls.xlim[1]} km/s")

    def test_spectrum_initialization(self):
        """
        Test Spectrum Initialization
        
        Validates:
        - Successful loading of spectral data
        - Presence of required attributes (wavelength, flux, error)
        - Basic data integrity checks
        """
        print("\n--- Testing Spectrum Initialization ---")
        
        # Check for essential spectral attributes
        self.assertTrue(hasattr(self.spec, 'wave'), 
            "Spectrum missing wavelength data")
        self.assertTrue(hasattr(self.spec, 'flux'), 
            "Spectrum missing flux data")
        self.assertTrue(hasattr(self.spec, 'error'), 
            "Spectrum missing error data")
        
        # Validate data lengths and basic properties
        self.assertGreater(len(self.spec.wave), 0, 
            "Wavelength array is empty")
        self.assertEqual(len(self.spec.wave), len(self.spec.flux), 
            "Wavelength and flux arrays have different lengths")
        self.assertEqual(len(self.spec.wave), len(self.spec.error), 
            "Wavelength and error arrays have different lengths")
        
        print("Spectrum initialization checks passed.")

    def test_spectrum_shift(self):
        """
        Test Spectrum Wavelength Shifting
        
        Validates:
        - Correct conversion to rest frame
        - Accurate redshift application
        """
        print("\n--- Testing Spectrum Rest Frame Shift ---")
        
        # Shift spectrum to rest frame
        wrest, shifted_zabs = self.spec.shift_spec(self.zabs)
        
        # Verify redshift application
        self.assertEqual(shifted_zabs, self.zabs, 
            "Incorrect redshift application")
        
        # Check wavelength transformation
        self.assertTrue(np.allclose(wrest, self.spec.wave / (1 + self.zabs)), 
            "Rest frame wavelength conversion incorrect")
        
        print(f"Spectrum successfully shifted to rest frame (z = {self.zabs})")

    def test_spectrum_slicing(self):
        """
        Test Spectrum Slicing
        
        Validates:
        - Ability to extract specific spectral segments
        - Correct velocity-based slicing
        """
        print("\n--- Testing Spectrum Slicing ---")
        
        # Shift spectrum first
        self.spec.shift_spec(self.zabs)
        
        # Slice spectrum around specific transition
        self.spec.slice_spec(
            lam_rest=self.transition,  # Transition wavelength
            lam_min=self.xlim[0],      # Minimum velocity 
            lam_max=self.xlim[1],      # Maximum velocity
            use_vel=True               # Use velocity space
        )
        
        # Verify slicing attributes
        self.assertTrue(hasattr(self.spec, 'wave_slice'), 
            "Sliced wavelength array missing")
        self.assertTrue(hasattr(self.spec, 'flux_slice'), 
            "Sliced flux array missing")
        self.assertTrue(hasattr(self.spec, 'velo'), 
            "Velocity array missing")
        
        # Check slice content
        self.assertGreater(len(self.spec.wave_slice), 0, 
            "Wavelength slice is empty")
        self.assertGreater(len(self.spec.flux_slice), 0, 
            "Flux slice is empty")
        self.assertGreater(len(self.spec.velo), 0, 
            "Velocity array is empty")
        
        print("Spectrum slicing successful")

    def test_continuum_fitting(self):
        """
        Test Continuum Fitting
        
        Validates:
        - Successful application of continuum fitting
        - Normalization quality
        - Handling of masked regions
        """
        print("\n--- Testing Continuum Fitting ---")
        
        # Prepare spectrum for fitting
        self.spec.shift_spec(self.zabs)
        self.spec.slice_spec(self.transition, 
                             lam_min=self.xlim[0], 
                             lam_max=self.xlim[1], 
                             use_vel=True)
        
        # Perform continuum fitting with specific parameters
        try:
            self.spec.fit_continuum(
                mask=self.mask,      # Regions to exclude
                domain=self.xlim,    # Velocity domain
                Legendre=3           # Polynomial order
            )
            
            # Verify continuum attributes
            self.assertTrue(hasattr(self.spec, 'cont'), 
                "Continuum fitting failed to produce continuum")
            self.assertTrue(hasattr(self.spec, 'fnorm'), 
                "Normalized flux missing")
            
            # Validate continuum properties
            self.assertTrue(np.all(self.spec.cont > 0), 
                "Continuum contains non-positive values")
            
            # Check normalization quality
            fnorm = self.spec.fnorm
            outliers = np.sum((fnorm < 0) | (fnorm > 3))
            total_points = len(fnorm)
            outlier_percentage = (outliers / total_points) * 100
            
            self.assertLess(outlier_percentage, 10, 
                f"Too many outliers in normalized spectrum: {outliers} out of {total_points} points")
            
            # Statistical check on normalization
            self.assertAlmostEqual(np.median(fnorm), 1.0, delta=0.5, 
                msg="Normalization deviates significantly from median=1")
            
            print("Continuum fitting successful")
            print(f"Outliers: {outliers}/{total_points} points ({outlier_percentage:.2f}%)")
            print(f"Median normalized flux: {np.median(fnorm):.3f}")
        
        except Exception as e:
            self.fail(f"Continuum fitting failed: {str(e)}")

    def test_equivalent_width_computation(self):
        """
        Test Equivalent Width Computation
        
        Validates:
        - Successful EW calculation
        - Physically meaningful column density
        - Proper error propagation from linear to log units
        """
        print("\n--- Testing Equivalent Width Computation ---")
        
        # Prepare spectrum
        self.spec.shift_spec(self.zabs)
        self.spec.slice_spec(self.transition, 
                             lam_min=self.xlim[0], 
                             lam_max=self.xlim[1], 
                             use_vel=True)
        
        # Perform continuum fitting
        self.spec.fit_continuum(
            mask=self.mask,
            domain=self.xlim,
            Legendre=3
        )
        
        # Compute equivalent width
        try:
            self.spec.compute_EW(
                lam_cen=self.transition, 
                vmin=-200, 
                vmax=200
            )
            
            # Verify EW attributes
            self.assertTrue(hasattr(self.spec, 'W'), 
                "Equivalent width not computed")
            self.assertTrue(hasattr(self.spec, 'W_e'), 
                "Equivalent width uncertainty not computed")
            
            # Verify column density attributes (both linear and logarithmic)
            self.assertTrue(hasattr(self.spec, 'N'), 
                "Linear column density not computed")
            self.assertTrue(hasattr(self.spec, 'N_e'), 
                "Linear column density uncertainty not computed")
            self.assertTrue(hasattr(self.spec, 'logN'), 
                "Log column density not computed")
            self.assertTrue(hasattr(self.spec, 'logN_e'), 
                "Log column density uncertainty not computed")
            
            # Validate physical constraints
            self.assertGreater(self.spec.W, 0, 
                "Equivalent width is non-positive")
            self.assertGreater(self.spec.N, 0, 
                "Linear column density is non-positive")
            self.assertGreater(self.spec.logN, 0, 
                "Log column density is non-positive")
            
            # Validate error propagation consistency
            # Check if logN and logN_e are consistent with N and N_e
            expected_logN = np.log10(self.spec.N) if self.spec.N > 0 else 0
            expected_logN_e = 0.434 * self.spec.N_e/self.spec.N if self.spec.N > 0 else 0
            
            self.assertAlmostEqual(self.spec.logN, expected_logN, places=5,
                msg="Log column density calculation inconsistent with linear value")
            self.assertAlmostEqual(self.spec.logN_e, expected_logN_e, places=5,
                msg="Log column density error propagation is inconsistent")
            
            print("Equivalent Width Computation Successful")
            print(f"Equivalent Width: {self.spec.W:.3f} Å")
            print(f"Column Density: N = {self.spec.N:.3e}, logN = {self.spec.logN:.3f}±{self.spec.logN_e:.3f}")
        
        except Exception as e:
            self.fail(f"Equivalent width computation failed: {str(e)}")
    def test_spectrum_saving(self):
        """
        Test Spectrum Saving Methods
        
        Validates:
        - JSON export functionality
        - Pickle export functionality
        """
        print("\n--- Testing Spectrum Saving Methods ---")
        
        # Prepare spectrum for saving
        self.spec.shift_spec(self.zabs)
        self.spec.slice_spec(self.transition, 
                             lam_min=self.xlim[0], 
                             lam_max=self.xlim[1], 
                             use_vel=True)
        
        # Perform continuum fitting
        self.spec.fit_continuum(
            mask=self.mask,
            domain=self.xlim,
            Legendre=3
        )
        self.spec.compute_EW(self.transition, vmin=-200, vmax=200)
        
        # Test saving methods
        import tempfile
        
        # Test JSON saving
        with tempfile.NamedTemporaryFile(delete=False, suffix='.json') as temp_json:
            try:
                self.spec.save_slice(temp_json.name, file_format='json')
                self.assertTrue(os.path.exists(temp_json.name), 
                    "JSON file not created")
                print(f"JSON saved successfully: {temp_json.name}")
            except Exception as e:
                self.fail(f"JSON saving failed: {str(e)}")
        
        # Test Pickle saving
        with tempfile.NamedTemporaryFile(delete=False, suffix='.p') as temp_pickle:
            try:
                self.spec.save_slice(temp_pickle.name, file_format='pickle')
                self.assertTrue(os.path.exists(temp_pickle.name), 
                    "Pickle file not created")
                print(f"Pickle saved successfully: {temp_pickle.name}")
            except Exception as e:
                self.fail(f"Pickle saving failed: {str(e)}")

def run_tests():
    """Run the unit tests with detailed output."""
    unittest.main(argv=[''], exit=False)

if __name__ == '__main__':
    run_tests()