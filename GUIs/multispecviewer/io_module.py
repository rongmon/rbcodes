import os
import numpy as np
from linetools.spectra.xspectrum1d import XSpectrum1D
from astropy.io import fits

def read_fits_files(file_paths):
    """
    Reads multiple FITS files using linetools.XSpectrum1D and stores data in a dictionary.
    
    :param file_paths: List of FITS file paths
    :return: Dictionary with filenames as keys and wave, flux, error arrays as values
    """
    spectra_dict = {}
    
    for file_path in file_paths:
        try:
            spec = XSpectrum1D.from_file(file_path)
            wave = spec.wavelength.value  # Extract wavelength array
            flux = spec.flux.value        # Extract flux array
            error = spec.sig.value if spec.sig is not None else np.zeros_like(flux)  # Extract error array
            
            spectra_dict[os.path.basename(file_path)] = {'wave': wave, 'flux': flux, 'error': error}
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    
    return spectra_dict
