""" Rebin 1D spectrum to new pixel scale."""

import math
import numpy as np

def rb_specbin(flux, nbin, **kwargs):
    """This function bins up 1D spectra in integer pixels. The routine returns a
       structure of flux and wavelength and variance that has been rebinned.
    
    Parameters
    -----------
     
        flux      - Input flux array
        nbin      - Number of pixels to bin on
        var=      - Input variance array [Optional]
        wave=     - Input wavelength array [Optional]
    
    Returns
    --------
        output    - Dictionary containing rebinned data with the following keys:
                    'flux': Rebinned flux array
                    'error': Rebinned error array (if var was provided)
                    'wave': Rebinned wavelength array (if wave was provided)
      
    Example
    --------
        output = rb_specbin(fx, 3)
        output = rb_specbin(fx, 3, var=var_array, wave=wave_array)
     
     
    REVISION HISTORY:
       Written by RB. June 2015
       Updated by RB. April 2025 - Improved error handling, documentation, and efficiency
    """
    # Input validation
    if not isinstance(flux, np.ndarray):
        flux = np.asarray(flux)
    
    if not isinstance(nbin, int) or nbin <= 0:
        raise ValueError("nbin must be a positive integer")
        
    if len(flux) == 0:
        raise ValueError("Input flux array is empty")
        
    if nbin > len(flux):
        raise ValueError(f"nbin ({nbin}) cannot be larger than the length of flux array ({len(flux)})")
    
    # Extract optional parameters
    has_variance = False
    has_wavelength = False
    variance = None
    wavelength = None
    
    if 'var' in kwargs:
        variance = kwargs['var']
        if not isinstance(variance, np.ndarray):
            variance = np.asarray(variance)
        if len(variance) != len(flux):
            raise ValueError("Length of variance array must match length of flux array")
        has_variance = True
        
    if 'wave' in kwargs:
        wavelength = kwargs['wave']
        if not isinstance(wavelength, np.ndarray):
            wavelength = np.asarray(wavelength)
        if len(wavelength) != len(flux):
            raise ValueError("Length of wavelength array must match length of flux array")
        has_wavelength = True

    # Calculate dimensions for the rebinned arrays
    wavePixels = len(flux)
    if (wavePixels % nbin) != 0:
        numPix = math.floor(wavePixels/nbin) + 1
        first = wavePixels % nbin
    else:
        numPix = math.floor(wavePixels/nbin)
        first = nbin
    
    # Pre-allocate arrays with correct size
    newFlux = np.zeros(numPix)
    newVar = np.zeros(numPix) if has_variance else None
    newWave = np.zeros(numPix) if has_wavelength else None
    
    # Binning
    for qq in range(numPix):
        start_idx = qq * nbin
        
        # Handle the last bin specially if it's smaller than nbin
        if qq == numPix - 1 and first != nbin:
            end_idx = start_idx + first
        else:
            end_idx = start_idx + nbin
        
        # Make sure we don't go out of bounds
        end_idx = min(end_idx, wavePixels)
        
        if start_idx >= wavePixels:
            break
            
        # Get the indices for this bin
        index = np.arange(start_idx, end_idx)
        
        # Compute mean values for this bin
        newFlux[qq] = np.mean(flux[index])
        
        if has_variance:
            newVar[qq] = np.mean(variance[index])
            
        if has_wavelength:
            newWave[qq] = np.mean(wavelength[index])
    
    # Prepare output dictionary
    output = {'flux': newFlux}
    
    # Add error (calculated from variance) if variance was provided
    if has_variance:
        output['error'] = np.sqrt(newVar/nbin)
        
    # Add wavelength if it was provided
    if has_wavelength:
        output['wave'] = newWave
    
    return output