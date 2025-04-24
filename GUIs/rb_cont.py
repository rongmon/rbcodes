#!/usr/bin/env python
"""
Interactive continuum fitter for 1D spectra.

This script loads a spectrum from a file and launches an interactive
continuum fitting interface. The user can select points by clicking on
the plot, and fit a cubic spline through these points to create a continuum.
The normalized spectrum and fitted continuum can be saved to a file.

Example usage:
    python rb_cont.py filename [filetype]

Where:
    filename : Path to the file containing the spectrum
    filetype : (Optional) Type of file - 'fits', 'ascii', 'p' (pickle), or 'xfits'
               If not provided, it will be inferred from the file extension.

Controls:
    Mouse:
        Left Click  : Select the median flux value within +/- 2.5 units from 
                      the x-coordinate for continuum fitting.
        Right Click : Delete the nearest continuum point.
    
    Keyboard:
        b     : Select a point for continuum fit at the exact (x,y) coordinate 
                of the cursor.
        enter : Perform a spline fit to create a continuum.
        n     : Show the normalized spectrum.
        w     : After pressing 'n', this will save the continuum.
        h     : Display the help screen.
        r     : Reset fit.
        q     : Quit the interactive session.
"""

import sys
import os
import argparse
import warnings
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# Import the interactive continuum fitter class
try:
    from rbcodes.GUIs.rb_fit_interactive_continuum import rb_fit_interactive_continuum
except ImportError:
    # If not available in rbcodes, try to import from local directory
    try:
        from rb_fit_interactive_continuum import rb_fit_interactive_continuum
    except ImportError:
        print("Error: rb_fit_interactive_continuum module not found.")
        print("Please make sure it's either in the current directory or in rbcodes.GUIs.")
        sys.exit(1)

def print_help():
    """Print help information about the interactive continuum fitter."""
    print("""
    ---------------------------------------------------------------------------
    Interactive Continuum Fitter for 1D Spectrum
    ---------------------------------------------------------------------------
    
    This program loads a spectrum from a file and allows interactive fitting
    of a continuum using cubic splines. The normalized spectrum and fitted
    continuum can be saved to a file.
    
    The program only works properly if none of the toolbar buttons in the figure
    is activated.
    
    Controls:
    
        Mouse Clicks:
        
            Left Click  : Select the median flux value within +/- 2.5 units from
                         the x-coordinate for continuum fitting.
            Right Click : Delete the nearest continuum point.
        
        Keyboard:
        
            b     : Select a point for continuum fit at the exact (x,y) coordinate
                   of the cursor.
            enter : Perform a spline fit to create a continuum.
            n     : Show the normalized spectrum.
            w     : After pressing 'n', this will save the continuum.
            h     : Display this help screen.
            r     : Reset fit.
            q     : Quit the interactive session.
    
    ---------------------------------------------------------------------------
    """)

def load_spectrum(filename, filetype=None):
    """
    Load a spectrum from a file.
    
    Parameters
    ----------
    filename : str
        Path to the file containing the spectrum.
    filetype : str, optional
        Type of file - 'fits', 'ascii', 'p' (pickle), or 'xfits'.
        If None, it will be inferred from the file extension.
    
    Returns
    -------
    tuple
        (wave, flux, error, metadata) where:
        - wave is the wavelength array
        - flux is the flux array
        - error is the error array
        - metadata is a dictionary with additional information
    
    Raises
    ------
    FileNotFoundError
        If the file doesn't exist.
    ValueError
        If the file type is not supported or the file doesn't contain valid data.
    ImportError
        If required modules are not available.
    """
    # Check if file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    # Infer file type if not provided
    if filetype is None:
        ext = os.path.splitext(filename)[1].lower()
        if ext in ['.txt', '.dat']:
            filetype = 'ascii'
        elif ext == '.fits':
            filetype = 'fits'
        elif ext in ['.p', '.pkl']:
            filetype = 'p'
        else:
            filetype = ext[1:] if ext.startswith('.') else ext
        print(f"Inferred file type: {filetype}")
    
    metadata = {'filename': filename, 'filetype': filetype}
    
    # Read the file based on its type
    try:
        if filetype == 'ascii':
            from astropy.io import ascii
            try:
                data = ascii.read(filename)
                if len(data) == 0:
                    raise ValueError(f"Empty data file: {filename}")
                
                keys = data.keys()
                if len(keys) < 2:
                    raise ValueError(f"Not enough columns in file: {filename}. Need at least wavelength and flux.")
                
                wave = np.array(data[keys[0]])
                flux = np.array(data[keys[1]])
                
                if len(keys) >= 3:
                    error = np.array(data[keys[2]])
                else:
                    print("Warning: No error column found. Assuming 10% error.")
                    error = 0.1 * np.abs(flux)
                
                metadata['format'] = 'ascii'
                metadata['columns'] = keys
            except Exception as e:
                raise ValueError(f"Error reading ASCII file: {str(e)}")
        
        elif filetype in ['fits', 'FITS']:
            # Try to use linetools if available
            try:
                from linetools.spectra.xspectrum1d import XSpectrum1D
                
                print(f"Loading {filename} with linetools.spectra.xspectrum1d...")
                sp = XSpectrum1D.from_file(filename)
                wave = sp.wavelength.value
                flux = sp.flux.value
                
                if sp.sig_is_set:
                    error = sp.sig.value
                else:
                    print("Warning: No error information in XSpectrum. Assuming 10% error.")
                    error = 0.1 * np.abs(flux)
                
                if sp.co_is_set:
                    print("Continuum found in file. Using it for normalization.")
                    continuum = sp.co.value
                    metadata['original_continuum'] = continuum
                
                metadata['format'] = 'linetools'
                
            except ImportError:
                # Fall back to astropy.io.fits if linetools is not available
                print("linetools not available. Using astropy.io.fits instead.")
                from astropy.io import fits
                
                try:
                    with fits.open(filename) as file:
                        if len(file) < 2:
                            raise ValueError(f"FITS file does not have expected HDU structure: {filename}")
                        
                        data = file[1].data
                        if data is None or len(data) == 0:
                            raise ValueError(f"Empty data in FITS file: {filename}")
                        
                        if 'WAVE' in data.names:
                            wave = np.array(data['WAVE'])
                        elif 'wave' in data.names:
                            wave = np.array(data['wave'])
                        else:
                            raise ValueError(f"Wavelength column not found in FITS file: {filename}")
                        
                        if 'FLUX' in data.names:
                            flux = np.array(data['FLUX'])
                        elif 'flux' in data.names:
                            flux = np.array(data['flux'])
                        else:
                            raise ValueError(f"Flux column not found in FITS file: {filename}")
                        
                        if 'ERROR' in data.names:
                            error = np.array(data['ERROR'])
                        elif 'error' in data.names:
                            error = np.array(data['error'])
                        elif 'ERR' in data.names:
                            error = np.array(data['ERR'])
                        elif 'err' in data.names:
                            error = np.array(data['err'])
                        else:
                            print("Warning: No error column found in FITS file. Assuming 10% error.")
                            error = 0.1 * np.abs(flux)
                        
                        metadata['format'] = 'fits'
                        metadata['columns'] = data.names
                except Exception as e:
                    raise ValueError(f"Error reading FITS file: {str(e)}")
        
        elif filetype == 'xfits':
            try:
                from linetools.spectra.xspectrum1d import XSpectrum1D
                
                print(f"Loading {filename} as XSpectrum1D...")
                sp = XSpectrum1D.from_file(filename)
                wave = sp.wavelength.value
                flux = sp.flux.value
                
                if sp.sig_is_set:
                    error = sp.sig.value
                else:
                    print("Warning: No error information in XSpectrum. Assuming 10% error.")
                    error = 0.1 * np.abs(flux)
                
                if sp.co_is_set:
                    print("Continuum found in file. Using it for normalization.")
                    continuum = sp.co.value
                    metadata['original_continuum'] = continuum
                
                metadata['format'] = 'xfits'
            except ImportError:
                raise ImportError("linetools is required to read xfits files. Please install linetools.")
            except Exception as e:
                raise ValueError(f"Error reading XFITS file: {str(e)}")
        
        elif filetype in ['p', 'pkl', 'pickle']:
            try:
                import pickle
                
                with open(filename, 'rb') as f:
                    data = pickle.load(f)
                
                if not isinstance(data, dict):
                    raise ValueError(f"Pickle file does not contain a dictionary: {filename}")
                
                if 'wave' not in data:
                    raise ValueError(f"'wave' key not found in pickle file: {filename}")
                if 'flux' not in data:
                    raise ValueError(f"'flux' key not found in pickle file: {filename}")
                
                wave = np.array(data['wave'])
                flux = np.array(data['flux'])
                
                if 'error' in data:
                    error = np.array(data['error'])
                else:
                    print("Warning: No 'error' key found in pickle file. Assuming 10% error.")
                    error = 0.1 * np.abs(flux)
                
                # Add any additional keys from the pickle file to metadata
                for key, value in data.items():
                    if key not in ['wave', 'flux', 'error']:
                        metadata[key] = value
                
                metadata['format'] = 'pickle'
            except Exception as e:
                raise ValueError(f"Error reading pickle file: {str(e)}")
        
        else:
            raise ValueError(f"Unsupported file type: {filetype}")
        
        # Validate the loaded arrays
        if len(wave) == 0 or len(flux) == 0:
            raise ValueError(f"Empty arrays loaded from file: {filename}")
        
        if len(wave) != len(flux) or len(wave) != len(error):
            raise ValueError(f"Mismatched array lengths: wave={len(wave)}, flux={len(flux)}, error={len(error)}")
        
        # Check for NaN or infinity values
        if np.any(np.isnan(wave)) or np.any(np.isinf(wave)):
            warnings.warn(f"Wavelength array contains {np.sum(np.isnan(wave))} NaN and {np.sum(np.isinf(wave))} infinite values")
        
        if np.any(np.isnan(flux)) or np.any(np.isinf(flux)):
            warnings.warn(f"Flux array contains {np.sum(np.isnan(flux))} NaN and {np.sum(np.isinf(flux))} infinite values")
            
            # Replace NaN/Inf values in flux with interpolated or median values
            bad_indices = np.isnan(flux) | np.isinf(flux)
            if np.any(bad_indices):
                good_indices = ~bad_indices
                if np.any(good_indices):
                    try:
                        from scipy.interpolate import interp1d
                        x_good = wave[good_indices]
                        y_good = flux[good_indices]
                        f = interp1d(x_good, y_good, bounds_error=False, fill_value=np.median(y_good))
                        flux[bad_indices] = f(wave[bad_indices])
                        print(f"Replaced {np.sum(bad_indices)} NaN/Inf flux values with interpolated values")
                    except ImportError:
                        # Fall back to using median if scipy is not available
                        flux[bad_indices] = np.median(flux[good_indices])
                        print(f"Replaced {np.sum(bad_indices)} NaN/Inf flux values with median value")
                else:
                    # All values are bad, use zeros
                    flux[:] = 0.0
                    print("All flux values are NaN/Inf. Replaced with zeros.")
        
        if np.any(np.isnan(error)) or np.any(np.isinf(error)) or np.any(error <= 0):
            warnings.warn(f"Error array contains {np.sum(np.isnan(error))} NaN and {np.sum(np.isinf(error))} infinite values")
            
            # Replace NaN/Inf/negative values in error with 10% of flux or median error
            bad_indices = np.isnan(error) | np.isinf(error) | (error <= 0)
            if np.any(bad_indices):
                good_indices = ~bad_indices
                if np.any(good_indices):
                    median_error = np.median(error[good_indices])
                    error[bad_indices] = median_error
                    print(f"Replaced {np.sum(bad_indices)} invalid error values with median error")
                else:
                    # All values are bad, use 10% of flux
                    error[:] = 0.1 * np.abs(flux)
                    print("All error values are invalid. Using 10% of flux as error.")
        
        return wave, flux, error, metadata
        
    except Exception as e:
        print(f"Error loading spectrum: {str(e)}")
        raise

def save_spectrum(wave, flux, error, continuum, filename, filetype, metadata=None):
    """
    Save a spectrum with its fitted continuum to a file.
    
    Parameters
    ----------
    wave : array-like
        Wavelength array.
    flux : array-like
        Flux array.
    error : array-like
        Error array.
    continuum : array-like
        Fitted continuum array.
    filename : str
        Path to the original input file.
    filetype : str
        Type of file to save - 'fits', 'ascii', or 'p' (pickle).
    metadata : dict, optional
        Additional metadata to include in the file.
    
    Returns
    -------
    str
        Path to the saved file.
    
    Raises
    ------
    ValueError
        If the file type is not supported or if there's an error saving the file.
    ImportError
        If required modules are not available.
    
    Notes
    -----
    The output file will be saved with the same name as the input file but with
    '_norm' appended before the extension. For example, if the input file is
    'spectrum.fits', the output file will be 'spectrum_norm.fits'.
    """
    if metadata is None:
        metadata = {}
    
    # Create output filename
    outfilename = os.path.splitext(filename)[0] + '_norm'
    if filetype == 'ascii':
        outfilename += '.txt'
    elif filetype in ['fits', 'FITS', 'xfits']:
        outfilename += '.fits'
    elif filetype in ['p', 'pkl', 'pickle']:
        outfilename += '.p'
    else:
        outfilename += '.' + filetype
    
    try:
        if filetype == 'ascii':
            from astropy.table import Table
            
            # Create table with wavelength, flux, error, and continuum
            table = Table([wave, flux, error, continuum, flux/continuum], 
                          names=['wave', 'flux', 'error', 'continuum', 'normalized_flux'])
            
            # Add metadata as table metadata
            for key, value in metadata.items():
                if isinstance(value, (str, int, float, bool)):
                    table.meta[key] = value
            
            # Save table to file
            from astropy.io import ascii as asc
            asc.write(table, outfilename, overwrite=True)
            
        elif filetype in ['fits', 'FITS', 'xfits']:
            # Try to use linetools if available
            try:
                from linetools.spectra.xspectrum1d import XSpectrum1D
                
                # Create XSpectrum1D object
                sp = XSpectrum1D.from_tuple((wave, flux, error, continuum))
                
                # Add metadata as header
                for key, value in metadata.items():
                    if isinstance(value, (str, int, float, bool)):
                        sp.meta[key] = value
                
                # Save to file
                sp.write_to_fits(outfilename)
                
            except ImportError:
                # Fall back to astropy.io.fits if linetools is not available
                from astropy.table import Table
                from astropy.io import fits
                
                # Create table with wavelength, flux, error, and continuum
                table = Table([wave, flux, error, continuum, flux/continuum], 
                              names=['wave', 'flux', 'error', 'continuum', 'normalized_flux'])
                
                # Add metadata as table metadata
                for key, value in metadata.items():
                    if isinstance(value, (str, int, float, bool)):
                        table.meta[key] = value
                
                # Save table to file
                table.write(outfilename, format='fits', overwrite=True)
                
        elif filetype in ['p', 'pkl', 'pickle']:
            import pickle
            
            # Create dictionary with data
            data = {
                'wave': wave,
                'flux': flux,
                'error': error,
                'continuum': continuum,
                'normalized_flux': flux/continuum
            }
            
            # Add metadata
            for key, value in metadata.items():
                data[key] = value
            
            # Save to file
            with open(outfilename, 'wb') as f:
                pickle.dump(data, f, protocol=2)
                
        else:
            raise ValueError(f"Unsupported file type for saving: {filetype}")
        
        print(f"Saved normalized spectrum to {outfilename}")
        return outfilename
        
    except Exception as e:
        print(f"Error saving spectrum: {str(e)}")
        return None

def main():
    """
    Main function to run the interactive continuum fitting script.
    
    This function parses command-line arguments, loads the spectrum,
    launches the interactive fitter, and saves the results.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Interactive continuum fitter for 1D spectra",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Controls:
            Mouse:
                Left Click  : Select median flux in window
                Right Click : Delete point
            
            Keyboard:
                b     : Select exact point
                enter : Fit continuum
                n     : Show normalized spectrum
                w     : Save continuum
                h     : Help screen
                r     : Reset fit
                q     : Quit
        """
    )
    
    parser.add_argument("filename", help="Path to the spectrum file")
    parser.add_argument("filetype", nargs="?", default=None,
                        help="Type of file (fits, ascii, p). If not provided, inferred from extension.")
    parser.add_argument("--help-controls", action="store_true",
                        help="Show detailed help about the controls")
    
    args = parser.parse_args()
    
    # Show help about controls if requested
    if args.help_controls:
        print_help()
        return
    
    try:
        # Load the spectrum
        print(f"Loading spectrum from {args.filename}...")
        wave, flux, error, metadata = load_spectrum(args.filename, args.filetype)
        
        print(f"Loaded spectrum with {len(wave)} points")
        print(f"Wavelength range: {wave.min():.2f} - {wave.max():.2f}")
        
        # Launch the interactive fitter
        print("\nLaunching interactive continuum fitter...")
        print("Press 'h' in the plot window for help on controls.")
        
        fitter = rb_fit_interactive_continuum(wave, flux, error)
        
        # Check if a continuum was fitted and saved
        if fitter.cont is not None:
            print("\nContinuum fitting successful!")
            
            # Save the result
            filetype = metadata.get('filetype', None)
            if filetype:
                save_spectrum(wave, flux, error, fitter.cont, args.filename, filetype, metadata)
            else:
                print("Warning: Could not determine file type for saving. Please specify.")
        else:
            print("\nNo continuum was fitted or saved.")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())