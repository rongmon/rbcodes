#!/usr/bin/env python
"""
Spectrum Viewer GUI: A tool for visualizing and analyzing astronomical spectra.

This script provides a command-line interface to load and display various 
spectrum file formats, with options to normalize and analyze spectral data.
"""

import argparse
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from pathlib import Path
import logging

from linetools.spectra.xspectrum1d import XSpectrum1D
from rbcodes.gui import PlotSpec_Integrated as r

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('rb_specgui')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Spectrum Viewer GUI for astronomical data analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('filename', help='Spectrum file to load')
    
    parser.add_argument(
        '-t', '--filetype', 
        choices=['ascii', 'linetools', 'fits', 'lt', 'lt_cont_norm', 'p'],
        help='File format type. If not specified, will be guessed from file extension'
    )
    
    parser.add_argument(
        '-e', '--error-file', 
        dest='efil',
        help='Optional separate error file'
    )
    
    parser.add_argument(
        '-n', '--normalize',
        action='store_true',
        help='Apply continuum normalization if available'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '-z', '--redshift',
        type=float,
        default=0.0,
        help='Initial redshift (zabs) value'
    )
    
    return parser.parse_args()


def guess_filetype(filename):
    """
    Guess the file type based on file extension.
    
    Parameters
    ----------
    filename : str
        Path to the spectrum file
    
    Returns
    -------
    str
        Guessed file type identifier
    """
    extension = os.path.splitext(filename)[1].lower()
    
    # Map file extensions to file types
    extension_map = {
        '.txt': 'ascii',
        '.dat': 'ascii',
        '.fits': 'linetools',
        '.fit': 'linetools',
        '.p': 'p',
        '.pkl': 'p'
    }
    
    return extension_map.get(extension, extension[1:])


def load_spectrum(filename, filetype, efil=None, normalize=False):
    """
    Load a spectrum file based on its format.
    
    Parameters
    ----------
    filename : str
        Path to the spectrum file
    filetype : str
        Format of the spectrum file
    efil : str, optional
        Path to error file if separate
    normalize : bool, optional
        Whether to apply continuum normalization
    
    Returns
    -------
    tuple
        (wavelength, flux, error) arrays
    """
    cwd = os.getcwd()
    full_path = os.path.join(cwd, filename)
    logger.info(f"Loading spectrum from {full_path}")
    
    # ASCII format
    if filetype == 'ascii':
        from astropy.io import ascii
        dat = ascii.read(full_path)
        columns = dat.keys()
        wave = np.array(dat[columns[0]])
        flux = np.array(dat[columns[1]])
        error = np.array(dat[columns[2]]) if len(columns) >= 3 else 0.1 * flux
    
    # LineTools format
    elif filetype in ('fits', 'lt', 'linetools', 'lt_cont_norm'):
        try:
            sp = XSpectrum1D.from_file(filename, efil=efil) if efil else XSpectrum1D.from_file(filename)
            wave = sp.wavelength.value
            
            # Handle continuum normalization if requested
            if filetype == 'lt_cont_norm' or normalize:
                if sp.co_is_set:
                    logger.info('Using provided continuum for normalization')
                    flux = sp.flux.value / sp.co.value
                    error = sp.sig.value / sp.co.value if sp.sig_is_set else 0.1 * flux
                else:
                    logger.warning('No continuum found, using raw flux')
                    flux = sp.flux.value
                    error = sp.sig.value if sp.sig_is_set else 0.1 * flux
            else:
                flux = sp.flux.value
                error = sp.sig.value if sp.sig_is_set else 0.1 * flux
                
            if not sp.sig_is_set:
                logger.warning('No error array found, assuming 10% of flux')
        except Exception as e:
            logger.error(f"Error loading LineTools spectrum: {e}")
            raise
    
    # Pickle format
    elif filetype == 'p':
        import pickle
        with open(full_path, "rb") as f:
            dat = pickle.load(f)
        columns = dat.keys()
        wave = np.array(dat['wave'])
        flux = np.array(dat['flux'])
        error = np.array(dat['error']) if 'error' in columns else 0.1 * flux
    
    else:
        raise ValueError(f"Unsupported file type: {filetype}")
    
    return wave, flux, error


def normalize_flux(flux, error):
    """
    Normalize flux by its median value.
    
    Parameters
    ----------
    flux : numpy.ndarray
        Flux values
    error : numpy.ndarray
        Error values
    
    Returns
    -------
    tuple
        (normalization_constant, normalized_flux, normalized_error)
    """
    norm_factor = np.nanmedian(flux)
    
    # Handle special cases
    if norm_factor == 0:
        norm_factor = np.nanmean(flux)
        
    if (norm_factor == 0) and (np.nanmean(flux) == 0):
        norm_factor = 1.0
        
    return norm_factor, flux / norm_factor, error / norm_factor


def main():
    """Main function to run the spectrum GUI."""
    args = parse_arguments()
    
    # Set logging level based on verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Determine file type
    filetype = args.filetype
    if not filetype:
        filetype = guess_filetype(args.filename)
        logger.info(f"File type not specified, guessed as: {filetype}")
    
    try:
        # Load spectrum data
        wave, flux, error = load_spectrum(
            args.filename, 
            filetype, 
            efil=args.efil, 
            normalize=args.normalize
        )
        
        # Normalize flux
        norm_factor, norm_flux, norm_error = normalize_flux(flux, error)
        
        # Start the GUI with the specified redshift
        logger.info(f"Starting GUI application with initial redshift: {args.redshift}")
        r.rb_plotspec(wave, norm_flux, norm_error, zabs=args.redshift)
        
    except Exception as e:
        logger.error(f"Error running spectrum GUI: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())