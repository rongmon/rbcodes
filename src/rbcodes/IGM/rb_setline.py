"""
Read in atomic line information for a given or approximate rest frame wavelength.

This module provides functionality to find atomic line information based on different
matching methods including exact wavelength matching, closest wavelength matching,
or matching by species name.
"""

from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from astropy.io import ascii
from pkg_resources import resource_filename
import os
import logging
from pathlib import Path
from typing import Dict, List, Union, Optional, Any
import warnings

# Set up logging
logger = logging.getLogger(__name__)

# Global cache for line lists to avoid repeated file loading
_LINE_LIST_CACHE = {}


def rb_setline(lambda_rest: float, 
               method: str, 
               linelist: str = 'atom', 
               target_name: Optional[str] = None) -> Dict[str, Any]:
    """
    Function to read in atomic line information for a given rest frame wavelength,
    for the line matching the closest wavelength, or by name.

    Parameters
    ----------
    lambda_rest : float
        Rest Frame wavelength (in Å) of the line to match
    method : str
        'closest' -> If set will match the closest line.
        'Exact'   -> If set will match the exact wavelength.
        'Name'    -> Match by name, USE WITH CARE. MUST INPUT OPTIONAL target_name
    linelist : str, optional
        The line list to use. Default is 'atom'.
        Available options: 'atom', 'LLS', 'LLS Small', 'DLA', 'LBG', 'Gal',
        'Eiger_Strong', 'Gal_Em', 'Gal_Abs', 'Gal_long', 'AGN', 'HI_recomb',
        'HI_recomb_light'
    target_name : str, optional
        Required when method='Name'. The name of the target species to match.

    Returns
    -------
    dict
        Dictionary with the following keys:
        - 'wave': Rest frame wavelength(s) of the matched line(s)
        - 'fval': Oscillator strength value(s)
        - 'name': Species name(s)
        - 'gamma': Radiation damping constant (if available in the line list)

    Examples
    --------
    Match the closest line to a given wavelength:
    
    >>> result = rb_setline(2796.3, 'closest')
    >>> print(f"Found: {result['name']} at {result['wave']} Å")
    
    Match a line by exact wavelength:
    
    >>> result = rb_setline(1215.67, 'Exact')
    
    Match a line by name:
    
    >>> result = rb_setline(0, 'Name', target_name='HI 1215')

    Raises
    ------
    ValueError
        If method is not one of 'closest', 'Exact', or 'Name'
        If method='Name' but target_name is not provided
        If method='Exact' but no lines match the wavelength
    FileNotFoundError
        If the requested line list file cannot be found

    Notes
    -----
    - For 'Exact' method, a match is considered if the wavelength is within 0.001 Å
    - The 'Name' method requires the target_name parameter to be set
    - The full list of available line lists can be found in the read_line_list function
    - Line lists are cached to improve performance when called multiple times

    History
    -------
    Written By: Rongmon Bordoloi                Jan 2018, Python 2.7
    Edit:       Rongmon Bordoloi                Sep 2018, Depreciated kwargs to be compatible with python 3
    Edit:       Improved version                Apr 2025, Various improvements while maintaining compatibility
    """
    # Validate inputs
    if method not in ['Exact', 'closest', 'Name']:
        raise ValueError(f"Method must be one of 'Exact', 'closest', or 'Name', got '{method}'")
    
    if method == 'Name' and target_name is None:
        raise ValueError("target_name must be provided when method='Name'")
    
    # Get the line list data
    try:
        line_str = read_line_list(linelist)
    except Exception as e:
        logger.error(f"Error loading line list '{linelist}': {str(e)}")
        raise
    
    if not line_str:
        logger.warning(f"Line list '{linelist}' is empty or could not be loaded")
        return {'wave': np.array([]), 'fval': np.array([]), 'name': np.array([])}
    
    # Extract data from the line list
    wavelist = np.zeros((len(line_str),))
    name = np.empty(len(line_str), dtype='object')
    fval = np.zeros((len(line_str),))
    
    if linelist == 'atom':
        gamma = np.zeros((len(line_str),))

    # Populate arrays from line list data
    for i in range(len(wavelist)):
        wavelist[i] = np.double(line_str[i]['wrest'])
        fval[i] = float(line_str[i]['fval'])
        name[i] = str(line_str[i]['ion'])
        if linelist == 'atom':
            gamma[i] = float(line_str[i]['gamma'])

    # Match based on the requested method
    if method == 'Exact':
        # Look for exact wavelength matches within a small tolerance
        q = np.where(np.abs(lambda_rest - wavelist) < 1e-3)
        
        if len(q[0]) == 0:
            # No matches found
            logger.warning(f"No lines found within 0.001 Å of {lambda_rest}")
            return {'wave': np.array([]), 'fval': np.array([]), 'name': np.array([])}
        
        if linelist == 'atom':
            outstr = {'wave': wavelist[q], 'fval': fval[q], 'name': name[q], 'gamma': gamma[q]}
        else:
            outstr = {'wave': wavelist[q], 'fval': fval[q], 'name': name[q]}

    elif method == 'Name':
        # Match by name
        q = np.where(name == target_name)
        
        if len(q[0]) == 0:
            # No matches found
            logger.warning(f"No lines found with name '{target_name}'")
            return {'wave': np.array([]), 'fval': np.array([]), 'name': np.array([])}
        
        if linelist == 'atom':
            outstr = {'wave': wavelist[q], 'fval': fval[q], 'name': name[q], 'gamma': gamma[q]}
        else:
            outstr = {'wave': wavelist[q], 'fval': fval[q], 'name': name[q]}

    elif method == 'closest':
        # Find the closest wavelength
        idx = (np.abs(lambda_rest - wavelist)).argmin()
        match_wavelength = wavelist[idx]
        
        # Log information about the match
        wavelength_diff = abs(lambda_rest - match_wavelength)
        logger.debug(f"Closest match to {lambda_rest} Å is {match_wavelength} Å (diff: {wavelength_diff:.4f} Å)")
        
        if wavelength_diff > 1.0:
            # Warn if the closest match is quite far from the requested wavelength
            logger.warning(f"Closest match is {wavelength_diff:.4f} Å away from requested wavelength")
        
        if linelist == 'atom':
            outstr = {'wave': wavelist[idx], 'fval': fval[idx], 'name': name[idx], 'gamma': gamma[idx]}
        else:
            outstr = {'wave': wavelist[idx], 'fval': fval[idx], 'name': name[idx]}
    
    return outstr


def read_line_list(label: str) -> List[Dict[str, Any]]:
    """
    Read a line list defined by the label.

    Parameters
    ----------
    label : str
        Label string identifying which line list to read
        Available options: 'atom', 'LLS', 'LLS Small', 'DLA', 'LBG', 'Gal',
        'Eiger_Strong', 'Gal_Em', 'Gal_Abs', 'Gal_long', 'AGN', 'HI_recomb',
        'HI_recomb_light'

    Returns
    -------
    list of dict
        A list of dictionaries containing line data with keys:
        - 'wrest': Rest wavelength
        - 'ion': Ion name
        - 'fval': f-value (oscillator strength)
        - 'gamma': Radiation damping constant (for 'atom' line list only)

    Raises
    ------
    FileNotFoundError
        If the line list file cannot be found
    ValueError
        If an invalid line list label is provided
    """
    # Check if we already have this line list cached
    if label in _LINE_LIST_CACHE:
        logger.debug(f"Using cached line list '{label}'")
        return _LINE_LIST_CACHE[label]

    # Define the mapping from label to filename
    filename_mapping = {
        'atom': resource_filename('rbcodes.IGM', 'lines/atom_full.dat'),
        'LLS': resource_filename('rbcodes.IGM', 'lines/lls.lst'),
        'LLS Small': resource_filename('rbcodes.IGM', 'lines/lls_sub.lst'),
        'DLA': resource_filename('rbcodes.IGM', 'lines/dla.lst'),
        'LBG': resource_filename('rbcodes.IGM', 'lines/lbg.lst'),
        'Gal': resource_filename('rbcodes.IGM', 'lines/gal_vac.lst'),
        'Eiger_Strong': resource_filename('rbcodes.IGM', 'lines/Eiger_Strong.lst'),
        'Gal_Em': resource_filename('rbcodes.IGM', 'lines/Galaxy_emission_Lines.lst'),
        'Gal_Abs': resource_filename('rbcodes.IGM', 'lines/Galaxy_absorption_Lines.lst'),
        'Gal_long': resource_filename('rbcodes.IGM', 'lines/Galaxy_Long_E_n_A.lst'),
        'AGN': resource_filename('rbcodes.IGM', 'lines/AGN.lst'),
        'HI_recomb': resource_filename('rbcodes.IGM', 'lines/HI_recombination.lst'),
        'HI_recomb_light': resource_filename('rbcodes.IGM', 'lines/HI_recombination_light.lst')
    }

    # Check if the label is valid
    if label not in filename_mapping:
        valid_labels = ", ".join(sorted(filename_mapping.keys()))
        raise ValueError(f"Invalid line list label: '{label}'. Valid options are: {valid_labels}")

    filename = filename_mapping[label]
    
    # Check if the file exists
    if not os.path.exists(filename):
        file_path = Path(filename)
        raise FileNotFoundError(f"Line list file not found: {file_path}")

    logger.info(f"Reading line list '{label}' from {filename}")
    data = []

    try:
        if label == 'atom':
            # Read atom format file
            s = ascii.read(filename)
            
            for line in range(len(s['col1'])):
                source = {
                    'wrest': float(s['col2'][line]),
                    'ion': f"{s['col1'][line]} {int(s['col2'][line])}",
                    'fval': float(s['col3'][line]),
                    'gamma': float(s['col4'][line])
                }
                data.append(source)

        elif label in ('LBG', 'Gal'):
            # Read LBG or Gal format file
            s = ascii.read(filename)
            
            for line in range(len(s['wrest'])):
                source = {
                    'wrest': float(s['wrest'][line]),
                    'ion': f"{s['name'][line]} {s['transition'][line]}",
                    'fval': float(s['ID'][line]),
                    'gamma': float(s['ID'][line])  # Use ID as gamma for consistency
                }
                data.append(source)

        elif label in ('Eiger_Strong', 'Gal_Em', 'Gal_Abs', 'Gal_long', 'AGN'):
            # Read simple format files
            s = ascii.read(filename)
            
            for line in range(len(s['wrest'])):
                source = {
                    'wrest': float(s['wrest'][line]),
                    'ion': s['name'][line],
                    'fval': 0.0,  # Default to 0
                    'gamma': 0.0  # Default to 0
                }
                data.append(source)

        elif label in ('HI_recomb', 'HI_recomb_light'):
            # Read recombination line lists 
            s = ascii.read(filename)
            
            for line in range(len(s['wrest'])):
                source = {
                    'wrest': float(s['wrest'][line]) * 1e4,  # Convert to Angstroms
                    'ion': s['name'][line],
                    'fval': 0.0,  # Default to 0
                    'gamma': 0.0  # Default to 0
                }
                data.append(source)

        else:
            # Generic format for other line lists
            with open(filename, 'r') as f:
                # Skip header line
                header1 = f.readline()
                
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                        
                    columns = line.split()
                    
                    # Basic validation
                    if len(columns) < 4:
                        logger.warning(f"Skipping invalid line in {filename}: {line}")
                        continue
                        
                    try:
                        source = {
                            'wrest': float(columns[0]),
                            'ion': f"{columns[1]} {columns[2]}",
                            'fval': float(columns[3])
                        }
                        data.append(source)
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing line in {filename}: {line} - {str(e)}")

        # Cache the result for future calls
        _LINE_LIST_CACHE[label] = data
        logger.debug(f"Cached '{label}' line list with {len(data)} entries")
        
        return data
    
    except Exception as e:
        logger.error(f"Error reading line list '{label}' from {filename}: {str(e)}")
        raise

# Setup default logging configuration if not already configured
if not logging.root.handlers:
    logging.basicConfig(
        level=logging.WARNING,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )