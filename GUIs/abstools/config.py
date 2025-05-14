"""
Configuration module for the absorption line analysis toolbox.
This module centralizes constants, settings, and utility functions
used throughout the application.
"""

import numpy as np
from rbcodes.utils import rb_utility as rt

# Physical constants
SPEED_OF_LIGHT = 2.9979e5  # km/s

# Color settings
COLORS = rt.rb_set_color()

# Default settings
DEFAULT_POLY_ORDER = 6
DEFAULT_MASK_RANGE = [-200, 200]
DEFAULT_WINDOW_LIMITS = [-1000, 1000]
DEFAULT_EW_LIMITS = [None, None]

# GUI settings
MAX_TABS = 5
MAX_IONS_PER_TAB = 6
MAX_TOTAL_IONS = MAX_TABS * MAX_IONS_PER_TAB

# Plot settings
PLOT_SETTINGS = {
    'linewidth': 0.9,
    'normalized_ylim': [0, 2.2],
    'text_location': {
        'fvals': [0.85, 0.85],
        'ew_text': [0.05, 0.85]
    }
}

# Help text
HELP_TEXT = '''
        ---------------------------------------------------------------------------
        This is an interactive 1D absorption line measurement toolbox.
        This allows for interactive continuum fitting and equivalent width measurement 
        of CGM/IGM/ISM absorption lines.

       Screen Layout:
            LHS/RHS = Left Hand Side/Right Hand Side
            LMB/RMB = Left Mouse Button/Right Mouse Button
            
            LHS shows raw spectrum with overlaid legendre poly for continuum fitting
            ---grayed regions indicate masked regions
            
            RHS shows normalized spectrum with velocity limits
        ------------------------------Mouse Clicks------------------------------------

            
        Useful Mouse Clicks:
            
            LHS LMB     : Add wavelengths within region set by two clicks to continuum fit.
            LHS RMB     : Remove wavelengths from continuum fit.
            RHS LMB     : Set lower velocity limit
            RHS RMB     : Set upper velocity limit
            
        Useful Keystrokes:            

            v           : place mouse on desired subplot
                                   LHS: manually enter regions to mask continuum 
                                   RHS: manually enter EW intergration limits
            
            V (RHS only): Use active subplot velocity limits for all RHS plots
            
            Up arrow    : Increase Maximum Polynomial Order [default 6]
            Down arrow  : Decrease Maximum Polynomial Order [default 6]

            m           : Measure EW/N for active subplot
            M           : Measure EW/N for ALL subplots
            
            1/2/0 (RHS only): flag absorber as
                              (0) positive detection
                              (1) upper limit 
                              (2) lower limit

            t (RHS only): Cycle text printed on absorbers. Displays logN, or EW
            ------------------------------------------------------------------------------
            Each tab displays up to 6 transitions. There are maximum 5 tabs allowed.
            This limits the total number of transitions that can be simultaneously analyized to 30.



 

            Written By: Sean Clark, Rongmon Bordoloi [2021]

                    '''

def shift2vel(z1, z2, rest_wavelength=1215.67):
    """
    Compute velocity shift given two redshifts.
    
    Parameters:
    -----------
    z1 : float
        Reference redshift (at rest)
    z2 : float
        Redshift for which relative velocity is computed
    rest_wavelength : float, optional
        Rest wavelength in Angstroms, default: 1215.67 (Lyman-alpha)
        
    Returns:
    --------
    vel : float
        Velocity in km/s
    """
    vel = ((rest_wavelength * (1. + z2)) - (rest_wavelength * (1.0 + z1))) * SPEED_OF_LIGHT / (rest_wavelength * (1.0 + z1))
    return vel
