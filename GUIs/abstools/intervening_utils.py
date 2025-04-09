"""
Utilities for handling intervening absorption lines.
This module provides functions for reading and plotting
intervening absorption line data.
"""

import numpy as np
import pandas as pd
from config import SPEED_OF_LIGHT, COLORS
from config import shift2vel

def grab_intervening_linelist(filename, z_gal, wrest_galaxy, wavelength):
    """
    Read intervening absorption line data from a file and filter based on wavelength range.
    
    Parameters:
    -----------
    filename : str
        Path to the file containing intervening line data
    z_gal : float
        Redshift of the host galaxy
    wrest_galaxy : float
        Rest wavelength of the line being analyzed
    wavelength : ndarray
        Array of wavelengths in the analysis window
        
    Returns:
    --------
    outlist : dict
        Dictionary containing filtered intervening line data
    """
    try:
        # Read the file
        s = pd.read_csv(filename, sep=' ')
        
        # Extract data
        ion = s['Name'].values
        wobs = s['Wave_obs'].values
        zobs = s['Zabs'].values
        wrest = wobs / (1. + zobs)
        
        # Compute relative velocity difference with the host galaxy
        delv = shift2vel(z_gal, zobs, wrest_galaxy)
        
        # Select only lines within the given wavelength range
        window_max = np.max(wavelength)
        window_min = np.min(wavelength)
        q = np.where((wobs >= window_min) & (wobs <= window_max))
        
        outlist = {}
        if len(q[0]) > 0:
            # If there are lines, create a new list
            outlist['ion'] = ion[q]
            outlist['wrest'] = wrest[q]
            outlist['wobs'] = wobs[q]
            outlist['zobs'] = zobs[q]
            wobs_small = wobs[q]
            outlist['number'] = len(wobs_small)
            outlist['delv'] = delv[q]
            
            # Computing velocity for each
            vel = np.zeros((len(wobs_small),))
            for i in range(0, len(wobs_small)):
                vel[i] = (outlist['wobs'][i] - wrest_galaxy * (1.0 + z_gal)) * SPEED_OF_LIGHT / (wrest_galaxy * (1.0 + z_gal))
            outlist['vel'] = vel
        else:
            outlist['number'] = 0
            
        return outlist
        
    except Exception as e:
        print(f"Error reading intervening line list: {e}")
        return {'number': 0}

def plot_intervening_lines(ax, outlist, delv):
    """
    Plot intervening absorption lines on the given axes.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to plot on
    outlist : dict
        Dictionary containing intervening line data
    delv : float
        Velocity difference threshold
        
    Returns:
    --------
    None
    """
    try:
        # Check if there are intervening lines to plot
        if outlist['number'] > 0:
            vellist = outlist['vel']
            relative_vel_z = outlist['delv']
            
            for index in range(0, outlist['number']):
                if (np.abs(vellist[index]) < delv):
                    # Determine color based on velocity difference
                    if np.abs(relative_vel_z[index]) > 200.:
                        color = COLORS['pale_red']
                    else:
                        color = 'b'
                        
                    # Add text labels for the lines
                    ax.text(
                        vellist[index], 
                        1.05, 
                        f"{outlist['ion'][index]} {outlist['wrest'][index]:.0f}",
                        fontsize=8, 
                        rotation=90, 
                        rotation_mode='anchor', 
                        color=color
                    )
                    
                    ax.text(
                        vellist[index] + 50., 
                        1.05, 
                        f"z = {outlist['zobs'][index]:.3f}",
                        fontsize=8, 
                        rotation=90, 
                        rotation_mode='anchor', 
                        color=color
                    )
    except Exception as e:
        print(f"Error plotting intervening lines: {e}")
