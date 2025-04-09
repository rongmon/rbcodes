"""
Equivalent Width module for the absorption line analysis toolbox.
This module handles the calculation of equivalent widths, column densities,
and other line measurements.
"""

from IGM import compute_EW
from plotting import Plotting

class EquivalentWidth:
    """
    Class for calculating equivalent widths and related properties
    for absorption lines.
    """
    
    @staticmethod
    def calculate(parent, page, ii, lims):
        """
        Calculate equivalent width, column density, and other properties
        for an absorption line.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all the necessary data and UI elements
        page : int
            The page index in the tabbed interface
        ii : int
            Index of the ion to calculate properties for
        lims : list
            The velocity limits [vmin, vmax] for the integration
        """
        try:
            # Determine which page is being accessed
            key_idx = ii + 6 * parent.page
                
            # Define variables for readability
            vel = parent.ions[parent.keys[key_idx]]['vel']
            wave = parent.ions[parent.keys[key_idx]]['wave']
            error = parent.ions[parent.keys[key_idx]]['error']
            flux = parent.ions[parent.keys[key_idx]]['flux']
            name = parent.ions[parent.keys[key_idx]]['name']
            zabs = parent.ions[parent.keys[key_idx]]['z']
            f0 = parent.ions[parent.keys[key_idx]]['f']
            lam_0 = parent.ions[parent.keys[key_idx]]['lam_0']
            cont = parent.ions[parent.keys[key_idx]]['cont']

            # Compute EW/N/med_vel
            output = compute_EW.compute_EW(
                wave, flux/cont, lam_0, lims, error/cont, 
                plot=False, zabs=zabs, f0=f0
            )
            
            # Save variables in ion's respective dictionary
            parent.ions[parent.keys[key_idx]]['N'] = output['col']
            parent.ions[parent.keys[key_idx]]['Nsig'] = output['colerr']
            parent.ions[parent.keys[key_idx]]['EW'] = output['ew_tot'] * 1000  # Convert to mÅ
            parent.ions[parent.keys[key_idx]]['EWsig'] = output['err_ew_tot'] * 1000  # Convert to mÅ
            parent.ions[parent.keys[key_idx]]['med_vel'] = output['med_vel']
            parent.ions[parent.keys[key_idx]]['EWlims'] = lims
            
            # Plot EW on page
            Plotting.plot(parent, ii, modify=False, Print=True)
            
            return True
            
        except Exception as e:
            print(f"Error calculating equivalent width: {e}")
            return False
