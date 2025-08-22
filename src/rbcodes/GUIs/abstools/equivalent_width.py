"""
Equivalent Width module for the absorption line analysis toolbox.
This module handles the calculation of equivalent widths, column densities,
and other line measurements.
Enhanced with improved error handling and signal-slot communication.
"""

import numpy as np
from rbcodes.IGM import compute_EW
from rbcodes.GUIs.abstools.plotting import Plotting

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
            
        Returns:
        --------
        bool
            True if calculation was successful, False otherwise
        """
        try:
            # Determine which page is being accessed
            key_idx = ii + 6 * parent.page
            
            # Signal start of calculation
            if hasattr(parent, 'signals'):
                parent.signals.data_updated.emit({
                    "action": "ew_calculation_started",
                    "ion_index": key_idx,
                    "ion_name": parent.keys[key_idx]
                })
                
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

            # Check for missing velocity limits
            if lims[0] is None or lims[1] is None:
                if hasattr(parent, 'signals'):
                    parent.signals.error_occurred.emit(f"Cannot calculate EW for {name}: Missing velocity limits")
                return False
                
            # Check for inverted velocity limits
            if lims[0] > lims[1]:
                lims = [lims[1], lims[0]]  # Swap limits to ensure correct order
                
                # Update the limits in the ion dictionary
                parent.ions[parent.keys[key_idx]]['EWlims'] = lims
                
                if hasattr(parent, 'signals'):
                    parent.signals.status_message.emit(f"Velocity limits were inverted, corrected to: {lims[0]:.1f}-{lims[1]:.1f}")

            # Compute EW/N/med_vel
            try:
                output = compute_EW.compute_EW(
                    wave, flux/cont, lam_0, lims, error/cont, 
                    plot=False, zabs=zabs, f0=f0
                )
            except Exception as e:
                if hasattr(parent, 'signals'):
                    parent.signals.error_occurred.emit(f"Error computing EW for {name}: {str(e)}")
                print(f"Error in compute_EW: {e}")
                return False
            
            # Save variables in ion's respective dictionary
            try:
                parent.ions[parent.keys[key_idx]]['N'] = output['col']
                parent.ions[parent.keys[key_idx]]['Nsig'] = output['colerr']
                parent.ions[parent.keys[key_idx]]['EW'] = output['ew_tot'] * 1000  # Convert to mÅ
                parent.ions[parent.keys[key_idx]]['EWsig'] = output['err_ew_tot'] * 1000  # Convert to mÅ
                parent.ions[parent.keys[key_idx]]['med_vel'] = output['med_vel']
                parent.ions[parent.keys[key_idx]]['EWlims'] = lims
                
                # Signal successful calculation
                if hasattr(parent, 'signals'):
                    ew_value = parent.ions[parent.keys[key_idx]]['EW']
                    ew_error = parent.ions[parent.keys[key_idx]]['EWsig']
                    col_value = parent.ions[parent.keys[key_idx]]['N']
                    col_error = parent.ions[parent.keys[key_idx]]['Nsig']
                    
                    parent.signals.data_updated.emit({
                        "action": "ew_calculation_complete",
                        "ion_index": key_idx,
                        "ion_name": name,
                        "ew_value": ew_value,
                        "ew_error": ew_error,
                        "col_value": col_value,
                        "col_error": col_error,
                        "med_vel": output['med_vel'],
                        "velocity_limits": lims
                    })
                    
                    # Format log10(N) for status message if available
                    if col_value > 0:
                        log_N = np.log10(col_value)
                        log_N_err = np.log10(col_value + col_error) - log_N
                        parent.signals.status_message.emit(
                            f"EW = {ew_value:.1f} ± {ew_error:.1f} mÅ, "
                            f"log(N) = {log_N:.2f} ± {log_N_err:.2f} for {name}"
                        )
                    else:
                        parent.signals.status_message.emit(
                            f"EW = {ew_value:.1f} ± {ew_error:.1f} mÅ for {name}"
                        )
                    
            except Exception as e:
                if hasattr(parent, 'signals'):
                    parent.signals.error_occurred.emit(f"Error saving EW calculation results for {name}: {str(e)}")
                print(f"Error saving calculation results: {e}")
                return False
            
            # Plot EW on page
            try:
                Plotting.plot(parent, ii, modify=False, Print=True)
            except Exception as e:
                if hasattr(parent, 'signals'):
                    parent.signals.error_occurred.emit(f"Error updating plot with EW results for {name}: {str(e)}")
                print(f"Error updating plot: {e}")
                # Don't return False here as the calculation itself was successful
            
            return True
            
        except Exception as e:
            print(f"Error calculating equivalent width: {e}")
            import traceback
            traceback.print_exc()
            
            if hasattr(parent, 'signals'):
                parent.signals.error_occurred.emit(f"Error in equivalent width calculation: {str(e)}")
            
            return False
    
    @staticmethod
    def calculate_all(parent):
        """
        Calculate equivalent width for all transitions in all pages.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all the necessary data and UI elements
            
        Returns:
        --------
        tuple
            (successful_count, total_count) Number of successful calculations and total ions
        """
        successful = 0
        total = 0
        
        # Store current page to restore later
        original_page = parent.page
        
        try:
            # Iterate through all pages
            for page_idx in range(len(parent.figs)):
                parent.page = page_idx
                
                # Iterate through all ions on this page
                for ii in range(parent.nions[page_idx]):
                    total += 1
                    key_idx = ii + 6 * page_idx
                    
                    # Get velocity limits
                    lims = parent.ions[parent.keys[key_idx]]['EWlims']
                    
                    # Calculate EW
                    if EquivalentWidth.calculate(parent, page_idx, ii, lims):
                        successful += 1
                        
            # Restore original page
            parent.page = original_page
            
            # Signal completion
            if hasattr(parent, 'signals'):
                parent.signals.data_updated.emit({
                    "action": "ew_calculation_batch_complete",
                    "successful": successful,
                    "total": total
                })
                
                if successful == total:
                    parent.signals.status_message.emit(f"Completed EW calculations for all {total} transitions")
                else:
                    parent.signals.status_message.emit(f"Completed {successful} of {total} EW calculations")
                    
        except Exception as e:
            print(f"Error in batch EW calculation: {e}")
            import traceback
            traceback.print_exc()
            
            if hasattr(parent, 'signals'):
                parent.signals.error_occurred.emit(f"Error in batch equivalent width calculation: {str(e)}")
            
            # Restore original page
            parent.page = original_page
            
        return (successful, total)