# spectrum_controller.py
import numpy as np
from PyQt5.QtCore import QObject, pyqtSignal
from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object
import matplotlib.pyplot as plt

class SpectrumController(QObject):
    """Controller class that manages the rb_spec instance and operations."""
    
    # Signals
    spectrum_changed = pyqtSignal()  # Emitted when spectrum data changes
    status_message = pyqtSignal(str)  # For status bar messages
    
    def __init__(self):
        super().__init__()
        self.spec = None  # rb_spec instance
        self.current_filename = None
        self.history = []  # Simple history for undo
    
    def reset_state(self):
        """Reset all analysis-related attributes."""
        # Keep only the basic spec object reference but clear its analysis attributes
        if hasattr(self, 'spec') and self.spec is not None:
            # Clear specific attributes while keeping the spec object itself
            for attr in ['cont', 'fnorm', 'enorm', 'velo', 'wave_slice', 'flux_slice', 
                        'error_slice', 'zabs', 'transition', 'transition_name', 
                        'vmin', 'vmax', 'W', 'W_e', 'N', 'N_e', 'logN', 'logN_e', 
                        'Tau', 'vel_centroid', 'vel_disp', 'SNR']:
                if hasattr(self.spec, attr):
                    delattr(self.spec, attr)
        
        # Reset controller state
        self.current_filename = None
        self.history = []  # Reset history
    
    def load_from_file(self, filename, filetype=None):
        """Load spectrum from file using rb_spec."""
        try:
            # Reset state before loading new file
            self.reset_state()
            
            self.spec = rb_spec.from_file(filename, filetype=filetype)
            self.current_filename = filename
            self.history = []  # Reset history
            self.spectrum_changed.emit()
            return True, False, False  # success, no redshift or transition info
        except Exception as e:
            print(f"Error loading file: {str(e)}")
            return False, False, False
    
    def load_from_data(self, wave, flux, error):
        """Load spectrum from numpy arrays."""
        try:
            # Reset state before loading new data
            self.reset_state()
            
            self.spec = rb_spec.from_data(wave, flux, error)
            self.current_filename = None
            self.history = []  # Reset history
            self.spectrum_changed.emit()
            return True, False, False  # success, no redshift or transition info
        except Exception as e:
            print(f"Error loading data: {str(e)}")
            return False, False, False
    
    def load_from_json(self, filename):
        """Load spectrum from saved JSON analysis."""
        try:
            # Reset state before loading new file
            self.reset_state()
            
            self.spec = load_rb_spec_object(filename)
            self.current_filename = filename
            self.history = []  # Reset history
            
            # Check if the loaded file has redshift and transition info
            has_redshift = hasattr(self.spec, 'zabs')
            has_transition = hasattr(self.spec, 'trans_wave') and hasattr(self.spec, 'trans')
            
            self.spectrum_changed.emit()
            return True, has_redshift, has_transition
        except Exception as e:
            print(f"Error loading JSON: {str(e)}")
            return False, False, False
    
    def get_json_info(self):
        """Get information from a loaded JSON file."""
        if not self.has_spectrum():
            return None, None, None, None
        
        zabs = getattr(self.spec, 'zabs', None)
        trans_wave = getattr(self.spec, 'trans_wave', None)
        trans_name = getattr(self.spec, 'trans', None)
        
        # For velocity limits, use the sliced spectrum range instead of EW values
        vmin = None
        vmax = None
        
        # If velo attribute exists, get min and max
        if hasattr(self.spec, 'velo') and self.spec.velo is not None and len(self.spec.velo) > 0:
            vmin = min(self.spec.velo)
            vmax = max(self.spec.velo)
    
        return zabs, trans_wave, trans_name, (vmin, vmax)

    def has_spectrum(self):
        """Check if a spectrum is loaded."""
        return hasattr(self, 'spec') and self.spec is not None

    def get_spectrum_data(self):
        """Get the current spectrum data for plotting."""
        if not self.has_spectrum():
            return None, None, None
        
        # Return the wavelength, flux, and error arrays
        return self.spec.wave, self.spec.flux, self.spec.error
    
    def apply_redshift(self, zabs):
        """Apply redshift to the spectrum using rb_spec."""
        if not self.has_spectrum():
            return False
            
        try:
            # Call rb_spec's shift_spec method
            self.spec.shift_spec(zabs)
            self.spectrum_changed.emit()
            return True
        except Exception as e:
            print(f"Error applying redshift: {str(e)}")
            return False

    def reset_after_redshift(self):
        """Reset analysis after redshift change, preserving the loaded spectrum."""
        if not self.has_spectrum():
            return False
        
        # Clear transition-specific attributes that depend on redshift
        for attr in ['wave_slice', 'flux_slice', 'error_slice', 'velo', 'cont',
                    'fnorm', 'enorm', 'trans', 'fval', 'trans_wave', 'vmin',
                    'vmax', 'W', 'W_e', 'N', 'N_e', 'logN', 'logN_e',
                    'vel_centroid', 'vel_disp', 'transition', 'linelist']:
            if hasattr(self.spec, attr):
                delattr(self.spec, attr)
        
        # Signal that spectrum data has changed
        self.spectrum_changed.emit()
        return True


    def reset_after_transition(self):
        """Reset analysis after transition change, preserving spectrum and redshift."""
        if not self.has_spectrum():
            return False
        
        # Clear attributes specific to analysis after slicing
        for attr in ['cont', 'fnorm', 'enorm', 'vmin', 'vmax', 
                    'W', 'W_e', 'N', 'N_e', 'logN', 'logN_e',
                    'vel_centroid', 'vel_disp']:
            if hasattr(self.spec, attr):
                delattr(self.spec, attr)
        
        # Signal that spectrum data has changed
        self.spectrum_changed.emit()
        return True
            

    def slice_spectrum(self, transition, lam_min, lam_max, use_vel=True, linelist="atom", method="closest"):
        """Slice the spectrum around a transition using rb_spec."""
        if not self.has_spectrum():
            return False
            
        try:
            # Call rb_spec's slice_spec method
            self.spec.slice_spec(transition, lam_min, lam_max, method=method, linelist=linelist, use_vel=use_vel)
            self.spectrum_changed.emit()
            return True
        except Exception as e:
            print(f"Error slicing spectrum: {str(e)}")
            return False
    
    def get_sliced_data(self):
        """Get the sliced spectrum data if available."""
        if not self.has_spectrum() or not hasattr(self.spec, 'velo'):
            return None, None, None
        
        try:
            return self.spec.velo, self.spec.flux_slice, self.spec.error_slice
        except AttributeError:
            return None, None, None

    def get_transition_info(self):
        """Get information about the currently selected transition."""
        if not self.has_spectrum() or not hasattr(self.spec, 'transition'):
            return None, None
        
        try:
            return self.spec.transition, self.spec.transition_name
        except AttributeError:
            return None, None




    def fit_continuum(self):
        """Launch interactive continuum fitter as a dialog."""
        if not self.has_spectrum() or not hasattr(self.spec, 'velo'):
            return False
            
        try:
            from rbcodes.GUIs.interactive_continuum_fit import launch_interactive_continuum_fit_dialog
            
            # Prepare the input parameters for the interactive fitter
            input_params = {
                'wave': self.spec.wave_slice,
                'flux': self.spec.flux_slice,
                'error': self.spec.error_slice,
                'velocity': self.spec.velo,
                'existing_masks': getattr(self.spec, 'continuum_masks', []),
                'order': 3,  # Default polynomial order
                'use_weights': False,
                'domain': [min(self.spec.velo), max(self.spec.velo)]
            }
            
            # Launch the interactive GUI as a dialog
            result = launch_interactive_continuum_fit_dialog(**input_params)
            
            # Check if fitting was cancelled
            if result is None or result.get('cancelled', True):
                print("Continuum fitting was cancelled.")
                return False
            
            # Update the rb_spec object with the results
            self.spec.cont = result.get('continuum')
            self.spec.fnorm = self.spec.flux_slice / self.spec.cont
            self.spec.enorm = self.spec.error_slice / self.spec.cont
            
            # Store mask information if available
            if 'masks' in result:
                self.spec.continuum_masks = result.get('masks', [])
            
            self.spectrum_changed.emit()
            return True
        except Exception as e:
            print(f"Error fitting continuum: {str(e)}")
            return False
        
    def get_continuum(self):
        """Get the fitted continuum if available."""
        if not self.has_spectrum() or not hasattr(self.spec, 'cont'):
            return None
        
        return self.spec.cont

    

    def compute_equivalent_width(self, vmin, vmax, snr=False, binsize=1, plot=False):
        """Compute equivalent width using rb_spec's compute_EW method."""
        if not self.has_spectrum() or not hasattr(self.spec, 'velo'):
            return False, None
            
        if not hasattr(self.spec, 'cont'):
            print("No continuum fitted.")
            return False, None
        
        try:
            # Call rb_spec's compute_EW method
            self.spec.compute_EW(
                self.spec.transition,
                vmin=vmin,
                vmax=vmax,
                SNR=snr,
                _binsize=binsize,
                plot=plot
            )
            
            # Collect results into a dictionary
            results = {
                'transition_name': getattr(self.spec, 'trans', 'Unknown'),
                'transition_wave': getattr(self.spec, 'trans_wave', 0),
                'vmin': vmin,
                'vmax': vmax,
                'W': getattr(self.spec, 'W', 0),
                'W_e': getattr(self.spec, 'W_e', 0),
                'N': getattr(self.spec, 'N', 0),
                'N_e': getattr(self.spec, 'N_e', 0),
                'logN': getattr(self.spec, 'logN', 0),
                'logN_e': getattr(self.spec, 'logN_e', 0),
                'vel_centroid': getattr(self.spec, 'vel_centroid', 0),
                'vel_disp': getattr(self.spec, 'vel_disp', 0),
                'vel50_err': getattr(self.spec, 'vel50_err', 0)
            }
            
            # Add SNR if calculated
            if snr and hasattr(self.spec, 'SNR'):
                results['SNR'] = self.spec.SNR
            
            # Add saturation information if available
            try:
                # rb_spec compute_EW returns a dictionary with saturation info
                line_sat = getattr(self.spec, '_EW_result', {}).get('line_saturation', False)
                sat_frac = getattr(self.spec, '_EW_result', {}).get('saturation_fraction', 0)
                
                results['line_saturation'] = line_sat
                results['saturation_fraction'] = sat_frac
            except:
                pass  # Optional, older rb_spec versions might not have this
            
            self.spectrum_changed.emit()
            return True, results
        except Exception as e:
            print(f"Error computing equivalent width: {str(e)}")
            return False, None
    
    def has_continuum(self):
        """Check if the spectrum has a fitted continuum."""
        return self.has_spectrum() and hasattr(self.spec, 'cont')
    
    def get_normalized_data(self):
        """Get the normalized spectrum data."""
        if not self.has_spectrum() or not hasattr(self.spec, 'velo'):
            return None, None, None, None
        
        if not hasattr(self.spec, 'cont'):
            return None, None, None, None
        
        try:
            return self.spec.velo, self.spec.flux_slice, self.spec.error_slice, self.spec.cont
        except AttributeError:
            return None, None, None, None
    
    def get_measurement_figure(self):
        """Get the figure from rb_spec's compute_EW for display."""
        # This requires modifying compute_EW to save the figure
        # as an attribute of the object
        # For now, we'll return None and handle it in the update_plot method
        return plt.gcf() if plt.get_fignums() else None



    def save_analysis(self, filepath, format='json'):
        """Save the analysis to a file."""
        if not self.has_spectrum():
            return False
        
        try:
            # Call rb_spec's save_slice method
            self.spec.save_slice(filepath, file_format=format)
            return True
        except Exception as e:
            print(f"Error saving analysis: {str(e)}")
            return False
    
    def export_continuum_plot(self, filepath):
        """Export the continuum fit plot to a file."""
        if not self.has_spectrum() or not hasattr(self.spec, 'cont'):
            return False
        
        try:
            # Use rb_spec's plot_continuum_fit method
            fig = self.spec.plot_continuum_fit(outfilename=filepath)
            plt.close(fig)  # Close the figure to avoid memory leak
            return True
        except Exception as e:
            print(f"Error exporting continuum plot: {str(e)}")
            return False
    
    
    def get_analysis_summary(self):
        """Get a summary of the analysis for display."""
        if not self.has_spectrum():
            return None, None, None, None, None
        
        # Get redshift
        zabs = getattr(self.spec, 'zabs', None)
        
        # Get transition info
        transition = getattr(self.spec, 'trans_wave', None)
        name = getattr(self.spec, 'trans', None)
        
        # Get EW measurements
        ew = None
        if hasattr(self.spec, 'W') and hasattr(self.spec, 'W_e'):
            ew = {
                'value': self.spec.W,
                'error': self.spec.W_e
            }
        
        # Get column density - handle negative N case
        col_dens = None
        if hasattr(self.spec, 'logN') and hasattr(self.spec, 'logN_e'):
            logN_display = self.spec.logN
            logN_e_display = self.spec.logN_e
            
            # Handle negative N case for display
            if hasattr(self.spec, 'N') and hasattr(self.spec, 'N_e'):
                if self.spec.N < 0 and self.spec.N_e > 0:
                    import numpy as np
                    logN_display = 0.0
                    logN_e_display = np.log10(self.spec.N_e)
            
            col_dens = {
                'value': logN_display,
                'error': logN_e_display
            }
        
        return zabs, transition, name, ew, col_dens

    def apply_flat_continuum(self):
        """Apply a flat continuum (value = 1.0) to the spectrum."""
        if not self.has_spectrum() or not hasattr(self.spec, 'velo'):
            return False
            
        try:
            # Create a flat continuum
            self.spec.cont = np.ones_like(self.spec.flux_slice)
            self.spec.fnorm = self.spec.flux_slice.copy()  # Normalized flux = original flux
            self.spec.enorm = self.spec.error_slice.copy()  # Normalized error = original error
            
            self.spectrum_changed.emit()
            return True
        except Exception as e:
            print(f"Error applying flat continuum: {str(e)}")
            return False

    def export_measurement_plot(self, filepath):
        """Export the equivalent width measurement plot directly from compute_EW."""
        if not self.has_spectrum() or not hasattr(self.spec, 'W'):
            return False
        
        try:
            # Create a fresh plot instead of trying to use an existing one
            # This avoids potential event loop issues
            
            # Create a new figure for the EW plot
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
            
            # Get data
            velocity = self.spec.velo
            flux = self.spec.flux_slice
            error = self.spec.error_slice
            continuum = self.spec.cont
            norm_flux = self.spec.fnorm
            norm_error = self.spec.enorm
            
            # Get EW range
            vmin = self.spec.vmin
            vmax = self.spec.vmax
            
            # Get EW measurement info
            ew = self.spec.W
            ew_err = self.spec.W_e
            logN = getattr(self.spec, 'logN', None)
            logN_err = getattr(self.spec, 'logN_e', None)
            
            # Plot flux and continuum in top panel
            ax1.step(velocity, flux, 'k-', where='mid', label='Flux')
            ax1.step(velocity, error, 'r-', where='mid', alpha=0.5, label='Error')
            ax1.plot(velocity, continuum, 'g-', linewidth=2, label='Continuum')
            ax1.axvspan(vmin, vmax, alpha=0.2, color='blue', label='EW Region')
            ax1.set_ylabel('Flux')
            ax1.legend()
            
            # Plot normalized flux in bottom panel
            ax2.step(velocity, norm_flux, 'k-', where='mid', label='Normalized Flux')
            ax2.step(velocity, norm_error, 'r-', where='mid', alpha=0.5, label='Error')
            ax2.axhline(y=1.0, color='g', linestyle='--', alpha=0.7)
            ax2.axvspan(vmin, vmax, alpha=0.2, color='blue')
            ax2.axvline(x=vmin, color='blue', linestyle='--', alpha=0.7)
            ax2.axvline(x=vmax, color='blue', linestyle='--', alpha=0.7)
            ax2.set_ylabel('Normalized Flux')
            ax2.set_xlabel('Velocity (km/s)')
            
            # Add measurement text
            info_text = f"W = {ew:.3f} ± {ew_err:.3f} Å"
            if logN is not None:
                info_text += f"\nlog N = {logN:.2f} ± {logN_err:.2f}"
            
            # Check for saturation
            if hasattr(self.spec, '_EW_result') and self.spec._EW_result.get('line_saturation', False):
                info_text += "\n(SATURATED)"
            
            ax2.text(0.98, 0.05, info_text, transform=ax2.transAxes, 
                     horizontalalignment='right', verticalalignment='bottom',
                     bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'))
            
            # Add title
            fig.suptitle(f"{self.spec.trans}: EW Measurement")
            
            # Adjust layout and save
            fig.tight_layout()
            fig.subplots_adjust(top=0.95)  # Make room for title
            
            # Save and close the figure
            fig.savefig(filepath, bbox_inches='tight')
            plt.close(fig)
            
            return True
        except Exception as e:
            print(f"Error exporting measurement plot: {str(e)}")
            return False 
