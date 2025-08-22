"""
Plotting module for the absorption line analysis toolbox.
This module handles the plotting functionality for both continuum fitting
and normalized spectra.
Enhanced with improved error handling and signal-slot communication.
"""

import numpy as np
#import numpy.polynomial.legendre as L
from rbcodes.IGM import rb_iter_contfit as rf
from rbcodes.GUIs.abstools.config import COLORS, PLOT_SETTINGS
from rbcodes.GUIs.abstools.intervening_utils import grab_intervening_linelist, plot_intervening_lines

class Plotting:
    """
    Class for handling the plotting functionality of absorption line data.
    """
    
    @staticmethod
    def plot(parent, ii, modify=False, Print=False):
        """
        Main plotting function for both continuum and normalized spectra.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all the necessary data and UI elements
        ii : int
            Index of the ion to plot
        modify : bool, optional
            Whether to recalculate the continuum fit, default: False
        Print : bool, optional
            Whether to print the equivalent width text, default: False
        """
        try:
            # Determine which canvas it is on for proper indexing (page1/page2)
            key_idx = ii + 6 * parent.page
            
            # Define variables for readability
            vel = parent.ions[parent.keys[key_idx]]['vel']
            wave = parent.ions[parent.keys[key_idx]]['wave']
            error = parent.ions[parent.keys[key_idx]]['error']
            flux = parent.ions[parent.keys[key_idx]]['flux']
            weight = parent.ions[parent.keys[key_idx]]['weight']
            name = parent.ions[parent.keys[key_idx]]['name']
            wc = np.array(parent.ions[parent.keys[key_idx]]['wc'] & (error != 0))  # error != 0 is a bad pixel mask 
            cont = parent.ions[parent.keys[key_idx]]['cont']
            window_lim = parent.ions[parent.keys[key_idx]]['window_lim']
            order = parent.ions[parent.keys[key_idx]]['order']
            EWlims = parent.ions[parent.keys[key_idx]]['EWlims']
            lam_0 = parent.ions[parent.keys[key_idx]]['lam_0']
            fvals = parent.ions[parent.keys[key_idx]]['f']

            if not Print:
                # If modify == True, adjustments on LHS(continuum) are needed otherwise, only replot the RHS
                if modify:
                    try:
                        # Re-evaluate continuum
                        try:
                            # Try the optimal polynomial method first
                            result = rf.fit_optimal_polynomial(wave[wc], flux[wc], min_order=1, max_order=order, 
                                   maxiter=20, sigma=3.0, include_model=True, silent=True)
                            parent.ions[parent.keys[key_idx]]['pco'] = result['model']
                        except Exception as e:
                            print(f"Optimal polynomial fitting failed: {e}, falling back to Legendre fit")
                            # Fallback to original Legendre method
                            import numpy.polynomial.legendre as L
                            parent.ions[parent.keys[key_idx]]['pco'] = L.Legendre.fit(
                                wave[wc], flux[wc], order, w=weight[wc]
                            )

                        parent.ions[parent.keys[key_idx]]['cont'] = parent.ions[parent.keys[key_idx]]['pco'](wave)
                        cont = parent.ions[parent.keys[key_idx]]['cont']
                        
                        # Signal that continuum fit has been updated
                        if hasattr(parent, 'signals'):
                            parent.signals.continuum_fitted.emit(parent.ions[parent.keys[key_idx]])
                            parent.signals.data_updated.emit({
                                "action": "continuum_updated",
                                "ion_index": key_idx,
                                "ion_name": parent.keys[key_idx]
                            })
                    except Exception as e:
                        print(f"Error fitting continuum: {e}")
                        if hasattr(parent, 'signals'):
                            parent.signals.error_occurred.emit(f"Error fitting continuum for {name}: {str(e)}")
                        return
                    
                    # Gray_idx is to avoid line plotted through spectra from discontinuity in flux/vel/err from wc mask.
                    # Uses np.diff to find where wc (true/false array) changes and plots the line segments instead of the full line with missing values
                    try:
                        if wc[0] == True:
                            gray_idx = np.where(np.diff(wc, prepend=np.nan))[0][1:]
                        else:
                            gray_idx = np.where(np.diff(wc, prepend=np.nan))[0]
                            
                        parent.axesL[parent.page][ii].clear()
                        parent.axesL[parent.page][ii].step(vel, flux, color='k', where='mid')
                        parent.axesL[parent.page][ii].step(vel, error, color='r', where='mid')
                        parent.axesL[parent.page][ii].step(vel, cont, color='b', where='mid')
                        
                        # Plot grayed out regions that are masked
                        for zz in range(int(len(gray_idx) / 2)):
                            parent.axesL[parent.page][ii].step(
                                vel[gray_idx[zz * 2]:gray_idx[2 * zz + 1]],
                                flux[gray_idx[2 * zz]:gray_idx[2 * zz + 1]],
                                where='mid',
                                color='lightgray',
                                linewidth=1.3,
                                alpha=1
                            )
                            parent.axesL[parent.page][ii].step(
                                vel[gray_idx[zz * 2]:gray_idx[2 * zz + 1]],
                                error[gray_idx[2 * zz]:gray_idx[2 * zz + 1]],
                                where='mid',
                                color='lightgray',
                                linewidth=1.3,
                                alpha=1
                            )
                    except Exception as e:
                        print(f"Error plotting masked regions: {e}")
                        if hasattr(parent, 'signals'):
                            parent.signals.error_occurred.emit(f"Error plotting masked regions for {name}: {str(e)}")

                # Clear axes to redraw modifications
                parent.axesR[parent.page][ii].clear()

                # Replot right (flux; error)
                try:
                    parent.axesR[parent.page][ii].step(vel, flux / cont, color='k', where='mid')
                    parent.axesR[parent.page][ii].step(vel, error / cont, color='r', where='mid')
                    parent.axesR[parent.page][ii].axhline(y=1, xmin=window_lim[0], xmax=window_lim[1], ls='--', c='b')
                except Exception as e:
                    print(f"Error plotting normalized spectrum: {e}")
                    if hasattr(parent, 'signals'):
                        parent.signals.error_occurred.emit(f"Error plotting normalized spectrum for {name}: {str(e)}")
                
                # Clear y ticks and label plots
                try:
                    parent.axesL[parent.page][ii].set_yticks([])
                    parent.axesL[parent.page][ii].set_ylabel(name)
                    parent.axesR[parent.page][ii].set_ylabel(name)

                    # Set axes bounds
                    parent.axesL[parent.page][ii].set_xlim(window_lim)
                    parent.axesR[parent.page][ii].set_xlim(window_lim)
                    parent.axesR[parent.page][ii].set_ylim(PLOT_SETTINGS['normalized_ylim'])

                    # Set x ticks only on bottom row
                    if ii != parent.nions[parent.page] - 1:
                        parent.axesL[parent.page][ii].set_xticks([])
                        parent.axesR[parent.page][ii].set_xticks([])
                    else:
                        if ii > 0:  # Ensures if a new subplot is added the xticks will be removed
                            parent.axesL[parent.page][ii - 1].set_xticks([])
                            parent.axesR[parent.page][ii - 1].set_xticks([])
                        parent.axesL[parent.page][ii].set_xlabel('Velocity (km/s)')
                        parent.axesR[parent.page][ii].set_xlabel('Velocity (km/s)')
                        
                    # Plot column titles
                    if ii == 0:
                        parent.axesL[parent.page][0].set_title('Continuum Fitter')
                        parent.axesR[parent.page][0].set_title('Normalized Spectra')

                    # Plot vel = 0 line
                    parent.axesL[parent.page][ii].axvline(0, ymin=0, ymax=parent.axesL[parent.page][ii].get_ylim()[1], ls='--', c='k')
                    parent.axesR[parent.page][ii].axvline(0, ymin=0, ymax=parent.axesR[parent.page][ii].get_ylim()[1], ls='--', c='k')
                except Exception as e:
                    print(f"Error setting up plot labels/limits: {e}")
                    if hasattr(parent, 'signals'):
                        parent.signals.error_occurred.emit(f"Error setting up plot for {name}: {str(e)}")

                # Plot EW velocity limit
                try:
                    if EWlims[0] is not None:
                        parent.axesR[parent.page][ii].axvline(EWlims[0], ymin=0, ymax=2.5, ls='--', c='b')
                    if EWlims[1] is not None:
                        parent.axesR[parent.page][ii].axvline(EWlims[1], ymin=0, ymax=2.5, ls='--', c='r')
                except Exception as e:
                    print(f"Error plotting EW limits: {e}")
                    if hasattr(parent, 'signals'):
                        parent.signals.error_occurred.emit(f"Error plotting EW limits for {name}: {str(e)}")

                # Plot intervening lines if provided
                try:
                    if parent.intervening != False:
                        outlist = grab_intervening_linelist(parent.intervening, np.double(parent.z), lam_0, wave)
                        plot_intervening_lines(parent.axesR[parent.page][ii], outlist, np.max(vel))
                except Exception as e:
                    print(f"Error plotting intervening lines: {e}")
                    if hasattr(parent, 'signals'):
                        parent.signals.error_occurred.emit(f"Error plotting intervening lines for {name}: {str(e)}")

                # Plot fvals
                try:
                    text_loc = PLOT_SETTINGS['text_location']['fvals']
                    parent.axesR[parent.page][ii].text(
                        text_loc[0], text_loc[1], 
                        'f: ' + str('%.4f' % fvals), 
                        transform=parent.axesR[parent.page][ii].transAxes, 
                        color=COLORS['teal']
                    )
                except Exception as e:
                    print(f"Error plotting f-values: {e}")
                    if hasattr(parent, 'signals'):
                        parent.signals.error_occurred.emit(f"Error plotting f-values for {name}: {str(e)}")
                
                # Redraw (must be last item in list)
                try:
                    parent.figs[parent.page].canvas.draw_idle()  # Use draw_idle for better performance
                except Exception as e:
                    print(f"Error during drawing: {e}")
                    # Fall back to regular draw if draw_idle fails
                    try:
                        parent.figs[parent.page].canvas.draw()
                    except Exception as e2:
                        print(f"Critical drawing error: {e2}")
                        if hasattr(parent, 'signals'):
                            parent.signals.error_occurred.emit(f"Error updating plot: {str(e2)}")

            # Plot EW measurements on frame
            if Print:
                try:
                    from rbcodes.GUIs.abstools.text_utils import plot_text
                    
                    if parent.ions[parent.keys[key_idx]]['EW_text'] is not None:
                        try:
                            parent.ions[parent.keys[key_idx]]['EW_text'].remove()
                        except:
                            # Handle case where text element might have been already removed
                            pass

                    # Generate the text to display
                    plot_text(parent, parent.ions[parent.keys[key_idx]])
                    text = parent.ions[parent.keys[key_idx]]['text']
                    
                    # Position and add the text to the plot
                    text_loc = PLOT_SETTINGS['text_location']['ew_text']
                    EWtoggle = parent.axesR[parent.page][ii].text(
                        text_loc[0], text_loc[1], 
                        text, 
                        transform=parent.axesR[parent.page][ii].transAxes
                    )
                    
                    parent.ions[parent.keys[key_idx]]['EW_text'] = EWtoggle
                    
                    # Signal that measurements have been displayed
                    if hasattr(parent, 'signals'):
                        parent.signals.data_updated.emit({
                            "action": "measurements_displayed",
                            "ion_index": key_idx,
                            "ion_name": parent.keys[key_idx],
                            "display_mode": parent.pFlag
                        })
                    
                    # Redraw the canvas
                    try:
                        parent.figs[parent.page].canvas.draw_idle()
                    except Exception as e:
                        print(f"Error during drawing: {e}")
                        parent.figs[parent.page].canvas.draw()
                        
                except Exception as e:
                    print(f"Error plotting measurement text: {e}")
                    if hasattr(parent, 'signals'):
                        parent.signals.error_occurred.emit(f"Error displaying measurements for {name}: {str(e)}")
                        
        except Exception as e:
            print(f"Error in plotting: {e}")
            import traceback
            traceback.print_exc()