"""
Plotting module for the absorption line analysis toolbox.
This module handles the plotting functionality for both continuum fitting
and normalized spectra.
"""

import numpy as np
import numpy.polynomial.legendre as L
from config import COLORS, PLOT_SETTINGS
from intervening_utils import grab_intervening_linelist, plot_intervening_lines

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
                # Re-evaluate continuum
                parent.ions[parent.keys[key_idx]]['pco'] = L.Legendre.fit(
                    wave[wc], flux[wc], order, w=weight[wc]
                )
                parent.ions[parent.keys[key_idx]]['cont'] = parent.ions[parent.keys[key_idx]]['pco'](wave)
                cont = parent.ions[parent.keys[key_idx]]['cont']
                
                # Gray_idx is to avoid line plotted through spectra from discontinuity in flux/vel/err from wc mask.
                # Uses np.diff to find where wc (true/false array) changes and plots the line segments instead of the full line with missing values
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

            # Clear axes to redraw modifications
            parent.axesR[parent.page][ii].clear()

            # Replot right (flux; error)
            parent.axesR[parent.page][ii].step(vel, flux / cont, color='k', where='mid')
            parent.axesR[parent.page][ii].step(vel, error / cont, color='r', where='mid')
            parent.axesR[parent.page][ii].axhline(y=1, xmin=window_lim[0], xmax=window_lim[1], ls='--', c='b')
            
            # Clear y ticks and label plots
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

            # Plot EW velocity limit
            if EWlims[0] is not None:
                parent.axesR[parent.page][ii].axvline(EWlims[0], ymin=0, ymax=2.5, ls='--', c='b')
            if EWlims[1] is not None:
                parent.axesR[parent.page][ii].axvline(EWlims[1], ymin=0, ymax=2.5, ls='--', c='r')

            # Plot intervening lines if provided
            if parent.intervening != False:
                outlist = grab_intervening_linelist(parent.intervening, np.double(parent.z), lam_0, wave)
                plot_intervening_lines(parent.axesR[parent.page][ii], outlist, np.max(vel))

            # Plot fvals
            text_loc = PLOT_SETTINGS['text_location']['fvals']
            parent.axesR[parent.page][ii].text(
                text_loc[0], text_loc[1], 
                'f: ' + str('%.4f' % fvals), 
                transform=parent.axesR[parent.page][ii].transAxes, 
                color=COLORS['teal']
            )
            
            # Redraw (must be last item in list)
            try:
                parent.figs[parent.page].canvas.draw_idle()  # Use draw_idle for better performance
            except Exception as e:
                print(f"Error during drawing: {e}")
                # Fall back to regular draw if draw_idle fails
                parent.figs[parent.page].canvas.draw()

        # Plot EW measurements on frame
        if Print:
            from text_utils import plot_text
            
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
            
            # Redraw the canvas
            try:
                parent.figs[parent.page].canvas.draw_idle()
            except Exception as e:
                print(f"Error during drawing: {e}")
                parent.figs[parent.page].canvas.draw()
