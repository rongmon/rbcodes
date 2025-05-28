import sys
import os
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt



# PyQt5 imports
from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, 
                             QVBoxLayout, QWidget, QFileDialog, QSplitter,
                             QLabel, QHBoxLayout, QSizePolicy, QInputDialog,
                             QDialog, QListWidget,QCheckBox)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPalette, QColor
from PyQt5 import QtCore, QtGui

# Matplotlib imports
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

# Spectral data handling
from linetools.spectra.xspectrum1d import XSpectrum1D
import astropy.units as u

# Convolution and data processing
from astropy.convolution import convolve, Box1DKernel

import pandas as pd

# Custom local imports
from rbcodes.GUIs.multispecviewer.RedshiftInputWidget import RedshiftInputWidget
from rbcodes.GUIs.multispecviewer.MessageBox import MessageBox
from rbcodes.GUIs.multispecviewer.AbsorberManager import AbsorberManager
from rbcodes.GUIs.multispecviewer.io_manager import IOManager
from rbcodes.IGM import rb_setline as rb_setline
from rbcodes.utils import rb_utility as rt
from rbcodes.GUIs.multispecviewer.LineSelectionDialog import LineSelectionDialog




# Matplotlib dark theme setup
plt.style.use('dark_background')
matplotlib.rcParams.update({
    'figure.facecolor': '#353535',     # Match the dark background
    'axes.facecolor': '#353535',       # Dark background for axes
    'axes.edgecolor': '#FFFFFF',       # White axes edges
    'axes.labelcolor': '#FFFFFF',      # White label color
    'xtick.color': '#FFFFFF',          # White x-tick labels
    'ytick.color': '#FFFFFF',          # White y-tick labels
    'grid.color': '#666666',           # Dark gray grid lines
    'text.color': '#FFFFFF',           # White text
    'lines.color': '#FFFFFF',          # White line color
})

class SpectralPlot(FigureCanvas):
    """
    A Matplotlib canvas for displaying spectral plots.
    Supports interactive key commands for adjusting axis limits and quitting the application.
    """
    def __init__(self, parent=None, width=12, height=10, dpi=100, message_box=None):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(SpectralPlot, self).__init__(self.fig)
        self.setParent(parent)
        self.spectra = []  # List to store spectra data
        self.axes = None  # Stores the Matplotlib subplot axes
        self.parent_window = parent  # Reference to the parent window for quitting
        self.scale = 1. # 1D spec convolution kernel size
        self.redshift_lines = []  # Store references to plotted lines
        self.quickid_lines = []  # New list to store quick line ID markers
        self.line_list = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs'])

        
        # Store reference to the message box
        self.message_box = message_box
        
        # Configure size policy
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        # Enable focus to receive key events
        self.setFocusPolicy(Qt.StrongFocus)
        
        # Connect Matplotlib's key press event
        self.mpl_connect('key_press_event', self.on_key_press)
        
        # Connect right-click event
        self.mpl_connect('button_press_event', self.on_mouse_press)
    

    def on_mouse_press(self, event):
        """
        Handles mouse press events. Right-click shows a menu of possible spectral lines.
        """
        if event.button == 3:  # Right mouse button (button 3)
            if not self.spectra:
                if self.message_box:
                    self.message_box.on_sent_message("Load spectra", "#FF0000")
                return
            if not hasattr(self, 'linelist') or self.linelist == "None":
                if self.message_box:
                    self.message_box.on_sent_message("No line list selected, defaulting to LLS", "#FFA500")
                self.linelist='LLS'
            observed_wavelength = event.xdata
            if observed_wavelength is None:
                if self.message_box:
                    self.message_box.on_sent_message("Click on a valid position on the plot", "#FF0000")
                return
                
            # Read line list
            line_list = rb_setline.read_line_list(self.linelist)
            
            
            # Create the dialog but don't show it yet
            dialog = LineSelectionDialog(self.parent_window, observed_wavelength, line_list)
            
            # Define the handler for line selection
            def on_line_selected(transition_name, rest_wavelength, new_redshift):
                if self.message_box:
                    self.message_box.on_sent_message(
                        f"Line identified as {transition_name}, λrest = {rest_wavelength:.2f} Å, "
                        f"implied z = {new_redshift:.6f}", "#008000")
                
                # Update redshift input and plot lines with new redshift
                if hasattr(self.parent_window, 'redshift_widget'):
                    # Get the current color from the redshift widget
                    current_color = self.parent_window.redshift_widget.color_combo.currentText()
                    current_linelist = self.parent_window.redshift_widget.linelist_combo.currentText()
                    
                    self.parent_window.redshift_widget.set_redshift(new_redshift)
                    # Update the redshift value and replot lines with the current color
                    self.redshift = new_redshift
                    self.linelist = current_linelist
                    self.line_color = current_color
                    self.plot_redshift_lines()
            
            # Connect the signal to our handler
            dialog.lineSelected.connect(on_line_selected)
            
            # Show the dialog non-modally
            dialog.show()
    
    def plot_spectra(self, spectra):
        """
        Plots the given list of spectra as subplots.
        
        :param spectra: List of spectra objects from linetools.XSpectrum1D
        """
        self.spectra = spectra
        self.fig.clear()
        
        num_spectra = len(spectra)
        self.num_spectra = num_spectra
        
        # Create subplots with reduced spacing
        self.axes = self.fig.subplots(num_spectra, 1, sharex=True)
        
        # Adjust the spacing between subplots
        self.fig.subplots_adjust(hspace=0.0)  # Reduce vertical space between subplots
        # Use tight_layout with minimal padding
        self.fig.tight_layout(pad=0.5, h_pad=0.0, w_pad=0.5)
        
        if num_spectra == 1:
            self.axes = [self.axes]  # Ensure axes is always a list
        
        # Reset the original view limits storage
        self.original_xlims = []
        self.original_ylims = []
        
        for i, spec in enumerate(spectra):
            wave = spec.wavelength.value
            flux = spec.flux.value
            error = spec.sig.value if hasattr(spec, 'sig') else None
            filename = os.path.basename(spec.filename)

            self.plot_one_spec(wave, flux, error, i, filename)
            
            # Store original view limits
            self.original_xlims.append(self.axes[i].get_xlim())
            self.original_ylims.append(self.axes[i].get_ylim())
        
        # Use tight_layout with padding parameters
        self.fig.tight_layout(pad=1.0, h_pad=0.5, w_pad=1.0)
        self.draw()
        
        # Send message if message box is available
        if self.message_box:
            self.message_box.on_sent_message(f"Plotted {num_spectra} spectra", "#008000")
    
    def on_key_press(self, event):
        """
        Handles key press events for interactive adjustments of the plots.
        """
        if not self.spectra or self.axes is None or len(self.axes) == 0:
            if event.key == 'q':
                print("Exiting application...")
                if self.message_box:
                    self.message_box.on_sent_message("Exiting application...", "#FF0000")
                self.parent_window.close()
            return
        
        x, y = event.xdata, event.ydata
        
        if event.key == 'q':
            print("Exiting application...")
            if self.message_box:
                self.message_box.on_sent_message("Exiting application...", "#FF0000")
            self.parent_window.close()
            return
        
        # Handle reset key 'r' to show original spectra and remove any transient lines
        if event.key == 'r':
            self.reset_view()
            return
        # Handle capital R key to just remove the line identifications
        if event.key == 'R':
            self.clear_quickid_lines()
            self.clear_redshift_lines()
            self.draw()
            if self.message_box:
                self.message_box.on_sent_message("Cleared all line identifications", "#008000")
            return

        if x is None or y is None:
            print("No valid coordinates at cursor position")
            if self.message_box:
                self.message_box.on_sent_message("No valid coordinates at cursor position", "#FF0000")
            return

        ax_index = 0
        for i, ax in enumerate(self.axes):
            if event.inaxes == ax:
                ax_index = i
                break
        
        if event.key == 'x':
            current_xlim = self.axes[0].get_xlim()
            self.axes[0].set_xlim(x, current_xlim[1])
            if self.message_box:
                self.message_box.on_sent_message(f"Set minimum x-limit to {x:.2f}", "#8AB4F8")
            self.draw()
        elif event.key == 'X':
            current_xlim = self.axes[0].get_xlim()
            self.axes[0].set_xlim(current_xlim[0], x)
            if self.message_box:
                self.message_box.on_sent_message(f"Set maximum x-limit to {x:.2f}", "#8AB4F8")
            self.draw()
        elif event.key == 't':
            current_ylim = self.axes[ax_index].get_ylim()
            self.axes[ax_index].set_ylim(current_ylim[0], y)
            if self.message_box:
                self.message_box.on_sent_message(f"Set maximum y-limit to {y:.2f} on panel {ax_index+1}", "#8AB4F8")
            self.draw()
        elif event.key == 'b':
            current_ylim = self.axes[ax_index].get_ylim()
            self.axes[ax_index].set_ylim(y, current_ylim[1])
            if self.message_box:
                self.message_box.on_sent_message(f"Set minimum y-limit to {y:.2f} on panel {ax_index+1}", "#8AB4F8")
            self.draw()
        elif event.key == 'S':
            self.scale += 2
            if self.scale % 2 == 0:  # ensure odd
                self.scale += 1
        
            self.scale = max(1, self.scale)  # ensure at least 1
        
            # Grab spectrum
            spec = self.spectra[ax_index]
            filename = os.path.basename(spec.filename)
            wave = spec.wavelength.value
            flux = spec.flux.value
            error = spec.sig.value if hasattr(spec, 'sig') else None
        
            new_flux = convolve(flux, Box1DKernel(self.scale))
            new_err = convolve(error, Box1DKernel(self.scale)) if error is not None else None
        
            self.replot(wave, new_flux, new_err, ax_index, filename)
            self.draw()
        
            self.message_box.on_sent_message(f'Convolutional kernel size = {int(self.scale)}.', "#8AB4F8")
        
        elif event.key == 'U':
            self.scale -= 2
            if self.scale % 2 == 0:  # ensure odd
                self.scale -= 1
        
            self.scale = max(1, self.scale)  # ensure at least 1
        
            # Grab spectrum
            spec = self.spectra[ax_index]
            filename = os.path.basename(spec.filename)
            wave = spec.wavelength.value
            flux = spec.flux.value
            error = spec.sig.value if hasattr(spec, 'sig') else None
        
            new_flux = convolve(flux, Box1DKernel(self.scale))
            new_err = convolve(error, Box1DKernel(self.scale)) if error is not None else None
        
            self.replot(wave, new_flux, new_err, ax_index, filename)
            self.draw()
        
            self.message_box.on_sent_message(f'Convolutional kernel size = {int(self.scale)}.', "#8AB4F8")
        elif event.key=='[':
            xlim = self.axes[ax_index].get_xlim()
            delx = (xlim[-1] - xlim[0])
            self.axes[ax_index].set_xlim([xlim[0] - delx, xlim[0]])
            if self.message_box:
                self.message_box.on_sent_message(f"Shifted view left on panel {ax_index+1}", "#8AB4F8")
            self.draw()
        elif event.key==']':
            xlim = self.axes[ax_index].get_xlim()
            delx = (xlim[-1] - xlim[0])
            self.axes[ax_index].set_xlim([xlim[1], xlim[1] + delx])
            if self.message_box:
                self.message_box.on_sent_message(f"Shifted view right on panel {ax_index+1}", "#8AB4F8")
            self.draw()
        elif event.key == 'Y':
            Windowname='Manual y-Limits'
            instruction='Input range (e.g. 0.,2.)'
            ylim, ok = QInputDialog.getText(self, Windowname, instruction)
            if ok:
                ylimit = ylim.split(',')
                ylimit = np.array(ylimit).astype('float32')
                self.axes[ax_index].set_ylim(ylimit)
                self.draw()
        # zoom out of xrange
        elif (event.key == 'o'):
            xlim = self.axes[ax_index].get_xlim()
            ylim = self.axes[ax_index].get_ylim()
            xcen = (xlim[0]+xlim[1])/2.0
            delx = xlim[1] - xcen
            self.axes[ax_index].set_xlim([xcen - 1.5*delx, xcen + 1.5*delx])
            self.draw()
        


        # Add to the on_key_press method in SpectralPlot class
        elif event.key == 'h' or event.key == 'H':
            # Show help dialog
            from rbcodes.GUIs.multispecviewer.utils import show_help_dialog
            show_help_dialog(self.parent_window)
            if self.message_box:
                self.message_box.on_sent_message("Displayed help window", "#8AB4F8")
        


        # New quick line identification keystroke handlers
        elif (event.key == 'C'):
            self.check_lineid(event.xdata, 'CIV', event.ydata, ax_index)
        elif (event.key == 'M'):
            self.check_lineid(event.xdata, 'MgII', event.ydata, ax_index)
        elif (event.key == 'F'):
            self.check_lineid(event.xdata, 'FeII', event.ydata, ax_index)
        elif (event.key == '6'):
            self.check_lineid(event.xdata, 'OVI', event.ydata, ax_index)
        elif (event.key == '4'):
            self.check_lineid(event.xdata, 'SiIV', event.ydata, ax_index)
        elif (event.key == '8'):
            self.check_lineid(event.xdata, 'NeVIII', event.ydata, ax_index)
        elif (event.key == '2'):
            self.check_lineid(event.xdata, 'Lyb', event.ydata, ax_index)
        elif (event.key == '1'):
            self.check_lineid(event.xdata, 'Lya', event.ydata, ax_index)
        

        # Intialize vstack
        elif event.key == 'v' or event.key == 'V':
            # Check if redshift is set
            if hasattr(self, 'redshift') and hasattr(self, 'linelist'):
                # Initialize line_list if it doesn't exist
                if not hasattr(self, 'line_list'):
                    self.line_list = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs'])
                
                # Get velocity limits
                if event.key == 'V':
                    # Show dialog to get custom velocity limits
                    Windowname = 'Manual velocity-Limits'
                    instruction = 'Input range (e.g. -1500,1500)'
                    vlim_input, ok = QInputDialog.getText(self, Windowname, instruction)
                    
                    if not ok:
                        return  # User cancelled
                        
                    try:
                        # Parse velocity limits
                        vlim = vlim_input.split(',')
                        vlim = np.array(vlim).astype('float32')
                    except Exception as e:
                        if self.message_box:
                            self.message_box.on_sent_message(f"Error parsing velocity limits: {str(e)}", "#FF0000")
                        return
                else:
                    # Default velocity limits
                    vlim = [-1000, 1000]
                
                # Determine which spectrum to use based on where the user clicked
                spec_index = 0  # Default to first spectrum
                
                if event.inaxes is not None:
                    # Find which axes was clicked
                    for i, ax in enumerate(self.axes):
                        if event.inaxes == ax:
                            spec_index = i
                            break
                    
                    if spec_index >= len(self.spectra):
                        spec_index = 0  # Fallback to first spectrum
                        
                # Get the selected spectrum data
                spec = self.spectra[spec_index]
                wave = spec.wavelength.value
                flux = spec.flux.value
                error = spec.sig.value if hasattr(spec, 'sig') else None
                
                # Launch vStack with the selected spectrum
                try:
                    from rbcodes.GUIs.multispecviewer.vStack import vStack
                    self.ion_selection = vStack(self, wave, flux, error, 
                                                self.linelist, zabs=self.redshift, vlim=vlim)
                    
                    if self.message_box:
                        filename = os.path.basename(spec.filename) if hasattr(spec, 'filename') else f"Spectrum {spec_index+1}"
                        self.message_box.on_sent_message(
                            f"Launched vStack for {filename} at z={self.redshift:.6f}", "#008000")
                except Exception as e:
                    if self.message_box:
                        self.message_box.on_sent_message(f"Error launching vStack: {str(e)}", "#FF0000")
            else:
                if self.message_box:
                    self.message_box.on_sent_message("Set a redshift and line list first before using vStack", "#FF0000")

    def reset_view(self):
        """
        Resets the view to original x and y limits and removes all lines.
        Called when 'r' key is pressed.
        """
        if not self.spectra or self.axes is None or len(self.axes) == 0:
            return
        
        # Store current redshift and linelist information to reapply later if needed
        had_redshift_lines = hasattr(self, 'redshift') and hasattr(self, 'linelist')
        if had_redshift_lines:
            current_redshift = self.redshift
            current_linelist = self.linelist
            
        # Clear any redshift and quick ID lines
        self.clear_redshift_lines()
        self.clear_quickid_lines()
        
        # Reset the scale to default
        self.scale = 1.0
        
        # Replot each spectrum with original data and limits
        for i, spec in enumerate(self.spectra):
            if i < len(self.original_xlims) and i < len(self.original_ylims):
                # Get original data
                wave = spec.wavelength.value
                flux = spec.flux.value
                error = spec.sig.value if hasattr(spec, 'sig') else None
                filename = os.path.basename(spec.filename)
                
                # Clear and replot
                self.axes[i].clear()
                self.plot_one_spec(wave, flux, error, i, filename)
                
                # Set back to original limits
                self.axes[i].set_xlim(self.original_xlims[i])
                self.axes[i].set_ylim(self.original_ylims[i])
        
        # Update the canvas
        self.draw()
        
        #
        ## Reapply the redshift lines if they existed before
        #if had_redshift_lines:
        #    # Re-plot the redshift lines with the previously used values
        #    self.plot_redshift_lines()
        
        if self.message_box:
            self.message_box.on_sent_message("Reset view to original state", "#008000")

    def replot(self, wave, new_spec, new_err, ax_index, filename):
        '''Re-plot smoothed/unsmoothed spectrum
        '''
        self.cur_xlims = self.axes[ax_index].get_xlim()
        self.cur_ylims = self.axes[ax_index].get_ylim()

        self.axes[ax_index].clear() 

        self.plot_one_spec(wave, new_spec, new_err, ax_index, filename)
        
        self.axes[ax_index].set_ylim(self.cur_ylims)
        self.axes[ax_index].set_xlim(self.cur_xlims)
        self.draw()

    def plot_one_spec(self, wave, flux, error, index, filename):
        """Plot a single spectrum panel in step mode with specific dark theme colors"""
        # Define colors for dark theme
        pale_red = '#FF6B6B'  # Error spectrum color
        white = '#FFFFFF'     # Flux spectrum color
        
        # Plot error spectrum in step mode
        if error is not None:
            self.axes[index].step(wave, error, '-', where='mid', color=pale_red, lw=0.5, alpha=0.5)
        
        # Plot flux spectrum in step mode
        flux_line, = self.axes[index].step(wave, flux, '-', where='mid', color=white, lw=0.5)
        
        # Configure axes with dark theme colors
        self.axes[index].set_facecolor('#353535')  # Dark background
        self.axes[index].grid(True, linestyle='--', color='#666666', alpha=0.3)
        
        if index == self.num_spectra - 1:
            self.axes[index].set_xlabel('Wavelength', color='#FFFFFF')
        self.axes[index].set_ylabel('Flux', color='#FFFFFF')
        
        # Add legend with filename
        self.axes[index].legend([flux_line], [filename], loc='upper right', 
                                 fontsize='small', 
                                 facecolor='#353535', 
                                 edgecolor='#FFFFFF', 
                                 framealpha=0.5, 
                                 labelcolor='#FFFFFF')
        

    def set_redshift_data(self, redshift, linelist, color='white'):
        """
        Receives redshift and linelist data from the main window.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        :param color: str, the color to use for line plotting (default 'white')
        """
        # Get color dictionary
        clr = rt.rb_set_color()
        
        self.redshift = redshift
        self.linelist = linelist
        self.line_color = color  # Store the color name, not the RGB value
        
        if self.linelist == "None":
            # Clear any redshift lines
            self.clear_redshift_lines()
            self.fig.canvas.draw_idle()  # Ensure the figure updates
            if self.message_box:
                self.message_box.on_sent_message(f"{linelist} line list. Clearing lines...", "#8AB4F8")
        
        # If spectra are already loaded, update the plot with the new redshift lines
        elif self.spectra and len(self.spectra) > 0:
            # Add debug message
            if self.message_box:
                self.message_box.on_sent_message(f"Setting redshift z={redshift} with {linelist} line list in {color}", "#8AB4F8")
            # Call plot_redshift_lines directly here to plot the lines immediately
            self.plot_redshift_lines()
            # Explicitly update the canvas
            self.fig.canvas.draw()

            if self.message_box:
                self.message_box.on_sent_message(f"Applied redshift z={redshift} with {linelist} line list", "#008000")
        else:
            if self.message_box:
                self.message_box.on_sent_message(f"No spectra loaded. Load spectra before applying redshift.", "#FF0000")
    
    def plot_redshift_lines(self):
        """
        Plots spectral lines based on the redshift and linelist.
        """
        if not hasattr(self, 'redshift') or not hasattr(self, 'linelist'):
            if self.message_box:
                self.message_box.on_sent_message("Redshift or linelist not set", "#FF0000")
            return
        
        # Clear any existing redshift lines
        self.clear_redshift_lines()
        
        # Use the color from the dictionary, with a fallback
        clr = rt.rb_set_color()
        line_color = clr.get(self.line_color, clr['white'])
        
        # Get the line list data
        try:
            line_data = rb_setline.read_line_list(self.linelist)
        except Exception as e:
            if self.message_box:
                self.message_box.on_sent_message(f"Error reading line list: {str(e)}", "#FF0000")
            return
        
        # Get the full wavelength range
        wave_arrays = [spec.wavelength.value for spec in self.spectra]
        all_wavelengths = np.concatenate(wave_arrays)
        global_min_wave = np.min(all_wavelengths)
        global_max_wave = np.max(all_wavelengths)
        
        # Draw vertical lines for each spectral feature
        for ix in range(0, len(line_data)):
            rest_wavelength = line_data[ix]['wrest']
            transition_name = line_data[ix]['ion']
            observed_wavelength = rest_wavelength * (1. + self.redshift)
            
            if global_min_wave <= observed_wavelength <= global_max_wave:
                for i, ax in enumerate(self.axes):
                    line = ax.axvline(x=observed_wavelength, color=line_color, linestyle='--', alpha=0.7)
                    self.redshift_lines.append(line)
                    
                    # Only add text labels in the top panel (index 0)
                    if i == 0:  # Only for the first/top panel
                        ylim = ax.get_ylim()
                        # Position the text within the plot bounds
                        y_pos = ylim[0] + 0.85 * (ylim[1] - ylim[0])  # Position at 85% of y-axis range
                        # Use transform=ax.get_xaxis_transform() to keep the text fixed in the y-direction
                        text = ax.text(observed_wavelength, y_pos, f"{transition_name} z={self.redshift:.4f}", 
                                    rotation=90, color=line_color, fontsize=8,
                                    horizontalalignment='right', verticalalignment='top')
                        self.redshift_lines.append(text)
        self.fig.canvas.draw()  # Change from draw_idle() to draw() for immediate update


    def clear_redshift_lines(self):
        """
        Removes any previously plotted redshift lines.
        """
        if hasattr(self, 'redshift_lines') and self.redshift_lines:
            for line_obj in self.redshift_lines:
                try:
                    line_obj.remove()
                except:
                    pass
            self.redshift_lines = []
            
    def clear_quickid_lines(self):
        """
        Removes any previously plotted quick ID lines.
        """
        if hasattr(self, 'quickid_lines') and self.quickid_lines:
            for line_obj in self.quickid_lines:
                try:
                    line_obj.remove()
                except:
                    pass
            self.quickid_lines = []
    
    def check_lineid(self, wave0, ionname, yval, ax_index):
        """
        This method quickly draws some doublet/multiplet lines on the canvas for a quicklook
        
        :param wave0: The observed wavelength where the user clicked
        :param ionname: The ion to identify (e.g., 'CIV', 'MgII')
        :param yval: The y-position where the user clicked
        :param ax_index: The index of the active axis
        """
        if wave0 is None or yval is None:
            if self.message_box:
                self.message_box.on_sent_message("Invalid cursor position", "#FF0000")
            return
            
        # Set the current axis
        ax = self.axes[ax_index]
        
        # Get current y-limits to calculate relative positions
        y_min, y_max = ax.get_ylim()
        y_range = y_max - y_min
        
        # Calculate relative positions (as a fraction of the y-axis range)
        # Position the base line at 70% of the y-axis height for better visibility
        rel_pos = 0.7
        
        # Line height will be 15% of the y-axis range
        line_height = 0.15
        # Text position will be 20% of the y-axis range above the base
        text_height = 0.2
        
        # Calculate the redshift and other lines based on the ion
        if (ionname == 'CIV'):
            wave1 = wave0 * 1550.77845 / 1548.2049
            z = wave0 / 1548.2049 - 1
            msg = f"CIV: z = {z:.6f}"
        elif (ionname == 'MgII'):
            wave1 = wave0 * 2803.5314853 / 2796.3542699
            z = wave0 / 2796.354 - 1
            msg = f"MgII: z = {z:.6f}"
        elif (ionname == 'FeII'):
            wave1 = wave0 * 2586.6495659 / 2600.1724835
            wave2 = wave0 * 2382.7641781 / 2600.1724835
            z = wave0 / 2600.1724835 - 1
            msg = f"FeII: z = {z:.6f}"
        elif (ionname == 'OVI'):
            wave1 = wave0 * 1037.6167 / 1031.9261
            z = wave0 / 1031.9261 - 1
            msg = f"OVI: z = {z:.6f}"
        elif (ionname == 'NeVIII'):
            wave1 = wave0 * 780.324 / 770.409
            z = wave0 / 770.409 - 1
            msg = f"NeVIII: z = {z:.6f}"
        elif (ionname == 'SiIV'):
            wave1 = wave0 * 1402.77291 / 1393.76018
            z = wave0 / 1393.76018 - 1
            msg = f"SiIV: z = {z:.6f}"
        elif (ionname == 'Lyb'):
            wave1 = wave0 * 1215.6701 / 1025.7223
            z = wave0 / 1025.7223 - 1
            msg = f"HI Lyb: z = {z:.6f}"
        elif (ionname == 'Lya'):
            wave1 = wave0 * 1025.7223 / 1215.6701
            z = wave0 / 1215.6701 - 1
            msg = f"HI Lya: z = {z:.6f}"
        else:
            if self.message_box:
                self.message_box.on_sent_message(f"Unknown ion: {ionname}", "#FF0000")
            return
            
        # Display message
        print(msg)
        if self.message_box:
            self.message_box.on_sent_message(msg, "#008000")
            
        ## Convert relative positions back to data coordinates
        #y_pos = y_min + rel_pos * y_range
        #y_line_top = y_pos + line_height * y_range
        #y_text_pos = y_pos + text_height * y_range
        #set y-values based on where it was clicked
        y_pos = yval + yval*0.15
        y_line_top = y_pos + line_height * y_range
        y_text_pos = y_pos + text_height * y_range

        
        # Draw the line pattern with relative positioning
        line, = ax.plot([wave0, wave0, wave1, wave1], [y_pos, y_line_top, y_line_top, y_pos], color='r')
        self.quickid_lines.append(line)
        
        if ionname == 'FeII':
            # Special case for FeII which has 3 lines
            text = ax.text(0.5*(wave0+wave2), y_text_pos, f"{ionname} z: {z:.4f}", 
                          rotation=0, verticalalignment='bottom', color='r')
            line2, = ax.plot([wave1, wave1, wave2, wave2], [y_pos, y_line_top, y_line_top, y_pos], color='r')
            self.quickid_lines.append(line2)
        else:
            text = ax.text(0.5*(wave0+wave1), y_text_pos, f"{ionname} z: {z:.4f}", 
                          rotation=0, verticalalignment='bottom', color='r')
            
        self.quickid_lines.append(text)
        self.draw()
        
        

    def plot_absorber_lines(self, absorber_id, z_abs, line_list, color,**kwargs):
        """
        Plot spectral lines for a specific absorber system.
        
        :param absorber_id: Unique identifier for this absorber system
        :param z_abs: Redshift of the absorber
        :param line_list: Name of the line list to use
        :param color: Color name to use for the lines
        :return: Success status
        """
        alpha = kwargs.get('alpha', 0.7) 
        linewidth=kwargs.get('linewidth',1)
        if not self.spectra or self.axes is None or len(self.axes) == 0:
            if self.message_box:
                self.message_box.on_sent_message("No spectra loaded - cannot plot absorber lines", "#FF0000")
            return False
        
        # Skip if None is selected as the line list
        if line_list == "None":
            if self.message_box:
                self.message_box.on_sent_message(f"'None' selected as line list - nothing to plot", "#FFA500")
            return True
        
        # Store current axis limits before plotting
        current_xlims = []
        current_ylims = []
        for ax in self.axes:
            current_xlims.append(ax.get_xlim())
            current_ylims.append(ax.get_ylim())
        
        # Initialize absorber_lines dictionary if it doesn't exist
        if not hasattr(self, 'absorber_lines'):
            self.absorber_lines = {}
        
        # Create entry for this absorber if needed
        if absorber_id not in self.absorber_lines:
            self.absorber_lines[absorber_id] = []
        # Clear any existing lines for this absorber
        else:
            for line_obj in self.absorber_lines[absorber_id]:
                try:
                    line_obj.remove()
                except:
                    pass
            self.absorber_lines[absorber_id] = []
        
        # Read the line list data
        try:
            line_data = rb_setline.read_line_list(line_list)
        except Exception as e:
            if self.message_box:
                self.message_box.on_sent_message(f"Error reading line list {line_list}: {str(e)}", "#FF0000")
            return False
        
        # Get the full wavelength range from all spectra
        wave_arrays = [spec.wavelength.value for spec in self.spectra]
        all_wavelengths = np.concatenate(wave_arrays)
        global_min_wave = np.min(all_wavelengths)
        global_max_wave = np.max(all_wavelengths)
        
        # Get the color values from the dictionary
        clr = rt.rb_set_color()
        # Default to white if color isn't in the dictionary
        color_value = clr.get(color, clr['white'])
        
        # Draw vertical lines for each spectral feature
        lines_plotted = 0
        
        for ix in range(0, len(line_data)):
            # Loop through each transition
            rest_wavelength = line_data[ix]['wrest']
            transition_name = line_data[ix]['ion']
            # Calculate observed wavelength using the redshift
            observed_wavelength = rest_wavelength * (1. + z_abs)
            
            # Only plot lines within the full wavelength range of our spectra
            if global_min_wave <= observed_wavelength <= global_max_wave:
                # Draw on all axes
                for i, ax in enumerate(self.axes):
                    # Get the y-limits for this axis
                    ylim = current_ylims[i]
                    
                    # Draw the line with color from the clr dictionary
                    line_obj = ax.axvline(x=observed_wavelength, color=color_value, linestyle='--', alpha=alpha,lw=linewidth)
                    self.absorber_lines[absorber_id].append(line_obj)
                    
                    # Only add text labels in the top panel (index 0)
                    if i == 0:  # Only for the first/top panel
                        y_range = ylim[1] - ylim[0]
                        y_pos = ylim[0] + 0.85 * y_range  # 85% of the way up
                        
                        text = ax.text(observed_wavelength, y_pos, transition_name, rotation=90, 
                               horizontalalignment='right', verticalalignment='top',
                               fontsize=8, color=color_value)
                        
                        # Add redshift indicator to the text
                        text_with_z = f"{transition_name} z={z_abs:.4f}"
                        text.set_text(text_with_z)
                        
                        self.absorber_lines[absorber_id].append(text)
                        
                    lines_plotted += 1
        
        # Restore the axis limits after drawing lines
        for i, ax in enumerate(self.axes):
            ax.set_xlim(current_xlims[i])
            ax.set_ylim(current_ylims[i])
        
        self.fig.canvas.draw_idle()  # Ensure the figure updates
        
        if self.message_box:
            if lines_plotted > 0:
                self.message_box.on_sent_message(f"Plotted {lines_plotted} lines for absorber at z={z_abs:.6f}", "#008000")
            else:
                self.message_box.on_sent_message(f"No lines found within wavelength range for z={z_abs:.6f}", "#FFA500")
        
        return True  
        
    def remove_absorber_lines(self, absorber_id):
            """
            Remove all spectral lines for a specific absorber system.
            
            :param absorber_id: Identifier for the absorber system to remove
            :return: Success status
            """
            if not hasattr(self, 'absorber_lines') or absorber_id not in self.absorber_lines:
                return False
            
            # Remove each line object
            for line_obj in self.absorber_lines[absorber_id]:
                try:
                    line_obj.remove()
                except:
                    pass
            
            # Clear the list
            self.absorber_lines[absorber_id] = []
            
            # Update the plot
            self.fig.canvas.draw_idle()
            
            if self.message_box:
                self.message_box.on_sent_message(f"Removed absorber system {absorber_id}", "#008000")
            
            return True


    
    
class MainWindow(QMainWindow):
    """
    Main application window for displaying spectral plots.
    Includes a file selection button and a Matplotlib canvas.
    """
    def __init__(self):
        super(MainWindow, self).__init__()
        
        self.setWindowTitle("FITS File Viewer")
        self.setGeometry(100, 100, 1700, 900)
        
        # Dark theme palette similar to PlotSpec_Integrated
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(53, 53, 53))          # Dark background
        palette.setColor(QPalette.WindowText, QtCore.Qt.white)         # White text
        palette.setColor(QPalette.Base, QColor(25, 25, 25))            # Darker base
        palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))   # Alternate dark background
        palette.setColor(QPalette.Button, QColor(53, 53, 53))          # Button background
        palette.setColor(QPalette.ButtonText, QtCore.Qt.white)         # Button text
        palette.setColor(QPalette.BrightText, QtCore.Qt.red)           # Bright text (e.g., for warnings)
        palette.setColor(QPalette.Link, QColor(42, 130, 218))          # Link color
        palette.setColor(QPalette.Highlight, QColor(42, 130, 218))     # Highlight color
        palette.setColor(QPalette.Text, QtCore.Qt.white)               # Text color
        
        # Apply the palette to the application
        self.setPalette(palette)
        self.setStyleSheet("QMainWindow { color: #FFFFFF; }")  # Add this line                
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Initialize IO Manager
        self.io_manager = IOManager()
        
        # Create a single container for both the left and right buttons
        top_buttons_container = QWidget()
        top_buttons_layout = QHBoxLayout(top_buttons_container)
        top_buttons_layout.setContentsMargins(5, 5, 5, 5)
        
        # Create the Clear All Lines button (left side)
        self.clear_lines_button = QPushButton("Clear All Lines")
        self.clear_lines_button.setStyleSheet("""
            QPushButton {
                background-color: #8B3A3A;
                color: #F2F2F2;
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #A04545;
            }
            QPushButton:pressed {
                background-color: #702E2E;
            }
        """)
        self.clear_lines_button.clicked.connect(self.clear_all_absorber_lines)
        
        # Add the clear button to the left side
        top_buttons_layout.addWidget(self.clear_lines_button, 0)
        
        # Add a stretch to push the next set of widgets to the right
        top_buttons_layout.addStretch(1)
        
        # Create the Select FITS Files button and label (right side)
        self.select_button = QPushButton("Select FITS Files")
        self.select_button.clicked.connect(self.select_fits_files)
        self.select_button.setStyleSheet("""
            QPushButton {
                background-color: #474747;
                color: #F2F2F2;
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050;
            }
            QPushButton:pressed {
                background-color: #2A2A2A;
            }
        """)
        
        self.file_label = QLabel("No files selected")
        self.file_label.setStyleSheet("color: white;")
        
        # Add the select button and label to the right side
        top_buttons_layout.addWidget(self.select_button)
        top_buttons_layout.addWidget(self.file_label)
        
        # Add the container to the main layout
        main_layout.addWidget(top_buttons_container, 0)




        
        # Create message box first so we can pass it to the SpectralPlot
        self.message_box = MessageBox()
        self.io_manager.set_message_box(self.message_box)
        
        # Create SpectralPlot with reference to the message box
        self.canvas = SpectralPlot(self, message_box=self.message_box)
        self.toolbar = NavigationToolbar(self.canvas, self)
    
        # Create a splitter for the main display area
        self.main_splitter = QSplitter(Qt.Horizontal)
        self.main_splitter.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Create the absorber manager
        #clr = rt.rb_set_color()
        self.absorber_manager = AbsorberManager(self)#, colors=list(clr.keys())[1:])
        self.absorber_manager.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.absorber_manager.setMinimumWidth(200)
        self.absorber_manager.setMaximumWidth(300)
        
        # Add absorber manager to splitter
        self.main_splitter.addWidget(self.absorber_manager)
        
        # Create a widget for the right side of the splitter (canvas and toolbar)
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas, 10)  # Canvas should expand to fill space
        
        # Set the toolbar and canvas to expand
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        
        # Add right widget to splitter
        self.main_splitter.addWidget(right_widget)
        
        # Set initial sizes (20% for absorber manager, 80% for canvas)
        self.main_splitter.setSizes([20, 80])
        
        # Add the splitter to the main layout
        main_layout.addWidget(self.main_splitter, 10)
        
        # Create bottom widget container with redshift inputs and message box
        bottom_widget = QWidget()
        bottom_layout = QHBoxLayout(bottom_widget)
        bottom_widget.setMaximumHeight(150)  # Limit the maximum height of the bottom widget
        
        # Create redshift input widget with updated features
        self.redshift_widget = RedshiftInputWidget()
        self.redshift_widget.setMaximumWidth(300)  # Limit the width of the redshift widget
        
        # Connect to signals from updated RedshiftInputWidget
        self.redshift_widget.submitted.connect(self.handle_redshift_submission)
        self.redshift_widget.linelist_changed.connect(self.handle_linelist_changed)
        self.redshift_widget.catalog_clicked.connect(self.handle_catalog_clicked)
        
        # Add widgets to the bottom layout
        bottom_layout.addWidget(self.redshift_widget)
        bottom_layout.addWidget(self.message_box)
        
        # Add the bottom widget to the main layout
        main_layout.addWidget(bottom_widget)

        # Create a container for the action buttons (Load, Save, Show)
        button_container = QWidget()
        button_layout = QVBoxLayout(button_container)
        button_layout.setSpacing(5)
        button_container.setMaximumWidth(100)  # Set fixed width for the buttons
        
        # Create the "Load" button
        self.load_button = QPushButton("Load")
        self.load_button.setStyleSheet("""
            QPushButton {
                background-color: #474747;  
                color: #F2F2F2;             
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050; 
            }
            QPushButton:pressed {
                background-color: #2A2A2A;  
            }
        """)
        self.load_button.clicked.connect(self.handle_load_clicked)
        
        # Create the "Save" button
        self.save_button = QPushButton("Save")
        self.save_button.setStyleSheet("""
            QPushButton {
                background-color: #474747;  
                color: #F2F2F2;             
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050; 
            }
            QPushButton:pressed {
                background-color: #2A2A2A;  
            }
        """)
        self.save_button.clicked.connect(self.handle_save_clicked)
        
        # Create the "Show" button
        self.show_button = QPushButton("Show")
        self.show_button.setStyleSheet("""
            QPushButton {
                background-color: #474747;  
                color: #F2F2F2;             
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050; 
            }
            QPushButton:pressed {
                background-color: #2A2A2A;  
            }
        """)
        #connect the show button to the function we want to use
        self.show_button.clicked.connect(self.handle_show_clicked)

        # Create the "List" button
        self.list_button = QPushButton("List")
        self.list_button.setStyleSheet("""
            QPushButton {
                background-color: #474747;  
                color: #F2F2F2;             
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050; 
            }
            QPushButton:pressed {
                background-color: #2A2A2A;  
            }
        """)
        self.list_button.clicked.connect(self.display_line_list)
        
        
        # Add buttons to the container layout
        button_layout.addWidget(self.load_button)
        button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.show_button)
        button_layout.addWidget(self.list_button)
        
        # Add the button container to the bottom layout
        bottom_layout.addWidget(button_container)
        
        self.spectra = []  # Stores loaded spectra
        
        status_text = "Ready - Click on plot then use keys: x = set min x-limit, X = set max x-limit, "
        status_text += "t = set max y-limit, b = set min y-limit, q = quit application, "
        status_text += "Right-click = identify spectral lines"
        self.statusBar().showMessage(status_text)
        self.message_box.on_sent_message("Application ready. Select FITS files to begin.", "#FFFFFF")
    

    def clear_all_absorber_lines(self):
        """Remove all absorber lines from the plot and reset all checkboxes in the absorber manager"""
        # Clear redshift lines from redshift widget
        self.canvas.clear_redshift_lines()
        
        # Clear quick ID lines
        self.canvas.clear_quickid_lines()
        
        # Clear all absorber lines
        if hasattr(self.canvas, 'absorber_lines'):
            for absorber_id in list(self.canvas.absorber_lines.keys()):
                self.canvas.remove_absorber_lines(absorber_id)
        
        # Clear line objects if they exist
        if hasattr(self.canvas, 'line_objects') and self.canvas.line_objects:
            for line_obj in self.canvas.line_objects:
                try:
                    line_obj.remove()
                except Exception:
                    pass
            self.canvas.line_objects = []
        
        # Update the canvas
        self.canvas.draw()
        
        # Update UI to reflect cleared state
        if hasattr(self, 'absorber_manager'):
            for row in range(self.absorber_manager.table.rowCount()):
                checkbox_container = self.absorber_manager.table.cellWidget(row, 3)
                if checkbox_container:
                    for child in checkbox_container.children():
                        if isinstance(child, QCheckBox):
                            child.setChecked(False)
                            # Update the absorbers_df to reflect unchecked state
                            if row < len(self.absorber_manager.absorbers_df):
                                self.absorber_manager.absorbers_df.at[row, 'Visible'] = False
                            break
        
        # Show success message
        self.message_box.on_sent_message("Cleared all spectral lines from display", "#FFA500")    
    
    def handle_redshift_submission(self, redshift, linelist, color):
        """
        Processes the redshift and linelist data submitted from the RedshiftInputWidget.
        Creates the appropriate lines on the plot.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        :param color: str, the selected color
        """
        try:
            # Convert redshift to float (in case it's received as string)
            redshift = float(redshift)
            
            # Send a message to the message box
            self.message_box.on_sent_message(f"Processing redshift z={redshift} with {linelist} line list in {color}", "#008000")
            
            # Set the redshift data on the canvas directly
            self.canvas.set_redshift_data(redshift, linelist, color)
            
            # Explicitly call the plot_redshift_lines method to ensure lines are drawn
            self.canvas.plot_redshift_lines()
            
        except ValueError:
            error_msg = f"Invalid redshift value: {redshift}. Please enter a number."
            print(error_msg)
            self.message_box.on_sent_message(error_msg, "#FF0000")
        except Exception as e:
            error_msg = f"Error processing redshift: {str(e)}"
            print(error_msg)
            self.message_box.on_sent_message(error_msg, "#FF0000")
    

    def handle_linelist_changed(self, redshift, linelist, color):
        """
        Handle linelist or color selection change in the redshift widget
        
        :param redshift: float, the current redshift value
        :param linelist: str, the newly selected line list
        :param color: str, the selected color
        """
        try:
            # Clear existing lines and replot with new settings
            self.canvas.clear_redshift_lines()
            
            if linelist == "None":
                # Don't plot anything for "None" linelist
                pass
            elif self.canvas.spectra and len(self.canvas.spectra) > 0:
                # Set the data and plot with new linelist/color
                self.canvas.set_redshift_data(redshift, linelist, color)
        except Exception as e:
            self.message_box.on_sent_message(f"Error updating line list or color: {str(e)}", "#FF0000")    
    
    
    def handle_catalog_clicked(self, redshift, linelist, color):
        """
        Handles the catalog button click. Adds the current redshift to the absorber manager.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        :param color: str, the selected color
        """
        try:
            # Add to absorber manager
            if self.absorber_manager.add_absorber(redshift, linelist, color):
                self.message_box.on_sent_message(f"Added absorber z={redshift:.6f} to manager", "#008000")
                
                # Get the row of the newly added absorber
                row = self.absorber_manager.get_absorber_count() - 1
                
                # Clear the redshift window lines
                self.canvas.clear_redshift_lines()
                # Automatically plot the absorber lines
                self.canvas.plot_absorber_lines(row, redshift, linelist, color)
            else:
                self.message_box.on_sent_message(f"Failed to add absorber z={redshift:.6f}", "#FF0000")
        except Exception as e:
            self.message_box.on_sent_message(f"Error adding to catalog: {str(e)}", "#FF0000")      
    
    
    
    
        
    # Add these methods to MainWindow:
    
    def plot_absorber_lines(self, row, z_abs, line_list, color, **kwargs):
        """Wrapper to call the canvas's plot_absorber_lines method"""
        return self.canvas.plot_absorber_lines(row, z_abs, line_list, color, **kwargs)
    
    def remove_absorber_lines(self, row):
        """Wrapper to call the canvas's remove_absorber_lines method"""
        return self.canvas.remove_absorber_lines(row)
    
    def update_absorber_redshift(self, row, z_abs):
        """
        Update the redshift value for an absorber and replot its lines.
        
        :param row: The row index in the absorber manager
        :param z_abs: The new redshift value
        """
        # Get the current line list and color for this absorber
        absorber_data = self.absorber_manager.get_absorber_data(row)
        if not absorber_data:
            return
        
        # Remove existing lines
        self.canvas.remove_absorber_lines(row)
        
        # Replot with the new redshift
        self.canvas.plot_absorber_lines(row, z_abs, absorber_data['LineList'], absorber_data['Color'])  
    def select_fits_files(self):
        """
        Opens a file dialog for selecting FITS files and loads them for plotting.
        """
        options = QFileDialog.Options()
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select FITS Files", "", 
                                                     "FITS Files (*.fits *.fit);;All Files (*)", 
                                                     options=options)
        
        if file_paths:
            self.statusBar().showMessage(f"Loading {len(file_paths)} files...")
            self.message_box.on_sent_message(f"Loading {len(file_paths)} files...", "#8AB4F8")
            self.load_fits_files(file_paths)
            self.file_label.setText(f"{len(file_paths)} files selected")
            self.statusBar().showMessage(f"Loaded {len(file_paths)} files - Use keyboard shortcuts for adjustments")
            self.message_box.on_sent_message(f"Successfully loaded {len(file_paths)} files", "#008000")
            self.canvas.setFocus()
    
    def load_fits_files(self, file_paths):
        """
        Loads and plots the selected FITS files using the IO Manager.
        """
        try:
            # Use IO Manager to load the FITS files
            self.spectra, error = self.io_manager.load_fits_files(file_paths)
            
            # Exit if error occurred or no spectra loaded
            if error or not self.spectra:
                if not self.spectra:
                    self.message_box.on_sent_message("No spectra could be loaded.", "#FF0000")
                return
            
            # Plot the loaded spectra
            self.canvas.plot_spectra(self.spectra)
            self.canvas.setFocus()
        except Exception as e:
            error_message = f"Error loading FITS files: {str(e)}"
            self.statusBar().showMessage(error_message)
            self.message_box.on_sent_message(error_message, "#FF0000")    


    def handle_load_clicked(self):
        """
        Handle the Load button click event.
        Uses IOManager to load line list and absorber data.
        Provides option to append or overwrite existing data.
        """
        try:
            # Check if we already have data that might be overwritten
            has_existing_data = (
                (hasattr(self.canvas, 'line_list') and not self.canvas.line_list.empty) or
                (self.absorber_manager.get_absorber_count() > 0)
            )
            
            append_mode = False
            
            # If we have existing data, ask if the user wants to append or overwrite
            if has_existing_data:
                from PyQt5.QtWidgets import QMessageBox
                
                # Create a custom message box
                msgBox = QMessageBox(self)
                msgBox.setWindowTitle("Append or Overwrite?")
                msgBox.setText("Do you want to append to existing data or overwrite it?")
                
                # Add custom buttons
                appendButton = msgBox.addButton("Append", QMessageBox.ActionRole)
                overwriteButton = msgBox.addButton("Overwrite", QMessageBox.ActionRole)
                cancelButton = msgBox.addButton(QMessageBox.Cancel)
                
                # Set default button
                msgBox.setDefaultButton(appendButton)
                
                # Show the dialog and get the result
                msgBox.exec_()
                
                # Check which button was clicked
                if msgBox.clickedButton() == cancelButton:
                    return
                elif msgBox.clickedButton() == appendButton:
                    append_mode = True
                # else: Overwrite mode (clicked overwriteButton)
            
            # Use IO Manager to load the data
            line_list, absorbers_df, error = self.io_manager.integrated_load_data(self)
            
            if error:
                return
            
            # Update line list if loaded
            if line_list is not None and not line_list.empty:
                # Initialize line_list if it doesn't exist
                if not hasattr(self.canvas, 'line_list'):
                    self.canvas.line_list = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs'])
                
                if append_mode:
                    # Append to existing line list
                    self.canvas.line_list = pd.concat([self.canvas.line_list, line_list], ignore_index=True)
                    self.message_box.on_sent_message(f"Appended {len(line_list)} line identifications", "#008000")
                else:
                    # Replace with loaded data
                    self.canvas.line_list = line_list
                    self.message_box.on_sent_message(f"Loaded {len(line_list)} line identifications", "#008000")
            
            # Update absorbers if loaded
            if absorbers_df is not None and not absorbers_df.empty:
                # Set all absorbers to not visible
                if 'Visible' in absorbers_df.columns:
                    absorbers_df['Visible'] = False
                else:
                    absorbers_df['Visible'] = False
                
                if append_mode:
                    # Get existing absorbers data
                    existing_absorbers = self.absorber_manager.get_all_absorber_data()
                    
                    # For append mode, only clear the plot, not the absorber manager
                    if hasattr(self.canvas, 'absorber_lines'):
                        for absorber_id in self.canvas.absorber_lines:
                            for line_obj in self.canvas.absorber_lines[absorber_id]:
                                try:
                                    line_obj.remove()
                                except:
                                    pass
                        self.canvas.absorber_lines = {}
                        self.canvas.draw()
                    
                    # Add new absorbers to manager
                    for _, row in absorbers_df.iterrows():
                        self.absorber_manager.add_absorber(
                            row['Zabs'], row['LineList'], row['Color'], visible=False
                        )
                    
                    self.message_box.on_sent_message(
                        f"Appended {len(absorbers_df)} absorber systems (not displayed)", "#008000"
                    )
                else:
                    # Clear existing absorber lines from plot
                    if hasattr(self.canvas, 'absorber_lines'):
                        for absorber_id in self.canvas.absorber_lines:
                            for line_obj in self.canvas.absorber_lines[absorber_id]:
                                try:
                                    line_obj.remove()
                                except:
                                    pass
                        self.canvas.absorber_lines = {}
                        self.canvas.draw()
                    
                    # Clear existing absorbers from manager
                    for i in range(self.absorber_manager.get_absorber_count()-1, -1, -1):
                        self.absorber_manager.remove_absorber(i)
                    
                    # Add loaded absorbers
                    for _, row in absorbers_df.iterrows():
                        self.absorber_manager.add_absorber(
                            row['Zabs'], row['LineList'], row['Color'], visible=False
                        )
                    
                    self.message_box.on_sent_message(
                        f"Loaded {len(absorbers_df)} absorber systems (not displayed)", "#008000"
                    )
                
        except Exception as e:
            self.message_box.on_sent_message(f"Error loading data: {str(e)}", "#FF0000")
            import traceback
            traceback.print_exc()
    
            
    # In MainWindow.handle_save_clicked method
    def handle_save_clicked(self):
        """
        Handle the Save button click event.
        Uses IOManager to save line list and absorber data.
        """
        try:
            # Get spectrum filenames if available
            spectrum_filenames = []
            if hasattr(self, 'spectra') and self.spectra:
                spectrum_filenames = [os.path.basename(spec.filename) for spec in self.spectra]
            
            # Get line list if available
            line_list = None
            if hasattr(self.canvas, 'line_list') and not self.canvas.line_list.empty:
                line_list = self.canvas.line_list
            
            # Get absorbers if available
            absorbers_df = None
            if hasattr(self, 'absorber_manager'):
                absorbers_df = self.absorber_manager.get_all_absorber_data()
            
            # Use IO Manager to save the data
            success = self.io_manager.integrated_save_data(
                self,  # parent window
                line_list,
                absorbers_df,
                spectrum_filenames
            )
            
            if not success:
                self.message_box.on_sent_message("Save operation was not completed.", "#FFA500")
                
        except Exception as e:
            self.message_box.on_sent_message(f"Error saving data: {str(e)}", "#FF0000")


    def handle_show_clicked(self):
        """
        Handle the Show button click event.
        Displays/hides all identified lines from line_list on the plot.
        Toggle behavior: first click shows lines, second click hides them.
        """
        try:
            # Initialize storage for line objects if not already present
            if not hasattr(self.canvas, 'line_objects'):
                self.canvas.line_objects = []
            
            # If lines are currently shown, remove them
            if self.canvas.line_objects:
                # Remove all line objects
                for line_obj in self.canvas.line_objects:
                    try:
                        line_obj.remove()
                    except Exception:
                        pass  # In case the line was already removed
                
                # Clear the list
                self.canvas.line_objects = []
                
                # Update the canvas without redrawing everything
                self.canvas.draw()
                
                self.message_box.on_sent_message("Removed all shown lines", "#008000")
                return
            
            # Check if there's a line_list with identifications
            if not hasattr(self.canvas, 'line_list') or self.canvas.line_list.empty:
                self.message_box.on_sent_message("No line identifications to show", "#FFA500")
                return
            
            # Get unique redshifts and assign colors
            unique_redshifts = self.canvas.line_list['Zabs'].unique()
            
            # Define a set of distinct colors
            colors = ['cyan', 'magenta', 'yellow', 'lime', 'orange', 'pink', 'white', 
                     'red', 'green', 'blue', 'purple', 'gold', 'silver', 'turquoise']
            
            # Assign a color to each redshift
            redshift_colors = {}
            for i, z in enumerate(unique_redshifts):
                redshift_colors[z] = colors[i % len(colors)]
            
            # Plot each identified line
            lines_shown = 0
            for _, row in self.canvas.line_list.iterrows():
                wave_obs = row['Wave_obs']
                name = row['Name']
                zabs = row['Zabs']
                
                # Get the color for this redshift
                color = redshift_colors[zabs]
                
                # Plot a vertical line on all axes
                for i, ax in enumerate(self.canvas.axes):
                    # Draw the line with specified properties
                    line = ax.axvline(x=wave_obs, color=color, linestyle='--', 
                                     alpha=0.5, linewidth=0.5)
                    self.canvas.line_objects.append(line)
                    
                    # Add a text label (only on the top panel)
                    if i == 0:  # Only for the first subplot
                        # Get the current y limits
                        ylim = ax.get_ylim()
                        y_range = ylim[1] - ylim[0]
                        
                        # Position the text at 90% of the way up the panel
                        y_pos = ylim[0] + 0.9 * y_range
                        
                        # Create vertical text with name and redshift
                        text = ax.text(wave_obs, y_pos, f"{name}\nz={zabs:.4f}", 
                                      rotation=90, color=color, fontsize=8,
                                      horizontalalignment='right', verticalalignment='top')
                        self.canvas.line_objects.append(text)
                
                lines_shown += 1
            
            # Update the canvas
            self.canvas.draw()
            
            self.message_box.on_sent_message(f"Displayed {lines_shown} identified lines from {len(unique_redshifts)} absorber systems", "#008000")
            
        except Exception as e:
            self.message_box.on_sent_message(f"Error showing lines: {str(e)}", "#FF0000")

    # Display linelist and edit

    def display_line_list(self):
        """
        Display a sortable table of all identified lines with options to select and delete entries.
        """
        try:
            # Check if there's a line_list with identifications
            if not hasattr(self.canvas, 'line_list') or self.canvas.line_list.empty:
                self.message_box.on_sent_message("No line identifications to display", "#FFA500")
                return
                
            # Create a dialog to show the line list
            from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, 
                                        QPushButton, QHBoxLayout, QCheckBox, QMessageBox,
                                        QHeaderView, QLabel)
            
            dialog = QDialog(self)
            dialog.setWindowTitle("Identified Lines")
            dialog.setMinimumSize(650, 400)
            
            layout = QVBoxLayout(dialog)
            
            # Add header label with instructions
            header_label = QLabel("Select lines to delete and click 'Delete Selected' button. Click column headers to sort.")
            layout.addWidget(header_label)
            
            # Create table
            table = QTableWidget()
            table.setColumnCount(4)  # Checkbox, Name, Wavelength, Redshift
            table.setHorizontalHeaderLabels(['Delete', 'Name', 'Wavelength (Å)', 'Redshift'])
            
            # Enable sorting on data columns only (not the checkbox column)
            table.setSortingEnabled(True)
            table.horizontalHeader().setSortIndicatorShown(True)
            
            # Reset DataFrame indices to ensure clean 0,1,2,3... sequence
            line_list = self.canvas.line_list.reset_index(drop=True)
            
            # Add data
            table.setRowCount(len(line_list))
            
            for i, (idx, row) in enumerate(line_list.iterrows()):
                # Add checkbox in first column
                checkbox_widget = QWidget()
                checkbox_layout = QHBoxLayout(checkbox_widget)
                checkbox_layout.setContentsMargins(5, 0, 5, 0)
                checkbox_layout.setAlignment(Qt.AlignCenter)
                checkbox = QCheckBox()
                
                checkbox_layout.addWidget(checkbox)
                table.setCellWidget(i, 0, checkbox_widget)
                
                # Store data tuple for robust identification (Name, Wave_obs, Zabs, DataFrame_index)
                data_tuple = (str(row['Name']), float(row['Wave_obs']), float(row['Zabs']), idx)
                checkbox.setProperty("data_identifier", data_tuple)
                
                # Add other data columns with proper data types for sorting
                table.setItem(i, 1, QTableWidgetItem(str(row['Name'])))
                
                # Create numeric items for proper sorting
                wavelength_item = QTableWidgetItem()
                wavelength_item.setData(Qt.DisplayRole, f"{row['Wave_obs']:.2f}")
                wavelength_item.setData(Qt.UserRole, float(row['Wave_obs']))  # Store numeric value for sorting
                table.setItem(i, 2, wavelength_item)
                
                redshift_item = QTableWidgetItem()
                redshift_item.setData(Qt.DisplayRole, f"{row['Zabs']:.6f}")
                redshift_item.setData(Qt.UserRole, float(row['Zabs']))  # Store numeric value for sorting
                table.setItem(i, 3, redshift_item)
            
            # Disable sorting on the checkbox column
            table.horizontalHeader().setSortIndicator(1, Qt.AscendingOrder)  # Default sort by Name
            table.horizontalHeaderItem(0).setFlags(table.horizontalHeaderItem(0).flags() & ~Qt.ItemIsEnabled)
            
            # Set column widths
            table.setColumnWidth(0, 60)  # Checkbox column
            table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)  # Name column
            table.setColumnWidth(2, 120)  # Wavelength column
            table.setColumnWidth(3, 120)  # Redshift column
            
            # Add to layout
            layout.addWidget(table)
            
            # Buttons layout
            buttons_layout = QHBoxLayout()
            
            # Add select all / deselect all checkboxes
            select_widget = QWidget()
            select_layout = QHBoxLayout(select_widget)
            select_layout.setContentsMargins(0, 0, 0, 0)
            
            select_all = QCheckBox("Select All")
            select_none = QCheckBox("Deselect All")
            
            select_layout.addWidget(select_all)
            select_layout.addWidget(select_none)
            buttons_layout.addWidget(select_widget)
            
            # Function to get all checkboxes from table
            def get_all_checkboxes():
                checkboxes = []
                for i in range(table.rowCount()):
                    checkbox_widget = table.cellWidget(i, 0)
                    if checkbox_widget:
                        checkbox = checkbox_widget.findChild(QCheckBox)
                        if checkbox:
                            checkboxes.append(checkbox)
                return checkboxes
            
            # Connect select all/none checkboxes
            def on_select_all_changed(state):
                checkboxes = get_all_checkboxes()
                for checkbox in checkboxes:
                    checkbox.setChecked(state == Qt.Checked)
                # Uncheck deselect_all if select_all is checked
                if state == Qt.Checked:
                    select_none.blockSignals(True)
                    select_none.setChecked(False)
                    select_none.blockSignals(False)
            
            def on_select_none_changed(state):
                if state == Qt.Checked:
                    checkboxes = get_all_checkboxes()
                    for checkbox in checkboxes:
                        checkbox.setChecked(False)
                    # Uncheck select_all
                    select_all.blockSignals(True)
                    select_all.setChecked(False)
                    select_all.blockSignals(False)
                    # Reset deselect_all checkbox
                    select_none.blockSignals(True)
                    select_none.setChecked(False)
                    select_none.blockSignals(False)
            
            select_all.stateChanged.connect(on_select_all_changed)
            select_none.stateChanged.connect(on_select_none_changed)
            
            # Add spacer to push buttons to the right
            buttons_layout.addStretch(1)
            
            # Add delete button
            delete_button = QPushButton("Delete Selected")
            delete_button.setStyleSheet("""
                QPushButton {
                    background-color: #FF453A;
                    color: white;
                    border: none;
                    border-radius: 6px;
                    padding: 8px 16px;
                    font-size: 14px;
                }
                QPushButton:hover {
                    background-color: #FF6961;
                }
                QPushButton:pressed {
                    background-color: #D91C0D;
                }
            """)
            
            # Add export button
            export_button = QPushButton("Export to CSV")
            export_button.setStyleSheet("""
                QPushButton {
                    background-color: #30D158;
                    color: white;
                    border: none;
                    border-radius: 6px;
                    padding: 8px 16px;
                    font-size: 14px;
                }
                QPushButton:hover {
                    background-color: #4CD964;
                }
                QPushButton:pressed {
                    background-color: #248A3D;
                }
            """)
            
            # Add close button
            close_button = QPushButton("Close")
            close_button.setStyleSheet("""
                QPushButton {
                    background-color: #3A3A3C;
                    color: white;
                    border: none;
                    border-radius: 6px;
                    padding: 8px 16px;
                    font-size: 14px;
                }
                QPushButton:hover {
                    background-color: #48484A;
                }
                QPushButton:pressed {
                    background-color: #303032;
                }
            """)
            
            buttons_layout.addWidget(export_button)
            buttons_layout.addWidget(delete_button)
            buttons_layout.addWidget(close_button)
            
            layout.addLayout(buttons_layout)
            
            # Function to delete selected lines using data tuples
            def delete_selected():
                # Get all checked rows and their data identifiers
                rows_to_delete_data = []
                
                checkboxes = get_all_checkboxes()
                for checkbox in checkboxes:
                    if checkbox.isChecked():
                        data_identifier = checkbox.property("data_identifier")
                        if data_identifier:
                            rows_to_delete_data.append(data_identifier)
                
                # If nothing selected, show message and return
                if not rows_to_delete_data:
                    QMessageBox.information(dialog, "No Selection", 
                                          "No lines selected for deletion.", 
                                          QMessageBox.Ok)
                    return
                
                # Confirm deletion
                confirm_msg = f"Delete {len(rows_to_delete_data)} selected lines? This cannot be undone."
                reply = QMessageBox.question(dialog, "Confirm Deletion", 
                                            confirm_msg, 
                                            QMessageBox.Yes | QMessageBox.No, 
                                            QMessageBox.No)
                
                if reply == QMessageBox.Yes:
                    # Find matching rows in the current DataFrame and collect indices to delete
                    indices_to_delete = []
                    
                    for data_tuple in rows_to_delete_data:
                        name, wave_obs, zabs, stored_idx = data_tuple
                        
                        # Find rows that match the data (handle potential duplicates)
                        matching_mask = (
                            (self.canvas.line_list['Name'] == name) & 
                            (abs(self.canvas.line_list['Wave_obs'] - wave_obs) < 0.01) &  # Small tolerance for float comparison
                            (abs(self.canvas.line_list['Zabs'] - zabs) < 1e-6)  # Small tolerance for redshift
                        )
                        
                        matching_indices = self.canvas.line_list[matching_mask].index.tolist()
                        
                        if matching_indices:
                            # If we have the exact stored index and it matches, use it
                            if stored_idx in matching_indices:
                                indices_to_delete.append(stored_idx)
                            else:
                                # Otherwise, take the first match (handles case where indices changed)
                                indices_to_delete.append(matching_indices[0])
                    
                    # Remove duplicates from indices_to_delete
                    indices_to_delete = list(set(indices_to_delete))
                    
                    # Remove from DataFrame and reset indices
                    self.canvas.line_list = self.canvas.line_list.drop(indices_to_delete).reset_index(drop=True)
                    
                    # If lines are currently displayed, remove them and redraw
                    if hasattr(self.canvas, 'line_objects') and self.canvas.line_objects:
                        # Remove all line objects and redraw
                        for line_obj in self.canvas.line_objects:
                            try:
                                line_obj.remove()
                            except Exception:
                                pass
                        self.canvas.line_objects = []
                        self.canvas.draw()
                    
                    # Notify success
                    self.message_box.on_sent_message(f"Deleted {len(indices_to_delete)} lines", "#008000")
                    
                    # Close dialog
                    dialog.accept()
            
            # Function to export to CSV
            def export_csv():
                options = QFileDialog.Options()
                file_path, _ = QFileDialog.getSaveFileName(
                    dialog, "Export Line List", "", 
                    "CSV Files (*.csv);;All Files (*)", options=options
                )
                
                if file_path:
                    # Add .csv extension if not present
                    if not file_path.lower().endswith('.csv'):
                        file_path += '.csv'
                    
                    # Export to CSV
                    line_list.to_csv(file_path, index=False)
                    self.message_box.on_sent_message(f"Exported line list to: {file_path}", "#008000")
            
            # Connect buttons
            delete_button.clicked.connect(delete_selected)
            export_button.clicked.connect(export_csv)
            close_button.clicked.connect(dialog.reject)
            
            # Show dialog
            dialog.exec_()
            
        except Exception as e:
            self.message_box.on_sent_message(f"Error displaying line list: {str(e)}", "#FF0000")    

    #def display_line_list(self):
    #    """
    #    Display a table of all identified lines with options to select and delete entries.
    #    """
    #    try:
    #        # Check if there's a line_list with identifications
    #        if not hasattr(self.canvas, 'line_list') or self.canvas.line_list.empty:
    #            self.message_box.on_sent_message("No line identifications to display", "#FFA500")
    #            return
    #            
    #        # Create a dialog to show the line list
    #        from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, 
    #                                    QPushButton, QHBoxLayout, QCheckBox, QMessageBox,
    #                                    QHeaderView, QLabel)
    #        
    #        dialog = QDialog(self)
    #        dialog.setWindowTitle("Identified Lines")
    #        dialog.setMinimumSize(650, 400)
    #        
    #        layout = QVBoxLayout(dialog)
    #        
    #        # Add header label with instructions
    #        header_label = QLabel("Select lines to delete and click 'Delete Selected' button")
    #        #header_label.setStyleSheet("color: #FFFFFF; font-weight: bold;")
    #        layout.addWidget(header_label)
    #        
    #        # Create table
    #        table = QTableWidget()
    #        table.setColumnCount(4)  # Added checkbox column
    #        table.setHorizontalHeaderLabels(['Delete', 'Name', 'Wavelength (Å)', 'Redshift'])
    #        
    #        ## Set table properties
    #        #table.setStyleSheet("""
    #        #    QTableWidget {
    #        #        background-color: #252525;
    #        #        color: #F2F2F7;
    #        #        gridline-color: #636366;
    #        #        border-radius: 6px;
    #        #    }
    #        #    QTableWidget::item {
    #        #        padding: 4px;
    #        #    }
    #        #    QHeaderView::section {
    #        #        background-color: #3A3A3C;
    #        #        color: #F2F2F7;
    #        #        padding: 6px;
    #        #        border: 1px solid #636366;
    #        #    }
    #        #""")
    #        # Add data
    #        line_list = self.canvas.line_list
    #        table.setRowCount(len(line_list))
    #        
    #        for i, (idx, row) in enumerate(line_list.iterrows()):
    #            # Add checkbox in first column
    #            checkbox_widget = QWidget()
    #            checkbox_layout = QHBoxLayout(checkbox_widget)
    #            checkbox_layout.setContentsMargins(5, 0, 5, 0)
    #            checkbox_layout.setAlignment(Qt.AlignCenter)
    #            checkbox = QCheckBox()
    #            
    #            #checkbox.setStyleSheet("""
    #            #    QCheckBox {
    #            #        background-color: transparent;
    #            #    }
    #            #    QCheckBox::indicator {
    #            #        width: 16px;
    #            #        height: 16px;
    #            #    }
    #            #    QCheckBox::indicator:unchecked {
    #            #        background-color: #3A3A3C;
    #            #        border: 1px solid #636366;
    #            #        border-radius: 3px;
    #            #    }
    #            #    QCheckBox::indicator:checked {
    #            #        background-color: #0A84FF;
    #            #        border: 1px solid #0A84FF;
    #            #        border-radius: 3px;
    #            #    }
    #            #""")
    #            
    #            checkbox_layout.addWidget(checkbox)
    #            table.setCellWidget(i, 0, checkbox_widget)
    #            
    #            # Store the DataFrame index with the checkbox for easier deletion
    #            checkbox.setProperty("row_index", idx)
    #            
    #            # Add other data columns
    #            table.setItem(i, 1, QTableWidgetItem(str(row['Name'])))
    #            table.setItem(i, 2, QTableWidgetItem(f"{row['Wave_obs']:.2f}"))
    #            table.setItem(i, 3, QTableWidgetItem(f"{row['Zabs']:.6f}"))
    #        
    #        # Set column widths
    #        table.setColumnWidth(0, 60)  # Checkbox column
    #        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)  # Name column
    #        table.setColumnWidth(2, 120)  # Wavelength column
    #        table.setColumnWidth(3, 120)  # Redshift column
    #        
    #        # Add to layout
    #        layout.addWidget(table)
    #        
    #        # Buttons layout
    #        buttons_layout = QHBoxLayout()
    #        
    #        # Add select all / deselect all checkboxes
    #        select_widget = QWidget()
    #        select_layout = QHBoxLayout(select_widget)
    #        select_layout.setContentsMargins(0, 0, 0, 0)
    #        
    #        select_all = QCheckBox("Select All")
    #        #select_all.setStyleSheet("color: #FFFFFF;")
    #        select_none = QCheckBox("Deselect All")
    #        #select_none.setStyleSheet("color: #FFFFFF;")
    #        
    #        select_layout.addWidget(select_all)
    #        select_layout.addWidget(select_none)
    #        buttons_layout.addWidget(select_widget)
    #        
    #        # Connect select all/none checkboxes
    #        def on_select_all_changed(state):
    #            for i in range(table.rowCount()):
    #                checkbox_widget = table.cellWidget(i, 0)
    #                if checkbox_widget:
    #                    checkbox = checkbox_widget.findChild(QCheckBox)
    #                    if checkbox:
    #                        checkbox.setChecked(state == Qt.Checked)
    #            # Uncheck deselect_all if select_all is checked
    #            if state == Qt.Checked:
    #                select_none.blockSignals(True)
    #                select_none.setChecked(False)
    #                select_none.blockSignals(False)
    #        
    #        def on_select_none_changed(state):
    #            if state == Qt.Checked:
    #                for i in range(table.rowCount()):
    #                    checkbox_widget = table.cellWidget(i, 0)
    #                    if checkbox_widget:
    #                        checkbox = checkbox_widget.findChild(QCheckBox)
    #                        if checkbox:
    #                            checkbox.setChecked(False)
    #                # Uncheck select_all
    #                select_all.blockSignals(True)
    #                select_all.setChecked(False)
    #                select_all.blockSignals(False)
    #                # Reset deselect_all checkbox
    #                select_none.blockSignals(True)
    #                select_none.setChecked(False)
    #                select_none.blockSignals(False)
    #        
    #        select_all.stateChanged.connect(on_select_all_changed)
    #        select_none.stateChanged.connect(on_select_none_changed)
    #        
    #        # Add spacer to push buttons to the right
    #        buttons_layout.addStretch(1)
    #        
    #        # Add delete button
    #        delete_button = QPushButton("Delete Selected")
    #        delete_button.setStyleSheet("""
    #            QPushButton {
    #                background-color: #FF453A;
    #                color: white;
    #                border: none;
    #                border-radius: 6px;
    #                padding: 8px 16px;
    #                font-size: 14px;
    #            }
    #            QPushButton:hover {
    #                background-color: #FF6961;
    #            }
    #            QPushButton:pressed {
    #                background-color: #D91C0D;
    #            }
    #        """)
    #        
    #        # Add export button
    #        export_button = QPushButton("Export to CSV")
    #        export_button.setStyleSheet("""
    #            QPushButton {
    #                background-color: #30D158;
    #                color: white;
    #                border: none;
    #                border-radius: 6px;
    #                padding: 8px 16px;
    #                font-size: 14px;
    #            }
    #            QPushButton:hover {
    #                background-color: #4CD964;
    #            }
    #            QPushButton:pressed {
    #                background-color: #248A3D;
    #            }
    #        """)
    #        
    #        # Add close button
    #        close_button = QPushButton("Close")
    #        close_button.setStyleSheet("""
    #            QPushButton {
    #                background-color: #3A3A3C;
    #                color: white;
    #                border: none;
    #                border-radius: 6px;
    #                padding: 8px 16px;
    #                font-size: 14px;
    #            }
    #            QPushButton:hover {
    #                background-color: #48484A;
    #            }
    #            QPushButton:pressed {
    #                background-color: #303032;
    #            }
    #        """)
    #        
    #        buttons_layout.addWidget(export_button)
    #        buttons_layout.addWidget(delete_button)
    #        buttons_layout.addWidget(close_button)
    #        
    #        layout.addLayout(buttons_layout)
    #        
    #        # Function to delete selected lines
    #        def delete_selected():
    #            # Get all checked rows
    #            rows_to_delete = []
    #            indices_to_delete = []
    #            
    #            for i in range(table.rowCount()):
    #                checkbox_widget = table.cellWidget(i, 0)
    #                if checkbox_widget:
    #                    checkbox = checkbox_widget.findChild(QCheckBox)
    #                    if checkbox and checkbox.isChecked():
    #                        # Get the stored DataFrame index
    #                        idx = checkbox.property("row_index")
    #                        if idx is not None:
    #                            indices_to_delete.append(idx)
    #                            rows_to_delete.append(i)
    #            
    #            # If nothing selected, show message and return
    #            if not rows_to_delete:
    #                QMessageBox.information(dialog, "No Selection", 
    #                                      "No lines selected for deletion.", 
    #                                      QMessageBox.Ok)
    #                return
    #            
    #            # Confirm deletion
    #            confirm_msg = f"Delete {len(rows_to_delete)} selected lines? This cannot be undone."
    #            reply = QMessageBox.question(dialog, "Confirm Deletion", 
    #                                        confirm_msg, 
    #                                        QMessageBox.Yes | QMessageBox.No, 
    #                                        QMessageBox.No)
    #            
    #            if reply == QMessageBox.Yes:
    #                # Remove from DataFrame - must convert indices to list to avoid issues
    #                # with changing indices during deletion
    #                self.canvas.line_list = self.canvas.line_list.drop(indices_to_delete).reset_index(drop=True)
    #                
    #                # If lines are currently displayed, remove them and redraw
    #                if hasattr(self.canvas, 'line_objects') and self.canvas.line_objects:
    #                    # Remove all line objects and redraw
    #                    for line_obj in self.canvas.line_objects:
    #                        try:
    #                            line_obj.remove()
    #                        except Exception:
    #                            pass
    #                    self.canvas.line_objects = []
    #                    self.canvas.draw()
    #                
    #                # Notify success
    #                self.message_box.on_sent_message(f"Deleted {len(rows_to_delete)} lines", "#008000")
    #                
    #                # Close dialog
    #                dialog.accept()
    #        
    #        # Function to export to CSV
    #        def export_csv():
    #            options = QFileDialog.Options()
    #            file_path, _ = QFileDialog.getSaveFileName(
    #                dialog, "Export Line List", "", 
    #                "CSV Files (*.csv);;All Files (*)", options=options
    #            )
    #            
    #            if file_path:
    #                # Add .csv extension if not present
    #                if not file_path.lower().endswith('.csv'):
    #                    file_path += '.csv'
    #                
    #                # Export to CSV
    #                line_list.to_csv(file_path, index=False)
    #                self.message_box.on_sent_message(f"Exported line list to: {file_path}", "#008000")
    #        
    #        # Connect buttons
    #        delete_button.clicked.connect(delete_selected)
    #        export_button.clicked.connect(export_csv)
    #        close_button.clicked.connect(dialog.reject)
    #        
    #        # Show dialog
    #        dialog.exec_()
    #        
    #    except Exception as e:
    #        self.message_box.on_sent_message(f"Error displaying line list: {str(e)}", "#FF0000")

    def handle_convert_clicked(self):
        """
        Handle conversion between file formats.
        """
        try:
            from PyQt5.QtWidgets import QFileDialog, QMessageBox
            
            # Ask user for source file
            options = QFileDialog.Options()
            source_path, _ = QFileDialog.getOpenFileName(
                self, "Select File to Convert", "", 
                "All Files (*.json *.txt *.csv);;JSON Files (*.json);;Text Files (*.txt);;CSV Files (*.csv)",
                options=options
            )
            
            if not source_path:
                return  # User cancelled
                
            # Determine source type from extension
            ext = os.path.splitext(source_path)[1].lower()
            
            if ext == '.json':
                # Converting from JSON to text/CSV
                target_dir = os.path.dirname(source_path)
                success = self.io_manager.convert_json_to_text(source_path, target_dir)
                if success:
                    QMessageBox.information(self, "Conversion Complete", 
                        "Successfully converted JSON to text/CSV files.")
            else:
                # Converting from text/CSV to JSON
                options = QFileDialog.Options()
                target_path, _ = QFileDialog.getSaveFileName(
                    self, "Save JSON File", os.path.splitext(source_path)[0] + ".json", 
                    "JSON Files (*.json)", options=options
                )
                
                if not target_path:
                    return  # User cancelled
                    
                # Determine the companion file if needed
                companion_path = None
                if ext == '.txt':
                    # Look for absorbers CSV
                    potential_csv = os.path.splitext(source_path)[0] + "_absorbers.csv"
                    if os.path.exists(potential_csv):
                        companion_path = potential_csv
                elif ext == '.csv':
                    # Look for line list TXT
                    potential_txt = os.path.splitext(source_path)[0] + "_lines.txt"
                    if os.path.exists(potential_txt):
                        companion_path = potential_txt
                
                # Determine which is which
                line_list_path = source_path if ext == '.txt' else companion_path
                absorbers_path = source_path if ext == '.csv' else companion_path
                
                # Convert
                success = self.io_manager.convert_text_to_json(
                    line_list_path, absorbers_path, target_path,
                    "Converted from text/CSV", None
                )
                
                if success:
                    QMessageBox.information(self, "Conversion Complete", 
                        "Successfully converted to JSON file.")
                    
        except Exception as e:
            self.message_box.on_sent_message(f"Error converting file: {str(e)}", "#FF0000")
    

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())