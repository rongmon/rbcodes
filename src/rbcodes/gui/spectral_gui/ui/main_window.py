#!/usr/bin/env python
"""
Main Window - The primary user interface for the spectrum viewer.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, 
    QPushButton, QComboBox, QLineEdit, QLabel, QMessageBox,
    QTableWidget, QTableWidgetItem, QSpacerItem, QSizePolicy
)
from PyQt5.QtGui import QPalette, QColor

from rbcodes.GUIs.spectral_gui.core.spectrum import Spectrum
from rbcodes.GUIs.spectral_gui.core.absorber import AbsorberManager
from rbcodes.GUIs.spectral_gui.core.line_manager import LineList
from rbcodes.GUIs.spectral_gui.utils.constants import COLORS, LINE_LISTS
from rbcodes.GUIs.spectral_gui.ui.absorber_manager import AbsorberManagerWidget
from rbcodes.GUIs.spectral_gui.ui.help_window import HelpWindow


class MainWindow(QMainWindow):
    """
    Main window for the spectrum viewer application.
    
    This class defines the primary user interface and coordinates communication
    between various components of the application.
    
    Parameters
    ----------
    wavelength : array-like
        Wavelength array
    flux : array-like
        Flux array
    error : array-like
        Error array
    redshift : float, optional
        Initial redshift value, by default 0.0
    parent : QWidget, optional
        Parent widget, by default None
    """
    
    def __init__(self, wavelength, flux, error, redshift=0.0, parent=None):
        """Initialize the main window."""
        super(MainWindow, self).__init__(parent)
        
        # Set window properties
        self.setWindowTitle('Spectrum Viewer')
        self.resize(1700, 900)
        
        # Initialize data models
        self.spectrum = Spectrum(wavelength, flux, error, redshift)
        self.absorber_manager = AbsorberManager()
        
        # Track state for various UI operations
        self.active_redshift = redshift
        self.active_line_list = 'LLS'
        self.active_color = 'white'
        self.ew_points = []
        self.ew_y_values = []
        self.gaussian_points = []
        self.gaussian_y_values = []
        self.identified_lines_active = False
        
        # Set up the UI
        self._create_ui()
        self._connect_signals()
        self._initialize_plot()
        
        # Center the window on screen
        self._center_window()
        
        # Initialize with focus on the canvas for keyboard events
        self.canvas.setFocus()
    
    def _create_ui(self):
        """Create the user interface components."""
        # Create central widget and main layout
        central_widget = QWidget()
        main_layout = QHBoxLayout(central_widget)
        
        # Create the absorber manager widget
        self.absorber_widget = AbsorberManagerWidget(self.absorber_manager)
        
        # Create save/load buttons
        save_layout = QHBoxLayout()
        self.save_button = QPushButton("Save")
        self.load_button = QPushButton("Load")
        save_layout.addWidget(self.save_button)
        save_layout.addWidget(self.load_button)
        
        # Create identified lines button
        self.identified_lines_button = QPushButton("Plot Identified Lines")
        
        # Combine left panel widgets
        left_panel = QVBoxLayout()
        left_panel.addWidget(self.absorber_widget)
        left_panel.addLayout(save_layout)
        left_panel.addWidget(self.identified_lines_button)
        
        # Create spectrum plot
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        
        # Create active values section (bottom of main panel)
        active_layout = QFormLayout()
        
        # Redshift input
        self.redshift_label = QLabel("Active Redshift")
        self.redshift_input = QLineEdit()
        self.redshift_input.setText(str(self.active_redshift))
        active_layout.addRow(self.redshift_label, self.redshift_input)
        
        # Line list combo
        self.line_list_label = QLabel("Active LineList")
        self.line_list_combo = QComboBox()
        for line_list in LINE_LISTS:
            self.line_list_combo.addItem(line_list)
        active_layout.addRow(self.line_list_label, self.line_list_combo)
        
        # Color combo
        self.color_label = QLabel("Color")
        self.color_combo = QComboBox()
        for color in COLORS:
            self.color_combo.addItem(color)
        active_layout.addRow(self.color_label, self.color_combo)
        
        # Message window
        self.message_window = QLabel("Message Window")
        self.message_window.setStyleSheet('background-color: black')
        
        # Action buttons
        action_layout = QVBoxLayout()
        self.plot_button = QPushButton("Plot")
        self.catalog_button = QPushButton("Catalog")
        self.refresh_button = QPushButton("Refresh")
        action_layout.addWidget(self.plot_button)
        action_layout.addWidget(self.catalog_button)
        action_layout.addWidget(self.refresh_button)
        
        # Create bottom bar layout
        bottom_layout = QHBoxLayout()
        bottom_layout.addWidget(self.message_window, stretch=1)
        bottom_layout.addLayout(active_layout, stretch=1)
        bottom_layout.addLayout(action_layout)
        
        # Add spacer to push components to the left
        spacer = QSpacerItem(100, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        bottom_layout.addItem(spacer)
        
        # Create main plot layout
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.toolbar, stretch=1)
        plot_layout.addWidget(self.canvas, stretch=5)
        plot_layout.addLayout(bottom_layout, stretch=1)
        
        # Add layouts to main layout
        main_layout.addLayout(left_panel, 28)
        main_layout.addLayout(plot_layout, 80)
        
        # Set central widget
        self.setCentralWidget(central_widget)
    
    def _connect_signals(self):
        """Connect UI signals to handlers."""
        # Keyboard and mouse events on canvas
        self.canvas.mpl_connect('key_press_event', self._on_key_press)
        self.canvas.mpl_connect('button_press_event', self._on_mouse_click)
        
        # Input changes
        self.redshift_input.textChanged.connect(self._on_redshift_changed)
        self.line_list_combo.currentIndexChanged.connect(self._on_line_list_changed)
        self.color_combo.currentIndexChanged.connect(self._on_color_changed)
        
        # Button clicks
        self.save_button.clicked.connect(self._on_save_catalog)
        self.load_button.clicked.connect(self._on_load_catalog)
        self.identified_lines_button.clicked.connect(self._on_toggle_identified_lines)
        self.plot_button.clicked.connect(self._on_plot_absorber)
        self.catalog_button.clicked.connect(self._on_catalog_absorber)
        self.refresh_button.clicked.connect(self._on_refresh)
        
        # Absorber widget signals
        self.absorber_widget.plotRequested.connect(self._on_absorber_plot)
        self.absorber_widget.hideRequested.connect(self._on_absorber_hide)
        self.absorber_widget.removeRequested.connect(self._on_absorber_remove)
    
    
    def _initialize_plot(self):
        """Initialize the spectrum plot."""
        # Create figure with dark background
        self.fig.patch.set_facecolor('#303030')
        self.ax.set_facecolor('#303030')
        
        # Plot error spectrum
        self.ax.step(self.spectrum.wavelength, self.spectrum.error, '-', 
                     lw=0.5, color='red', alpha=0.7, zorder=2)
        
        # Plot flux spectrum
        self.ax.step(self.spectrum.wavelength, self.spectrum.flux, '-', 
                     lw=0.5, color='white')
        
        # Configure axes
        self.ax.set_xlabel('Wavelength', color='white')
        self.ax.set_ylabel('Flux', color='white')
        self.ax.set_xlim(self.spectrum.initial_xlim)
        self.ax.set_ylim(self.spectrum.initial_ylim)
        
        # Set tick colors
        self.ax.tick_params(axis='both', colors='white')
        
        # Set spines color
        for spine in self.ax.spines.values():
            spine.set_color('white')
        
        # Draw the canvas
        self.canvas.draw()
    
    
    def _center_window(self):
        """Center the window on the screen."""
        frame_geometry = self.frameGeometry()
        center_point = QtWidgets.QDesktopWidget().availableGeometry().center()
        frame_geometry.moveCenter(center_point)
        self.move(frame_geometry.topLeft())
    
    def _on_key_press(self, event):
        """
        Handle keyboard events on the canvas.
        
        Parameters
        ----------
        event : matplotlib.backend_bases.KeyEvent
            Matplotlib key event
        """
        # Get current axes limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        
        # Reset view
        if event.key == 'r':
            self._reset_view(full=True)
        
        # Clear lines but keep view
        elif event.key == 'R':
            self._clear_lines()
        
        # Set y-axis maximum
        elif event.key == 't':
            if event.ydata:
                self.ax.set_ylim(ylim[0], event.ydata)
                self.canvas.draw()
        
        # Set y-axis minimum
        elif event.key == 'b':
            if event.ydata:
                self.ax.set_ylim(event.ydata, ylim[1])
                self.canvas.draw()
        
        # Smooth spectrum (increase)
        elif event.key == 'S':
            self._smooth_spectrum(increment=True)
        
        # Smooth spectrum (decrease)
        elif event.key == 'U':
            self._smooth_spectrum(increment=False)
        
        # Set x-axis maximum
        elif event.key == 'X':
            if event.xdata:
                self.ax.set_xlim(xlim[0], event.xdata)
                self.canvas.draw()
        
        # Set x-axis minimum
        elif event.key == 'x':
            if event.xdata:
                self.ax.set_xlim(event.xdata, xlim[1])
                self.canvas.draw()
        
        # Pan right
        elif event.key == ']':
            delta_x = xlim[1] - xlim[0]
            self.ax.set_xlim(xlim[1], xlim[1] + delta_x)
            self.canvas.draw()
        
        # Pan left
        elif event.key == '[':
            delta_x = xlim[1] - xlim[0]
            self.ax.set_xlim(xlim[0] - delta_x, xlim[0])
            self.canvas.draw()
        
        # Zoom out
        elif event.key == 'o':
            center = (xlim[0] + xlim[1]) / 2.0
            half_width = (xlim[1] - xlim[0]) / 2.0
            self.ax.set_xlim(center - 1.5 * half_width, center + 1.5 * half_width)
            self.canvas.draw()
        
        # Set manual y limits
        elif event.key == 'Y':
            self._set_manual_ylimits()
        
        # Mark CIV doublet
        elif event.key == 'C':
            if event.xdata and event.ydata:
                self._mark_doublet('CIV', event.xdata, event.ydata)
        
        # Mark MgII doublet
        elif event.key == 'M':
            if event.xdata and event.ydata:
                self._mark_doublet('MgII', event.xdata, event.ydata)
        
        # Mark FeII triplet
        elif event.key == 'F':
            if event.xdata and event.ydata:
                self._mark_doublet('FeII', event.xdata, event.ydata)
        
        # Mark OVI doublet
        elif event.key == '6':
            if event.xdata and event.ydata:
                self._mark_doublet('OVI', event.xdata, event.ydata)
        
        # Mark SiIV doublet
        elif event.key == '4':
            if event.xdata and event.ydata:
                self._mark_doublet('SiIV', event.xdata, event.ydata)
        
        # Mark NeVIII doublet
        elif event.key == '8':
            if event.xdata and event.ydata:
                self._mark_doublet('NeVIII', event.xdata, event.ydata)
        
        # Mark Lyb/Lya pair
        elif event.key == '2':
            if event.xdata and event.ydata:
                self._mark_doublet('Lyb', event.xdata, event.ydata)
        
        # Mark Lya/Lyb pair
        elif event.key == '1':
            if event.xdata and event.ydata:
                self._mark_doublet('Lya', event.xdata, event.ydata)
        
        # Measure equivalent width
        elif event.key == 'E':
            if event.xdata and event.ydata:
                self._measure_equivalent_width(event.xdata, event.ydata)
        
        # Fit Gaussian
        elif event.key == 'G':
            if event.xdata and event.ydata:
                self._fit_gaussian(event.xdata, event.ydata)
        
        # Show help window
        elif event.key == 'H' or event.key == 'h':
            self._show_help()
        
        # Open vstack display
        elif event.key == 'v':
            self._open_vstack()
        
        # Open vstack with custom velocity limits
        elif event.key == 'V':
            self._open_vstack_custom_velocity()
        
        # Open manual transition selector
        elif event.key == 'j':
            if event.xdata:
                self._open_manual_transition(event.xdata)
    
    def _on_mouse_click(self, event):
        """
        Handle mouse clicks on the canvas.
        
        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Matplotlib mouse event
        """
        # Right-click opens manual transition selector
        if event.button == 3 and event.xdata:
            self._open_manual_transition(event.xdata)
    
    def _on_redshift_changed(self):
        """Handle changes to the redshift input."""
        try:
            self.active_redshift = float(self.redshift_input.text())
            self.spectrum.set_redshift(self.active_redshift)
        except ValueError:
            self.redshift_input.setText("Invalid input")
    
    def _on_line_list_changed(self):
        """Handle changes to the line list combo box."""
        self.active_line_list = self.line_list_combo.currentText()
    
    def _on_color_changed(self):
        """Handle changes to the color combo box."""
        self.active_color = self.color_combo.currentText()
    
    def _on_save_catalog(self):
        """Handle save catalog button click."""
        # Create a save dialog to get file path
        from PyQt5.QtWidgets import QFileDialog
        from pathlib import Path
        
        # Default save directory
        default_dir = Path.cwd() / 'SpecPlot_Projects'
        default_dir.mkdir(exist_ok=True)
        
        # Get file path
        file_path, _ = QFileDialog.getSaveFileName(
            self, 'Save Catalog', 
            str(default_dir / 'Absorber_Catalog.csv'),
            'CSV Files (*.csv)'
        )
        
        if file_path:
            # Save the catalog
            success = self.absorber_manager.save(file_path)
            if success:
                self.message_window.setText(f"Catalog saved to {file_path}")
            else:
                self.message_window.setText("Error saving catalog")
    
    def _on_load_catalog(self):
        """Handle load catalog button click."""
        # Create a load dialog to get file path
        from PyQt5.QtWidgets import QFileDialog
        
        # Get file path
        file_path, _ = QFileDialog.getOpenFileName(
            self, 'Load Catalog', 
            str(Path.cwd()),
            'CSV Files (*.csv)'
        )
        
        if file_path:
            # Load the catalog
            success = self.absorber_manager.load(file_path)
            if success:
                # Update the absorber widget
                self.absorber_widget.refresh()
                
                # Update the plot with loaded absorbers
                self._clear_lines()
                self.message_window.setText(f"Catalog loaded from {file_path}")
            else:
                self.message_window.setText("Error loading catalog")
    
    def _on_toggle_identified_lines(self):
        """Toggle display of identified lines."""
        if self.identified_lines_active:
            self._hide_identified_lines()
        else:
            self._show_identified_lines()
    
    def _show_identified_lines(self):
        """Show identified lines on the plot."""
        # Set flag
        self.identified_lines_active = True
        
        # Update button style
        self.identified_lines_button.setStyleSheet('background-color: green')
        
        # Get line data from all absorbers
        lines_df = self.absorber_manager.export_lines_to_dataframe()
        
        if lines_df.empty:
            self.message_window.setText("No identified lines to display")
            return
        
        # Get current y limits
        ylim = self.ax.get_ylim()
        
        # Plot each line
        for _, line in lines_df.iterrows():
            xdata = [line['Wave_obs'], line['Wave_obs']]
            ydata = [ylim[0] * 1.5, ylim[1] * 0.75]
            
            # Choose line style based on flag
            linestyle = '--'
            if '[b]' in line['Name']:  # Blended detection
                linestyle = '-.'
            elif '[p]' in line['Name']:  # Low confidence detection
                linestyle = ':'
            
            # Get color from the absorber
            abs_idx, absorber = self.absorber_manager.get_absorber_by_redshift(line['Zabs'])
            color = 'yellow'  # Default if absorber not found
            if absorber:
                color = absorber.color
            
            # Plot the line
            line_plot = self.ax.plot(xdata, ydata, linestyle, color=COLORS[color])[0]
            text_obj = self.ax.text(
                xdata[0], ydata[1], 
                f"{line['Name']}  z={line['Zabs']:.4f}", 
                rotation=90
            )
            
            # Store these objects for later removal
            self.identified_lines.append(line_plot)
            self.identified_texts.append(text_obj)
        
        # Redraw the canvas
        self.canvas.draw()
        
        # Update message
        self.message_window.setText(f"Displaying {len(lines_df)} identified lines")
    
    def _hide_identified_lines(self):
        """Hide identified lines from the plot."""
        # Set flag
        self.identified_lines_active = False
        
        # Update button style
        self.identified_lines_button.setStyleSheet('')
        
        # Remove all identified lines and texts
        for line in self.identified_lines:
            try:
                line.remove()
            except:
                pass
        
        for text in self.identified_texts:
            try:
                text.remove()
            except:
                pass
        
        # Clear lists
        self.identified_lines = []
        self.identified_texts = []
        
        # Redraw the canvas
        self.canvas.draw()
        
        # Update message
        self.message_window.setText("Identified lines hidden")
    
    def _on_plot_absorber(self):
        """Plot lines for current redshift and line list."""
        # Draw line list at current redshift with current color
        self._draw_line_list(
            self.active_line_list, 
            self.active_redshift, 
            self.active_color
        )
    
    def _on_catalog_absorber(self):
        """Add current redshift and line list to catalog."""
        # Create a new absorber
        from ..core.absorber import Absorber
        absorber = Absorber(
            self.active_redshift,
            color=self.active_color,
            line_list=self.active_line_list
        )
        
        # Add to manager
        idx = self.absorber_manager.add_absorber(absorber)
        
        # Update the widget
        self.absorber_widget.refresh()
        
        # Store visualization objects
        line_plots, text_objects = self._get_current_visualization()
        self.absorber_manager.store_visualization(idx, line_plots, text_objects)
        
        # Update message
        self.message_window.setText(
            f"Added absorber at z={self.active_redshift:.4f} with {self.active_line_list} line list"
        )
    
    def _on_refresh(self):
        """Refresh the plot to initial view."""
        self._reset_view()
    
    def _on_absorber_plot(self, index):
        """
        Plot lines for an absorber.
        
        Parameters
        ----------
        index : int
            Index of the absorber
        """
        # Get the absorber
        absorber = self.absorber_manager.get_absorber(index)
        if not absorber:
            return
        
        # Update active values
        self.active_redshift = absorber.redshift
        self.active_line_list = absorber.line_list
        self.active_color = absorber.color
        
        # Update UI
        self.redshift_input.setText(str(self.active_redshift))
        self.line_list_combo.setCurrentText(self.active_line_list)
        self.color_combo.setCurrentText(self.active_color)
        
        # Draw the line list
        self._draw_line_list(
            absorber.line_list,
            absorber.redshift,
            absorber.color
        )
    
    def _on_absorber_hide(self, index):
        """
        Hide lines for an absorber.
        
        Parameters
        ----------
        index : int
            Index of the absorber
        """
        # Get visualization objects
        line_plots, text_objects = self.absorber_manager.get_visualization(index)
        
        # Remove from plot
        if line_plots:
            for line in line_plots:
                try:
                    line.remove()
                except:
                    pass
        
        if text_objects:
            for text in text_objects:
                try:
                    text.remove()
                except:
                    pass
        
        # Clear stored visualization
        self.absorber_manager.store_visualization(index, [], [])
        
        # Update the plot
        self.canvas.draw()
        
        # Update message
        absorber = self.absorber_manager.get_absorber(index)
        if absorber:
            self.message_window.setText(f"Hid absorber at z={absorber.redshift:.4f}")
    
    def _on_absorber_remove(self, index):
        """
        Remove an absorber.
        
        Parameters
        ----------
        index : int
            Index of the absorber
        """
        # Get the absorber
        absorber = self.absorber_manager.get_absorber(index)
        if not absorber:
            return
        
        # Hide first
        self._on_absorber_hide(index)
        
        # Then remove
        redshift = absorber.redshift
        self.absorber_manager.remove_absorber(index)
        
        # Update the widget
        self.absorber_widget.refresh()
        
        # Update message
        self.message_window.setText(f"Removed absorber at z={redshift:.4f}")
    
    def _reset_view(self, full=False):
        """
        Reset the plot view.
        
        Parameters
        ----------
        full : bool, optional
            If True, also clear all lines, by default False
        """
        # Reset axes limits
        self.ax.set_xlim(self.spectrum.initial_xlim)
        self.ax.set_ylim(self.spectrum.initial_ylim)
        
        # Clear temporary markers
        self.ew_points = []
        self.ew_y_values = []
        self.gaussian_points = []
        self.gaussian_y_values = []
        
        # Clear all lines if requested
        if full:
            self._clear_lines()
        
        # Redraw the canvas
        self.canvas.draw()
    
    def _clear_lines(self):
        """Clear all lines from the plot but keep the spectrum."""
        # Store references to flux and error lines
        flux_line = self.ax.lines[0]
        error_line = self.ax.lines[1]
        
        # Clear all lines and texts
        self.ax.lines = []
        self.ax.texts = []
        self.ax.collections = []
        
        # Add back the spectrum lines
        self.ax.lines.append(flux_line)
        self.ax.lines.append(error_line)
        
        # Clear stored visualization objects
        for i in range(self.absorber_manager.count):
            self.absorber_manager.store_visualization(i, [], [])
        
        # Hide identified lines
        if self.identified_lines_active:
            self.identified_lines_active = False
            self.identified_lines_button.setStyleSheet('')
            self.identified_lines = []
            self.identified_texts = []
        
        # Redraw the canvas
        self.canvas.draw()
    
    def _smooth_spectrum(self, increment=True):
        """
        Smooth the spectrum.
        
        Parameters
        ----------
        increment : bool, optional
            If True, increase smoothing; if False, decrease, by default True
        """
        # Store original colors
        flux_color = self.ax.lines[0].get_color()
        error_color = self.ax.lines[1].get_color()
        
        # Apply smoothing
        smoothed_flux, smoothed_error = self.spectrum.smooth(increment)
        
        # Update the plot with preserved colors
        self.ax.lines[0].set_ydata(smoothed_flux)
        self.ax.lines[1].set_ydata(smoothed_error)
        
        # Ensure colors are preserved
        self.ax.lines[0].set_color(flux_color)
        self.ax.lines[1].set_color(error_color)
        
        # Redraw the canvas
        self.canvas.draw()
        
        # Update message
        self.message_window.setText(f"Smoothing level: {self.spectrum.smooth_level}")    
    def _set_manual_ylimits(self):
        """Set manual y-axis limits."""
        # Show input dialog
        from PyQt5.QtWidgets import QInputDialog
        ylim, ok = QInputDialog.getText(
            self, 'Manual Y-Limits', 
            'Enter y-axis limits (min,max):'
        )
        
        if ok:
            try:
                # Parse input
                ymin, ymax = map(float, ylim.split(','))
                
                # Set new limits
                self.ax.set_ylim(ymin, ymax)
                
                # Redraw the canvas
                self.canvas.draw()
                
                # Update message
                self.message_window.setText(f"Y-axis limits set to [{ymin}, {ymax}]")
            except Exception as e:
                self.message_window.setText(f"Invalid input: {e}")
    
    def _mark_doublet(self, doublet_type, wave, y_value):
        """
        Mark a doublet feature on the plot.
        
        Parameters
        ----------
        doublet_type : str
            Type of doublet to mark
        wave : float
            Observed wavelength
        y_value : float
            Y-value for plotting
        """
        from ..core.line_manager import LineIdentifier
        
        # Identify the doublet
        doublet = LineIdentifier.DOUBLETS.get(doublet_type)
        if not doublet:
            self.message_window.setText(f"Unknown doublet type: {doublet_type}")
            return
        
        # Calculate redshift based on first wavelength
        redshift = wave / doublet[0] - 1.0
        
        # Calculate wavelengths for all components
        wavelengths = [w * (1.0 + redshift) for w in doublet]
        
        # Plot the doublet lines
        self.ax.plot([wave, wave, wavelengths[1], wavelengths[1]], 
                     [y_value, y_value + 0.5, y_value + 0.5, y_value], 
                     color='red')
        
        # For triplets like FeII
        if len(wavelengths) > 2:
            for i in range(2, len(wavelengths)):
                self.ax.plot([wavelengths[i-1], wavelengths[i-1], wavelengths[i], wavelengths[i]], 
                             [y_value, y_value + 0.5, y_value + 0.5, y_value], 
                             color='red')
        
        # Add label
        self.ax.text(0.5 * (wavelengths[0] + wavelengths[-1]), y_value + 0.7,
                     f"{doublet_type} z: {redshift:.4f}", 
                     rotation=90, verticalalignment='bottom')
        
        # Redraw the canvas
        self.canvas.draw()
        
        # Update message and redshift
        self.message_window.setText(f"{doublet_type}: z = {redshift:.4f}")
        self.redshift_input.setText(f"{redshift:.4f}")
        self.active_redshift = redshift
    
    def _measure_equivalent_width(self, wave, y_value):
        """
        Measure equivalent width between two points.
        
        Parameters
        ----------
        wave : float
            Wavelength
        y_value : float
            Y-value
        """
        # Add point to list
        self.ew_points.append(wave)
        self.ew_y_values.append(y_value)
        
        # Mark the point
        self.ax.plot(wave, y_value, 'rs', ms=5, markeredgecolor='k')
        self.canvas.draw()
        
        # If we have two points, calculate EW
        if len(self.ew_points) == 2:
            # Sort points
            idx = np.argsort(self.ew_points)
            sorted_points = [self.ew_points[i] for i in idx]
            sorted_y_values = [self.ew_y_values[i] for i in idx]
            
            # Calculate EW
            ew_result = self.spectrum.compute_ew(sorted_points, sorted_y_values)
            
            # Plot continuum
            self.ax.plot(ew_result['wave_slice'], ew_result['continuum'], 'r--')
            
            # Convert to mÅ
            ew_ma = ew_result['ew'] * 1000.0
            err_ma = ew_result['error'] * 1000.0
            
            # Format result message
            message = f"EW [mÅ]: {ew_ma:.1f} ± {err_ma:.1f}"
            
            # Check if there's a line at this position
            line_data = self._find_line_at_position(sorted_points)
            if line_data:
                # Calculate column density
                col_result = self.spectrum.compute_column_density(
                    sorted_points, sorted_y_values, 
                    line_data['f'], line_data['wrest']
                )
                
                log_n = np.log10(col_result['column_density'])
                log_n_err = np.log10((col_result['column_density'] + col_result['error']) / 
                                       col_result['column_density'])
                
                # Add column density to message
                message += f"\nlog N({line_data['ion']}) [cm⁻²]: {log_n:.2f} ± {log_n_err:.2f}"
            
            # Display result
            midpoint = np.mean(sorted_points)
            max_y = max(sorted_y_values)
            self.ax.text(midpoint, max_y + 0.2, message, 
                         rotation=90, verticalalignment='bottom')
            
            # Update message window
            self.message_window.setText(message)
            
            # Clear points for next measurement
            self.ew_points = []
            self.ew_y_values = []
            
            # Redraw canvas
            self.canvas.draw()
    
    def _find_line_at_position(self, wave_range):
        """
        Find a line at a given wavelength range.
        
        Parameters
        ----------
        wave_range : list
            Wavelength range [min, max]
            
        Returns
        -------
        dict or None
            Line data if found, None otherwise
        """
        # Get current line list
        try:
            line_list = LineList.get_instance(self.active_line_list)
            data = line_list.data
            
            # Convert to rest frame
            rest_range = [w / (1.0 + self.active_redshift) for w in wave_range]
            
            # Find lines in this range
            matches = []
            for line in data:
                if rest_range[0] <= line['wrest'] <= rest_range[1]:
                    matches.append(line)
            
            # Return the first match (could be improved to handle multiple matches)
            if matches:
                return matches[0]
        except Exception as e:
            print(f"Error finding line: {e}")
        
        return None
    
    def _fit_gaussian(self, wave, y_value):
        """
        Fit a Gaussian to spectral feature.
        
        Parameters
        ----------
        wave : float
            Wavelength
        y_value : float
            Y-value
        """
        # Add point to list
        self.gaussian_points.append(wave)
        self.gaussian_y_values.append(y_value)
        
        # Mark the point
        self.ax.plot(wave, y_value, 'rs', ms=5, markeredgecolor='k')
        self.canvas.draw()
        
        # If we have three points, fit a Gaussian
        if len(self.gaussian_points) == 3:
            from astropy.modeling import models, fitting
            from scipy.interpolate import splrep, splev
            
            # Sort points by x value
            idx = np.argsort(self.gaussian_points)
            x_points = [self.gaussian_points[i] for i in idx]
            y_points = [self.gaussian_y_values[i] for i in idx]
            
            # Select data in range
            mask = ((self.spectrum.wavelength >= x_points[0]) & 
                    (self.spectrum.wavelength <= x_points[2]))
            x_data = self.spectrum.wavelength[mask]
            y_data = self.spectrum.flux[mask]
            
            # Fit a continuum through the first and last points
            spline = splrep([x_points[0], x_points[2]], [y_points[0], y_points[2]], k=1)
            continuum = splev(x_data, spline)
            
            # Initialize a Gaussian model
            amplitude = y_points[1] - np.interp(x_points[1], [x_points[0], x_points[2]], 
                                             [y_points[0], y_points[2]])
            mean = x_points[1]
            stddev = 0.5 * (x_points[2] - x_points[0])
            
            # Check if absorption or emission
            is_absorption = (y_points[1] < y_points[0]) and (y_points[1] < y_points[2])
            
            # Prepare data for fitting
            if is_absorption:
                norm_data = 1.0 - (y_data / continuum)
                g_init = models.Gaussian1D(amplitude=abs(amplitude), mean=mean, stddev=stddev)
            else:
                norm_data = (y_data / continuum) - 1.0
                g_init = models.Gaussian1D(amplitude=amplitude, mean=mean, stddev=stddev)
            
            # Fit the model
            fitter = fitting.LevMarLSQFitter()
            g_fit = fitter(g_init, x_data, norm_data)
            
            # Generate the final model
            if is_absorption:
                model_flux = (1.0 - g_fit(x_data)) * continuum
            else:
                model_flux = (1.0 + g_fit(x_data)) * continuum
            
            # Plot the result
            self.ax.plot(x_data, model_flux, 'r-')
            
            # Display fit parameters
            mid_x = np.mean(x_data)
            max_y = max(y_points) + 0.2
            self.ax.text(mid_x, max_y, f"Amp: {g_fit.amplitude.value:.3f}", 
                         rotation=90, verticalalignment='bottom')
            self.ax.text(mid_x + stddev, max_y, f"Center: {g_fit.mean.value:.3f}", 
                         rotation=90, verticalalignment='bottom')
            self.ax.text(mid_x - stddev, max_y, f"Sigma: {g_fit.stddev.value:.3f}", 
                         rotation=90, verticalalignment='bottom')
            
            # Update message window
            self.message_window.setText(
                f"Gaussian fit:\n"
                f"Amplitude: {g_fit.amplitude.value:.3f}\n"
                f"Center: {g_fit.mean.value:.3f}\n"
                f"Sigma: {g_fit.stddev.value:.3f}"
            )
            
            # Clear points for next fit
            self.gaussian_points = []
            self.gaussian_y_values = []
            
            # Redraw canvas
            self.canvas.draw()
    
    def _show_help(self):
        """Show the help window."""
        self.help_window = HelpWindow(self)
        self.help_window.show()
    
    def _open_vstack(self):
        """Open the vstack display."""
        from .vstack_display import VStackDisplay
        
        # Check if a vstack exists for this redshift
        message = self._check_existing_linelist()
        if message:
            reply = QMessageBox.question(
                self, "Reevaluate",
                f"{message}\nReevaluate and overwrite?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        # Open vstack
        self.vstack = VStackDisplay(
            self.spectrum, 
            self.active_redshift,
            self.active_line_list
        )
        
        # Connect signals
        self.vstack.completed.connect(self._on_vstack_completed)
        
        # Show the vstack
        self.vstack.show()
    
    def _open_vstack_custom_velocity(self):
        """Open vstack with custom velocity limits."""
        from PyQt5.QtWidgets import QInputDialog
        
        # Get velocity limits
        vlim, ok = QInputDialog.getText(
            self, 'Manual Velocity Limits', 
            'Enter velocity limits (min,max) in km/s:'
        )
        
        if ok:
            try:
                # Parse input
                vmin, vmax = map(float, vlim.split(','))
                
                # Open vstack with custom limits
                self._open_vstack()
                if hasattr(self, 'vstack'):
                    self.vstack.set_velocity_limits([vmin, vmax])
            except Exception as e:
                self.message_window.setText(f"Invalid input: {e}")
    
    def _open_manual_transition(self, wavelength):
        """
        Open manual transition selector.
        
        Parameters
        ----------
        wavelength : float
            Observed wavelength
        """
        from .line_selector import LineSelector
        
        # Create and show the selector
        self.line_selector = LineSelector(
            wavelength, 
            self.active_line_list,
            self.active_redshift
        )
        
        # Connect signals
        self.line_selector.lineSelected.connect(self._on_line_selected)
        
        # Show the selector
        self.line_selector.show()
    
    def _on_line_selected(self, data):
        """
        Handle line selection.
        
        Parameters
        ----------
        data : dict
            Information about the selected line
        """
        # Update redshift
        self.active_redshift = data['redshift']
        self.redshift_input.setText(f"{self.active_redshift:.4f}")
        
        # Update line list if changed
        if data.get('line_list'):
            self.active_line_list = data['line_list']
            self.line_list_combo.setCurrentText(self.active_line_list)
        
        # Draw the line list
        self._draw_line_list(
            self.active_line_list,
            self.active_redshift,
            self.active_color
        )
        
        # Update message
        self.message_window.setText(
            f"Selected {data['ion']} at z={self.active_redshift:.4f}"
        )
    
    def _on_vstack_completed(self, results):
        """
        Handle vstack completion.
        
        Parameters
        ----------
        results : dict
            Results from the vstack
        """
        # Update the absorber manager
        if 'absorber' in results:
            # Find existing absorber or create new one
            match = self.absorber_manager.get_absorber_by_redshift(results['absorber'].redshift)
            if match:
                idx, _ = match
                self.absorber_manager.remove_absorber(idx)
            
            # Add the new absorber
            self.absorber_manager.add_absorber(results['absorber'])
            
            # Update the widget
            self.absorber_widget.refresh()
            
            # Update message
            self.message_window.setText(
                f"Updated absorber at z={results['absorber'].redshift:.4f} "
                f"with {results['absorber'].line_count} lines"
            )
        
        # Re-enable identified lines if they were active
        if self.identified_lines_active:
            self._show_identified_lines()
    
    def _check_existing_linelist(self):
        """
        Check if a line list exists for the current redshift.
        
        Returns
        -------
        str or None
            Message if line list exists, None otherwise
        """
        # Get line data
        lines_df = self.absorber_manager.export_lines_to_dataframe()
        
        # Check if redshift exists
        if not lines_df.empty:
            redshifts = lines_df['Zabs'].unique()
            if self.active_redshift in redshifts:
                count = lines_df[lines_df['Zabs'] == self.active_redshift].shape[0]
                return f"Current redshift z={self.active_redshift:.4f} already has {count} identified lines."
        
        return None
    
    def _draw_line_list(self, line_list_name, redshift, color_name):
        """
        Draw absorption lines for a line list at a given redshift.
        
        Parameters
        ----------
        line_list_name : str
            Name of the line list
        redshift : float
            Redshift
        color_name : str
            Name of the color to use
        
        Returns
        -------
        tuple
            Lists of line plots and text objects
        """
        # Load the line list
        line_list = LineList.get_instance(line_list_name)
        if not line_list:
            self.message_window.setText(f"Error loading line list: {line_list_name}")
            return [], []
        
        # Get data
        data = line_list.data
        if not data:
            self.message_window.setText(f"No data in line list: {line_list_name}")
            return [], []
        
        # Get plot limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        color = COLORS.get(color_name, 'white')
        
        # Store visualization objects
        line_plots = []
        text_objects = []
        
        # Draw each line in the wavelength range
        for line in data:
            # Calculate observed wavelength
            obs_wave = line['wrest'] * (1.0 + redshift)
            
            # Skip if outside current view
            if obs_wave < xlim[0] or obs_wave > xlim[1]:
                continue
            
            # Draw the line
            y_values = [0, ylim[1]]
            line_plot = self.ax.plot([obs_wave, obs_wave], y_values, '--', color=color)[0]
            
            # Add label
            text_obj = self.ax.text(
                obs_wave, 0.75 * ylim[1], 
                f"{line['ion']}  z={redshift:.4f}", 
                rotation=90, picker=True
            )
            
            # Store objects
            line_plots.append(line_plot)
            text_objects.append(text_obj)
        
        # Redraw the canvas
        self.canvas.draw()
        
        # Update message
        self.message_window.setText(
            f"Plotted {len(line_plots)} lines from {line_list_name} at z={redshift:.4f}"
        )
        
        return line_plots, text_objects
    
    def _get_current_visualization(self):
        """
        Get current visualization objects.
        
        Returns
        -------
        tuple
            Lists of line plots and text objects
        """
        # Filter line plots and texts that belong to the current redshift
        line_plots = [line for line in self.ax.lines[2:]]
        text_objects = [text for text in self.ax.texts]
        
        return line_plots, text_objects