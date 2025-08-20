"""
Interactive masking and continuum fitting tool for spectral analysis.

This module provides a PyQt5-based GUI for interactive selection of mask regions,
continuum fitting with both polynomial and spline methods.
"""

import sys
import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, 
                           QWidget, QPushButton, QLabel, QSpinBox, QCheckBox, 
                           QTabWidget, QGroupBox, QFormLayout, QDoubleSpinBox,
                           QMessageBox, QStatusBar, QDialog, QDialogButtonBox,
                           QComboBox,QScrollArea,QGridLayout,QTextBrowser)
from PyQt5.QtCore import Qt
from scipy.interpolate import splrep, splev
import warnings
import datetime
import os
import pkg_resources
from PyQt5.QtCore import QUrl
from PyQt5.QtGui import QIcon
# Import the minimum required functions
try:
    from rbcodes.IGM.rb_iter_contfit import fit_optimal_polynomial, rb_iter_contfit
except ImportError:
    try:
        from IGM.rb_iter_contfit import fit_optimal_polynomial, rb_iter_contfit
    except ImportError:
        print("Warning: rb_iter_contfit not found. Some functionality may be limited.")

class MplCanvas(FigureCanvasQTAgg):
    """Matplotlib canvas for embedding in PyQt5."""
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        # Two panels: flux and normalized flux
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        
        super(MplCanvas, self).__init__(self.fig)
        self.setMinimumSize(600, 400)  # Set minimum size for better visibility

class InteractiveContinuumFitWindow(QMainWindow):
    """Main window for interactive masking and continuum fitting."""
    
    def __init__(self, wave, flux, error, velocity=None, existing_masks=None,
                 order=3, use_weights=False, domain=None):
        super(InteractiveContinuumFitWindow, self).__init__()
        
        # Store data
        self.wave = np.array(wave)
        self.flux = np.array(flux)
        self.error = np.array(error)
        
        # Handle velocity/wavelength x-axis setup
        if velocity is not None:
            self.velocity = np.array(velocity)
            self.x_axis = self.velocity
            self.x_axis_type = 'velocity'
            self.x_axis_label = 'Velocity (km/s)'
            self.x_axis_unit = 'km/s'
        else:
            self.velocity = None
            self.x_axis = self.wave
            self.x_axis_type = 'wavelength'
            self.x_axis_label = 'Wavelength (Å)'
            self.x_axis_unit = 'Å'
        
        self.masks = list(existing_masks) if existing_masks is not None else []
        self.domain = domain if domain is not None else [min(self.x_axis), max(self.x_axis)]
        self.velocity_domain = domain if domain is not None else [min(self.velocity), max(self.velocity)] if self.velocity is not None else None
        self.wave_domain = [min(self.wave), max(self.wave)]
        
        # Working variables
        self.current_click = None  # Store first click of a pair
        self.continuum = None
        self.normalized_flux = None
        self.spline_points = []   # (x, y) coordinates for spline points
        self.fit_params = {
            'method': 'polynomial',  # 'polynomial' or 'spline'
            'order': order,
            'use_weights': use_weights,
            'optimize_cont': True,
            'min_order': 1,
            'max_order': 6,
            'sigma': 3.0,
            'spline_degree': 3,      # Cubic spline by default
            'spline_smooth': 0.0,    # Smoothing factor for spline
            'median_window': 5.0,    # Window for median calculation (pixels or units)
            'window_type': 'pixels'  # 'pixels' or 'units'
        }
        self.fit_result = None
        
        # Setup UI
        self.setup_ui()

        # Ensure canvas has focus capabilities
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        
        # Initialize plots
        self.update_plots()

        # Give focus to the canvas
        self.canvas.setFocus()
    
    def setup_ui(self):
        """Set up the user interface."""
        self.setWindowTitle("Interactive Continuum Fitting")
        self.setGeometry(100, 100, 1000, 700)  # Slightly smaller default size
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
         
        # Create plotting area
        self.canvas = MplCanvas(self, width=8, height=6)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        #Add a help window
        self.add_help_button_to_toolbar()
        
        # Create plot layout
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        
        # Add plot area to main layout
        main_layout.addLayout(plot_layout, 3)  # Plot takes 3/4 of width
        
        # Create sidebar with tabs
        sidebar = QTabWidget()
        sidebar.setMaximumWidth(300)  # Limit width of sidebar
        
        # Tab 1: Basic controls - make scrollable
        basic_tab = QWidget()
        basic_scroll = QScrollArea()
        basic_scroll.setWidgetResizable(True)
        basic_scroll.setFrameShape(QScrollArea.NoFrame)
        basic_scroll_content = QWidget()
        basic_layout = QVBoxLayout(basic_scroll_content)
        basic_layout.setSpacing(8)  # Reduce spacing between widgets
        basic_layout.setContentsMargins(8, 8, 8, 8)  # Reduce margins
        basic_scroll.setWidget(basic_scroll_content)
        basic_tab_layout = QVBoxLayout(basic_tab)
        basic_tab_layout.addWidget(basic_scroll)
        basic_tab_layout.setContentsMargins(0, 0, 0, 0)
    
        # Method selection group
        method_group = QGroupBox("Fitting Method")
        method_layout = QVBoxLayout()
        method_layout.setSpacing(4)  # Reduce spacing
        method_layout.setContentsMargins(8, 8, 8, 8)  # Reduce margins
        
        self.method_buttons = QtWidgets.QButtonGroup(self)
        self.poly_radio = QtWidgets.QRadioButton("Polynomial")
        self.spline_radio = QtWidgets.QRadioButton("Spline")
        self.poly_radio.setChecked(True)  # Default to polynomial
        
        self.method_buttons.addButton(self.poly_radio, 1)
        self.method_buttons.addButton(self.spline_radio, 2)
        
        method_layout.addWidget(self.poly_radio)
        method_layout.addWidget(self.spline_radio)
        
        # Connect signals
        self.poly_radio.toggled.connect(self.toggle_method)
        
        method_group.setLayout(method_layout)
        basic_layout.addWidget(method_group)
        
        # X-axis selection (only show if both wavelength and velocity are available)
        if self.velocity is not None:
            x_axis_group = QGroupBox("X-Axis Type")
            x_axis_layout = QVBoxLayout()
            x_axis_layout.setSpacing(4)  # Reduce spacing
            x_axis_layout.setContentsMargins(8, 8, 8, 8)  # Reduce margins
            
            self.x_axis_buttons = QtWidgets.QButtonGroup(self)
            self.velocity_radio = QtWidgets.QRadioButton("Velocity")
            self.wavelength_radio = QtWidgets.QRadioButton("Wavelength")
            self.velocity_radio.setChecked(True)
            
            self.x_axis_buttons.addButton(self.velocity_radio, 1)
            self.x_axis_buttons.addButton(self.wavelength_radio, 2)
            
            x_axis_layout.addWidget(self.velocity_radio)
            x_axis_layout.addWidget(self.wavelength_radio)
            
            # Connect signals
            self.velocity_radio.toggled.connect(self.toggle_x_axis)
            
            x_axis_group.setLayout(x_axis_layout)
            basic_layout.addWidget(x_axis_group)
        
        # Mask controls group
        mask_group = QGroupBox("Mask Controls")
        mask_layout = QVBoxLayout()
        mask_layout.setSpacing(4)  # Reduce spacing
        mask_layout.setContentsMargins(8, 8, 8, 8)  # Reduce margins
        
        # Add instructions label
        mask_instructions = QLabel("Left-click pairs to add masks\nRight-click pairs to remove masks")
        mask_layout.addWidget(mask_instructions)
        
        # Create a horizontal layout for mask buttons (2x2 grid)
        mask_buttons_layout = QGridLayout()
        mask_buttons_layout.setSpacing(4)  # Reduced spacing between buttons
        
        # Mask buttons
        self.reset_btn = QPushButton("Reset Masks")
        self.reset_btn.clicked.connect(self.reset_masks)
        
        self.undo_btn = QPushButton("Undo Last")
        self.undo_btn.clicked.connect(self.undo_last_mask)
        
        self.auto_btn = QPushButton("Auto-Mask")
        self.auto_btn.clicked.connect(self.auto_mask)
        
        self.manual_btn = QPushButton("Manual Entry")
        self.manual_btn.clicked.connect(self.manual_mask_entry)
        
        # Add buttons to grid layout, 2 buttons per row
        mask_buttons_layout.addWidget(self.reset_btn, 0, 0)
        mask_buttons_layout.addWidget(self.undo_btn, 0, 1)
        mask_buttons_layout.addWidget(self.auto_btn, 1, 0)
        mask_buttons_layout.addWidget(self.manual_btn, 1, 1)
        
        mask_layout.addLayout(mask_buttons_layout)
        mask_group.setLayout(mask_layout)
        basic_layout.addWidget(mask_group)
        
        # Create a stacked widget to hold method-specific controls
        self.method_stack = QtWidgets.QStackedWidget()
    
        # Create polynomial fitting panel (with scroll area)
        poly_panel = QWidget()
        poly_scroll = QScrollArea()
        poly_scroll.setWidgetResizable(True)
        poly_scroll.setFrameShape(QScrollArea.NoFrame)
        poly_scroll_content = QWidget()
        poly_layout = QVBoxLayout(poly_scroll_content)
        poly_layout.setSpacing(4)  # Reduce spacing
        poly_layout.setContentsMargins(4, 4, 4, 4)  # Reduce margins
        poly_scroll.setWidget(poly_scroll_content)
        poly_panel_layout = QVBoxLayout(poly_panel)
        poly_panel_layout.addWidget(poly_scroll)
        poly_panel_layout.setContentsMargins(0, 0, 0, 0)
    
        # Polynomial fitting controls group
        poly_group = QGroupBox("Polynomial Fitting Options")
        poly_form_layout = QFormLayout()
        poly_form_layout.setSpacing(4)  # Reduce spacing
        poly_form_layout.setContentsMargins(8, 8, 8, 8)  # Reduce margins
    
        self.order_spin = QSpinBox()
        self.order_spin.setRange(0, 10)
        self.order_spin.setValue(self.fit_params['order'])
    
        self.use_weights_check = QCheckBox()
        self.use_weights_check.setChecked(self.fit_params['use_weights'])
    
        self.optimize_check = QCheckBox()
        self.optimize_check.setChecked(True)  
        self.optimize_check.toggled.connect(self.toggle_optimization)
    
        poly_form_layout.addRow("Polynomial Order:", self.order_spin)
        poly_form_layout.addRow("Use Weights:", self.use_weights_check)
        poly_form_layout.addRow("Auto Optimize:", self.optimize_check)
        poly_group.setLayout(poly_form_layout)
        poly_layout.addWidget(poly_group)
    
        # Create spline fitting panel (with scroll area)
        spline_panel = QWidget()
        spline_scroll = QScrollArea()
        spline_scroll.setWidgetResizable(True)
        spline_scroll.setFrameShape(QScrollArea.NoFrame)
        spline_scroll_content = QWidget()
        spline_layout = QVBoxLayout(spline_scroll_content)
        spline_layout.setSpacing(4)  # Reduce spacing
        spline_layout.setContentsMargins(4, 4, 4, 4)  # Reduce margins
        spline_scroll.setWidget(spline_scroll_content)
        spline_panel_layout = QVBoxLayout(spline_panel)
        spline_panel_layout.addWidget(spline_scroll)
        spline_panel_layout.setContentsMargins(0, 0, 0, 0)
    
        # Spline fitting controls group
        spline_group = QGroupBox("Spline Fitting Options")
        spline_form_layout = QFormLayout()
        spline_form_layout.setSpacing(4)  # Reduce spacing
        spline_form_layout.setContentsMargins(8, 8, 8, 8)  # Reduce margins
    
        self.spline_degree_spin = QSpinBox()
        self.spline_degree_spin.setRange(1, 5)
        self.spline_degree_spin.setValue(self.fit_params['spline_degree'])
    
        self.spline_smooth_spin = QDoubleSpinBox()
        self.spline_smooth_spin.setRange(0.0, 100.0)
        self.spline_smooth_spin.setSingleStep(0.1)
        self.spline_smooth_spin.setValue(self.fit_params['spline_smooth'])
    
        # Add median window control
        self.median_window_spin = QDoubleSpinBox()
        self.median_window_spin.setRange(1.0, 100.0)
        self.median_window_spin.setSingleStep(1.0)
        self.median_window_spin.setValue(self.fit_params['median_window'])
    
        # Add window type selector
        self.window_type_combo = QComboBox()
        self.window_type_combo.addItems(["Pixels", "Units"])
        self.window_type_combo.setCurrentText(self.fit_params['window_type'].capitalize())
        self.window_type_combo.currentTextChanged.connect(self.update_window_type)
    
        self.spline_clear_btn = QPushButton("Clear Points")
        self.spline_clear_btn.clicked.connect(self.clear_spline_points)
    
        spline_form_layout.addRow("Spline Degree:", self.spline_degree_spin)
        spline_form_layout.addRow("Smoothing:", self.spline_smooth_spin)
        spline_form_layout.addRow("Median Window:", self.median_window_spin)
        spline_form_layout.addRow("Window Type:", self.window_type_combo)
        spline_form_layout.addRow("", self.spline_clear_btn)
    
        # Add spline instructions
        spline_instructions = QLabel("For spline fitting:\n- Left-click adds median point\n- Right-click removes closest point\n- Press 'b' for exact point")
        spline_form_layout.addRow(spline_instructions)
    
        spline_group.setLayout(spline_form_layout)
        spline_layout.addWidget(spline_group)
    
        # Add panels to the stacked widget
        self.method_stack.addWidget(poly_panel)   # Index 0
        self.method_stack.addWidget(spline_panel)  # Index 1
    
        # Add stacked widget to the main layout
        basic_layout.addWidget(self.method_stack)
    
        # Connect radio buttons to switch stacked widget
        self.poly_radio.toggled.connect(lambda checked: self.method_stack.setCurrentIndex(0) if checked else None)
        self.spline_radio.toggled.connect(lambda checked: self.method_stack.setCurrentIndex(1) if checked else None)

        # Add this right after them:
        self.poly_radio.toggled.connect(lambda checked: self.canvas.setFocus() if checked else None)
        self.spline_radio.toggled.connect(lambda checked: self.canvas.setFocus() if checked else None)

        # For velocity/wavelength radio buttons (if they exist)
        if hasattr(self, 'velocity_radio') and hasattr(self, 'wavelength_radio'):
            self.velocity_radio.toggled.connect(lambda checked: self.canvas.setFocus() if checked else None)
            self.wavelength_radio.toggled.connect(lambda checked: self.canvas.setFocus() if checked else None)

    
        # Set initial panel based on default method
        self.method_stack.setCurrentIndex(0 if self.fit_params['method'] == 'polynomial' else 1)
        
        # Create grid layout for action buttons
        action_buttons_layout = QGridLayout()
        action_buttons_layout.setSpacing(4)  # Reduced spacing
        
        # Fit button (applies to both polynomial and spline)
        self.fit_btn = QPushButton("Fit Continuum")
        self.fit_btn.clicked.connect(self.fit_continuum)
        
        # Reset All button
        self.reset_all_btn = QPushButton("Reset All")
        self.reset_all_btn.clicked.connect(self.reset_all)
        
        # Action buttons
        self.accept_btn = QPushButton("Accept & Return")
        self.accept_btn.clicked.connect(self.accept_results)
        self.accept_btn.setEnabled(False)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.cancel)
        
        # Add buttons to grid layout, 2 buttons per row
        action_buttons_layout.addWidget(self.fit_btn, 0, 0)
        action_buttons_layout.addWidget(self.reset_all_btn, 0, 1)
        action_buttons_layout.addWidget(self.accept_btn, 1, 0)
        action_buttons_layout.addWidget(self.cancel_btn, 1, 1)
        
        # After button signal connections in setup_ui, add these:
        self.fit_btn.clicked.connect(lambda: self.canvas.setFocus())
        self.reset_btn.clicked.connect(lambda: self.canvas.setFocus())
        self.reset_all_btn.clicked.connect(lambda: self.canvas.setFocus())
        self.auto_btn.clicked.connect(lambda: self.canvas.setFocus())
        self.spline_clear_btn.clicked.connect(lambda: self.canvas.setFocus())


        basic_layout.addLayout(action_buttons_layout)
        basic_layout.addStretch()  # Add stretch to push everything to the top
        
        # Tab 2: Advanced options - make scrollable
        advanced_tab = QWidget()
        advanced_scroll = QScrollArea()
        advanced_scroll.setWidgetResizable(True)
        advanced_scroll.setFrameShape(QScrollArea.NoFrame)
        advanced_scroll_content = QWidget()
        advanced_layout = QFormLayout(advanced_scroll_content)
        advanced_layout.setSpacing(8)  # Reduce spacing
        advanced_layout.setContentsMargins(8, 8, 8, 8)  # Reduce margins
        advanced_scroll.setWidget(advanced_scroll_content)
        advanced_tab_layout = QVBoxLayout(advanced_tab)
        advanced_tab_layout.addWidget(advanced_scroll)
        advanced_tab_layout.setContentsMargins(0, 0, 0, 0)
        
        self.min_order_spin = QSpinBox()
        self.min_order_spin.setRange(0, 10)
        self.min_order_spin.setValue(self.fit_params['min_order'])
        
        self.max_order_spin = QSpinBox()
        self.max_order_spin.setRange(1, 15)
        self.max_order_spin.setValue(self.fit_params['max_order'])
        
        self.sigma_spin = QDoubleSpinBox()
        self.sigma_spin.setRange(0.5, 10.0)
        self.sigma_spin.setSingleStep(0.5)
        self.sigma_spin.setValue(self.fit_params['sigma'])
        
        # Add new spinboxes for auto-masking parameters
        self.auto_mask_sigma_spin = QDoubleSpinBox()
        self.auto_mask_sigma_spin.setRange(0.5, 10.0)
        self.auto_mask_sigma_spin.setSingleStep(0.5)
        self.auto_mask_sigma_spin.setValue(2.0)  # Separate default for auto-masking
        
        self.min_mask_width_spin = QDoubleSpinBox()
        self.min_mask_width_spin.setRange(1.0, 500.0)
        self.min_mask_width_spin.setSingleStep(10.0)
        self.min_mask_width_spin.setValue(20.0)  # Default 20 km/s minimum mask width
        # Set suffix based on current x-axis type
        self.min_mask_width_spin.setSuffix(f" {self.x_axis_unit}")


        self.mask_gap_spin = QDoubleSpinBox()
        self.mask_gap_spin.setRange(1.0, 500.0)
        self.mask_gap_spin.setSingleStep(10.0)
        self.mask_gap_spin.setValue(100.0)  # Default 100 km/s grouping window
        self.mask_gap_spin.setSuffix(f" {self.x_axis_unit}")


        advanced_layout.addRow("Min Polynomial Order:", self.min_order_spin)
        advanced_layout.addRow("Max Polynomial Order:", self.max_order_spin)
        advanced_layout.addRow("Fitting Sigma:", self.sigma_spin)
        advanced_layout.addRow("Auto-Mask Sigma:", self.auto_mask_sigma_spin)
        advanced_layout.addRow("Min Mask Width:", self.min_mask_width_spin)
        advanced_layout.addRow("Grouping Window:", self.mask_gap_spin)
        
        # Add tabs to sidebar
        sidebar.addTab(basic_tab, "Basic")
        sidebar.addTab(advanced_tab, "Advanced")
        
        # Add sidebar to main layout
        main_layout.addWidget(sidebar, 1)  # Sidebar takes 1/4 of width
        
        # Status bar
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.showMessage("Ready. Left-click pairs to add masks, right-click pairs to remove masks.")
        
        # Connect canvas events
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
        
        # Focus on canvas whenever clicked
        self.canvas.mpl_connect('button_press_event', lambda event: self.canvas.setFocus())
        
        # Set initial UI state based on default method
        self.toggle_optimization(True)  # Initial state matches the default


    def show_help_dialog(self):
        """Open the help dialog."""
        dialog = HelpDialog(self)
        dialog.exec_()

    # After creating the toolbar, add a help button:
    def add_help_button_to_toolbar(self):
        """Add a help button to the matplotlib toolbar."""
        spacer = QWidget()
        spacer.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.toolbar.addWidget(spacer)
        
        # Create help button with icon if available
        help_button = QPushButton("?")
        help_button.setToolTip("Show Help")
        help_button.setMaximumWidth(30)
        help_button.clicked.connect(self.show_help_dialog)
        
        # Try to set an icon if available
        try:
            icon_path = pkg_resources.resource_filename(
                'rbcodes.GUIs.interactive_continuum_fit_help', 'images/help_icon.png'
            )
            if os.path.exists(icon_path):
                help_button.setIcon(QIcon(icon_path))
                help_button.setText("")  # Remove text if icon is used
        except (ImportError, pkg_resources.DistributionNotFound):
            pass  # Use text-only button if icon not found
        
        self.toolbar.addWidget(help_button)    


    def toggle_method(self, checked):
        """Toggle between polynomial and spline fitting methods."""
        if self.poly_radio.isChecked():
            self.fit_params['method'] = 'polynomial'
            self.method_stack.setCurrentIndex(0)
            self.statusBar.showMessage(f"Polynomial fitting mode. Mask absorption regions and set polynomial order. X-axis: {self.x_axis_type}.")
        else:
            self.fit_params['method'] = 'spline'
            self.method_stack.setCurrentIndex(1)
            self.statusBar.showMessage(f"Spline fitting mode. Left-click for median point, right-click to remove, 'b' for exact point. X-axis: {self.x_axis_type}.")
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        
    def toggle_x_axis(self, checked):
        """Toggle between velocity and wavelength x-axis."""
        # Get current view limits before switching
        ax_flux = self.canvas.axes[0]
        current_xlim = ax_flux.get_xlim()
        
        # Calculate equivalent range in the other coordinate system
        if self.velocity_radio.isChecked():
            # Switching to velocity
            self.x_axis = self.velocity
            self.x_axis_type = 'velocity'
            self.x_axis_label = 'Velocity (km/s)'
            self.x_axis_unit = 'km/s'
            
            # Set velocity-specific default values
            self.min_mask_width_spin.setRange(1.0, 500.0)
            self.min_mask_width_spin.setSingleStep(10.0)
            self.min_mask_width_spin.setValue(20.0)  # Default 20 km/s minimum mask width
            self.mask_gap_spin.setRange(1.0, 500.0)
            self.mask_gap_spin.setSingleStep(10.0)
            self.mask_gap_spin.setValue(100.0)  # Default 100 km/s grouping window

            
            # Convert current wavelength view to velocity if needed
            if not hasattr(self, 'prev_was_velocity') or not self.prev_was_velocity:
                # If we were in wavelength, convert the current view to equivalent velocity
                # Find indices of closest wavelength values
                left_idx = np.argmin(np.abs(self.wave - current_xlim[0]))
                right_idx = np.argmin(np.abs(self.wave - current_xlim[1]))
                # Get corresponding velocity values
                new_xlim = [self.velocity[left_idx], self.velocity[right_idx]]
                # Update domain
                self.domain = new_xlim
                self.velocity_domain = new_xlim
            else:
                # Use stored velocity domain
                self.domain = self.velocity_domain
            
            self.prev_was_velocity = True
        else:
            # Switching to wavelength
            self.x_axis = self.wave
            self.x_axis_type = 'wavelength'
            self.x_axis_label = 'Wavelength (Å)'
            self.x_axis_unit = 'Å'
            
            # Set wavelength-specific default values
            self.min_mask_width_spin.setRange(0.001, 500.0)
            self.min_mask_width_spin.setSingleStep(10.0)
            self.min_mask_width_spin.setValue(1.0)  # 1 Å for wavelength
            
            self.mask_gap_spin.setRange(0.001, 500.0)
            self.mask_gap_spin.setSingleStep(10.0)
            self.mask_gap_spin.setValue(5.0)  # 5 Å for wavelength


            
            # Convert current velocity view to wavelength if needed
            if hasattr(self, 'prev_was_velocity') and self.prev_was_velocity:
                # If we were in velocity, convert the current view to equivalent wavelength
                # Find indices of closest velocity values
                left_idx = np.argmin(np.abs(self.velocity - current_xlim[0]))
                right_idx = np.argmin(np.abs(self.velocity - current_xlim[1]))
                # Get corresponding wavelength values
                new_xlim = [self.wave[left_idx], self.wave[right_idx]]
                # Update domain
                self.domain = new_xlim
                self.wave_domain = new_xlim
            else:
                # Use stored wavelength domain
                self.domain = self.wave_domain
            
            self.prev_was_velocity = False
        
        # Update UI elements with correct units
        self.min_mask_width_spin.setSuffix(f" {self.x_axis_unit}")
        self.mask_gap_spin.setSuffix(f" {self.x_axis_unit}")
        
        # Update plots with new x-axis using the calculated domain
        self.update_plots()
        
        # Ensure canvas has focus capabilities
        self.canvas.setFocus()

        
        # Update status bar
        if self.fit_params['method'] == 'polynomial':
            self.statusBar.showMessage(f"Polynomial fitting mode. Mask absorption regions and set polynomial order. X-axis: {self.x_axis_type}.")
        else:
            self.statusBar.showMessage(f"Spline fitting mode. Left-click for median point, right-click to remove, 'b' for exact point. X-axis: {self.x_axis_type}.")    




    def update_window_type(self, text):
        """Update the window type for median calculation."""
        self.fit_params['window_type'] = text.lower()
        self.statusBar.showMessage(f"Window type set to {text}.")
        
    
    def toggle_optimization(self, checked):
        """Toggle optimization options."""
        if self.fit_params['method'] == 'polynomial':
            self.order_spin.setEnabled(not checked)
            self.min_order_spin.setEnabled(checked)
            self.max_order_spin.setEnabled(checked)
    
    def update_plots(self,**kwargs):
        """Update the plots with current data and masks."""

        # Get the input range or use domain
        input_xrange = kwargs.get('input_xrange', self.domain)

        # If we don't have a specific range, use domain
        if input_xrange is None:
            if self.x_axis_type == 'velocity' and self.velocity_domain is not None:
                input_xrange = self.velocity_domain
            else:
                input_xrange = self.wave_domain
        
        # Keep track of current domain for the active coordinate system
        if self.x_axis_type == 'velocity':
            self.velocity_domain = input_xrange
        else:
            self.wave_domain = input_xrange
        
        # Update the main domain property
        self.domain = input_xrange

        # Store current y limits before clearing
        ax_flux, ax_norm = self.canvas.axes
        current_flux_ylim = ax_flux.get_ylim() if ax_flux.get_ylim() != (0, 1.5) else None
        current_norm_ylim = ax_norm.get_ylim() if ax_norm.get_ylim() != (0, 1.5) else None

        # Clear axes
        for ax in self.canvas.axes:
            ax.clear()
        
        # Get axes
        ax_flux, ax_norm = self.canvas.axes
        
        # Plot original spectrum and error as separate step plots
        ax_flux.step(self.x_axis, self.flux, 'k-', where='mid', lw=0.5, label='Flux')
        ax_flux.step(self.x_axis, self.error, 'r-', where='mid', lw=0.5, alpha=0.7, label='Error')
        
        # Plot continuum if available
        if self.continuum is not None:
            ax_flux.plot(self.x_axis, self.continuum, 'r-', linewidth=2, label='Continuum')
        
        # Mark masked regions
        self.plot_masks(ax_flux)
        
        # Plot spline points if in spline mode
        if self.fit_params['method'] == 'spline' and self.spline_points:
            x_points, y_points = zip(*self.spline_points)
            ax_flux.plot(x_points, y_points, 'go', ms=6, alpha=0.8, label='Spline Points')
        
        # Plot normalized spectrum if available
        if self.normalized_flux is not None:
            ax_norm.step(self.x_axis, self.normalized_flux, 'k-', where='mid', lw=0.5)
            ax_norm.axhline(y=1.0, color='r', linestyle='--')
            
            # Mark masked regions in normalized plot
            self.plot_masks(ax_norm)
        else:
            ax_norm.text(0.5, 0.5, "Normalized spectrum will appear here after fitting", 
                       ha='center', va='center', transform=ax_norm.transAxes)
        
        # Set labels and limits
        ax_flux.set_ylabel('Flux')
        ax_flux.set_title('Interactive Continuum Fitting')
        #ax_flux.legend(loc='upper right')
        
        ax_norm.set_xlabel(self.x_axis_label)
        ax_norm.set_ylabel('Normalized Flux')
        ax_norm.set_xlim(input_xrange)
        
        
        # Set x-axis limits
        ax_flux.set_xlim(input_xrange)
        
        # Restore y limits if they existed, otherwise let matplotlib auto-scale
        if current_flux_ylim is not None:
            ax_flux.set_ylim(current_flux_ylim)
        if current_norm_ylim is not None and self.normalized_flux is not None:
            ax_norm.set_ylim(current_norm_ylim)
        
        # Update canvas
        self.canvas.draw()
    
    def plot_masks(self, ax):
        """Plot mask regions on the given axis."""
        ylim = ax.get_ylim()
        y_range = ylim[1] - ylim[0]
        y_min, y_max = ylim[0] - 0.05 * y_range, ylim[1] + 0.05 * y_range
        
        # Plot each mask pair
        for i in range(0, len(self.masks), 2):
            if i+1 < len(self.masks):
                vmin, vmax = self.masks[i], self.masks[i+1]
                ax.axvspan(vmin, vmax, alpha=0.2, color='gray')
                ax.axvline(x=vmin, color='red', linestyle=':')
                ax.axvline(x=vmax, color='red', linestyle=':')


    def on_click(self, event):
        """Handle mouse click events."""
        if event.inaxes is None or self.toolbar.mode != "":
            return  # Ignore clicks outside axes or when toolbar is active
    
        # Set focus to canvas to ensure key events work
        self.canvas.setFocus()
        
        # Handle spline point selection in spline mode
        if self.fit_params['method'] == 'spline' and event.inaxes == self.canvas.axes[0]:
            if event.button == 1:  # Left click
                self.add_median_spline_point(event.xdata, event.ydata)
                return
            elif event.button == 3:  # Right click
                self.remove_closest_spline_point(event.xdata, event.ydata)
                return
        
        # Handle mask operations (in both modes)
        if event.button == 1:  # Left click - always for adding masks
            self.handle_add_mask(event)
        elif event.button == 3:  # Right click - always for removing masks
            self.handle_remove_mask(event)
    
    def add_exact_spline_point(self, x, y):
        """Add a spline anchor point at the exact coordinates."""
        # Add the point to the list
        self.spline_points.append((x, y))
        
        # Sort points by x-coordinate for proper spline fitting
        self.spline_points.sort(key=lambda point: point[0])
        
        self.statusBar.showMessage(f"Exact spline point added at ({x:.1f}, {y:.2f}). Total: {len(self.spline_points)}")

        # Update plots
        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()
        self.update_plots()

    
    def add_median_spline_point(self, x, y):
        """Add a spline anchor point using the median flux in a window."""
        # Get window size from UI
        window_size = self.median_window_spin.value()
        window_type = self.fit_params['window_type']
        
        # Determine the window in data points
        if window_type == 'pixels':
            # Window size is in number of pixels
            half_window = int(window_size / 2)
            # Find the closest index to the click point
            closest_idx = np.argmin(np.abs(self.x_axis - x))
            # Define window indices with boundary checks
            start_idx = max(0, closest_idx - half_window)
            end_idx = min(len(self.x_axis) - 1, closest_idx + half_window)
        else:
            # Window size is in x-axis units
            # Define window bounds
            x_min, x_max = x - window_size/2, x + window_size/2
            # Find indices within this window
            window_indices = np.where((self.x_axis >= x_min) & (self.x_axis <= x_max))[0]
            if len(window_indices) == 0:
                # No points in window
                self.statusBar.showMessage(f"No data points found in window of size {window_size} {self.x_axis_unit}")
                return
            start_idx, end_idx = window_indices[0], window_indices[-1]
        
        # Calculate median flux in the window
        window_flux = self.flux[start_idx:end_idx+1]
        median_flux = np.median(window_flux)
        
        # Add the point at the clicked x-coordinate and the median flux
        self.spline_points.append((x, median_flux))
        
        # Sort points by x-coordinate for proper spline fitting
        self.spline_points.sort(key=lambda point: point[0])
        
        # Provide information about the window
        window_info = f"{window_size} {window_type}"
        if window_type != 'pixels':
            window_info += f", {end_idx - start_idx + 1} pixels"
        
        self.statusBar.showMessage(f"Median spline point added at ({x:.1f}, {median_flux:.2f}) using window: {window_info}. Total: {len(self.spline_points)}")
        
        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()
        # Update plots
        self.update_plots()

    
    def remove_closest_spline_point(self, x, y):
        """Remove the spline point closest to the given coordinates."""
        if not self.spline_points:
            self.statusBar.showMessage("No spline points to remove")
            return
        
        # Find the closest point
        distances = [(i, (point[0]-x)**2 + (point[1]-y)**2) 
                    for i, point in enumerate(self.spline_points)]
        closest_idx, _ = min(distances, key=lambda item: item[1])
        
        # Remove it
        removed_point = self.spline_points.pop(closest_idx)
        
        self.statusBar.showMessage(f"Removed spline point at ({removed_point[0]:.1f}, {removed_point[1]:.2f})")

        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()
        # Update plots
        self.update_plots()
    
    def clear_spline_points(self):
        """Clear all spline points."""
        if self.spline_points:
            self.spline_points = []
            self.statusBar.showMessage("All spline points cleared")

            ax_flux = self.canvas.axes[0]
            xlim = ax_flux.get_xlim()
            # Update plots
            self.update_plots()

        else:
            self.statusBar.showMessage("No spline points to clear")
        self.canvas.setFocus()

    
    def handle_add_mask(self, event):
        """Handle adding mask regions."""
        if self.current_click is None:
            # First click of a pair
            self.current_click = event.xdata
            self.statusBar.showMessage(f"First point selected at {self.x_axis_type}={event.xdata:.1f}. Click again to complete mask region.")
            
            # Mark the point
            self.canvas.axes[0].plot(event.xdata, event.ydata, 'ro', ms=5)
            self.canvas.draw()
        else:
            # Second click - create mask
            v_min = min(self.current_click, event.xdata)
            v_max = max(self.current_click, event.xdata)
            
            # Add the mask
            self.masks.extend([v_min, v_max])
            self.statusBar.showMessage(f"Mask added: {v_min:.1f} to {v_max:.1f} {self.x_axis_unit}")
    
            # Merge overlapping masks
            self.merge_masks()
            
            # Reset current click and update plot
            self.current_click = None

            ax_flux = self.canvas.axes[0]
            xlim = ax_flux.get_xlim()
            # Update plots
            self.update_plots()

    
    def handle_remove_mask(self, event):
        """Handle removing mask regions with right clicks."""
        if self.current_click is None:
            # First click of a pair
            self.current_click = event.xdata
            self.statusBar.showMessage(f"First removal point selected at {self.x_axis_type}={event.xdata:.1f}. Right-click again to complete removal region.")
            
            # Mark the point with a different color (blue) to distinguish from adding
            self.canvas.axes[0].plot(event.xdata, event.ydata, 'bo', ms=5)
            self.canvas.draw()
        else:
            # Second click - define removal range
            v_min = min(self.current_click, event.xdata)
            v_max = max(self.current_click, event.xdata)
            
            # Create a new list with only the masks we want to keep
            new_masks = []
            removed_count = 0
            
            # Check each mask pair
            for i in range(0, len(self.masks), 2):
                if i+1 < len(self.masks):
                    mask_min, mask_max = self.masks[i], self.masks[i+1]
                    
                    # Check for overlap - if no overlap, keep the mask
                    if not (mask_min <= v_max and mask_max >= v_min):
                        new_masks.extend([mask_min, mask_max])
                    else:
                        removed_count += 1
            
            # Replace the original masks list with our filtered list
            self.masks = new_masks
            
            if removed_count > 0:
                self.statusBar.showMessage(f"Removed {removed_count} mask regions in range [{v_min:.1f}, {v_max:.1f}] {self.x_axis_unit}")
            else:
                self.statusBar.showMessage(f"No masks found in range [{v_min:.1f}, {v_max:.1f}] {self.x_axis_unit}")
            
            # Reset current click and update plot
            self.current_click = None

            ax_flux = self.canvas.axes[0]
            xlim = ax_flux.get_xlim()
            # Update plots
            self.update_plots()

    
    def merge_masks(self):
        """Merge overlapping mask regions."""
        if len(self.masks) < 4:  # Need at least two pairs to merge
            return
        
        # Convert to list of pairs for easier manipulation
        mask_pairs = [(self.masks[i], self.masks[i+1]) for i in range(0, len(self.masks), 2)]
        
        # Sort by starting point
        mask_pairs.sort()
        
        # Merge overlapping intervals
        merged = []
        current_start, current_end = mask_pairs[0]
        
        for start, end in mask_pairs[1:]:
            if start <= current_end:  # Overlapping
                current_end = max(current_end, end)  # Extend current mask
            else:  # Not overlapping
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        
        # Add the last mask
        merged.append((current_start, current_end))
        
        # Convert back to flat list
        self.masks = [x for pair in merged for x in pair]
        
        if len(merged) < len(mask_pairs):
            self.statusBar.showMessage(f"Merged overlapping masks: {len(mask_pairs)} regions → {len(merged)} regions")

    def on_key_press(self, event):
        """Handle keyboard events."""
        if event.key == 'r':
            self.reset_masks()
        elif event.key == 'R':  # Capital R for reset all
            self.reset_all()
        elif event.key == 'z':
            self.undo_last_mask()
        elif event.key == 'c':
            self.clear_spline_points()
        elif event.key == 'f':
            self.fit_continuum()
        elif event.key == 'a':
            self.accept_results()
        elif event.key == 'escape':
            self.cancel()
        elif event.key == 'm':
            self.manual_mask_entry()
        elif event.key == 'p':
            self.poly_radio.setChecked(True)
        elif event.key == 's':
            self.spline_radio.setChecked(True)
        elif event.key == 'b' and self.fit_params['method'] == 'spline':
            # Add exact spline point at cursor position (if available)
            if hasattr(event, 'xdata') and hasattr(event, 'ydata') and event.inaxes:
                self.add_exact_spline_point(event.xdata, event.ydata)
            else:
                self.statusBar.showMessage("Position cursor over plot before pressing 'b'")
        elif event.key in ['+', '=']:  # Zoom in
            self.zoom_in()
        elif event.key in ['-', '_']:  # Zoom out
            self.zoom_out()
        elif event.key == '0':  # Reset zoom
            self.reset_zoom()

        elif event.key == 'h':
            self.show_help_dialog()
            
    def zoom_in(self):
        """Zoom in on the plot."""
        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()
        ylim = ax_flux.get_ylim()
        
        # Calculate new limits (zoom in by 25%)
        x_center = (xlim[0] + xlim[1]) / 2
        y_center = (ylim[0] + ylim[1]) / 2
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]
        new_xlim = (x_center - 0.375*x_range, x_center + 0.375*x_range)
        new_ylim = (y_center - 0.375*y_range, y_center + 0.375*y_range)
        
        # Apply new limits
        ax_flux.set_xlim(new_xlim)
        ax_flux.set_ylim(new_ylim)


        # Update the appropriate domain based on current x-axis type
        if self.x_axis_type == 'velocity':
            self.velocity_domain = new_xlim
        else:
            self.wave_domain = new_xlim
        
        self.domain = new_xlim
        self.canvas.draw()
        self.canvas.draw()
    
    def zoom_out(self):
        """Zoom out on the plot."""
        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()
        ylim = ax_flux.get_ylim()
        
        # Calculate new limits (zoom out by 25%)
        x_center = (xlim[0] + xlim[1]) / 2
        y_center = (ylim[0] + ylim[1]) / 2
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]
        new_xlim = (x_center - 0.625*x_range, x_center + 0.625*x_range)
        new_ylim = (y_center - 0.625*y_range, y_center + 0.625*y_range)
        
        # Apply new limits
        ax_flux.set_xlim(new_xlim)
        ax_flux.set_ylim(new_ylim)
        
        # Update the appropriate domain based on current x-axis type
        if self.x_axis_type == 'velocity':
            self.velocity_domain = new_xlim
        else:
            self.wave_domain = new_xlim
        
        self.domain = new_xlim

        self.canvas.draw()
    
    def reset_zoom(self):
        """Reset the plot zoom to the domain limits."""
        ax_flux = self.canvas.axes[0]

        # Get the appropriate domain limits based on x-axis type
        if self.x_axis_type == 'velocity':
            self.velocity_domain = self.domain if self.domain else [min(self.velocity), max(self.velocity)]
            ax_flux.set_xlim(self.velocity_domain)
        else:
            self.wave_domain = self.domain if self.domain else [min(self.wave), max(self.wave)]
            ax_flux.set_xlim(self.wave_domain)
        
        # Update the main domain property
        self.domain = ax_flux.get_xlim()

        # Reset y-axis scaling
        ax_flux.autoscale(axis='y')
        self.canvas.draw()
    
    def reset_masks(self):
        """Reset all masks."""
        if not self.masks:
            self.statusBar.showMessage("No masks to reset")
            return
            
        self.masks = []
        self.current_click = None
        self.statusBar.showMessage("All masks reset")
        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()

        # Update plots
        self.update_plots(input_xrange=xlim)

        # Ensure canvas has focus capabilities
        self.canvas.setFocus()


    
    def reset_all(self):
        """Reset everything: masks, spline points, and fitted continuum."""
        # Ask for confirmation
        if (self.masks or self.spline_points or self.continuum is not None):
            reply = QMessageBox.question(self, 'Reset All', 
                                        'Reset all masks, spline points, and fitted continuum?',
                                        QMessageBox.Yes | QMessageBox.No, 
                                        QMessageBox.No)
            if reply != QMessageBox.Yes:
                return
        
        # Reset all components
        self.masks = []
        self.spline_points = []
        self.continuum = None
        self.normalized_flux = None
        self.current_click = None
        self.fit_result = None
        
        # Disable accept button (no results to accept)
        self.accept_btn.setEnabled(False)
        
        self.statusBar.showMessage("Reset complete. All masks, spline points, and fitted continuum cleared.")
        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()

        # Update plots
        self.update_plots(input_xrange=xlim)
    
        # Ensure canvas has focus capabilities
        self.canvas.setFocus()



    
    def undo_last_mask(self):
        """Remove the last mask pair."""
        if len(self.masks) >= 2:
            self.masks.pop()
            self.masks.pop()
            self.statusBar.showMessage("Last mask removed")
            ax_flux = self.canvas.axes[0]
            xlim = ax_flux.get_xlim()

            # Update plots
            self.update_plots(input_xrange=xlim)

        else:
            self.statusBar.showMessage("No masks to remove")

        # Ensure canvas has focus capabilities
        self.canvas.setFocusPolicy(Qt.StrongFocus)


    def auto_mask(self):
        """
        Automatically generate masks using robust detection of spectral features,
        applicable in either wavelength or velocity space.
    
        Features:
        - Sigma clipping using MAD.
        - Error threshold masking.
        - Continuum fitting with fallback.
        - Minimum mask width and gap merging.
        - Applies in selected coordinate system: 'velocity' or 'wavelength'.
        """
        # Parameters
        sigma = self.auto_mask_sigma_spin.value()
        min_width = self.min_mask_width_spin.value()
        gap_threshold = self.mask_gap_spin.value() if hasattr(self, 'mask_gap_spin') else 10.0
        error_threshold = 0.5  # 50% relative error threshold
    
        coord_type = self.fit_params.get('window_type', 'velocity')  # 'velocity' or 'wavelength'
        x_data = self.velocity if coord_type == 'velocity' else self.x_axis
    
        if self.masks:
            reply = QMessageBox.question(
                self, 'Auto-Mask',
                'Replace existing masks?',
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            replace = reply == QMessageBox.Yes
        else:
            replace = True
    
        self.statusBar.showMessage("Generating masks automatically...")
    
        try:
            # Step 1: Build current unmasked index array
            current_mask = np.ones_like(x_data, dtype=bool)
            if not replace and self.masks:
                for i in range(0, len(self.masks), 2):
                    if i + 1 < len(self.masks):
                        x_min, x_max = self.masks[i], self.masks[i + 1]
                        current_mask[(x_data >= x_min) & (x_data <= x_max)] = False
    
            unmasked_x = x_data[current_mask]
            unmasked_flux = self.flux[current_mask]
            unmasked_error = self.error[current_mask]
    
            # Step 2: Continuum fit
            try:
                from rbcodes.IGM.rb_iter_contfit import rb_iter_contfit
                temp_order = min(2, max(0, len(unmasked_x) // 10 - 1))
                prelim_fit = rb_iter_contfit(
                    unmasked_x,
                    unmasked_flux,
                    error=unmasked_error,
                    order=temp_order,
                    sigma=3.0,
                    use_weights=False,
                    return_model=True
                )
                prelim_continuum = prelim_fit['model'](x_data)
            except Exception as e:
                self.statusBar.showMessage(f"Continuum fit failed: {str(e)}. Using median fallback.")
                prelim_continuum = np.full_like(x_data, np.nanmedian(unmasked_flux))
    
            # Step 3: Normalize flux, compute residuals
            normalized_flux = self.flux / prelim_continuum
            residual = normalized_flux - 1.0
    
            from astropy.stats import mad_std
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                robust_std = mad_std(residual, ignore_nan=True)

    
            # Step 4: Identify potential masked points
            potential_masks = (np.abs(residual) > sigma * robust_std) | \
                              (self.error / prelim_continuum > error_threshold)

            # Try wavelet method
            #from rbcodes.IGM.find_line_features_wavelet import find_spectral_lines_wavelet, create_feature_mask 
            #try:
            #    # Apply wavelet line finder to normalized residuals
            #    results = find_spectral_lines_wavelet(x_data, residual,
            #    line_type='absorption', 
            #    min_snr=sigma,
            #    min_scale=1,
            #    max_scale=5,
            #    num_scales=15,
            #    fit_lines=True,
            #    edge_buffer=10
            #)
            #   
            #    # Create a boolean mask (True = keep, False = mask out)
            #    wavelet_mask = create_feature_mask(x_data, results, width_scale_factor=5.0)
            #    #potential_masks = ~wavelet_mask  # invert because we want to mask the detected lines
            #    potential_masks = (( ~wavelet_mask) |
            #        (np.abs(residual) > sigma * robust_std) |(self.error / prelim_continuum > error_threshold)
            #        )
            #except Exception as e:
            #    self.statusBar.showMessage(f"Wavelet-based detection failed: {str(e)}. Falling back to MAD.")
            #    potential_masks = (np.abs(residual) > sigma * robust_std) | \
            #                    (self.error / prelim_continuum > error_threshold)
    
    
            # Step 5: Find contiguous mask regions
            mask_regions = []
            in_region = False
            start_idx = None
            for i, flag in enumerate(potential_masks):
                if flag and not in_region:
                    start_idx = i
                    in_region = True
                elif not flag and in_region:
                    end_idx = i - 1
                    in_region = False
                    mask_regions.append((start_idx, end_idx))
            if in_region:
                mask_regions.append((start_idx, len(potential_masks) - 1))
    
            # Step 6: Filter by min width in appropriate units
            filtered_regions = []
            for start_idx, end_idx in mask_regions:
                x_start = x_data[start_idx]
                x_end = x_data[end_idx]
                if x_end - x_start >= min_width:
                    filtered_regions.append((start_idx, end_idx))
    
            # Step 7: Merge regions closer than gap_threshold
            if len(filtered_regions) > 1:
                merged = [filtered_regions[0]]
                for s, e in filtered_regions[1:]:
                    last_s, last_e = merged[-1]
                    if x_data[s] - x_data[last_e] <= gap_threshold:
                        merged[-1] = (last_s, e)
                    else:
                        merged.append((s, e))
                filtered_regions = merged
    
            # Step 8: Finalize new masks
            new_masks = []
            for start_idx, end_idx in filtered_regions:
                new_masks.extend([x_data[start_idx], x_data[end_idx]])
    
            if replace:
                self.masks = new_masks
            else:
                self.masks.extend(new_masks)
                self.merge_masks()
    
            region_count = len(new_masks) // 2
            units = "km/s" if coord_type == "velocity" else "Å"
            self.statusBar.showMessage(
                f"Added {region_count} mask regions (σ={sigma}, MAD={robust_std:.3f}, "
                f"min_width={min_width} {units}, error>{error_threshold*100:.0f}%)"
            )
            xlim = self.canvas.figure.gca().get_xlim()

            self.update_plots(input_xrange=xlim)
            # Ensure canvas has focus capabilities
            self.canvas.setFocusPolicy(Qt.StrongFocus)
    
        except Exception as e:
            self.statusBar.showMessage(f"Error in auto-masking: {str(e)}")
            QMessageBox.warning(
                self, "Auto-Masking Error",
                f"Error generating masks:\n{str(e)}\n\n"
                "Try adjusting parameters or masking manually."
            )


    def manual_mask_entry(self):
        """Open a dialog for manual entry of mask values."""
        # Create dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Manual Mask Entry")
        layout = QVBoxLayout(dialog)
        
        # Create form for min/max values
        form_layout = QFormLayout()
        
        min_spin = QDoubleSpinBox()
        min_spin.setRange(-1e6, 1e6)  # Large range to accommodate both velocity and wavelength
        min_spin.setValue(0)
        
        max_spin = QDoubleSpinBox()
        max_spin.setRange(-1e6, 1e6)
        max_spin.setValue(0)
        
        # Add unit label
        unit_label = QLabel(f"Values in {self.x_axis_unit}")
        
        form_layout.addRow("Minimum value:", min_spin)
        form_layout.addRow("Maximum value:", max_spin)
        form_layout.addRow("", unit_label)
        
        layout.addLayout(form_layout)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            v_min = min(min_spin.value(), max_spin.value())
            v_max = max(min_spin.value(), max_spin.value())
            
            if v_min != v_max:
                self.masks.extend([v_min, v_max])
                self.statusBar.showMessage(f"Mask added: {v_min:.1f} to {v_max:.1f} {self.x_axis_unit}")
                
                ax_flux = self.canvas.axes[0]
                xlim = ax_flux.get_xlim()

                # Update plots
                self.update_plots(input_xrange=xlim)

            else:
                QMessageBox.warning(self, "Invalid Input", "Minimum and maximum values must be different")
        # Ensure canvas has focus capabilities
        self.canvas.setFocusPolicy(Qt.StrongFocus)

                
    def fit_continuum(self):
        """Fit continuum using current parameters and masks."""
        # Update parameters from UI
        self.fit_params['method'] = 'polynomial' if self.poly_radio.isChecked() else 'spline'
        self.fit_params['order'] = self.order_spin.value()
        self.fit_params['use_weights'] = self.use_weights_check.isChecked()
        self.fit_params['optimize_cont'] = self.optimize_check.isChecked()
        self.fit_params['min_order'] = self.min_order_spin.value()
        self.fit_params['max_order'] = self.max_order_spin.value()
        self.fit_params['sigma'] = self.sigma_spin.value()
        self.fit_params['spline_degree'] = self.spline_degree_spin.value()
        self.fit_params['spline_smooth'] = self.spline_smooth_spin.value()
        self.fit_params['median_window'] = self.median_window_spin.value()
        self.fit_params['window_type'] = self.window_type_combo.currentText().lower()
        
        # Check if we have enough data points
        if len(self.x_axis) < 3:
            QMessageBox.warning(self, "Insufficient Data", 
                              "Not enough data points for fitting.")
            return
            
        try:
            # Use different fitting approach based on method
            if self.fit_params['method'] == 'polynomial':
                self.fit_polynomial_continuum()
            else:
                self.fit_spline_continuum()
                
            # Calculate normalized flux
            self.normalized_flux = self.flux / self.continuum
            
            # Enable accept button
            self.accept_btn.setEnabled(True)
            
            ax_flux = self.canvas.axes[0]
            xlim = ax_flux.get_xlim()
            # Update plots
            self.update_plots(input_xrange=xlim)
            
        except Exception as e:
            self.statusBar.showMessage(f"Error fitting continuum: {str(e)}")
            QMessageBox.warning(self, "Fitting Error", 
                              f"Error fitting continuum: {str(e)}\n\n"
                              f"Try adjusting parameters or mask regions.")
        
        # Ensure canvas has focus capabilities
        self.canvas.setFocus()


    def fit_polynomial_continuum(self):
        """Fit continuum using polynomial method with masked regions excluded."""
        # Create mask array (True for unmasked, False for masked)
        mask_array = np.ones_like(self.x_axis, dtype=bool)
        
        for i in range(0, len(self.masks), 2):
            if i+1 < len(self.masks):
                v_min, v_max = self.masks[i], self.masks[i+1]
                mask_indices = (self.x_axis >= v_min) & (self.x_axis <= v_max)
                mask_array[mask_indices] = False
        
        # Extract unmasked points
        unmasked_x = self.x_axis[mask_array]
        unmasked_flux = self.flux[mask_array]
        unmasked_error = self.error[mask_array] if self.fit_params['use_weights'] else None
        
        if len(unmasked_x) < 3:
            raise ValueError("Not enough unmasked points for fitting")
        
        # If auto-optimizing, use BIC to find the best polynomial order
        if self.fit_params['optimize_cont']:
            # Check if we have enough points for highest order
            min_required = self.fit_params['max_order'] + 2
            
            if len(unmasked_x) < min_required:
                old_max = self.fit_params['max_order']
                self.fit_params['max_order'] = len(unmasked_x) - 2
                self.statusBar.showMessage(f"Reduced max order from {old_max} to {self.fit_params['max_order']} due to limited data points")
            
            try:
                result = fit_optimal_polynomial(
                    unmasked_x,
                    unmasked_flux,
                    error=unmasked_error,
                    min_order=self.fit_params['min_order'],
                    max_order=self.fit_params['max_order'],
                    use_weights=self.fit_params['use_weights'],
                    sigma=self.fit_params['sigma'],
                    plot=False,
                    include_model=True,
                    silent=True
                )
                
                self.fit_result = result
                best_order = result['best_order']
                self.continuum = result['model'](self.x_axis)
                
                self.statusBar.showMessage(f"Best polynomial order: {best_order} (BIC: {min([b for _, b in result['bic_results']]):.2f})")
                
            except Exception as e:
                raise ValueError(f"Optimization failed: {str(e)}. Try fixed order instead.")
        
        else:
            # Use specified order
            try:
                result = rb_iter_contfit(
                    unmasked_x,
                    unmasked_flux,
                    error=unmasked_error,
                    order=self.fit_params['order'],
                    sigma=self.fit_params['sigma'],
                    use_weights=self.fit_params['use_weights'],
                    return_model=True,
                    silent=True
                )
                
                self.fit_result = result
                self.continuum = result['model'](self.x_axis)
                
                self.statusBar.showMessage(f"Polynomial fit with order {self.fit_params['order']} (fit error: {result['fit_error']:.5f})")
                
            except Exception as e:
                raise ValueError(f"Fitting failed: {str(e)}")
        # Ensure canvas has focus capabilities
        self.canvas.setFocusPolicy(Qt.StrongFocus)

    
    def fit_spline_continuum(self):
        """Fit continuum using spline method with provided anchor points."""
        # Check if we have enough spline points
        if len(self.spline_points) < 3:
            raise ValueError(f"At least 3 spline points are needed, got {len(self.spline_points)}")
        
        try:
            # Extract x and y coordinates
            x_points, y_points = zip(*self.spline_points)
            x_points = np.array(x_points)
            y_points = np.array(y_points)
            
            # Sort by x-coordinate (just to be safe)
            sort_indices = np.argsort(x_points)
            x_points = x_points[sort_indices]
            y_points = y_points[sort_indices]
            
            # Fit spline
            spline_degree = self.fit_params['spline_degree']
            spline_smooth = self.fit_params['spline_smooth']
            
            # Add weight parameter if needed (uniform weight for simplicity)
            weights = np.ones_like(x_points)
            
            # Remember that splrep returns a tuple of (knots, coefficients, degree)
            spline_tck = splrep(x_points, y_points, w=weights, k=spline_degree, s=spline_smooth)
            
            # Evaluate the spline at velocity points
            self.continuum = splev(self.x_axis, spline_tck)
            
            # Store basic spline info in fit_result
            self.fit_result = {
                'method': 'spline',
                'spline_degree': spline_degree,
                'spline_smooth': spline_smooth,
                'spline_points': self.spline_points,
                'spline_tck': spline_tck
            }
            
            self.statusBar.showMessage(f"Spline fit with {len(self.spline_points)} points, degree {spline_degree}, smoothing {spline_smooth}")
            
        except Exception as e:
            raise ValueError(f"Spline fitting failed: {str(e)}")
        # Ensure canvas has focus capabilities
        self.canvas.setFocusPolicy(Qt.StrongFocus)

    
    def check_unmasked_points(self):
        """Check if there are enough unmasked points for fitting."""
        # Count unmasked points
        mask_array = np.ones_like(self.x_axis, dtype=bool)
        
        for i in range(0, len(self.masks), 2):
            if i+1 < len(self.masks):
                v_min, v_max = self.masks[i], self.masks[i+1]
                mask_indices = (self.x_axis >= v_min) & (self.x_axis <= v_max)
                mask_array[mask_indices] = False
        
        unmasked_count = np.sum(mask_array)
        
        # For optimal polynomial, need enough points for highest order + 1
        if self.fit_params['method'] == 'polynomial':
            if self.fit_params['optimize_cont']:
                min_required = self.fit_params['max_order'] + 2
            else:
                min_required = self.fit_params['order'] + 2
        else:
            # For spline, we just need 3 points minimum
            min_required = 3
        
        if unmasked_count < min_required:
            QMessageBox.warning(self, "Insufficient Data", 
                              f"Not enough unmasked points ({unmasked_count}) "
                              f"for fitting.\n\n"
                              f"At least {min_required} points are required.")
            self.statusBar.showMessage(f"Insufficient unmasked points: {unmasked_count}/{min_required}")
            return False
        
        return True
    
    def accept_results(self):
        """Accept results and close window."""
        if self.continuum is None:
            reply = QMessageBox.question(self, 'No Fit', 
                                         'No continuum has been fitted. Fit now?',
                                         QMessageBox.Yes | QMessageBox.No, 
                                         QMessageBox.Yes)
            
            if reply == QMessageBox.Yes:
                self.fit_continuum()
                if self.continuum is None:
                    return  # Fitting failed or was canceled
            else:
                return  # User canceled
        
        # Convert mask regions to wavelength if we're in velocity mode
        mask_wavelengths = []
        if self.velocity is not None and self.x_axis_type == 'velocity':
            for i in range(0, len(self.masks), 2):
                if i+1 < len(self.masks):
                    v_min, v_max = self.masks[i], self.masks[i+1]
                    # Find corresponding wavelength ranges
                    wmin = self.wave[np.abs(self.velocity - v_min).argmin()]
                    wmax = self.wave[np.abs(self.velocity - v_max).argmin()]
                    mask_wavelengths.extend([wmin, wmax])
        else:
            # We're already in wavelength mode, just copy the masks
            mask_wavelengths = self.masks.copy()
        
        # Prepare result dictionary
        self.result = {
            'masks': self.masks,
            'mask_wavelengths': mask_wavelengths,
            'continuum': self.continuum,
            'fit_params': self.fit_params,
            'normalized_flux': self.normalized_flux,
            'fit_method': self.fit_params['method'],
            'x_axis_type': self.x_axis_type,
            'timestamp': datetime.datetime.now().isoformat(),
            'cancelled': False
        }
        
        # Add method-specific results
        if self.fit_params['method'] == 'polynomial':
            if self.fit_result is not None:
                if 'fit_error' in self.fit_result:
                    self.result['fit_error'] = self.fit_result['fit_error']
                if 'best_order' in self.fit_result:
                    self.result['best_order'] = self.fit_result['best_order']
        elif self.fit_params['method'] == 'spline':
            self.result['spline_points'] = self.spline_points
            
        self.close()
    
    def cancel(self):
        """Cancel and close without saving."""
        self.result = {'cancelled': True}
        self.close()
    
    def closeEvent(self, event):
        """Handle window close event."""
        if not hasattr(self, 'result'):
            self.result = {'cancelled': True}
        event.accept()

def launch_interactive_continuum_fit(**kwargs):
    """
    Launch the interactive continuum fitting GUI.
    
    Parameters
    ----------
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    velocity : array, optional
        Velocity array. If not provided, wavelength will be used as x-axis.
    existing_masks : list, optional
        Existing mask regions
    order : int, optional
        Initial polynomial order
    use_weights : bool, optional
        Whether to use weights
    domain : list, optional
        Velocity or wavelength domain limits
        
    Returns
    -------
    dict
        Results including masks, continuum, and fit parameters,
        or {'cancelled': True} if the user cancelled.
    """
    # Create application if it doesn't exist
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    # Create and show window
    window = InteractiveContinuumFitWindow(**kwargs)
    window.show()
    
    # Run the application
    app.exec_()
    
    # Return the result
    return window.result

def launch_interactive_continuum_fit_dialog(**kwargs):
    """
    Launch the interactive continuum fitter as a modal dialog.
    
    This version uses QDialog.exec_() which works properly when called from
    an application that already has an event loop running.
    
    Parameters
    ----------
    Same parameters as launch_interactive_continuum_fit
    
    Returns
    -------
    dict
        Results including masks, continuum, and fit parameters,
        or {'cancelled': True} if the user cancelled.
    """
    from PyQt5.QtWidgets import QDialog, QVBoxLayout
    from PyQt5.QtCore import Qt
    
    # Create dialog to host the window
    dialog = QDialog()
    dialog.setWindowTitle("Interactive Continuum Fitting")
    dialog.resize(1000, 700)
    dialog.setWindowFlags(dialog.windowFlags() | Qt.WindowMaximizeButtonHint)
    
    # Create the continuum fit window as a widget inside the dialog
    window = InteractiveContinuumFitWindow(**kwargs)
    
    # Connect the window's close event to accept the dialog
    window.closeEvent = lambda event: dialog.accept()
    
    # Set up layout
    layout = QVBoxLayout(dialog)
    layout.addWidget(window)
    layout.setContentsMargins(0, 0, 0, 0)
    
    # Execute the dialog
    dialog.exec_()
    
    # Return the result from the window
    return getattr(window, 'result', {'cancelled': True})

#-- Help Window
class HelpDialog(QDialog):
    """Modal dialog to display HTML-formatted help content."""
    def __init__(self, parent=None):
        super(HelpDialog, self).__init__(parent)
        self.setWindowTitle("Interactive Continuum Fitting Help")
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowContextHelpButtonHint)
        self.resize(900, 700)  # Set a reasonable default size
        
        # Create layout
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # Create text browser for HTML content
        self.text_browser = QTextBrowser()
        self.text_browser.setOpenExternalLinks(True)  # Allow clicking on links
        layout.addWidget(self.text_browser)
        
        # Add close button at bottom
        button_box = QDialogButtonBox(QDialogButtonBox.Close)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        # Load the help content
        self.load_help_content()
    
    def load_help_content(self):
        """Load the HTML help content from file."""
        try:
            # Try to find the help file
            module_dir = os.path.dirname(os.path.abspath(__file__))
            help_file = os.path.join(module_dir, 'interactive_continuum_fit_help', 'help.html')
            
            # Check if file exists
            if os.path.exists(help_file):
                with open(help_file, 'r', encoding='utf-8') as f:
                    html_content = f.read()
                self.text_browser.setHtml(html_content)
                
                # Set base URL for relative resources (like images)
                base_url = QUrl.fromLocalFile(os.path.dirname(help_file) + '/')
                self.text_browser.setSearchPaths([os.path.dirname(help_file)])
                self.text_browser.setSource(QUrl())  # Reset source
                self.text_browser.document().setBaseUrl(base_url)
            else:
                # Fallback to basic help content
                self.load_fallback_help()
                
        except Exception as e:
            # If any error occurs, show fallback help
            print(f"Error loading help content: {str(e)}")
            self.load_fallback_help()
    
    def load_fallback_help(self):
        """Load a basic fallback help content if the file can't be found."""
        fallback_html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Interactive Continuum Fitting Help</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                h1 { color: #2c3e50; }
                h2 { color: #3498db; }
                .section { margin-bottom: 20px; }
                .key { display: inline-block; background: #f1f1f1; border: 1px solid #ddd; 
                       border-radius: 3px; padding: 2px 5px; font-family: monospace; }
            </style>
        </head>
        <body>
            <h1>Interactive Continuum Fitting Help</h1>
            
            <div class="section">
                <h2>Overview</h2>
                <p>This tool allows you to interactively fit a continuum to spectral data using either 
                polynomial fitting or spline interpolation methods.</p>
            </div>
            
            <div class="section">
                <h2>Basic Controls</h2>
                <p><b>Left-click</b>: Add mask region (first click = start, second click = end)</p>
                <p><b>Right-click</b>: Remove mask or spline point</p>
                <p><b>Key <span class="key">b</span></b>: Add exact spline point at cursor position</p>
                <p><b>Key <span class="key">f</span></b>: Fit continuum</p>
                <p><b>Key <span class="key">a</span></b>: Accept results</p>
                <p><b>Key <span class="key">r</span></b>: Reset masks</p>
                <p><b>Key <span class="key">R</span></b>: Reset all (masks, spline points, fit)</p>
            </div>
            
            <p><i>Note: This is a fallback help page. The full documentation should be in the 
            interactive_continuum_fit_help folder.</i></p>
        </body>
        </html>
        """
        self.text_browser.setHtml(fallback_html)



if __name__ == "__main__":
    # Example usage with synthetic data
    import numpy as np
    
    # Generate sample data
    n_points = 1000
    wave = np.linspace(4000, 4500, n_points)  # Wavelength in Angstroms
    
    # Optionally create velocity array (but can work without it)
    # Using arbitrary redshift for example
    z = 0.0
    rest_wave = np.mean(wave)  # Arbitrary rest wavelength
    velocity = (wave - rest_wave * (1.0 + z)) * 299792.458 / (rest_wave * (1.0 + z))
    
    # Create a realistic continuum with some complexity
    continuum = 1.0 + 0.0002 * (wave - 4000)
    
    # Add absorption features
    absorption = (1.0 - 0.5 * np.exp(-((wave-4100)/5)**2) - 
                 0.8 * np.exp(-((wave-4200)/10)**2) -
                 0.3 * np.exp(-((wave-4350)/8)**2))
    
    # Create flux with some noise
    flux = continuum * absorption + 0.05 * np.random.randn(n_points)
    
    # Error array
    error = 0.05 * np.ones_like(flux)
    
    # Launch the tool with some pre-existing masks
    # Note that we don't pass velocity for this example to show wavelength-only operation
    #result = launch_interactive_continuum_fit(
    #    wave=wave,
    #    flux=flux,
    #    error=error,
    #    velocity=None,  # Test wavelength-only mode
    #    existing_masks=[4090, 4110, 4190, 4210],
    #    order=3,
    #    use_weights=False,
    #    domain=[4050, 4450]
    #)
    
    # Test also with velocity
    # Uncomment to run with velocity x-axis
    result = launch_interactive_continuum_fit(
        wave=wave,
        flux=flux,
        error=error,
        velocity=velocity,  # Include velocity
        existing_masks=[-300, -100, 100, 200],
        order=3,
        use_weights=False,
        domain=[-15000, 4000]
    )
    
    # Print the result
    if not result.get('cancelled', False):
        print("Results:")
        print(f"X-axis type: {result['x_axis_type']}")
        print(f"Fit method: {result['fit_method']}")
        print(f"Masks: {result['masks']}")
        if result['fit_method'] == 'polynomial':
            if 'best_order' in result:
                print(f"Best Order: {result['best_order']}")
            else:
                print(f"Polynomial Order: {result['fit_params']['order']}")
        else:
            print(f"Spline degree: {result['fit_params']['spline_degree']}")
            print(f"Number of spline points: {len(result['spline_points'])}")
    else:
        print("Fitting was cancelled.")
