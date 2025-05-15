"""
Interactive masking and continuum fitting tool for spectral analysis.

This module provides a PyQt5-based GUI for interactive selection of mask regions
and continuum fitting for spectroscopic data.
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
                           QTextBrowser)
from PyQt5.QtCore import Qt

from astropy.stats import sigma_clip
import time
import warnings
import datetime

# Try to import from rbcodes, fall back to local imports if not available
try:
    from rbcodes.IGM.rb_iter_contfit import fit_optimal_polynomial
except ImportError:
    try:
        from IGM.rb_iter_contfit import fit_optimal_polynomial
    except ImportError:
        print("Warning: rb_iter_contfit not found. Some functionality may be limited.")

class MplCanvas(FigureCanvasQTAgg):
    """Matplotlib canvas for embedding in PyQt5."""
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        super(MplCanvas, self).__init__(self.fig)

class InteractiveMaskWindow(QMainWindow):
    """Main window for interactive masking and continuum fitting."""
    
    def __init__(self, wave, flux, error, velocity, existing_masks=None, 
                 order=3, use_weights=False, domain=None):
        super(InteractiveMaskWindow, self).__init__()
        
        # Store data
        self.wave = np.array(wave)
        self.flux = np.array(flux)
        self.error = np.array(error)
        self.velocity = np.array(velocity)
        self.masks = list(existing_masks) if existing_masks is not None else []
        self.domain = domain if domain is not None else [min(velocity), max(velocity)]
        
        # Working variables
        self.current_click = None  # Store first click of a pair
        self.continuum = None
        self.normalized_flux = None
        self.fit_params = {
            'order': order,
            'use_weights': use_weights,
            'optimize_cont': True,
            'min_order': 1,
            'max_order': 6,
            'sigma': 3.0,
            'min_mask_width': 20.0,  # Default 20 km/s
            'mask_gap': 100.0  # Default 100 km/s
        }
        self.fit_result = None
        #self.mode = 'add'  # 'add' or 'remove' mask mode
        
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
        self.setWindowTitle("Interactive Masking and Continuum Fitting")
        self.setGeometry(100, 100, 1000, 700)
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
        # Create plotting area
        self.canvas = MplCanvas(self, width=8, height=6)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        
        # Create plot layout
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        
        # Add plot area to main layout
        main_layout.addLayout(plot_layout, 3)  # Plot takes 3/4 of width
        
        # Create sidebar with tabs
        sidebar = QTabWidget()
        
        # Tab 1: Basic controls
        basic_tab = QWidget()
        basic_layout = QVBoxLayout(basic_tab)
        
        # Mask controls group
        mask_group = QGroupBox("Mask Controls")
        mask_layout = QVBoxLayout()
        
        # Add instructions label
        mask_instructions = QLabel("Left-click pairs to add masks\nRight-click pairs to remove masks")
        mask_layout.addWidget(mask_instructions)
        
        # Mask buttons
        self.reset_btn = QPushButton("Reset Masks")
        self.reset_btn.clicked.connect(self.reset_masks)
        
        self.undo_btn = QPushButton("Undo Last")
        self.undo_btn.clicked.connect(self.undo_last_mask)
        
        self.auto_btn = QPushButton("Auto-Mask")
        self.auto_btn.clicked.connect(self.auto_mask)
        
        self.manual_btn = QPushButton("Manual Entry")
        self.manual_btn.clicked.connect(self.manual_mask_entry)
        
        mask_layout.addWidget(self.reset_btn)
        mask_layout.addWidget(self.undo_btn)
        mask_layout.addWidget(self.auto_btn)
        mask_layout.addWidget(self.manual_btn)
        mask_group.setLayout(mask_layout)
        basic_layout.addWidget(mask_group)
        
        # Fitting controls group
        fit_group = QGroupBox("Fitting Controls")
        fit_layout = QFormLayout()
        
        self.order_spin = QSpinBox()
        self.order_spin.setRange(0, 10)
        self.order_spin.setValue(self.fit_params['order'])
        
        self.use_weights_check = QCheckBox()
        self.use_weights_check.setChecked(self.fit_params['use_weights'])
        
        self.optimize_check = QCheckBox()
        self.optimize_check.setChecked(True)  # Changed to True as default
        self.optimize_check.toggled.connect(self.toggle_optimization)
        
        self.fit_btn = QPushButton("Fit Continuum")
        self.fit_btn.clicked.connect(self.fit_continuum)
        
        fit_layout.addRow("Polynomial Order:", self.order_spin)
        fit_layout.addRow("Use Weights:", self.use_weights_check)
        fit_layout.addRow("Auto Optimize:", self.optimize_check)
        fit_layout.addWidget(self.fit_btn)
        fit_group.setLayout(fit_layout)
        basic_layout.addWidget(fit_group)
        
        # Action buttons
        self.accept_btn = QPushButton("Accept & Return")
        self.accept_btn.clicked.connect(self.accept_results)
        self.accept_btn.setEnabled(False)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.cancel)
        
        basic_layout.addWidget(self.accept_btn)
        basic_layout.addWidget(self.cancel_btn)
        basic_layout.addStretch()
        
        # Tab 2: Advanced options
        advanced_tab = QWidget()
        advanced_layout = QFormLayout(advanced_tab)
        
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
        self.min_mask_width_spin.setValue(2.0)  # Default 20 km/s minimum mask width
        self.min_mask_width_spin.setSuffix(" km/s")
        
        self.mask_gap_spin = QDoubleSpinBox()
        self.mask_gap_spin.setRange(1.0, 500.0)
        self.mask_gap_spin.setSingleStep(10.0)
        self.mask_gap_spin.setValue(100.0)  # Default 100 km/s grouping window
        self.mask_gap_spin.setSuffix(" km/s")
        
        advanced_layout.addRow("Min Polynomial Order:", self.min_order_spin)
        advanced_layout.addRow("Max Polynomial Order:", self.max_order_spin)
        advanced_layout.addRow("Fitting Sigma:", self.sigma_spin)
        advanced_layout.addRow("Auto-Mask Sigma:", self.auto_mask_sigma_spin)
        advanced_layout.addRow("Min Mask Width:", self.min_mask_width_spin)
        advanced_layout.addRow("Grouping Window:", self.mask_gap_spin)
        
        # Tab 3: Help tab with comprehensive information
        help_tab = QWidget()
        help_layout = QVBoxLayout(help_tab)
        
        # Create a text browser for help content
        help_text = QTextBrowser()
        help_text.setOpenExternalLinks(True)
        help_text.setHtml("""
            <h2>Interactive Masking and Continuum Fitting</h2>
            
            <h3>Overview</h3>
            <p>This tool allows you to interactively select regions of a spectrum to mask when fitting the continuum, 
            which is useful for excluding absorption or emission lines.</p>
            
            <h3>Basic Usage</h3>
            <ol>
                <li>Use <b>left-click pairs</b> to define mask regions (first click = start, second click = end)</li>
                <li>Use <b>right-click pairs</b> to remove mask regions within the selected range</li>
                <li>Click <b>Fit Continuum</b> to fit a continuum to the unmasked regions</li>
                <li>Click <b>Accept & Return</b> to save the continuum and masks</li>
            </ol>
            
            <h3>Advanced Features</h3>
            <h4>Auto-Masking</h4>
            <p>The Auto-Mask button automatically identifies potential absorption or emission features
            based on statistical deviations from an estimated continuum. Parameters include:</p>
            <ul>
                <li><b>Auto-Mask Sigma</b>: Detection threshold in terms of standard deviations</li>
                <li><b>Min Mask Width</b>: Minimum width of features to include (in km/s)</li>
                <li><b>Grouping Window</b>: Maximum gap between features to merge (in km/s)</li>
            </ul>
            
            <h4>Continuum Fitting Options</h4>
            <ul>
                <li><b>Auto Optimize</b>: When checked, automatically selects the best polynomial order</li>
                <li><b>Min/Max Order</b>: Range of polynomial orders to consider when auto-optimizing</li>
                <li><b>Fitting Sigma</b>: Threshold for outlier rejection during fitting</li>
                <li><b>Use Weights</b>: When checked, uses error values to weight the fit</li>
            </ul>
            
            <h3>Keyboard Shortcuts</h3>
            <table border="1" cellpadding="4">
                <tr><th>Key</th><th>Action</th></tr>
                <tr><td>r</td><td>Reset all masks</td></tr>
                <tr><td>z</td><td>Undo last mask</td></tr>
                <tr><td>f</td><td>Fit continuum</td></tr>
                <tr><td>a</td><td>Accept results</td></tr>
                <tr><td>c/Esc</td><td>Cancel</td></tr>
                <tr><td>m</td><td>Manual mask entry</td></tr>
            </table>
            
            <h3>Zoom Controls</h3>
            <p>In addition to the toolbar zoom tools, you can use these keyboard shortcuts:</p>
            <table border="1" cellpadding="4">
                <tr><th>Key</th><th>Action</th></tr>
                <tr><td>+</td><td>Zoom in</td></tr>
                <tr><td>-</td><td>Zoom out</td></tr>
                <tr><td>0</td><td>Reset view</td></tr>
            </table>
            
            <h3>Masking Tips</h3>
            <ul>
                <li><b>Absorption Lines</b>: Mask regions where flux dips below the continuum</li>
                <li><b>Emission Lines</b>: Mask regions where flux rises above the continuum</li>
                <li><b>Noisy Regions</b>: Mask regions with large error bars or data quality issues</li>
                <li><b>Overlapping Masks</b>: Will be automatically merged into a single mask region</li>
            </ul>
            
            <h3>Continuum Fitting Tips</h3>
            <ul>
                <li><b>Polynomial Order</b>: Use lower orders for smoother continua, higher for more complex shapes</li>
                <li><b>Weighting</b>: Enable weighting when error values are reliable and vary significantly</li>
                <li><b>Auto-Optimization</b>: Uses Bayesian Information Criterion to select the optimal polynomial order</li>
            </ul>
            """)
        
        help_layout.addWidget(help_text)
        
        # Add tabs to sidebar
        sidebar.addTab(basic_tab, "Basic")
        sidebar.addTab(advanced_tab, "Advanced")
        sidebar.addTab(help_tab, "Help")
        
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
        
        # Initial UI state
        self.toggle_optimization(True)  # Initial state matches the default
        
    def toggle_mode(self, checked):
        """Toggle between add and remove mask modes."""
        sender = self.sender()
        if sender == self.add_mode_btn and checked:
            self.mode = 'add'
            self.remove_mode_btn.setChecked(False)
            self.statusBar.showMessage("Add mask mode: Left-click pairs to define mask regions")
        elif sender == self.remove_mode_btn and checked:
            self.mode = 'remove'
            self.add_mode_btn.setChecked(False)
            self.statusBar.showMessage("Remove mask mode: Click near a mask region to remove it")
    
    def toggle_optimization(self, checked):
        """Toggle optimization options."""
        self.order_spin.setEnabled(not checked)
        self.min_order_spin.setEnabled(checked)
        self.max_order_spin.setEnabled(checked)
    
    def update_plots(self):
        """Update the plots with current data and masks."""
        # Clear axes
        for ax in self.canvas.axes:
            ax.clear()
        
        # Get axes
        ax_flux, ax_norm = self.canvas.axes
        
        # Plot original spectrum and error as separate step plots
        ax_flux.step(self.velocity, self.flux, 'k-', where='mid', lw=0.5, label='Flux')
        ax_flux.step(self.velocity, self.error, 'r-', where='mid', lw=0.5, label='Error')
        
        # Plot continuum if available
        if self.continuum is not None:
            ax_flux.plot(self.velocity, self.continuum, 'r-', linewidth=2, label='Continuum')
        
        # Mark masked regions
        self.plot_masks(ax_flux)
        
        # Plot normalized spectrum if available
        if self.normalized_flux is not None:
            ax_norm.step(self.velocity, self.normalized_flux, 'k-', where='mid', lw=0.5)
            ax_norm.axhline(y=1.0, color='r', linestyle='--')
            self.plot_masks(ax_norm)
        else:
            ax_norm.text(0.5, 0.5, "Normalized spectrum will appear here after fitting", 
                       ha='center', va='center', transform=ax_norm.transAxes)
        
        # Set labels and limits
        ax_flux.set_ylabel('Flux')
        ax_flux.set_title('Original Spectrum with Masks')
        ax_flux.legend()
        
        ax_norm.set_xlabel('Velocity (km/s)')
        ax_norm.set_ylabel('Normalized Flux')
        
        # Set velocity limits
        ax_flux.set_xlim(self.domain)
        
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
        
        if event.button == 1:  # Left click - always for adding masks
            self.handle_add_mask(event)
        elif event.button == 3:  # Right click - always for removing masks
            self.handle_remove_mask(event)
    
    def handle_add_mask(self, event):
        """Handle adding mask regions."""
        if self.current_click is None:
            # First click of a pair
            self.current_click = event.xdata
            self.statusBar.showMessage(f"First point selected at v={event.xdata:.1f}. Click again to complete mask region.")
            
            # Mark the point
            self.canvas.axes[0].plot(event.xdata, event.ydata, 'ro', ms=5)
            self.canvas.draw()
        else:
            # Second click - create mask
            v_min = min(self.current_click, event.xdata)
            v_max = max(self.current_click, event.xdata)
            
            # Add the mask
            self.masks.extend([v_min, v_max])
            self.statusBar.showMessage(f"Mask added: {v_min:.1f} to {v_max:.1f}")
    
            # Merge overlapping masks
            self.merge_masks()
            
            # Reset current click and update plot
            self.current_click = None
            self.update_plots()
    
    def handle_remove_mask(self, event):
        """Handle removing mask regions with right clicks."""
        if self.current_click is None:
            # First click of a pair
            self.current_click = event.xdata
            self.statusBar.showMessage(f"First removal point selected at v={event.xdata:.1f}. Right-click again to complete removal region.")
            
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
                self.statusBar.showMessage(f"Removed {removed_count} mask regions in range [{v_min:.1f}, {v_max:.1f}]")
            else:
                self.statusBar.showMessage(f"No masks found in range [{v_min:.1f}, {v_max:.1f}]")
            
            # Reset current click and update plot
            self.current_click = None
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
        elif event.key == 'z':
            self.undo_last_mask()
        elif event.key == 'f':
            self.fit_continuum()
        elif event.key == 'a':
            self.accept_results()
        elif event.key == 'c' or event.key == 'escape':
            self.cancel()
        elif event.key == 'm':
            self.manual_mask_entry()
        # Add zoom shortcuts
        elif event.key == '+' or event.key == '=':
            self.zoom_in()
        elif event.key == '-' or event.key == '_':
            self.zoom_out()
        elif event.key == '0':
            self.reset_zoom()
    def zoom_in(self):
        # Get current axis limits
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
        self.canvas.draw()
    
    def zoom_out(self):
        # Similar to zoom_in but with expansion
        ax_flux = self.canvas.axes[0]
        xlim = ax_flux.get_xlim()
        ylim = ax_flux.get_ylim()
        
        x_center = (xlim[0] + xlim[1]) / 2
        y_center = (ylim[0] + ylim[1]) / 2
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]
        new_xlim = (x_center - 0.625*x_range, x_center + 0.625*x_range)
        new_ylim = (y_center - 0.625*y_range, y_center + 0.625*y_range)
        
        ax_flux.set_xlim(new_xlim)
        ax_flux.set_ylim(new_ylim)
        self.canvas.draw()
    
    def reset_zoom(self):
        # Reset to initial view
        ax_flux = self.canvas.axes[0]
        ax_flux.set_xlim(self.domain)
        ax_flux.autoscale(axis='y')
        self.canvas.draw()
        
    def reset_masks(self):
        """Reset all masks."""
        self.masks = []
        self.current_click = None
        self.statusBar.showMessage("All masks reset")
        self.update_plots()
    
    def undo_last_mask(self):
        """Remove the last mask pair."""
        if len(self.masks) >= 2:
            self.masks.pop()
            self.masks.pop()
            self.statusBar.showMessage("Last mask removed")
            self.update_plots()
        else:
            self.statusBar.showMessage("No masks to remove")
    

    def auto_mask(self):
        """
        Automatically generate masks using improved detection of spectral features.
        
        This method uses a robust approach to detect both absorption and emission features,
        as well as regions with unreliable error estimates.
        """
        # Get current parameters
        sigma = self.auto_mask_sigma_spin.value()  # Use specific auto-mask sigma
        min_width = self.min_mask_width_spin.value()  # Minimum mask width in km/s
        gap_threshold = self.mask_gap_spin.value()  # Maximum gap to group masks in km/s
        error_threshold = 0.5  # Mask points where error > 50% of continuum level
        
        # Ask if user wants to replace or add to existing masks
        if self.masks:
            reply = QMessageBox.question(self, 'Auto-Mask', 
                                         'Replace existing masks?',
                                         QMessageBox.Yes | QMessageBox.No, 
                                         QMessageBox.No)
            replace = reply == QMessageBox.Yes
        else:
            replace = True
        
        self.statusBar.showMessage("Generating masks automatically...")
        
        try:
            # Step 1: Create a mask of currently masked regions
            current_mask = np.ones_like(self.velocity, dtype=bool)
            
            if not replace and self.masks:
                for i in range(0, len(self.masks), 2):
                    if i+1 < len(self.masks):
                        v_min, v_max = self.masks[i], self.masks[i+1]
                        mask_indices = (self.velocity >= v_min) & (self.velocity <= v_max)
                        current_mask[mask_indices] = False
            
            # Step 2: Fit a preliminary continuum to unmasked data
            unmasked_vel = self.velocity[current_mask]
            unmasked_flux = self.flux[current_mask]
            unmasked_error = self.error[current_mask]
            
            # Use a low-order polynomial (order 2 or 3) for robustness
            temp_order = min(2, max(0, len(unmasked_vel) // 10 - 1))
            
            try:
                from rbcodes.IGM.rb_iter_contfit import rb_iter_contfit
                
                # Fit with sigma clipping to further exclude potential absorption features
                prelim_fit = rb_iter_contfit(
                    unmasked_vel,
                    unmasked_flux,
                    error=unmasked_error,
                    order=temp_order,
                    sigma=3.0,  # Conservative sigma value for initial fit
                    use_weights=False,
                    return_model=True
                )
                
                # Get the continuum model function
                prelim_continuum_model = prelim_fit['model']
                
                # Evaluate the continuum on the full velocity grid
                prelim_continuum = prelim_continuum_model(self.velocity)
                
            except Exception as e:
                self.statusBar.showMessage(f"Error in preliminary fit: {str(e)}. Using median flux instead.")
                # Fallback to a simple median if continuum fitting fails
                prelim_continuum = np.ones_like(self.velocity) * np.median(unmasked_flux)
            
            # Step 3: Calculate residuals and robust standard deviation
            # Normalize the flux
            normalized_flux = self.flux / prelim_continuum
            
            # Calculate residuals (deviations from the continuum)
            residual = normalized_flux - 1.0
            
            # Calculate robust standard deviation using MAD
            from astropy.stats import mad_std
            # Handle possible NaN values
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                robust_std = mad_std(residual, ignore_nan=True)
            
            # Step 4: Create mask using combined criteria
            # Points that deviate significantly from the continuum OR have large relative errors
            potential_masks = (np.abs(residual) > sigma * robust_std) | \
                              (self.error/prelim_continuum > error_threshold)
            
            # Find contiguous regions
            mask_regions = []
            in_region = False
            start_idx = None
            
            for i in range(len(potential_masks)):
                if potential_masks[i] and not in_region:
                    # Start of a new region
                    in_region = True
                    start_idx = i
                elif not potential_masks[i] and in_region:
                    # End of a region
                    in_region = False
                    end_idx = i - 1
                    mask_regions.append((start_idx, end_idx))
            
            # Add the last region if we were still in one
            if in_region:
                mask_regions.append((start_idx, len(potential_masks) - 1))
            
            # Step 5: Filter regions by minimum width
            filtered_regions = []
            
            for start_idx, end_idx in mask_regions:
                start_vel = self.velocity[start_idx]
                end_vel = self.velocity[end_idx]
                
                # Check if the region is wide enough
                if end_vel - start_vel >= min_width:
                    filtered_regions.append((start_idx, end_idx))
            
            # Step 6: Group nearby regions (merge regions with small gaps)
            if len(filtered_regions) > 1:
                merged_regions = [filtered_regions[0]]
                
                for current_start, current_end in filtered_regions[1:]:
                    # Get the last merged region
                    last_start, last_end = merged_regions[-1]
                    
                    # Check if the gap is small enough to merge
                    if self.velocity[current_start] - self.velocity[last_end] <= gap_threshold:
                        # Merge with the previous region
                        merged_regions[-1] = (last_start, current_end)
                    else:
                        # Add as a new region
                        merged_regions.append((current_start, current_end))
                
                filtered_regions = merged_regions
            
            # Step 7: Convert to velocity pairs for masking
            new_masks = []
            
            for start_idx, end_idx in filtered_regions:
                v_min = self.velocity[start_idx]
                v_max = self.velocity[end_idx]
                new_masks.extend([v_min, v_max])
            
            # Step 8: Merge with existing masks or replace
            if replace:
                self.masks = new_masks
            else:
                self.masks.extend(new_masks)
                
                # Merge overlapping masks
                self.merge_masks()
            
            # Update plots
            count = len(new_masks) // 2
            detection_info = f"MAD = {robust_std:.3f}, error threshold = {error_threshold:.2f}"
            self.statusBar.showMessage(f"Added {count} mask regions using sigma={sigma}, {detection_info}")
            self.update_plots()
            
        except Exception as e:
            self.statusBar.showMessage(f"Error in auto-masking: {str(e)}")
            QMessageBox.warning(self, "Auto-Masking Error", 
                              f"Error generating masks: {str(e)}\n\n"
                              f"Try adjusting parameters or masking manually.")
            
    def manual_mask_entry(self):
        """Open a dialog for manual entry of mask values."""
        # Create dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Manual Mask Entry")
        layout = QVBoxLayout(dialog)
        
        # Create form for min/max values
        form_layout = QFormLayout()
        
        min_spin = QDoubleSpinBox()
        min_spin.setRange(-10000, 10000)
        min_spin.setValue(0)
        
        max_spin = QDoubleSpinBox()
        max_spin.setRange(-10000, 10000)
        max_spin.setValue(0)
        
        form_layout.addRow("Minimum Velocity:", min_spin)
        form_layout.addRow("Maximum Velocity:", max_spin)
        
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
                self.statusBar.showMessage(f"Mask added: {v_min:.1f} to {v_max:.1f}")
                self.update_plots()
            else:
                QMessageBox.warning(self, "Invalid Input", "Minimum and maximum values must be different")
    
    def fit_continuum(self):
        """Fit continuum using current parameters and masks."""
        # Get current parameters
        self.fit_params['order'] = self.order_spin.value()
        self.fit_params['use_weights'] = self.use_weights_check.isChecked()
        self.fit_params['optimize_cont'] = self.optimize_check.isChecked()
        self.fit_params['min_order'] = self.min_order_spin.value()
        self.fit_params['max_order'] = self.max_order_spin.value()
        self.fit_params['sigma'] = self.sigma_spin.value()
        
        # Check if we have enough unmasked points
        if not self.check_unmasked_points():
            return
        
        # Create a mask array (0 for masked, 1 for unmasked)
        mask_array = np.ones_like(self.velocity, dtype=bool)
        
        for i in range(0, len(self.masks), 2):
            if i+1 < len(self.masks):
                v_min, v_max = self.masks[i], self.masks[i+1]
                mask_indices = (self.velocity >= v_min) & (self.velocity <= v_max)
                mask_array[mask_indices] = False
        
        # Extract unmasked points
        unmasked_vel = self.velocity[mask_array]
        unmasked_flux = self.flux[mask_array]
        unmasked_error = self.error[mask_array] if self.fit_params['use_weights'] else None
        
        # Fit using optimal polynomial if requested
        self.statusBar.showMessage("Fitting continuum...")
        
        try:
            if self.fit_params['optimize_cont']:
                result = fit_optimal_polynomial(
                    unmasked_vel, 
                    unmasked_flux, 
                    error=unmasked_error,
                    min_order=self.fit_params['min_order'],
                    max_order=self.fit_params['max_order'],
                    use_weights=self.fit_params['use_weights'],
                    sigma=self.fit_params['sigma'],
                    plot=False,
                    include_model=True
                )
                
                self.fit_result = result
                best_order = result['best_order']
                self.continuum = result['model'](self.velocity)
                
                message = f"Best polynomial order: {best_order}, BIC: {min([b for _, b in result['bic_results']]):.2f}"
            else:
                # Use fixed polynomial order
                from rbcodes.IGM.rb_iter_contfit import rb_iter_contfit
                
                result = rb_iter_contfit(
                    unmasked_vel,
                    unmasked_flux,
                    error=unmasked_error,
                    order=self.fit_params['order'],
                    sigma=self.fit_params['sigma'],
                    use_weights=self.fit_params['use_weights'],
                    return_model=True
                )
                
                self.fit_result = result
                self.continuum = result['model'](self.velocity)
                
                message = f"Fixed polynomial order: {self.fit_params['order']}, Fit error: {result['fit_error']:.5f}"
            
            # Calculate normalized flux
            self.normalized_flux = self.flux / self.continuum
            
            # Update status and enable accept button
            self.statusBar.showMessage(f"Continuum fitted. {message}")
            self.accept_btn.setEnabled(True)
            
            # Update plots
            self.update_plots()
            
        except Exception as e:
            self.statusBar.showMessage(f"Error fitting continuum: {str(e)}")
            QMessageBox.warning(self, "Fitting Error", 
                              f"Error fitting continuum: {str(e)}\n\n"
                              f"Try adjusting mask regions or fitting parameters.")
    
    def check_unmasked_points(self):
        """Check if there are enough unmasked points for fitting."""
        # Count unmasked points
        mask_array = np.ones_like(self.velocity, dtype=bool)
        
        for i in range(0, len(self.masks), 2):
            if i+1 < len(self.masks):
                v_min, v_max = self.masks[i], self.masks[i+1]
                mask_indices = (self.velocity >= v_min) & (self.velocity <= v_max)
                mask_array[mask_indices] = False
        
        unmasked_count = np.sum(mask_array)
        
        # For optimal polynomial, need enough points for highest order + 1
        if self.fit_params['optimize_cont']:
            min_required = self.fit_params['max_order'] + 2
        else:
            min_required = self.fit_params['order'] + 2
        
        if unmasked_count < min_required:
            QMessageBox.warning(self, "Insufficient Data", 
                              f"Not enough unmasked points ({unmasked_count}) "
                              f"for fitting polynomial of order {min_required-2}.\n\n"
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
                    return
            else:
                return
        
        # Convert mask regions to wavelength
        mask_wavelengths = []
        for i in range(0, len(self.masks), 2):
            if i+1 < len(self.masks):
                v_min, v_max = self.masks[i], self.masks[i+1]
                # Find corresponding wavelength ranges
                wmin = self.wave[np.abs(self.velocity - v_min).argmin()]
                wmax = self.wave[np.abs(self.velocity - v_max).argmin()]
                mask_wavelengths.extend([wmin, wmax])
        
        # Prepare result
        self.result = {
            'masks': self.masks,
            'mask_wavelengths': mask_wavelengths,
            'continuum': self.continuum,
            'fit_params': self.fit_params,
            'cancelled': False
        }
        
        # Add additional info if available
        if self.fit_result is not None:
            if 'fit_error' in self.fit_result:
                self.result['fit_error'] = self.fit_result['fit_error']
            if 'best_order' in self.fit_result:
                self.result['best_order'] = self.fit_result['best_order']
        
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

def launch_interactive_mask(**kwargs):
    """
    Launch the interactive masking GUI.
    
    Parameters
    ----------
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    velocity : array
        Velocity array
    existing_masks : list, optional
        Existing mask regions
    order : int, optional
        Initial polynomial order
    use_weights : bool, optional
        Whether to use weights
    domain : list, optional
        Velocity domain limits
        
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
    window = InteractiveMaskWindow(**kwargs)
    window.show()
    
    # Run the application
    app.exec_()
    
    # Return the result
    return window.result

if __name__ == "__main__":
    # Example usage
    import numpy as np
    
    # Generate sample data
    n_points = 1000
    velocity = np.linspace(-1000, 1000, n_points)
    continuum = 1.0 + 0.001 * velocity + 0.0000002 * velocity**2
    absorption = 1.0 - 0.5 * np.exp(-(velocity+200)**2/50**2) - 0.8 * np.exp(-(velocity-150)**2/30**2)
    flux = continuum * absorption + 0.05 * np.random.randn(n_points)
    error = 0.05 * np.ones_like(flux)
    
    # Generate wavelength
    wave = 1215.67 * (1.0 + velocity/299792.458)
    
    # Launch the tool
    result = launch_interactive_mask(
        wave=wave,
        flux=flux,
        error=error,
        velocity=velocity,
        existing_masks=[-300, -100],
        order=3,
        use_weights=False,
        domain=[-500, 500]
    )
    
    # Print the result
    if not result.get('cancelled', False):
        print("Results:")
        print(f"Masks: {result['masks']}")
        print(f"Fit Parameters: {result['fit_params']}")
        if 'best_order' in result:
            print(f"Best Order: {result['best_order']}")                              