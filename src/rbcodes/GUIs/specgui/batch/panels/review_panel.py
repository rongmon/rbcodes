# rbcodes/GUIs/specgui/batch/panels/review_panel.py - COMPLETE IMPLEMENTATION
import os
import numpy as np
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QComboBox, QGroupBox, QSplitter,
                           QMessageBox, QFileDialog, QDialog, QDialogButtonBox, 
                           QFormLayout, QSpinBox, QRadioButton, QProgressBar,
                           QCheckBox, QDoubleSpinBox)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from datetime import datetime    
from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object
import copy

# Import the plotting module
from rbcodes.GUIs.specgui.batch.panels.spectrum_plotter import plot_spectrum_overview

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots."""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)

class AxisLimitsDialog(QDialog):
    """Dialog for setting custom X/Y axis limits and analysis settings."""
    
    def __init__(self, current_xlim, current_ylim_upper, current_ylim_lower, current_snr_settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Plot & Analysis Settings")
        self.setMinimumSize(450, 300)
        
        layout = QVBoxLayout(self)
        
        # Form layout
        form = QFormLayout()
        
        # X-axis range (applies to both panels)
        x_layout = QHBoxLayout()
        self.x_min = QDoubleSpinBox()
        self.x_min.setRange(-10000, 10000)
        self.x_min.setValue(current_xlim[0])
        self.x_min.setDecimals(1)
        self.x_min.setSuffix(" km/s")
        
        self.x_max = QDoubleSpinBox()
        self.x_max.setRange(-10000, 10000)
        self.x_max.setValue(current_xlim[1])
        self.x_max.setDecimals(1)
        self.x_max.setSuffix(" km/s")
        
        self.x_auto = QCheckBox("Auto")
        self.x_auto.toggled.connect(self.toggle_x_auto)
        
        x_layout.addWidget(QLabel("Min:"))
        x_layout.addWidget(self.x_min)
        x_layout.addWidget(QLabel("Max:"))
        x_layout.addWidget(self.x_max)
        x_layout.addWidget(self.x_auto)
        
        form.addRow("X-axis (Velocity):", x_layout)
        
        # Y-axis range for upper panel (flux + continuum)
        y_upper_layout = QHBoxLayout()
        self.y_upper_min = QDoubleSpinBox()
        self.y_upper_min.setRange(-10.0, 20.0)
        self.y_upper_min.setValue(current_ylim_upper[0])
        self.y_upper_min.setDecimals(2)
        
        self.y_upper_max = QDoubleSpinBox()
        self.y_upper_max.setRange(-10.0, 20.0)
        self.y_upper_max.setValue(current_ylim_upper[1])
        self.y_upper_max.setDecimals(2)
        
        self.y_upper_auto = QCheckBox("Auto")
        self.y_upper_auto.toggled.connect(self.toggle_y_upper_auto)
        
        y_upper_layout.addWidget(QLabel("Min:"))
        y_upper_layout.addWidget(self.y_upper_min)
        y_upper_layout.addWidget(QLabel("Max:"))
        y_upper_layout.addWidget(self.y_upper_max)
        y_upper_layout.addWidget(self.y_upper_auto)
        
        form.addRow("Y-axis Upper (Flux):", y_upper_layout)
        
        # Y-axis range for lower panel (normalized flux)
        y_lower_layout = QHBoxLayout()
        self.y_lower_min = QDoubleSpinBox()
        self.y_lower_min.setRange(-2.0, 5.0)
        self.y_lower_min.setValue(current_ylim_lower[0])
        self.y_lower_min.setDecimals(2)
        
        self.y_lower_max = QDoubleSpinBox()
        self.y_lower_max.setRange(-2.0, 5.0)
        self.y_lower_max.setValue(current_ylim_lower[1])
        self.y_lower_max.setDecimals(2)
        
        self.y_lower_auto = QCheckBox("Auto")
        self.y_lower_auto.toggled.connect(self.toggle_y_lower_auto)
        
        y_lower_layout.addWidget(QLabel("Min:"))
        y_lower_layout.addWidget(self.y_lower_min)
        y_lower_layout.addWidget(QLabel("Max:"))
        y_lower_layout.addWidget(self.y_lower_max)
        y_lower_layout.addWidget(self.y_lower_auto)
        
        form.addRow("Y-axis Lower (Normalized):", y_lower_layout)
        
        # Separator
        separator = QLabel()
        separator.setStyleSheet("QLabel { border-bottom: 1px solid gray; margin: 10px 0; }")
        form.addRow(separator)
        
        # SNR Analysis Settings
        snr_layout = QHBoxLayout()
        
        self.snr_enabled = QCheckBox("Calculate SNR")
        self.snr_enabled.setChecked(current_snr_settings.get('calculate_snr', True))
        self.snr_enabled.toggled.connect(self.toggle_snr_controls)
        
        self.binsize_spin = QSpinBox()
        self.binsize_spin.setRange(1, 10)
        self.binsize_spin.setValue(current_snr_settings.get('binsize', 3))
        self.binsize_spin.setEnabled(self.snr_enabled.isChecked())
        
        snr_layout.addWidget(self.snr_enabled)
        snr_layout.addWidget(QLabel("Bin size:"))
        snr_layout.addWidget(self.binsize_spin)
        snr_layout.addStretch()
        
        form.addRow("SNR Settings:", snr_layout)
        
        layout.addLayout(form)
        
        # Buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
    
    def toggle_x_auto(self, enabled):
        """Toggle X-axis manual controls."""
        self.x_min.setEnabled(not enabled)
        self.x_max.setEnabled(not enabled)
    
    def toggle_y_upper_auto(self, enabled):
        """Toggle upper Y-axis manual controls."""
        self.y_upper_min.setEnabled(not enabled)
        self.y_upper_max.setEnabled(not enabled)
    
    def toggle_y_lower_auto(self, enabled):
        """Toggle lower Y-axis manual controls."""
        self.y_lower_min.setEnabled(not enabled)
        self.y_lower_max.setEnabled(not enabled)
    
    def toggle_snr_controls(self, enabled):
        """Toggle SNR bin size control."""
        self.binsize_spin.setEnabled(enabled)
    
    def get_limits_and_settings(self):
        """Get the selected limits and analysis settings."""
        if self.x_auto.isChecked():
            xlim = None  # Auto
        else:
            xlim = (self.x_min.value(), self.x_max.value())
        
        if self.y_upper_auto.isChecked():
            ylim_upper = None  # Auto
        else:
            ylim_upper = (self.y_upper_min.value(), self.y_upper_max.value())
        
        if self.y_lower_auto.isChecked():
            ylim_lower = None  # Auto
        else:
            ylim_lower = (self.y_lower_min.value(), self.y_lower_max.value())
        
        snr_settings = {
            'calculate_snr': self.snr_enabled.isChecked(),
            'binsize': self.binsize_spin.value()
        }
        
        return xlim, ylim_upper, ylim_lower, snr_settings

class ReviewPanel(QWidget):
    """Panel for reviewing batch processing results - with navigation and filtering."""
    
    items_selected = pyqtSignal(list)
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.current_item = None
        self.current_spec = None
        self.current_index = -1  # Track current system index
        self.current_filter = "all"  # Track current filter
        self.filtered_items = []     # Cache filtered items
        self.importing_in_progress = False  # Flag to suppress warnings during import
        
        # Interactive EW selection state
        self.interactive_ew_mode = False
        self.ew_click_connections = []
        self.ew_click_count = 0  # Track two-click sequence
        
        # Custom axis limits
        self.custom_xlim = None
        self.custom_ylim_upper = None
        self.custom_ylim_lower = None

        self.init_ui()
        
        self.controller.table_changed.connect(self.refresh_results_table)
    
    def init_ui(self):
        """Initialize UI with navigation controls and filtering."""
        main_layout = QVBoxLayout(self)
        
        # Top navigation bar
        nav_group = QGroupBox("System Navigation")
        nav_layout = QVBoxLayout()
        
        # System selector row
        selector_layout = QHBoxLayout()
        
        # System dropdown
        selector_layout.addWidget(QLabel("System:"))
        self.system_combo = QComboBox()
        self.system_combo.setMinimumWidth(300)
        self.system_combo.currentIndexChanged.connect(self.on_system_selected)
        selector_layout.addWidget(self.system_combo)
        
        selector_layout.addStretch()
        
        # Position indicator
        self.position_label = QLabel("No systems loaded")
        selector_layout.addWidget(self.position_label)
        
        # Selection actions
        self.select_for_processing_btn = QPushButton("Select for Processing")
        self.select_for_processing_btn.clicked.connect(self.open_batch_selection_dialog)
        self.select_for_processing_btn.setEnabled(True)
        selector_layout.addWidget(self.select_for_processing_btn)
        
        nav_layout.addLayout(selector_layout)
        
        # Navigation buttons row
        button_layout = QHBoxLayout()
        
        self.first_btn = QPushButton("⏮ First")
        self.first_btn.clicked.connect(self.navigate_to_first)
        self.first_btn.setEnabled(False)
        
        self.prev_btn = QPushButton("⏪ Previous")
        self.prev_btn.clicked.connect(self.navigate_to_previous)
        self.prev_btn.setEnabled(False)
        
        self.next_btn = QPushButton("Next ⏩")
        self.next_btn.clicked.connect(self.navigate_to_next)
        self.next_btn.setEnabled(False)
        
        self.last_btn = QPushButton("Last ⏭")
        self.last_btn.clicked.connect(self.navigate_to_last)
        self.last_btn.setEnabled(False)
        
        button_layout.addWidget(self.first_btn)
        button_layout.addWidget(self.prev_btn)
        button_layout.addWidget(self.next_btn)
        button_layout.addWidget(self.last_btn)
        button_layout.addStretch()
        
        # Edit Continuum button in navigation row for more space
        self.edit_continuum_btn = QPushButton("Edit Continuum")
        self.edit_continuum_btn.clicked.connect(self.edit_continuum)
        self.edit_continuum_btn.setEnabled(False)
        button_layout.addWidget(self.edit_continuum_btn)
        
        nav_layout.addLayout(button_layout)
        
        # Status filter radio buttons + EW controls (Edit Continuum moved to navigation row)
        filter_layout = QHBoxLayout()
        
        filter_layout.addWidget(QLabel("Filter:"))
        
        self.filter_all = QRadioButton("All")
        self.filter_complete = QRadioButton("Complete") 
        self.filter_failed = QRadioButton("Failed")
        self.filter_pending = QRadioButton("Pending")
        
        # Set "All" as default
        self.filter_all.setChecked(True)
        
        # Connect radio button changes
        self.filter_all.toggled.connect(lambda: self.on_filter_changed("all"))
        self.filter_complete.toggled.connect(lambda: self.on_filter_changed("complete"))
        self.filter_failed.toggled.connect(lambda: self.on_filter_changed("failed"))
        self.filter_pending.toggled.connect(lambda: self.on_filter_changed("pending"))
        
        filter_layout.addWidget(self.filter_all)
        filter_layout.addWidget(self.filter_complete)
        filter_layout.addWidget(self.filter_failed)
        filter_layout.addWidget(self.filter_pending)
        
        filter_layout.addStretch()
        
        # EW Range controls - now with more space for larger Update button
        filter_layout.addWidget(QLabel("EW:"))
        
        filter_layout.addWidget(QLabel("vmin"))
        self.ew_vmin_spin = QSpinBox()
        self.ew_vmin_spin.setRange(-15000, 150000)
        self.ew_vmin_spin.setValue(-200)
        self.ew_vmin_spin.setSuffix(" km/s")
        self.ew_vmin_spin.setMaximumWidth(70)
        filter_layout.addWidget(self.ew_vmin_spin)
        
        filter_layout.addWidget(QLabel("vmax"))
        self.ew_vmax_spin = QSpinBox()
        self.ew_vmax_spin.setRange(-15000, 150000)
        self.ew_vmax_spin.setValue(200)
        self.ew_vmax_spin.setSuffix(" km/s")
        self.ew_vmax_spin.setMaximumWidth(70)
        filter_layout.addWidget(self.ew_vmax_spin)
        
        self.update_ew_btn = QPushButton("UPDATE")
        self.update_ew_btn.clicked.connect(self.update_ew_range)
        self.update_ew_btn.setEnabled(False)
        self.update_ew_btn.setMaximumWidth(80)  # Larger button
        self.update_ew_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        filter_layout.addWidget(self.update_ew_btn)
        
        self.interactive_ew_checkbox = QCheckBox("Interactive")
        self.interactive_ew_checkbox.toggled.connect(self.toggle_interactive_ew)
        self.interactive_ew_checkbox.setEnabled(False)
        filter_layout.addWidget(self.interactive_ew_checkbox)
        
        filter_layout.addStretch()
        
        nav_layout.addLayout(filter_layout)
        
        nav_group.setLayout(nav_layout)
        main_layout.addWidget(nav_group)
        
        # Large preview area (85% of space)
        preview_group = QGroupBox("System Preview")
        preview_layout = QVBoxLayout()
        
        # Large canvas
        self.spectrum_canvas = MatplotlibCanvas(self, width=10, height=8, dpi=100)
        self.spectrum_toolbar = NavigationToolbar2QT(self.spectrum_canvas, self)
        
        preview_layout.addWidget(self.spectrum_toolbar)
        preview_layout.addWidget(self.spectrum_canvas)
        
        preview_group.setLayout(preview_layout)
        main_layout.addWidget(preview_group)
        
        # Compact details at bottom (single line)
        self.details_label = QLabel("Select a system to view details.")
        self.details_label.setWordWrap(True)
        self.details_label.setMaximumHeight(25)  # Single line height
        self.details_label.setStyleSheet("QLabel { font-size: 8pt; padding: 3px; background-color: #f0f0f0; border: 1px solid #ccc; }")
        main_layout.addWidget(self.details_label)
        
        # Set up keyboard shortcuts
        self.setFocusPolicy(Qt.StrongFocus)
    
    def keyPressEvent(self, event):
        """Handle keyboard shortcuts."""
        if event.key() == Qt.Key_R and event.modifiers() == Qt.ShiftModifier:
            self.open_axis_limits_dialog()
        else:
            super().keyPressEvent(event)
    
    def open_axis_limits_dialog(self):
        """Open dialog for setting custom axis limits and analysis settings."""
        if not self.current_spec:
            self.controller.status_updated.emit("Please select a system first.")
            return
        
        # Get current axis limits from both panels
        axes = self.spectrum_canvas.fig.axes
        if len(axes) >= 2:
            # Two-panel plot
            upper_ax = axes[0]  # Flux panel
            lower_ax = axes[1]  # Normalized panel
            current_xlim = upper_ax.get_xlim()
            current_ylim_upper = upper_ax.get_ylim()
            current_ylim_lower = lower_ax.get_ylim()
        elif len(axes) == 1:
            # Single panel
            ax = axes[0]
            current_xlim = ax.get_xlim()
            current_ylim_upper = ax.get_ylim()
            current_ylim_lower = (0, 2)  # Default for normalized
        else:
            # No axes
            current_xlim = (-1500, 1500)
            current_ylim_upper = (0, 2)
            current_ylim_lower = (0, 2)
        
        # Get current SNR settings from current item or defaults
        if self.current_item:
            current_snr_settings = {
                'calculate_snr': self.current_item.analysis.calculate_snr,
                'binsize': self.current_item.analysis.binsize
            }
        else:
            current_snr_settings = {
                'calculate_snr': True,
                'binsize': 3
            }
        
        # Open dialog
        dialog = AxisLimitsDialog(current_xlim, current_ylim_upper, current_ylim_lower, current_snr_settings, self)
        
        if dialog.exec_() == QDialog.Accepted:
            xlim, ylim_upper, ylim_lower, snr_settings = dialog.get_limits_and_settings()
            
            # Store custom limits
            self.custom_xlim = xlim
            self.custom_ylim_upper = ylim_upper
            self.custom_ylim_lower = ylim_lower
            
            # Apply to current plot
            self.apply_custom_limits()
            
            # Update SNR settings in master table if we have a current item
            if self.current_item and self.current_index >= 0:
                self.controller.master_table.update_analysis(
                    self.current_index,
                    calculate_snr=snr_settings['calculate_snr'],
                    binsize=snr_settings['binsize']
                )
            
            self.controller.status_updated.emit("Settings applied. Use Shift+R to modify.")
    
    def apply_custom_limits(self):
        """Apply custom axis limits to the current plot."""
        axes = self.spectrum_canvas.fig.axes
        if not axes:
            return
        
        # Apply limits based on number of panels
        if len(axes) >= 2:
            # Two-panel plot
            upper_ax = axes[0]  # Flux + continuum panel
            lower_ax = axes[1]  # Normalized flux panel
            
            # Apply X limits to both panels
            if self.custom_xlim:
                upper_ax.set_xlim(self.custom_xlim)
                lower_ax.set_xlim(self.custom_xlim)
            
            # Apply Y limits to respective panels
            if self.custom_ylim_upper:
                upper_ax.set_ylim(self.custom_ylim_upper)
            
            if self.custom_ylim_lower:
                lower_ax.set_ylim(self.custom_ylim_lower)
                
        elif len(axes) == 1:
            # Single panel
            ax = axes[0]
            if self.custom_xlim:
                ax.set_xlim(self.custom_xlim)
            if self.custom_ylim_upper:
                ax.set_ylim(self.custom_ylim_upper)
        
        self.spectrum_canvas.draw()
    
    def toggle_interactive_ew(self):
        """Toggle interactive EW range selection mode."""
        if not self.current_spec:
            self.controller.status_updated.emit("Please select a system first.")
            self.interactive_ew_checkbox.setChecked(False)
            return
        
        self.interactive_ew_mode = self.interactive_ew_checkbox.isChecked()
        
        if self.interactive_ew_mode:
            # Enable interactive mode
            self.connect_ew_interactions()
            self.spectrum_canvas.setFocus()
            self.controller.status_updated.emit("Interactive EW mode enabled. Click on plot to set velocity limits. Press ESC to exit.")
        else:
            # Disable interactive mode
            self.disconnect_ew_interactions()
            self.controller.status_updated.emit("Interactive EW mode disabled.")
    
    def connect_ew_interactions(self):
        """Connect mouse and key events for interactive EW selection."""
        # Disconnect any existing connections
        self.disconnect_ew_interactions()
        
        # Connect new events
        click_cid = self.spectrum_canvas.mpl_connect('button_press_event', self.on_ew_plot_click)
        key_cid = self.spectrum_canvas.mpl_connect('key_press_event', self.on_ew_plot_key)
        
        self.ew_click_connections = [click_cid, key_cid]
    
    def disconnect_ew_interactions(self):
        """Disconnect interactive EW events."""
        for cid in self.ew_click_connections:
            self.spectrum_canvas.mpl_disconnect(cid)
        self.ew_click_connections = []
    
    def on_ew_plot_click(self, event):
        """Handle mouse clicks for EW range selection."""
        if not self.interactive_ew_mode or not event.inaxes:
            return
        
        # Get clicked velocity
        velocity = int(event.xdata)
        
        # Get current values
        vmin = self.ew_vmin_spin.value()
        vmax = self.ew_vmax_spin.value()
        
        # Case 1: Click outside current range → set nearest bound
        if velocity <= vmin:
            vmin = velocity
        elif velocity >= vmax:
            vmax = velocity
        else:
            # Case 2: Inside range → adjust whichever side is closer
            if abs(velocity - vmin) < abs(velocity - vmax):
                vmin = velocity
            else:
                vmax = velocity
        
        # Ensure ordering (vmin < vmax)
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        
        # Block signals during updates to avoid redundant redraws
        self.ew_vmin_spin.blockSignals(True)
        self.ew_vmax_spin.blockSignals(True)
        self.ew_vmin_spin.setValue(vmin)
        self.ew_vmax_spin.setValue(vmax)
        self.ew_vmin_spin.blockSignals(False)
        self.ew_vmax_spin.blockSignals(False)
        
        # Optimized redraw - only update EW region visualization
        self.update_ew_region_overlay()
        
        # Feedback
        self.controller.status_updated.emit(
            f"EW boundaries set to {vmin}–{vmax} km/s (clicked {velocity}). Click Update to apply changes."
        )
    
    def update_ew_region_overlay(self):
        """Optimized update of just the EW region visualization."""
        try:
            axes = self.spectrum_canvas.fig.axes
            if not axes:
                return
            
            # Get current EW range from spinboxes
            ew_vmin = self.ew_vmin_spin.value()
            ew_vmax = self.ew_vmax_spin.value()
            
            # Update EW region in all axes (both upper and lower panels)
            for ax in axes:
                # Remove existing EW region spans
                for collection in ax.collections[:]:
                    if hasattr(collection, '_ew_region_marker'):
                        collection.remove()
                
                # Remove existing EW boundary lines
                for line in ax.lines[:]:
                    if hasattr(line, '_ew_boundary_marker'):
                        line.remove()
                
                # Add new EW region
                span = ax.axvspan(ew_vmin, ew_vmax, alpha=0.15, color='blue', label='EW Region')
                span._ew_region_marker = True  # Mark for removal
                
                # Add new boundary lines
                vmin_line = ax.axvline(x=ew_vmin, color='b', linestyle=':', alpha=0.8, linewidth=0.8)
                vmax_line = ax.axvline(x=ew_vmax, color='b', linestyle=':', alpha=0.8, linewidth=0.8)
                vmin_line._ew_boundary_marker = True
                vmax_line._ew_boundary_marker = True
            
            # Quick redraw without full recomputation
            self.spectrum_canvas.draw_idle()
            
        except Exception as e:
            print(f"Error updating EW region overlay: {e}")
            # Fallback to full plot update if optimized version fails
            self.update_spectrum_plot(self.current_item, self.current_spec)
    
    def on_ew_plot_key(self, event):
        """Handle key presses in interactive EW mode."""
        if not self.interactive_ew_mode:
            return
        
        if event.key == 'escape':
            # Exit interactive mode
            self.interactive_ew_checkbox.setChecked(False)
            self.toggle_interactive_ew()
        elif event.key == 'r':
            # Reset to defaults
            self.ew_vmin_spin.setValue(-200)
            self.ew_vmax_spin.setValue(200)
            self.update_ew_region_overlay()
            self.controller.status_updated.emit("EW range reset to defaults.")

    def update_status_bar_info(self):
        """Send comprehensive info to main window status bar."""
        total_count, complete_count, failed_count = self.get_overall_statistics()
        filtered_count = len(self.filtered_items)
        current_pos = self.get_current_position_in_filter()
        
        # Build comprehensive status message
        if filtered_count > 0:
            filter_name = self.current_filter.title()
            progress_percent = int((complete_count / total_count) * 100) if total_count > 0 else 0
            
            status_msg = (f"Progress: {progress_percent}% | "
                         f"Overall: {complete_count}/{total_count} complete, {failed_count} failed | "
                         f"{filter_name} Filter: System {current_pos} of {filtered_count}")
        else:
            progress_percent = int((complete_count / total_count) * 100) if total_count > 0 else 0
            status_msg = f"Progress: {progress_percent}% | Overall: {complete_count}/{total_count} complete, {failed_count} failed | No systems match filter"
        
        # Send to main window status bar
        self.controller.status_updated.emit(status_msg)

    def on_filter_changed(self, filter_type):
        """Handle filter radio button changes."""
        if not self.sender().isChecked():  # Only respond to the selected radio button
            return
            
        self.current_filter = filter_type
        self.refresh_results_table()  # This will repopulate with filtered items
    
    def get_filtered_items(self):
        """Get items based on current filter."""
        all_items = self.controller.master_table.get_all_items()
        
        if self.current_filter == "all":
            return list(enumerate(all_items))  # Return (index, item) pairs
        
        filtered = []
        for i, item in enumerate(all_items):
            if self.current_filter == "complete":
                # Complete: status complete AND valid W/W_e
                if (item.analysis.processing_status == "complete" and 
                    hasattr(item.results, 'W') and hasattr(item.results, 'W_e') and
                    not np.isnan(item.results.W) and not np.isnan(item.results.W_e)):
                    filtered.append((i, item))
                    
            elif self.current_filter == "failed":
                # Failed: status error OR NaN W/W_e
                is_error_status = item.analysis.processing_status == "error"
                has_nan_results = (hasattr(item.results, 'W') and hasattr(item.results, 'W_e') and
                                 (np.isnan(item.results.W) or np.isnan(item.results.W_e)))
                if is_error_status or has_nan_results:
                    filtered.append((i, item))
                    
            elif self.current_filter == "pending":
                # Pending: ready, needs_processing, etc.
                if item.analysis.processing_status in ["ready", "needs_processing", "needs_ew_recalc"]:
                    filtered.append((i, item))
        
        return filtered
    
    def get_overall_statistics(self):
        """Calculate overall statistics for all items."""
        all_items = self.controller.master_table.get_all_items()
        
        total_count = len(all_items)
        complete_count = 0
        failed_count = 0
        
        for item in all_items:
            # Complete: status complete AND valid W/W_e
            if (item.analysis.processing_status == "complete" and 
                hasattr(item.results, 'W') and hasattr(item.results, 'W_e') and
                not np.isnan(item.results.W) and not np.isnan(item.results.W_e)):
                complete_count += 1
            
            # Failed: status error OR NaN W/W_e
            is_error_status = item.analysis.processing_status == "error"
            has_nan_results = (hasattr(item.results, 'W') and hasattr(item.results, 'W_e') and
                             (np.isnan(item.results.W) or np.isnan(item.results.W_e)))
            if is_error_status or has_nan_results:
                failed_count += 1
        
        return total_count, complete_count, failed_count
    
    def get_current_position_in_filter(self):
        """Get the current position within the filtered items."""
        if self.current_index < 0 or not self.filtered_items:
            return 0
        
        # Find position of current_index in filtered items
        for pos, (original_index, item) in enumerate(self.filtered_items):
            if original_index == self.current_index:
                return pos + 1  # 1-based position
        
        return 0
    
    def refresh_results_table(self):
        """Populate system dropdown from master table with current filter."""
        # Get filtered items
        self.filtered_items = self.get_filtered_items()
        
        # Preserve current selection if possible
        current_selection = self.current_index
        
        # Clear current state
        self.system_combo.clear()
        self.current_item = None
        self.current_spec = None
        self.current_index = -1
        self.update_navigation_state()
        self.clear_preview()
        
        if not self.filtered_items:
            self.position_label.setText("No systems match filter")
            self.update_status_bar_info()
            return
        
        # Populate dropdown with filtered items
        for i, (original_index, item) in enumerate(self.filtered_items):
            # Create descriptive label
            filename = os.path.basename(item.template.filename)
            transition_name = item.template.transition_name
            redshift = item.template.redshift
            status = item.analysis.processing_status
            
            # Status indicator
            if status == 'complete':
                # Check if results are valid
                if (hasattr(item.results, 'W') and hasattr(item.results, 'W_e') and
                    not np.isnan(item.results.W) and not np.isnan(item.results.W_e)):
                    status_icon = "✓"
                else:
                    status_icon = "✗"  # Complete but invalid results
            elif status == 'error':
                status_icon = "✗"
            elif 'processing' in status:
                status_icon = "⏳"
            else:
                status_icon = "⚠️"
            
            label = f"{status_icon} System {original_index+1}: {filename} | {transition_name} | z={redshift:.4f}"
            self.system_combo.addItem(label)
        
        # Update position indicator
        filter_name = self.current_filter.title()
        filtered_count = len(self.filtered_items)
        self.position_label.setText(f"{filter_name}: {filtered_count} systems")
        
        # Enable navigation and selection if we have items
        if filtered_count > 0:
            self.update_navigation_state()
            self.select_for_processing_btn.setEnabled(True)
            
            # Try to restore selection or select first
            restored = False
            if current_selection >= 0:
                # Look for the same original index in filtered items
                for combo_index, (original_index, item) in enumerate(self.filtered_items):
                    if original_index == current_selection:
                        self.system_combo.setCurrentIndex(combo_index)
                        restored = True
                        break
            
            if not restored:
                self.system_combo.setCurrentIndex(0)
        else:
            self.select_for_processing_btn.setEnabled(False)
        
        # Update progress display
        self.update_status_bar_info()
    
    def on_system_selected(self, combo_index):
        """Handle system selection from dropdown."""
        if combo_index < 0 or combo_index >= len(self.filtered_items):
            return
            
        # Get the original index and item from filtered list
        original_index, selected_item = self.filtered_items[combo_index]
        
        # Update current tracking
        self.current_index = original_index
        self.current_item = selected_item
        
        self.display_item_spectrum(selected_item)
        self.edit_continuum_btn.setEnabled(True)
        self.update_ew_btn.setEnabled(True)
        self.interactive_ew_checkbox.setEnabled(True)
        self.select_for_processing_btn.setEnabled(True)
        
        # Update EW spinboxes with current values
        self.ew_vmin_spin.setValue(selected_item.template.ew_vmin)
        self.ew_vmax_spin.setValue(selected_item.template.ew_vmax)
        
        # Update navigation state
        self.update_navigation_state()
        self.update_status_bar_info()
    
    def update_navigation_state(self):
        """Update navigation button enabled states."""
        filtered_count = len(self.filtered_items)
        current_combo_index = self.system_combo.currentIndex()
        
        # Enable/disable buttons based on position in filtered list
        self.first_btn.setEnabled(filtered_count > 0 and current_combo_index > 0)
        self.prev_btn.setEnabled(filtered_count > 0 and current_combo_index > 0)
        self.next_btn.setEnabled(filtered_count > 0 and current_combo_index < filtered_count - 1)
        self.last_btn.setEnabled(filtered_count > 0 and current_combo_index < filtered_count - 1)
    
    def navigate_to_first(self):
        """Navigate to first system in current filter."""
        if len(self.filtered_items) > 0:
            self.system_combo.setCurrentIndex(0)
    
    def navigate_to_previous(self):
        """Navigate to previous system in current filter."""
        current_index = self.system_combo.currentIndex()
        if current_index > 0:
            self.system_combo.setCurrentIndex(current_index - 1)
    
    def navigate_to_next(self):
        """Navigate to next system in current filter."""
        current_index = self.system_combo.currentIndex()
        if current_index < len(self.filtered_items) - 1:
            self.system_combo.setCurrentIndex(current_index + 1)
    
    def navigate_to_last(self):
        """Navigate to last system in current filter."""
        if len(self.filtered_items) > 0:
            self.system_combo.setCurrentIndex(len(self.filtered_items) - 1)
    
    def open_batch_selection_dialog(self):
        """Open the batch selection dialog for choosing systems to process."""
        from rbcodes.GUIs.specgui.batch.panels.batch_selection_dialog import show_batch_selection_dialog
        
        # Open dialog and get selection
        selected_indices = show_batch_selection_dialog(self.controller, parent=self)
        
        if selected_indices:
            # Emit signal with selected indices
            self.items_selected.emit(selected_indices)
            
            # Switch to processing tab
            parent = self.parentWidget()
            if parent and hasattr(parent, 'parent') and hasattr(parent.parent(), 'tabs'):
                tabs = parent.parent().tabs
                for i in range(tabs.count()):
                    if "Processing" in tabs.tabText(i):
                        tabs.setCurrentIndex(i)
                        break
        else:
            # User cancelled or selected nothing
            self.controller.status_updated.emit("No systems selected for processing.")
    
    def display_item_spectrum(self, item):
        """Get existing spectrum from master table and copy it for editing."""
        try:
            # Skip if import is in progress to avoid warnings
            if self.importing_in_progress:
                return
            # Get existing rb_spec object from master table
            master_spec = self.controller.master_table.get_rb_spec_object(self.current_index)
            
            if master_spec:
                # Deep copy for safe editing
                self.current_spec = copy.deepcopy(master_spec)
                self.update_spectrum_plot(item, self.current_spec)
                self.update_details(item, self.current_spec)
            else:
                self.clear_preview()
                print(f"No rb_spec object found for {item.template.transition_name}")
        except Exception as e:
            self.clear_preview()
            print(f"Error getting spectrum: {e}")

    
    def update_spectrum_plot(self, item, spec):
        """Updated plotting method using the new plotting module."""
        try:
            # Use the new plotting module to generate the complete figure
            plot_spectrum_overview(spec, item, figure=self.spectrum_canvas.fig, clear_figure=True)
            
            # Apply custom axis limits if set
            if self.custom_xlim or self.custom_ylim_upper or self.custom_ylim_lower:
                self.apply_custom_limits()
            
            # Redraw the canvas
            self.spectrum_canvas.draw()
            
        except Exception as e:
            print(f"Error updating spectrum plot: {e}")
            self.clear_preview()

    
    def update_details(self, item, spec):
        """Update details display with comprehensive system information."""
        # Build details string with all key information
        details_parts = []
        
        # 1. Filename (bold)
        filename = os.path.basename(item.template.filename)
        details_parts.append(f"<b>{filename}</b>")
        
        # 2. Transition info
        transition_name = item.template.transition_name
        transition_wave = item.template.transition
        details_parts.append(f"{transition_name} ({transition_wave:.2f} Å)")
        
        # 3. Redshift
        redshift = item.template.redshift
        details_parts.append(f"z={redshift:.6f}")
        
        # 4. EW Range
        ew_vmin = item.template.ew_vmin
        ew_vmax = item.template.ew_vmax
        details_parts.append(f"EW Range: [{ew_vmin}, {ew_vmax}] km/s")
        
        # 5. Status with color coding
        status = item.analysis.processing_status.replace('_', ' ').title()
        if status == 'Complete':
            status_html = f'<span style="color: green; font-weight: bold;">Status: {status}</span>'
        elif status == 'Error':
            status_html = f'<span style="color: red; font-weight: bold;">Status: {status}</span>'
        elif 'Processing' in status:
            status_html = f'<span style="color: orange; font-weight: bold;">Status: {status}</span>'
        else:
            status_html = f'Status: {status}'
        details_parts.append(status_html)
        
        # 6. Masks count
        if hasattr(item.analysis, 'continuum_masks') and item.analysis.continuum_masks:
            mask_count = len(item.analysis.continuum_masks)
            details_parts.append(f"Masks: {mask_count}")
        
        # 7. Additional useful info if available
        if hasattr(item.results, 'W') and not np.isnan(item.results.W):  # Fixed condition
            ew_value = item.results.W
            ew_error = item.results.W_e
            details_parts.append(f"EW: {ew_value:.3f}±{ew_error:.3f} Å")
            
            # Handle negative N case for display
            logn_value = item.results.logN
            logn_error = item.results.logN_e
            if hasattr(item.results, 'N') and item.results.N < 0 and item.results.N_e > 0:
                logn_value = 0.0
                logn_error = np.log10(item.results.N_e)
            
            details_parts.append(f"log N: {logn_value:.2f}±{logn_error:.2f}")
        
        # Join all parts with separators
        details = " | ".join(details_parts)
        self.details_label.setText(details)
    
    def clear_preview(self):
        """Clear preview area."""
        if hasattr(self, 'spectrum_canvas'):
            self.spectrum_canvas.fig.clear()
            ax = self.spectrum_canvas.fig.add_subplot(111)
            ax.set_title("No system selected")
            ax.text(0.5, 0.5, 'Select a system from the dropdown above', 
                   transform=ax.transAxes,
                   ha='center', va='center', fontsize=12, 
                   bbox=dict(facecolor='lightgray', alpha=0.8))
            self.spectrum_canvas.draw()
        
        if hasattr(self, 'details_label'):
            self.details_label.setText("Select a system to view details.")
    
    def update_ew_range(self):
        """Update EW range using the same logic as the existing EW range editor."""
        if not self.current_item or self.current_index < 0:
            self.controller.status_updated.emit("Please select a system first.")
            return
        
        try:
            # Get new values from spinboxes
            new_ew_vmin = self.ew_vmin_spin.value()
            new_ew_vmax = self.ew_vmax_spin.value()
            
            # Validate range
            if new_ew_vmin >= new_ew_vmax:
                QMessageBox.warning(self, "Invalid Range", "EW min must be less than EW max.")
                return
            
            # Update template parameters in master table
            self.controller.master_table.update_template(self.current_index, 
                                                       ew_vmin=new_ew_vmin, 
                                                       ew_vmax=new_ew_vmax)
            
            # Get the updated spectrum and recalculate EW
            master_spec = self.controller.master_table.get_rb_spec_object(self.current_index)
            if master_spec:
                updated_spec = copy.deepcopy(master_spec)
                
                # Recalculate EW with new range
                updated_spec.compute_EW(
                    self.current_item.template.transition,
                    vmin=new_ew_vmin,
                    vmax=new_ew_vmax,
                    SNR=self.current_item.analysis.calculate_snr,
                    _binsize=self.current_item.analysis.binsize
                )
                
                # Update results in master table
                self.controller.master_table.update_results(
                    self.current_index,
                    W=updated_spec.W,
                    W_e=updated_spec.W_e,
                    N=updated_spec.N,
                    N_e=updated_spec.N_e,
                    logN=updated_spec.logN,
                    logN_e=updated_spec.logN_e,
                    vel_centroid=getattr(updated_spec, 'vel_centroid', 0),
                    vel_disp=getattr(updated_spec, 'vel_disp', 0),
                    SNR=getattr(updated_spec, 'SNR', 0),
                    calculation_timestamp=datetime.now().isoformat()
                )
                
                # Mark as complete
                self.controller.master_table.update_analysis(self.current_index, processing_status="complete")
                
                # Store the updated rb_spec object
                self.controller.master_table.set_rb_spec_object(self.current_index, updated_spec)
                
                # Refresh display
                self.current_spec = updated_spec
                updated_item = self.controller.master_table.get_item(self.current_index)
                if updated_item:
                    self.current_item = updated_item
                    self.update_spectrum_plot(self.current_item, self.current_spec)
                    self.update_details(self.current_item, self.current_spec)
                
                self.controller.status_updated.emit("EW range updated successfully!")
                
        except Exception as e:
            self.controller.status_updated.emit(f"Error updating EW range: {str(e)}")
    
    def edit_continuum(self):
        """Launch continuum editor - using original working logic."""
        if not self.current_item or not self.current_spec or self.current_index < 0:
            self.controller.status_updated.emit("Please select a system first.")
            return
        
        try:
            from rbcodes.GUIs.interactive_continuum_fit import launch_interactive_continuum_fit_dialog
            
            # Prepare masks - ORIGINAL WORKING LOGIC
            existing_masks = []
            for mask_tuple in self.current_item.analysis.continuum_masks:
                if isinstance(mask_tuple, (list, tuple)) and len(mask_tuple) == 2:
                    existing_masks.extend([mask_tuple[0], mask_tuple[1]])
            
            # Launch dialog
            input_params = {
                'wave': self.current_spec.wave_slice,
                'flux': self.current_spec.flux_slice,
                'error': self.current_spec.error_slice,
                'velocity': self.current_spec.velo,
                'existing_masks': existing_masks,
                'order': self.current_item.analysis.continuum_order,
                'use_weights': self.current_item.analysis.use_weights,
                'domain': [min(self.current_spec.velo), max(self.current_spec.velo)]
            }
            
            result = launch_interactive_continuum_fit_dialog(**input_params)
            
            if result is None or result.get('cancelled', True):
                self.controller.status_updated.emit("Continuum fitting cancelled.")
                return
            
            # Process results - ORIGINAL WORKING LOGIC
            if 'masks' in result:
                masks = result['masks']
                mask_tuples = []
                for i in range(0, len(masks), 2):
                    if i + 1 < len(masks):
                        mask_tuples.append((masks[i], masks[i+1]))
                
                working_spec = self.current_spec
                working_item = self.current_item
                
                # Apply continuum
                if result.get('fit_method') == 'polynomial':
                    mask = []
                    for vmin, vmax in mask_tuples:
                        mask.extend([vmin, vmax])
                    
                    best_order = result.get('best_order', working_item.analysis.continuum_order)
                    working_spec.fit_continuum(
                        mask=mask if mask else False,
                        Legendre=best_order,
                        use_weights=working_item.analysis.use_weights
                    )
                    
                    continuum_method = "polynomial"
                    continuum_order = best_order

                elif result.get('fit_method') == 'spline':
                    cont = result.get('continuum')
                    working_spec.fit_continuum(prefit_cont=cont)
                    continuum_method = "polynomial"
                    continuum_order = result.get('spline_order', working_item.analysis.continuum_order)
                
                # Recalculate EW
                if hasattr(working_item, 'results'):
                    working_spec.compute_EW(
                        working_item.template.transition,
                        vmin=working_item.template.ew_vmin,
                        vmax=working_item.template.ew_vmax,
                        SNR=working_item.analysis.calculate_snr,
                        _binsize=working_item.analysis.binsize
                    )
                
                # Update master table
                self.controller.master_table.update_analysis(
                    self.current_index,
                    continuum_method=continuum_method,
                    continuum_order=continuum_order,
                    continuum_masks=mask_tuples,
                    processing_status="complete",
                    last_modified=datetime.now().isoformat()
                )
                
                if hasattr(working_item, 'results') and working_item.results.W > 0:
                    self.controller.master_table.update_results(
                        self.current_index,
                        W=working_spec.W,
                        W_e=working_spec.W_e,
                        N=working_spec.N,
                        N_e=working_spec.N_e,
                        logN=working_spec.logN,
                        logN_e=working_spec.logN_e,
                        vel_centroid=getattr(working_spec, 'vel_centroid', 0),
                        vel_disp=getattr(working_spec, 'vel_disp', 0),
                        SNR=getattr(working_spec, 'SNR', 0)
                    )
                
                # Store updated rb_spec object
                self.controller.master_table.set_rb_spec_object(self.current_index, working_spec)
                
                # ORIGINAL WORKING REFRESH LOGIC
                master_spec = self.controller.master_table.get_rb_spec_object(self.current_index)
                if master_spec:
                    self.current_spec = copy.deepcopy(master_spec)
                
                self.refresh_results_table()
                
                updated_item = self.controller.master_table.get_item(self.current_index)
                if updated_item:
                    self.current_item = updated_item
                    self.update_spectrum_plot(self.current_item, self.current_spec)
                    self.update_details(self.current_item, self.current_spec)
                    # Update EW spinboxes with potentially new values
                    self.ew_vmin_spin.setValue(updated_item.template.ew_vmin)
                    self.ew_vmax_spin.setValue(updated_item.template.ew_vmax)
                
                self.controller.status_updated.emit("Continuum updated successfully!")
                
        except Exception as e:
            self.controller.status_updated.emit(f"Error editing continuum: {str(e)}")