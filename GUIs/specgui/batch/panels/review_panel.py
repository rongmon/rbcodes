# rbcodes/GUIs/specgui/batch/panels/review_panel.py - PHASE 2: NAVIGATION CONTROLS
import os
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QComboBox, QGroupBox, QSplitter,
                           QMessageBox, QFileDialog, QDialog, QDialogButtonBox, 
                           QFormLayout, QSpinBox)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
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

class ReviewPanel(QWidget):
    """Panel for reviewing batch processing results - with navigation controls."""
    
    items_selected = pyqtSignal(list)
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.current_item = None
        self.current_spec = None
        self.current_index = -1  # Track current system index
        self.init_ui()
        
        self.controller.table_changed.connect(self.refresh_results_table)
    
    def init_ui(self):
        """Initialize UI with navigation controls."""
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
        
        # Selection actions
        self.select_for_processing_btn = QPushButton("Select for Processing")
        self.select_for_processing_btn.clicked.connect(self.open_batch_selection_dialog)
        self.select_for_processing_btn.setEnabled(True)  # Always enabled if we have systems
        
        button_layout.addWidget(self.select_for_processing_btn)
        
        nav_layout.addLayout(button_layout)
        nav_group.setLayout(nav_layout)
        main_layout.addWidget(nav_group)
        
        # Preview area (large)
        preview_group = QGroupBox("System Preview")
        preview_layout = QVBoxLayout()
        
        # Large canvas
        self.spectrum_canvas = MatplotlibCanvas(self, width=10, height=8, dpi=100)
        self.spectrum_toolbar = NavigationToolbar2QT(self.spectrum_canvas, self)
        
        preview_layout.addWidget(self.spectrum_toolbar)
        preview_layout.addWidget(self.spectrum_canvas)
        
        # Details display
        self.details_label = QLabel("Select a system to view details.")
        self.details_label.setWordWrap(True)
        self.details_label.setMaximumHeight(60)
        self.details_label.setStyleSheet("QLabel { font-size: 9pt; padding: 5px; background-color: #f0f0f0; border: 1px solid #ccc; }")
        preview_layout.addWidget(self.details_label)
        
        # Edit buttons
        edit_buttons = QHBoxLayout()
        
        self.edit_continuum_btn = QPushButton("Edit Continuum")
        self.edit_continuum_btn.clicked.connect(self.edit_continuum)
        self.edit_continuum_btn.setEnabled(False)
        
        self.edit_ew_btn = QPushButton("Edit EW Range")
        self.edit_ew_btn.clicked.connect(self.edit_ew_range)
        self.edit_ew_btn.setEnabled(False)
        
        edit_buttons.addWidget(self.edit_continuum_btn)
        edit_buttons.addWidget(self.edit_ew_btn)
        edit_buttons.addStretch()
        preview_layout.addLayout(edit_buttons)
        
        preview_group.setLayout(preview_layout)
        main_layout.addWidget(preview_group)
    
    def refresh_results_table(self):
        """
        Populate system dropdown from master table.
        
        Note: Despite the name 'refresh_results_table', this method now populates
        the navigation dropdown instead of a table (Phase 2 enhancement). The name
        is preserved for interface compatibility with existing code.
        """
        # Preserve current selection to avoid jumping back to first item
        current_selection = self.current_index
        # Clear current state
        self.system_combo.clear()
        self.current_item = None
        self.current_spec = None
        self.current_index = -1
        self.update_navigation_state()
        self.clear_preview()
        
        # Get all items from master table
        items = self.controller.master_table.get_all_items()
        
        if not items:
            self.position_label.setText("No systems loaded")
            return
        
        # Populate dropdown
        for i, item in enumerate(items):
            # Create descriptive label
            filename = os.path.basename(item.template.filename)
            transition_name = item.template.transition_name
            redshift = item.template.redshift
            status = item.analysis.processing_status
            
            # Status indicator
            if status == 'complete':
                status_icon = "✓"
            elif status == 'error':
                status_icon = "✗"
            elif 'processing' in status:
                status_icon = "⏳"
            else:
                status_icon = "⚠️"
            
            label = f"{status_icon} System {i+1}: {filename} | {transition_name} | z={redshift:.4f}"
            self.system_combo.addItem(label)
        
        # Update position indicator
        total_count = len(items)
        self.position_label.setText(f"Total: {total_count} systems")
        
        # Enable navigation and selection if we have items
        if total_count > 0:
            self.update_navigation_state()
            self.select_for_processing_btn.setEnabled(True)
            
            # Restore previous selection if valid, otherwise select first
            if current_selection >= 0 and current_selection < total_count:
                self.system_combo.setCurrentIndex(current_selection)
            else:
                self.system_combo.setCurrentIndex(0)
        else:
            self.select_for_processing_btn.setEnabled(False)
    
    def on_system_selected(self, index):
        """Handle system selection from dropdown."""
        if index < 0:
            return
            
        item_count = self.controller.master_table.get_item_count()
        if index >= item_count:
            return
            
        # Update current tracking
        self.current_index = index
        
        # Get selected item
        selected_item = self.controller.master_table.get_item(index)
        if selected_item:
            self.current_item = selected_item
            self.display_item_spectrum(selected_item)
            self.edit_continuum_btn.setEnabled(True)
            self.edit_ew_btn.setEnabled(True)
            self.select_for_processing_btn.setEnabled(True)
        
        # Update navigation state
        self.update_navigation_state()
        
        # Update position indicator
        total_count = self.controller.master_table.get_item_count()
        self.position_label.setText(f"System {index + 1} of {total_count}")
    
    def update_navigation_state(self):
        """Update navigation button enabled states."""
        total_count = self.controller.master_table.get_item_count()
        current_index = self.current_index
        
        # Enable/disable buttons based on position
        self.first_btn.setEnabled(total_count > 0 and current_index > 0)
        self.prev_btn.setEnabled(total_count > 0 and current_index > 0)
        self.next_btn.setEnabled(total_count > 0 and current_index < total_count - 1)
        self.last_btn.setEnabled(total_count > 0 and current_index < total_count - 1)
    
    def navigate_to_first(self):
        """Navigate to first system."""
        if self.controller.master_table.get_item_count() > 0:
            self.system_combo.setCurrentIndex(0)
    
    def navigate_to_previous(self):
        """Navigate to previous system."""
        if self.current_index > 0:
            self.system_combo.setCurrentIndex(self.current_index - 1)
    
    def navigate_to_next(self):
        """Navigate to next system."""
        total_count = self.controller.master_table.get_item_count()
        if self.current_index < total_count - 1:
            self.system_combo.setCurrentIndex(self.current_index + 1)
    
    def navigate_to_last(self):
        """Navigate to last system."""
        total_count = self.controller.master_table.get_item_count()
        if total_count > 0:
            self.system_combo.setCurrentIndex(total_count - 1)
    
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
        """
        Updated plotting method using the new plotting module.
        This completely redraws the figure from the rb_spec object.
        """
        try:
            # Use the new plotting module to generate the complete figure
            plot_spectrum_overview(spec, item, figure=self.spectrum_canvas.fig, clear_figure=True)
            
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
        if hasattr(item.results, 'W') and item.results.W > 0:
            ew_value = item.results.W
            ew_error = item.results.W_e
            details_parts.append(f"EW: {ew_value:.3f}±{ew_error:.3f} Å")
            
            logn_value = item.results.logN
            logn_error = item.results.logN_e
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
                
                self.controller.status_updated.emit("Continuum updated successfully!")
                
        except Exception as e:
            self.controller.status_updated.emit(f"Error editing continuum: {str(e)}")
    
    def edit_ew_range(self):
        """Launch EW range editor - using original working logic."""
        if not self.current_item or self.current_index < 0:
            self.controller.status_updated.emit("Please select a system first.")
            return
        
        try:
            from rbcodes.GUIs.specgui.batch.panels.ew_range_editor import edit_ew_range_dialog
            
            editor_item = type('EditorItem', (), {
                'template': self.current_item.template,
                'analysis': self.current_item.analysis,
                'results': self.current_item.results,
                'row_index': self.current_index
            })()
            
            success = edit_ew_range_dialog(
                editor_item, 
                self.current_spec, 
                self.controller,
                parent=self
            )
            
            if success:
                # ORIGINAL WORKING REFRESH LOGIC
                updated_item = self.controller.master_table.get_item(self.current_index)
                
                if updated_item:
                    master_spec = self.controller.master_table.get_rb_spec_object(self.current_index)
                    if master_spec:
                        updated_spec = copy.deepcopy(master_spec)
                        self.current_spec = updated_spec
                    else:
                        self.controller.status_updated.emit("Failed to get rb_spec object from master table.")
                        return
                    
                    self.current_item = updated_item
                    self.refresh_results_table()
                    self.update_spectrum_plot(self.current_item, self.current_spec)
                    self.update_details(self.current_item, self.current_spec)
                    self.controller.status_updated.emit("EW range updated successfully!")
                else:
                    self.controller.status_updated.emit("Failed to get updated item from master table.")
            
        except ImportError as e:
            self.controller.status_updated.emit(f"EW range editor not available: {str(e)}")
        except Exception as e:
            self.controller.status_updated.emit(f"Error editing EW range: {str(e)}")