# rbcodes/GUIs/specgui/batch/panels/review_panel.py - SIMPLIFIED CLEAN VERSION
import os
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QTableWidget, QTableWidgetItem,
                           QHeaderView, QCheckBox, QGroupBox, QSplitter,
                           QTabWidget, QScrollArea, QMessageBox, QFileDialog,
                           QDialog, QDialogButtonBox, QFormLayout, QSpinBox)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
from datetime import datetime    
from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object
import copy
class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots."""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)

class ReviewPanel(QWidget):
    """Panel for reviewing batch processing results - simplified approach."""
    
    items_selected = pyqtSignal(list)
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.current_item = None
        self.current_spec = None
        self.init_ui()
        
        self.controller.table_changed.connect(self.refresh_results_table)
    
    def init_ui(self):
        """Initialize UI with larger preview (30:70 ratio)."""
        main_layout = QVBoxLayout(self)
        splitter = QSplitter(Qt.Vertical)
        
        # Results table (30%)
        results_group = QGroupBox("Batch Results")
        results_layout = QVBoxLayout()
        
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(9)
        self.results_table.setHorizontalHeaderLabels([
            "", "Filename", "Redshift", "Transition", "EW (Å)", "logN", "SNR", "Status", "Actions"
        ])
        
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Fixed)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        for i in range(2, 9):
            header.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        
        self.results_table.setRowCount(0)
        self.results_table.itemClicked.connect(self.on_table_item_clicked)
        
        # Selection buttons
        selection_buttons = QHBoxLayout()
        
        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.clicked.connect(self.select_all_items)
        
        self.select_none_btn = QPushButton("Select None")
        self.select_none_btn.clicked.connect(self.select_no_items)
        
        self.select_failed_btn = QPushButton("Select Failed Items")
        self.select_failed_btn.clicked.connect(self.select_failed_items)
        
        self.process_selected_btn = QPushButton("Process Selected")
        self.process_selected_btn.clicked.connect(self.process_selected_items)
        
        selection_buttons.addWidget(self.select_all_btn)
        selection_buttons.addWidget(self.select_none_btn)
        selection_buttons.addWidget(self.select_failed_btn)
        selection_buttons.addWidget(self.process_selected_btn)
        
        results_layout.addWidget(self.results_table)
        results_layout.addLayout(selection_buttons)
        results_group.setLayout(results_layout)
        splitter.addWidget(results_group)
        
        # Preview area (70%)
        preview_group = QGroupBox("Item Preview")
        preview_layout = QVBoxLayout()
        
        # Large canvas
        self.spectrum_canvas = MatplotlibCanvas(self, width=8, height=6, dpi=100)
        self.spectrum_toolbar = NavigationToolbar2QT(self.spectrum_canvas, self)
        
        preview_layout.addWidget(self.spectrum_toolbar)
        preview_layout.addWidget(self.spectrum_canvas)
        
        # Compact details
        self.details_label = QLabel("Select an item to view details.")
        self.details_label.setWordWrap(True)
        self.details_label.setMaximumHeight(80)
        self.details_label.setStyleSheet("QLabel { font-size: 9pt; }")
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
        preview_layout.addLayout(edit_buttons)
        
        preview_group.setLayout(preview_layout)
        splitter.addWidget(preview_group)
        
        # 30:70 ratio
        splitter.setSizes([250, 750])
        main_layout.addWidget(splitter)
    
    def refresh_results_table(self):
        """Repopulate entire table from master table."""
        self.results_table.setRowCount(0)
        
        # Clear selection
        self.current_item = None
        self.current_spec = None
        self.edit_continuum_btn.setEnabled(False)
        self.edit_ew_btn.setEnabled(False)
        self.clear_preview()
        
        items = self.controller.master_table.get_all_items()
        
        for i, item in enumerate(items):
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)
            
            # Checkbox
            checkbox_item = QTableWidgetItem()
            checkbox_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            checkbox_item.setCheckState(Qt.Unchecked)
            self.results_table.setItem(row, 0, checkbox_item)
            
            # Data
            filename = os.path.basename(item.template.filename)
            self.results_table.setItem(row, 1, QTableWidgetItem(filename))
            self.results_table.setItem(row, 2, QTableWidgetItem(f"{item.template.redshift:.6f}"))
            
            transition_name = item.template.transition_name
            self.results_table.setItem(row, 3, QTableWidgetItem(f"{transition_name} ({item.template.transition:.2f} Å)"))
            
            # Results
            ew_text = f"{item.results.W:.3f} ± {item.results.W_e:.3f}"
            logn_text = f"{item.results.logN:.2f} ± {item.results.logN_e:.2f}"
            snr_text = f"{item.results.SNR:.1f}" if item.results.SNR > 0 else "0.0"
            
            self.results_table.setItem(row, 4, QTableWidgetItem(ew_text))
            self.results_table.setItem(row, 5, QTableWidgetItem(logn_text))
            self.results_table.setItem(row, 6, QTableWidgetItem(snr_text))
            
            # Status with color
            status = item.analysis.processing_status
            status_item = QTableWidgetItem(status.replace('_', ' ').title())
            
            if status == 'complete':
                status_item.setBackground(QColor(200, 255, 200))
            elif status == 'processing':
                status_item.setBackground(QColor(255, 255, 200))
            elif status == 'error':
                status_item.setBackground(QColor(255, 200, 200))
            elif status in ['needs_processing', 'needs_ew_recalc']:
                status_item.setBackground(QColor(255, 240, 200))
            
            self.results_table.setItem(row, 7, status_item)
            self.results_table.setItem(row, 8, QTableWidgetItem(""))    
    
    def on_table_item_clicked(self, item):
        """Handle table clicks."""
        row = item.row()
        item_count = self.controller.master_table.get_item_count()
        
        if row < item_count:
            selected_item = self.controller.master_table.get_item(row)
            if selected_item:
                # Store the row index for later use
                self.current_row_index = row
                self.display_item_spectrum(selected_item)
                self.edit_continuum_btn.setEnabled(True)
                self.edit_ew_btn.setEnabled(True)
                self.current_item = selected_item
    
    def display_item_spectrum(self, item):
        """Get existing spectrum from master table and copy it for editing."""
        try:
            # Get existing rb_spec object from master table
            master_spec = self.controller.master_table.get_rb_spec_object(self.current_row_index)
            
            if master_spec:
                # Deep copy for safe editing
                import copy
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
        """Update spectrum plot."""
        self.spectrum_canvas.axes.clear()
        
        # Plot data
        self.spectrum_canvas.axes.step(spec.velo, spec.fnorm, 'k-', where='mid', 
                                     linewidth=0.50, label='Normalized Flux')
        
        if hasattr(spec, 'enorm'):
            self.spectrum_canvas.axes.step(spec.velo, spec.enorm, 'r-', where='mid', 
                                         alpha=0.6, linewidth=0.5, label='Error')
        
        # Reference lines
        self.spectrum_canvas.axes.axhline(y=1.0, color='g', linestyle='--', alpha=0.7, linewidth=0.5)
        
        # EW region
        ew_vmin = item.template.ew_vmin
        ew_vmax = item.template.ew_vmax
        
        self.spectrum_canvas.axes.axvspan(ew_vmin, ew_vmax, alpha=0.15, color='blue', label='EW Region')
        self.spectrum_canvas.axes.axvline(x=ew_vmin, color='b', linestyle=':', alpha=0.8, linewidth=0.5)
        self.spectrum_canvas.axes.axvline(x=ew_vmax, color='b', linestyle=':', alpha=0.8, linewidth=0.5)
        
        # Continuum masks
        mask_count = 0
        for vmin, vmax in item.analysis.continuum_masks:
            label = 'Masked' if mask_count == 0 else ""
            self.spectrum_canvas.axes.axvspan(vmin, vmax, alpha=0.2, color='red', label=label)
            mask_count += 1
        
        # Limits
        x_buffer = (max(spec.velo) - min(spec.velo)) * 0.05
        self.spectrum_canvas.axes.set_xlim(min(spec.velo) - x_buffer, max(spec.velo) + x_buffer)
        
        y_data = spec.fnorm
        y_min = max(-0.1, min(y_data) - 0.1)
        y_max = min(2.0, max(y_data) + 0.3)
        self.spectrum_canvas.axes.set_ylim(y_min, y_max)
        
        # Labels
        transition_name = item.template.transition_name
        transition = item.template.transition
        self.spectrum_canvas.axes.set_title(f"{transition_name} ({transition:.2f} Å)", fontsize=12, fontweight='bold')
        self.spectrum_canvas.axes.set_xlabel('Velocity (km/s)', fontsize=11)
        self.spectrum_canvas.axes.set_ylabel('Normalized Flux', fontsize=11)
        
        # Results text
        ew = item.results.W
        ew_err = item.results.W_e
        logn = item.results.logN
        logn_err = item.results.logN_e
        
        info_text = f"EW = {ew:.3f} ± {ew_err:.3f} Å\nlog N = {logn:.2f} ± {logn_err:.2f}\nSNR = {item.results.SNR:.1f}"
        
        self.spectrum_canvas.axes.text(0.98, 0.05, info_text, 
                                     transform=self.spectrum_canvas.axes.transAxes, 
                                     horizontalalignment='right', verticalalignment='bottom',
                                     bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'),
                                     fontsize=10)
        
        self.spectrum_canvas.axes.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        self.spectrum_canvas.fig.tight_layout()
        self.spectrum_canvas.draw()
    
    def update_details(self, item, spec):
        """Update details display."""
        details = f"<b>{os.path.basename(item.template.filename)}</b> | "
        details += f"z={item.template.redshift:.4f} | "
        details += f"EW Range: [{item.template.ew_vmin}, {item.template.ew_vmax}] km/s | "
        details += f"Status: {item.analysis.processing_status.replace('_', ' ').title()}"
        
        if item.analysis.continuum_masks:
            mask_count = len(item.analysis.continuum_masks)
            details += f" | Masks: {mask_count}"
        
        self.details_label.setText(details)
    
    def clear_preview(self):
        """Clear preview area."""
        if hasattr(self, 'spectrum_canvas'):
            self.spectrum_canvas.axes.clear()
            self.spectrum_canvas.axes.set_title("No item selected")
            self.spectrum_canvas.axes.text(0.5, 0.5, 'Select an item from the table above', 
                                          transform=self.spectrum_canvas.axes.transAxes,
                                          ha='center', va='center', fontsize=12, 
                                          bbox=dict(facecolor='lightgray', alpha=0.8))
            self.spectrum_canvas.draw()
        
        if hasattr(self, 'details_label'):
            self.details_label.setText("Select an item to view details.")
    
    def edit_continuum(self):
        """Launch continuum editor - SIMPLE APPROACH."""
        if not self.current_item or not self.current_spec or not hasattr(self, 'current_row_index'):
            self.controller.status_updated.emit("Please select an item first.")
            return
        
        try:
            from rbcodes.GUIs.interactive_continuum_fit import launch_interactive_continuum_fit_dialog
            
            # Prepare masks
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
                        
            # Process masks
            if 'masks' in result:
                masks = result['masks']
                mask_tuples = []
                for i in range(0, len(masks), 2):
                    if i + 1 < len(masks):
                        mask_tuples.append((masks[i], masks[i+1]))
                
                # Save variables before any master table operations
                working_spec = self.current_spec
                working_item = self.current_item
                
                # Apply continuum to working spectrum
                if result.get('fit_method') == 'polynomial':
                    # Use polynomial fitting with updated order
                    mask = []
                    for vmin, vmax in mask_tuples:
                        mask.extend([vmin, vmax])
                    
                    best_order = result.get('best_order', working_item.analysis.continuum_order)
                    working_spec.fit_continuum(
                        mask=mask if mask else False,
                        Legendre=best_order,
                        use_weights=working_item.analysis.use_weights
                    )
                    
                    # Store method info for later master table update
                    continuum_method = "polynomial"
                    continuum_order = best_order

                
                
                elif result.get('fit_method') == 'spline':
                    # Use pre-fitted spline continuum
                    cont = result.get('continuum')
                    working_spec.fit_continuum(prefit_cont=cont)
                    
                    # Store method info for later master table update
                    continuum_method = "polynomial"
                    continuum_order = result.get('spline_order', working_item.analysis.continuum_order)  # Indicator for spline
                
                # Recalculate EW on working spectrum (if we have existing results)
                if hasattr(working_item, 'results'):
                    working_spec.compute_EW(
                        working_item.template.transition,
                        vmin=working_item.template.ew_vmin,
                        vmax=working_item.template.ew_vmax,
                        SNR=working_item.analysis.calculate_snr,
                        _binsize=working_item.analysis.binsize
                    )
                
                # Now update master table with all the results (do this AFTER all rb_spec work is done)
                
                # Update method and order
                self.controller.master_table.update_analysis(
                    self.current_row_index,
                    continuum_method=continuum_method,
                    continuum_order=continuum_order
                )
                
                # Update results (if EW was calculated)
                if hasattr(working_item, 'results') and working_item.results.W > 0:
                    self.controller.master_table.update_results(
                        self.current_row_index,
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
                
                # Update masks and status
                self.controller.master_table.update_analysis(
                    self.current_row_index,
                    continuum_masks=mask_tuples,
                    processing_status="complete",
                    last_modified=datetime.now().isoformat()
                )
                
                # Store the working rb_spec object in master table
                self.controller.master_table.set_rb_spec_object(self.current_row_index, working_spec)
                
                # Get fresh copy from master table for continued editing
                master_spec = self.controller.master_table.get_rb_spec_object(self.current_row_index)
                if master_spec:
                    self.current_spec = copy.deepcopy(master_spec)
                
                # Refresh UI
                self.refresh_results_table()
                
                # Get fresh item and update preview
                updated_item = self.controller.master_table.get_item(self.current_row_index)
                if updated_item:
                    self.current_item = updated_item
                    self.update_spectrum_plot(self.current_item, self.current_spec)
                    self.update_details(self.current_item, self.current_spec)
                
                self.controller.status_updated.emit("Continuum updated successfully!")
        except Exception as e:
            self.controller.status_updated.emit(f"Error editing continuum: {str(e)}")    
    
    def edit_ew_range(self):
        """Launch EW range editor - SIMPLE APPROACH."""
        if not self.current_item or not hasattr(self, 'current_row_index'):
            self.controller.status_updated.emit("Please select an item first.")
            return
        
        try:
            from rbcodes.GUIs.specgui.batch.panels.ew_range_editor import edit_ew_range_dialog
            
            # Create a simple object that mimics the old BatchItem structure for the editor
            editor_item = type('EditorItem', (), {
                'template': self.current_item.template,
                'analysis': self.current_item.analysis,
                'results': self.current_item.results,
                'row_index': self.current_row_index  # Pass row index for updates
            })()
            
            success = edit_ew_range_dialog(
                editor_item, 
                self.current_spec, 
                self.controller,
                parent=self
            )
            
            if success:
                # Get fresh item from master table using row index
                updated_item = self.controller.master_table.get_item(self.current_row_index)
                
                if updated_item:
                    
                    master_spec = self.controller.master_table.get_rb_spec_object(self.current_row_index)
                    if master_spec:
                        updated_spec = copy.deepcopy(master_spec)
                        self.current_spec = updated_spec
                    else:
                        updated_spec = None
                        self.controller.status_updated.emit("Failed to get rb_spec object from master table.")
                    
                    if updated_spec:
                        # Update current references to the fresh data
                        self.current_item = updated_item
                        self.current_spec = updated_spec
                        
                        # Repopulate table and update preview with fresh data
                        self.refresh_results_table()
                        self.update_spectrum_plot(self.current_item, self.current_spec)
                        self.update_details(self.current_item, self.current_spec)
                        self.controller.status_updated.emit("EW range updated successfully!")
                    else:
                        self.controller.status_updated.emit("Failed to regenerate spectrum.")
                else:
                    self.controller.status_updated.emit("Failed to get updated item from master table.")
            
        except ImportError as e:
            self.controller.status_updated.emit(f"EW range editor not available: {str(e)}")
        except Exception as e:
            self.controller.status_updated.emit(f"Error editing EW range: {str(e)}")    

    
    def select_all_items(self):
        """Select all items."""
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item:
                checkbox_item.setCheckState(Qt.Checked)
        self.emit_selected_indices()
    
    def select_no_items(self):
        """Deselect all items."""
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item:
                checkbox_item.setCheckState(Qt.Unchecked)
        self.emit_selected_indices()
    
    def select_failed_items(self):
        """Select failed items."""
        for row in range(self.results_table.rowCount()):
            status_item = self.results_table.item(row, 7)
            checkbox_item = self.results_table.item(row, 0)
            
            if status_item and checkbox_item:
                status = status_item.text().lower()
                if any(word in status for word in ['error', 'needs', 'processing']):
                    checkbox_item.setCheckState(Qt.Checked)
                else:
                    checkbox_item.setCheckState(Qt.Unchecked)
        self.emit_selected_indices()
    
    def emit_selected_indices(self):
        """Emit selected row indices."""
        selected_indices = []
        
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                selected_indices.append(row)  # Use row index directly
        
        self.items_selected.emit(selected_indices)
    
    def process_selected_items(self):
        """Process selected items."""
        selected_indices = []
        
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                selected_indices.append(row)  # Use row index directly
        
        if not selected_indices:
            self.controller.status_updated.emit("Please select items to process.")
            return
        
        self.items_selected.emit(selected_indices)
        
        # Switch to processing tab
        parent = self.parentWidget()
        if parent and isinstance(parent, QTabWidget):
            for i in range(parent.count()):
                if "Processing" in parent.tabText(i):
                    parent.setCurrentIndex(i)
                    break    