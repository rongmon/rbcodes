# rbcodes/GUIs/specgui/batch/panels/review_panel.py
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
# Import rb_spec for on-demand generation
from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots."""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)

class ReviewPanel(QWidget):
    """Simplified panel for reviewing batch processing results."""
    
    # Signals
    items_selected = pyqtSignal(list)  # Emitted when items are selected
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.current_item = None  # Currently selected item (direct reference)
        self.current_spec = None  # Currently generated spectrum
        self.init_ui()
        
        # Connect to controller signals
        self.controller.table_changed.connect(self.refresh_results_table)
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add a splitter to divide results table and preview
        splitter = QSplitter(Qt.Vertical)
        
        # Results table
        results_group = QGroupBox("Batch Results")
        results_layout = QVBoxLayout()
        
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(9)
        self.results_table.setHorizontalHeaderLabels([
            "", "Filename", "Redshift", "Transition", "EW (Å)", "logN", "SNR", "Status", "Actions"
        ])
        
        # Set column widths
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Fixed)  # Checkbox column
        header.setSectionResizeMode(1, QHeaderView.Stretch)  # Filename
        for i in range(2, 9):
            header.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        
        self.results_table.setRowCount(0)
        self.results_table.itemClicked.connect(self.on_table_item_clicked)
        
        # Buttons for managing selection
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
        
        # Preview area with tabs for different views
        preview_group = QGroupBox("Item Preview")
        preview_layout = QVBoxLayout()
        
        self.preview_tabs = QTabWidget()
        
        # Tab for spectrum plot
        self.spectrum_tab = QWidget()
        spectrum_layout = QVBoxLayout(self.spectrum_tab)
        
        self.spectrum_canvas = MatplotlibCanvas(self, width=5, height=4, dpi=100)
        self.spectrum_toolbar = NavigationToolbar2QT(self.spectrum_canvas, self)
        
        spectrum_layout.addWidget(self.spectrum_toolbar)
        spectrum_layout.addWidget(self.spectrum_canvas)
        
        # Tab for details
        self.details_tab = QWidget()
        details_layout = QVBoxLayout(self.details_tab)
        
        self.details_label = QLabel("Select an item to view details.")
        self.details_label.setWordWrap(True)
        details_layout.addWidget(self.details_label)
        
        # Add tabs to tab widget
        self.preview_tabs.addTab(self.spectrum_tab, "Spectrum")
        self.preview_tabs.addTab(self.details_tab, "Details")
        
        preview_layout.addWidget(self.preview_tabs)
        
        # Interactive editing buttons
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
        
        # Set initial splitter sizes
        splitter.setSizes([300, 400])  # 40% for results, 60% for preview
        
        main_layout.addWidget(splitter)
    
    def set_results(self, results):
        """Set the batch processing results (for backward compatibility)."""
        # The results are now managed by the master table
        # This method just triggers a refresh
        self.refresh_results_table()
        
        # Clear selection
        self.current_item = None
        self.current_spec = None
        self.edit_continuum_btn.setEnabled(False)
        self.edit_ew_btn.setEnabled(False)
        
        # Clear preview
        self.clear_preview()
    
    def refresh_results_table(self):
        """Update the results table from the master table."""
        # Clear existing rows
        self.results_table.setRowCount(0)
        
        # Clear current selection
        self.current_item = None
        self.current_spec = None
        self.edit_continuum_btn.setEnabled(False)
        self.edit_ew_btn.setEnabled(False)
        self.clear_preview()
        
        # Get all items from master table
        items = self.controller.master_table.get_all_items()
        
        # Add results
        for i, item in enumerate(items):
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)
            
            # Checkbox for selection
            checkbox_item = QTableWidgetItem()
            checkbox_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            checkbox_item.setCheckState(Qt.Unchecked)
            self.results_table.setItem(row, 0, checkbox_item)
            
            # Item data
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
            
            # Status with color coding
            status = item.analysis.processing_status
            status_item = QTableWidgetItem(status.replace('_', ' ').title())
            
            if status == 'complete':
                status_item.setBackground(QColor(200, 255, 200))  # Light green
            elif status == 'processing':
                status_item.setBackground(QColor(255, 255, 200))  # Light yellow
            elif status == 'error':
                status_item.setBackground(QColor(255, 200, 200))  # Light red
            elif status in ['needs_processing', 'needs_ew_recalc']:
                status_item.setBackground(QColor(255, 240, 200))  # Light orange
            
            self.results_table.setItem(row, 7, status_item)
            
            # Actions column - could add buttons here later
            self.results_table.setItem(row, 8, QTableWidgetItem(""))
    
    def on_table_item_clicked(self, item):
        """Handle clicks on table items - SIMPLIFIED VERSION."""
        row = item.row()
        
        # Get items from master table
        items = self.controller.master_table.get_all_items()
        
        if row < len(items):
            selected_item = items[row]  # Use direct indexing, no UUID
            print(f"Selected item: {selected_item.template.transition_name}")
            
            # Generate spectrum on-demand
            self.display_item_spectrum(selected_item)
            
            # Enable editing buttons (always enabled if item has parameters)
            self.edit_continuum_btn.setEnabled(True)
            self.edit_ew_btn.setEnabled(True)
            
            # Store for editing
            self.current_item = selected_item
    
    def display_item_spectrum(self, item):
        """Generate and display spectrum on-demand from parameters."""
        try:
            print(f"Generating spectrum for {item.template.transition_name}")
            
            # Generate rb_spec object fresh
            spec = self._generate_spectrum_from_item(item)
            
            if spec:
                self.current_spec = spec
                # Update displays
                self.update_spectrum_plot_simple(item, spec)
                self.update_details_simple(item, spec)
            else:
                self.clear_preview()
                
        except Exception as e:
            print(f"Error generating spectrum: {e}")
            self.clear_preview()

    def _generate_spectrum_from_item(self, item):
        """Generate rb_spec object from item parameters."""
        try:
            # Load spectrum file
            if item.template.filename.lower().endswith('.json'):
                spec = load_rb_spec_object(item.template.filename)
            else:
                # Determine file type
                ext = os.path.splitext(item.template.filename)[1].lower()
                filetype = 'linetools' if ext == '.fits' else ('ascii' if ext in ['.txt', '.dat'] else None)
                spec = rb_spec.from_file(item.template.filename, filetype=filetype)
            
            # Apply saved parameters
            spec.shift_spec(item.template.redshift)
            spec.slice_spec(
                item.template.transition,
                item.template.slice_vmin, item.template.slice_vmax,
                use_vel=True,
                linelist=item.template.linelist
            )
            
            # Apply continuum (from saved parameters)
            if item.analysis.continuum_method == 'polynomial':
                mask = []
                for vmin, vmax in item.analysis.continuum_masks:
                    mask.extend([vmin, vmax])
                spec.fit_continuum(
                    mask=mask if mask else False,
                    Legendre=item.analysis.continuum_order,
                    use_weights=item.analysis.use_weights
                )
            else:  # flat
                spec.cont = np.ones_like(spec.flux_slice)
                spec.fnorm = spec.flux_slice.copy()
                spec.enorm = spec.error_slice.copy()
            
            # Apply EW calculation results if we have them
            if item.results.W > 0:
                spec.W = item.results.W
                spec.W_e = item.results.W_e
                spec.N = item.results.N
                spec.N_e = item.results.N_e
                spec.logN = item.results.logN
                spec.logN_e = item.results.logN_e
                spec.vel_centroid = item.results.vel_centroid
                spec.vel_disp = item.results.vel_disp
                spec.SNR = item.results.SNR
                spec.vmin = item.template.ew_vmin
                spec.vmax = item.template.ew_vmax
                spec.trans = item.template.transition_name
                spec.trans_wave = item.template.transition
            
            return spec
            
        except Exception as e:
            print(f"Error generating spectrum: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def update_spectrum_plot_simple(self, item, spec):
        """Simple spectrum plotting without complex state checks."""
        # Clear existing plot
        self.spectrum_canvas.axes.clear()
        
        # Plot normalized spectrum
        self.spectrum_canvas.axes.step(spec.velo, spec.fnorm, 'k-', where='mid', label='Normalized Flux')
        
        # Plot normalized error
        if hasattr(spec, 'enorm'):
            self.spectrum_canvas.axes.step(spec.velo, spec.enorm, 'r-', where='mid', alpha=0.5, label='Error')
        
        # Reference line at y=1.0
        self.spectrum_canvas.axes.axhline(y=1.0, color='g', linestyle='--', alpha=0.7)
        
        # Add EW measurement region
        ew_vmin = item.template.ew_vmin
        ew_vmax = item.template.ew_vmax
        
        self.spectrum_canvas.axes.axvspan(ew_vmin, ew_vmax, alpha=0.1, color='blue', label='EW Region')
        self.spectrum_canvas.axes.axvline(x=ew_vmin, color='b', linestyle=':', alpha=0.7)
        self.spectrum_canvas.axes.axvline(x=ew_vmax, color='b', linestyle=':', alpha=0.7)
        
        # Show continuum masks if available
        for vmin, vmax in item.analysis.continuum_masks:
            self.spectrum_canvas.axes.axvspan(vmin, vmax, alpha=0.2, color='red', label='Masked')
        
        # Set axis limits based on the data
        x_buffer = (max(spec.velo) - min(spec.velo)) * 0.1
        self.spectrum_canvas.axes.set_xlim(min(spec.velo) - x_buffer, max(spec.velo) + x_buffer)
        
        # Try to set y limits nicely
        y_data = spec.fnorm
        y_min = max(0, min(y_data) - 0.1)
        y_max = min(2.0, max(y_data) + 0.3)
        self.spectrum_canvas.axes.set_ylim(y_min, y_max)
        
        # Add labels and measurement info
        transition_name = item.template.transition_name
        transition = item.template.transition
        self.spectrum_canvas.axes.set_title(f"{transition_name} ({transition:.2f} Å)")
        self.spectrum_canvas.axes.set_xlabel('Velocity (km/s)')
        self.spectrum_canvas.axes.set_ylabel('Normalized Flux')
        
        # Add EW and column density info
        ew = item.results.W
        ew_err = item.results.W_e
        logn = item.results.logN
        logn_err = item.results.logN_e
        
        
        info_text = f"EW = {ew:.3f} ± {ew_err:.3f} Å\nlog N = {logn:.2f} ± {logn_err:.2f}"
        info_text += f"\nSNR = {item.results.SNR:.1f}"    
        
        self.spectrum_canvas.axes.text(0.98, 0.05, info_text, 
                                     transform=self.spectrum_canvas.axes.transAxes, 
                                     horizontalalignment='right', verticalalignment='bottom',
                                     bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'))
        
        # Apply tight layout
        self.spectrum_canvas.fig.tight_layout()
        
        # Draw the updated plot
        self.spectrum_canvas.draw()
    
    def update_details_simple(self, item, spec):
        """Update the details view with the selected item's data."""
        # Format details text
        details = f"<b>Filename:</b> {item.template.filename}<br>"
        details += f"<b>Redshift:</b> {item.template.redshift:.6f}<br>"
        details += f"<b>Transition:</b> {item.template.transition_name} ({item.template.transition:.2f} Å)<br>"
        details += f"<b>Slice Range:</b> {item.template.slice_vmin} to {item.template.slice_vmax} km/s<br>"
        details += f"<b>EW Range:</b> {item.template.ew_vmin} to {item.template.ew_vmax} km/s<br>"
        details += f"<b>Line List:</b> {item.template.linelist}<br>"
        

        details += f"<b>Equivalent Width:</b> {item.results.W:.3f} ± {item.results.W_e:.3f} Å<br>"
        details += f"<b>Column Density (log):</b> {item.results.logN:.2f} ± {item.results.logN_e:.2f}<br>"
        details += f"<b>SNR:</b> {item.results.SNR:.1f}<br>"
        details += f"<b>Velocity Centroid:</b> {item.results.vel_centroid:.1f} km/s<br>"
        details += f"<b>Velocity Dispersion:</b> {item.results.vel_disp:.1f} km/s<br>"
        
        details += f"<b>Status:</b> {item.analysis.processing_status.replace('_', ' ').title()}<br>"
        
        # Add continuum info
        if item.analysis.continuum_masks:
            mask_text = ", ".join([f"[{vmin}, {vmax}]" for vmin, vmax in item.analysis.continuum_masks])
            details += f"<b>Continuum Masks:</b> {mask_text}<br>"
        
        details += f"<b>Last Modified:</b> {item.analysis.last_modified}<br>"
        
        self.details_label.setText(details)
    
    def clear_preview(self):
        """Clear the preview area."""
        # Clear spectrum plot
        if hasattr(self, 'spectrum_canvas'):
            self.spectrum_canvas.axes.clear()
            self.spectrum_canvas.axes.set_title("No item selected")
            self.spectrum_canvas.draw()
        
        # Clear details
        if hasattr(self, 'details_label'):
            self.details_label.setText("Select an item to view details.")
    
    def edit_continuum(self):
        """Launch the interactive continuum editor - UUID-FREE VERSION."""
        if not self.current_item or not self.current_spec:
            QMessageBox.warning(
                self, "Cannot Edit", 
                "Please select an item first."
            )
            return
        
        try:
            # Launch the interactive continuum fitter
            from rbcodes.GUIs.interactive_continuum_fit import launch_interactive_continuum_fit_dialog
            
            spec = self.current_spec
            
            # Fix continuum masks format - flatten tuples to list
            existing_masks = []
            for mask_tuple in self.current_item.analysis.continuum_masks:
                if isinstance(mask_tuple, (list, tuple)) and len(mask_tuple) == 2:
                    existing_masks.extend([mask_tuple[0], mask_tuple[1]])
                else:
                    print(f"Warning: Skipping malformed mask: {mask_tuple}")
            
            print(f"DEBUG: Continuum masks converted to: {existing_masks}")
            
            # Prepare input parameters
            input_params = {
                'wave': spec.wave_slice,
                'flux': spec.flux_slice,
                'error': spec.error_slice,
                'velocity': spec.velo,
                'existing_masks': existing_masks,  # Flattened list
                'order': self.current_item.analysis.continuum_order,
                'use_weights': self.current_item.analysis.use_weights,
                'domain': [min(spec.velo), max(spec.velo)]
            }
            
            # Launch the interactive GUI
            result = launch_interactive_continuum_fit_dialog(**input_params)
            
            # Check if fitting was cancelled
            if result is None or result.get('cancelled', True):
                QMessageBox.information(
                    self, "Cancelled", 
                    "Continuum fitting was cancelled. No changes made."
                )
                return
            
            # Update continuum masks directly on the item object
            if 'masks' in result:
                # Convert flat list back to tuples
                masks = result['masks']
                mask_tuples = []
                for i in range(0, len(masks), 2):
                    if i + 1 < len(masks):
                        mask_tuples.append((masks[i], masks[i+1]))
                
                print(f"DEBUG: Updating continuum masks to: {mask_tuples}")
                
                # Update directly on the item object - NO UUID LOOKUP
                self.current_item.analysis.continuum_masks = mask_tuples
                self.current_item.analysis.last_modified = datetime.now().isoformat()
                self.current_item.analysis.processing_status = "complete"
                
                # Regenerate spectrum with new continuum
                self.display_item_spectrum(self.current_item)
                self.refresh_results_table()
                
                QMessageBox.information(
                    self, "Continuum Updated", 
                    "Continuum has been updated successfully."
                )
            
        except Exception as e:
            QMessageBox.warning(
                self, "Error", 
                f"Error during continuum editing: {str(e)}"
            )
            import traceback
            traceback.print_exc()
    
    def edit_ew_range(self):
        """Launch the EW range editor for the selected item - UUID-FREE VERSION."""
        if not self.current_item:
            QMessageBox.warning(
                self, "Cannot Edit", 
                "Please select an item first."
            )
            return
        
        try:
            # Create a dialog to edit the velocity range
            dialog = QDialog(self)
            dialog.setWindowTitle("Edit EW Measurement Range")
            layout = QVBoxLayout(dialog)
            
            # Get current values
            current_ew_vmin = self.current_item.template.ew_vmin
            current_ew_vmax = self.current_item.template.ew_vmax
            slice_vmin = self.current_item.template.slice_vmin
            slice_vmax = self.current_item.template.slice_vmax
            
            # Form for velocity range
            form = QFormLayout()
            
            # Show current slice range for reference
            slice_info = QLabel(f"Slice Range: [{slice_vmin}, {slice_vmax}] km/s")
            slice_info.setStyleSheet("QLabel { color: gray; font-style: italic; }")
            form.addRow("Reference:", slice_info)
            
            # Velocity min/max inputs
            ew_vmin_spin = QSpinBox()
            ew_vmin_spin.setRange(-5000, 0)
            ew_vmin_spin.setValue(current_ew_vmin)
            ew_vmin_spin.setSingleStep(10)
            ew_vmin_spin.setSuffix(" km/s")
            
            ew_vmax_spin = QSpinBox()
            ew_vmax_spin.setRange(0, 5000)
            ew_vmax_spin.setValue(current_ew_vmax)
            ew_vmax_spin.setSingleStep(10)
            ew_vmax_spin.setSuffix(" km/s")
            
            form.addRow("EW Min:", ew_vmin_spin)
            form.addRow("EW Max:", ew_vmax_spin)
            
            # Validation label
            validation_label = QLabel("")
            form.addRow("", validation_label)
            
            # Validation function
            def validate_ew_range():
                ew_min = ew_vmin_spin.value()
                ew_max = ew_vmax_spin.value()
                
                # Check if EW range is valid
                if ew_min >= ew_max:
                    validation_label.setText("⚠️ EW min must be less than EW max")
                    validation_label.setStyleSheet("QLabel { color: red; }")
                    return False
                
                # Check if EW range is within slice range (with tolerance)
                tolerance = 50
                if ew_min < (slice_vmin - tolerance) or ew_max > (slice_vmax + tolerance):
                    validation_label.setText(f"⚠️ EW range should be within slice range ± {tolerance} km/s")
                    validation_label.setStyleSheet("QLabel { color: red; }")
                    return False
                
                validation_label.setText("✓ EW range is valid")
                validation_label.setStyleSheet("QLabel { color: green; }")
                return True
            
            # Connect validation to value changes
            ew_vmin_spin.valueChanged.connect(validate_ew_range)
            ew_vmax_spin.valueChanged.connect(validate_ew_range)
            
            # Initial validation
            validate_ew_range()
            
            # SNR calculation option
            snr_check = QCheckBox("Calculate SNR")
            snr_check.setChecked(self.current_item.analysis.calculate_snr)
            
            binsize_spin = QSpinBox()
            binsize_spin.setRange(1, 10)
            binsize_spin.setValue(self.current_item.analysis.binsize)
            binsize_spin.setEnabled(snr_check.isChecked())
            
            snr_check.toggled.connect(binsize_spin.setEnabled)
            
            snr_layout = QHBoxLayout()
            snr_layout.addWidget(snr_check)
            snr_layout.addWidget(QLabel("Bin Size:"))
            snr_layout.addWidget(binsize_spin)
            
            form.addRow("SNR Options:", snr_layout)
            
            # Interactive selection option
            interactive_selection = QCheckBox("Enable Interactive Selection")
            interactive_selection.setChecked(False)
            interactive_selection.setToolTip("Enable clicking on the plot to set velocity limits")
            form.addRow("", interactive_selection)
            
            # Tooltip help
            help_text = QLabel("When enabled, click on plot to set integration limits. Press 'r' to reset.")
            help_text.setWordWrap(True)
            help_text.setStyleSheet("color: gray; font-style: italic;")
            form.addRow("", help_text)
            
            layout.addLayout(form)
            
            # Add a small preview plot if spectrum is available
            if self.current_spec and hasattr(self.current_spec, 'velo'):
                preview_label = QLabel("Current measurement range:")
                layout.addWidget(preview_label)
                
                canvas = MatplotlibCanvas(dialog, width=4, height=3, dpi=100)
                spec_obj = self.current_spec
                canvas.axes.step(spec_obj.velo, spec_obj.fnorm, 'k-', where='mid')
                canvas.axes.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
                canvas.axes.axvspan(current_ew_vmin, current_ew_vmax, alpha=0.2, color='blue')
                canvas.axes.set_xlabel('Velocity (km/s)')
                canvas.axes.set_ylabel('Normalized Flux')
                canvas.axes.set_title('Current EW Range')
                
                # Function to update the preview plot
                def update_preview():
                    canvas.axes.clear()
                    canvas.axes.step(spec_obj.velo, spec_obj.fnorm, 'k-', where='mid')
                    canvas.axes.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
                    canvas.axes.axvspan(ew_vmin_spin.value(), ew_vmax_spin.value(), alpha=0.2, color='blue')
                    canvas.axes.set_xlabel('Velocity (km/s)')
                    canvas.axes.set_ylabel('Normalized Flux')
                    canvas.axes.set_title('Updated EW Range')
                    
                    # Add markers at vmin and vmax if interactive selection is enabled
                    if interactive_selection.isChecked():
                        # Add circle markers at vmin and vmax with labels
                        vmin_val = ew_vmin_spin.value()
                        vmax_val = ew_vmax_spin.value()
                        canvas.axes.plot([vmin_val], [1.0], 'bo', ms=6, zorder=10)
                        canvas.axes.plot([vmax_val], [1.0], 'bo', ms=6, zorder=10)
                        
                        # Add labels
                        canvas.axes.text(vmin_val, 1.05, f'vmin: {vmin_val}', 
                                        ha='center', va='bottom', fontsize=8, color='blue')
                        canvas.axes.text(vmax_val, 1.05, f'vmax: {vmax_val}', 
                                        ha='center', va='bottom', fontsize=8, color='blue')
                    
                    canvas.draw()
                
                # Interactive selection event handlers
                def on_canvas_click(event):
                    """Handle mouse clicks on the plot for setting vmin/vmax."""
                    if not event.inaxes or not interactive_selection.isChecked():
                        return
                        
                    # Get the x-coordinate (velocity)
                    velocity = event.xdata
                    
                    # Determine whether to set vmin or vmax based on which is closer
                    current_vmin = ew_vmin_spin.value()
                    current_vmax = ew_vmax_spin.value()
                    
                    # If velocity is outside current range, set the closest endpoint
                    if velocity < current_vmin:
                        ew_vmin_spin.setValue(int(velocity))
                    elif velocity > current_vmax:
                        ew_vmax_spin.setValue(int(velocity))
                    else:
                        # If inside current range, set whichever endpoint is closer
                        if abs(velocity - current_vmin) < abs(velocity - current_vmax):
                            ew_vmin_spin.setValue(int(velocity))
                        else:
                            ew_vmax_spin.setValue(int(velocity))
                
                def on_canvas_key_press(event):
                    """Handle key presses on the plot."""
                    if not interactive_selection.isChecked():
                        return
                        
                    if event.key == 'r':
                        # Reset velocity limits to current values
                        ew_vmin_spin.setValue(current_ew_vmin)
                        ew_vmax_spin.setValue(current_ew_vmax)
                
                def toggle_interactive_mode(enabled):
                    """Enable or disable interactive selection."""
                    if enabled:
                        # Connect events
                        canvas.mpl_connect('button_press_event', on_canvas_click)
                        canvas.mpl_connect('key_press_event', on_canvas_key_press)
                        canvas.setFocus()
                        update_preview()  # Refresh to show markers
                    else:
                        # Events will be automatically disconnected when canvas is destroyed
                        update_preview()  # Refresh to hide markers
                
                # Connect interactive selection toggle
                interactive_selection.toggled.connect(toggle_interactive_mode)
                
                # Connect value changes to update preview
                ew_vmin_spin.valueChanged.connect(update_preview)
                ew_vmax_spin.valueChanged.connect(update_preview)
                
                layout.addWidget(canvas)
            
            # Add buttons
            button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            
            def accept_with_validation():
                if not validate_ew_range():
                    QMessageBox.warning(dialog, "Validation Error", 
                                      "Please fix the EW range validation errors before proceeding.")
                    return
                dialog.accept()
            
            button_box.accepted.connect(accept_with_validation)
            button_box.rejected.connect(dialog.reject)
            layout.addWidget(button_box)
            
            # Show dialog
            if dialog.exec_() == QDialog.Accepted:
                # Get new values
                new_ew_vmin = ew_vmin_spin.value()
                new_ew_vmax = ew_vmax_spin.value()
                new_calculate_snr = snr_check.isChecked()
                new_binsize = binsize_spin.value()
                
                # UUID-FREE UPDATE: Work directly with the item object
                success = self._update_ew_range_direct(
                    self.current_item, 
                    new_ew_vmin, 
                    new_ew_vmax, 
                    new_calculate_snr, 
                    new_binsize
                )
                
                if success:
                    # Refresh displays
                    self.display_item_spectrum(self.current_item)
                    self.refresh_results_table()
                    
                    QMessageBox.information(
                        self, "EW Range Updated", 
                        "EW measurement range has been updated and measurements recalculated."
                    )
                else:
                    QMessageBox.warning(
                        self, "Update Failed", 
                        "Failed to update EW range."
                    )
                
        except Exception as e:
            QMessageBox.warning(
                self, "Error", 
                f"Error during EW range editing: {str(e)}"
            )
            import traceback
            traceback.print_exc()
    
    


    def _update_ew_range_direct(self, item, new_ew_vmin, new_ew_vmax, calculate_snr, binsize):
        """Update EW range directly on the item object - NO UUID LOOKUP."""
        try:
            print(f"DEBUG: Direct update for {item.template.transition_name}")
            print(f"DEBUG: Old EW range: [{item.template.ew_vmin}, {item.template.ew_vmax}]")
            print(f"DEBUG: New EW range: [{new_ew_vmin}, {new_ew_vmax}]")
            
            # Update template parameters directly
            item.template.ew_vmin = new_ew_vmin
            item.template.ew_vmax = new_ew_vmax
            
            # Update analysis settings directly
            item.analysis.calculate_snr = calculate_snr
            item.analysis.binsize = binsize
            item.analysis.last_modified = datetime.now().isoformat()
            
            # Generate fresh spectrum with new parameters
            updated_spec = self._generate_spectrum_from_item(item)
            
            if not updated_spec:
                print("ERROR: Failed to generate spectrum")
                return False
            
            # Recalculate EW with new range
            updated_spec.compute_EW(
                item.template.transition,
                vmin=new_ew_vmin,
                vmax=new_ew_vmax,
                SNR=calculate_snr,
                _binsize=binsize
            )
            
            print(f"DEBUG: New EW: {updated_spec.W:.3f} ± {updated_spec.W_e:.3f} Å")
            
            # Update results directly on the item
            item.results.W = updated_spec.W
            item.results.W_e = updated_spec.W_e
            item.results.N = updated_spec.N
            item.results.N_e = updated_spec.N_e
            item.results.logN = updated_spec.logN
            item.results.logN_e = updated_spec.logN_e
            item.results.vel_centroid = getattr(updated_spec, 'vel_centroid', 0)
            item.results.vel_disp = getattr(updated_spec, 'vel_disp', 0)
            item.results.SNR = getattr(updated_spec, 'SNR', 0)
            item.results.calculation_timestamp = datetime.now().isoformat()
            
            # Mark as complete
            item.analysis.processing_status = "complete"
            
            # Update the current spectrum for display
            self.current_spec = updated_spec
            
            print(f"DEBUG: Direct update completed successfully")
            return True
            
        except Exception as e:
            print(f"ERROR: Direct update failed: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def select_all_items(self):
        """Select all items in the table."""
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item:
                checkbox_item.setCheckState(Qt.Checked)
        
        self.emit_selected_indices()
    
    def select_no_items(self):
        """Deselect all items in the table."""
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item:
                checkbox_item.setCheckState(Qt.Unchecked)
        
        self.emit_selected_indices()
    
    def select_failed_items(self):
        """Select items with error or needs processing status."""
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
        """Emit signal with currently selected item IDs."""
        selected_ids = []
        items = self.controller.master_table.get_all_items()
        
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                if row < len(items):
                    selected_ids.append(items[row].id)
        
        self.items_selected.emit(selected_ids)
    
    def process_selected_items(self):
        """Process the currently selected items."""
        selected_ids = []
        items = self.controller.master_table.get_all_items()
        
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                if row < len(items):
                    selected_ids.append(items[row].id)
        
        if not selected_ids:
            QMessageBox.warning(
                self, "No Selection", 
                "Please select items to process."
            )
            return
        
        # Emit signal with selected IDs
        self.items_selected.emit(selected_ids)
        
        # Switch to processing tab
        parent = self.parentWidget()
        if parent and isinstance(parent, QTabWidget):
            # Find the processing tab
            for i in range(parent.count()):
                if "Processing" in parent.tabText(i):
                    parent.setCurrentIndex(i)
                    break