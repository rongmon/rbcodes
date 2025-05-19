# rbcodes/GUIs/specgui/batch/panels/review_panel.py
import os
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QTableWidget, QTableWidgetItem,
                           QHeaderView, QCheckBox, QGroupBox, QSplitter,
                           QTabWidget, QScrollArea, QMessageBox, QFileDialog)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots."""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)

class ReviewPanel(QWidget):
    """Panel for reviewing batch processing results."""
    
    # Signals
    items_selected = pyqtSignal(list)  # Emitted when items are selected
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.results = []  # Results from batch processing
        self.selected_item = None  # Currently selected item for review
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add a splitter to divide results table and preview
        splitter = QSplitter(Qt.Vertical)
        
        # Results table
        results_group = QGroupBox("Batch Results")
        results_layout = QVBoxLayout()
        
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(8)
        self.results_table.setHorizontalHeaderLabels([
            "", "Filename", "Redshift", "Transition", "λ (Å)", "EW (Å)", "logN", "Status"
        ])
        
        # Set column widths
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Fixed)  # Checkbox column
        header.setSectionResizeMode(1, QHeaderView.Stretch)  # Filename
        for i in range(2, 8):
            header.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        
        self.results_table.setRowCount(0)
        self.results_table.itemClicked.connect(self.on_table_item_clicked)
        
        # Buttons for managing selection
        selection_buttons = QHBoxLayout()
        
# rbcodes/GUIs/specgui/batch/panels/review_panel.py (continued)
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
        """Set the batch processing results and update the display."""
        self.results = results
        self.update_results_table()
        
        # Clear selection
        self.selected_item = None
        self.edit_continuum_btn.setEnabled(False)
        self.edit_ew_btn.setEnabled(False)
        
        # Clear preview
        self.clear_preview()
    
    def update_results_table(self):
        """Update the results table with the current results."""
        # Clear existing rows
        self.results_table.setRowCount(0)
        
        # Add results
        for i, result in enumerate(self.results):
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)
            
            # Checkbox for selection
            checkbox_item = QTableWidgetItem()
            checkbox_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            checkbox_item.setCheckState(Qt.Unchecked)
            self.results_table.setItem(row, 0, checkbox_item)
            
            # Item data
            filename = os.path.basename(result.get('filename', 'Unknown'))
            self.results_table.setItem(row, 1, QTableWidgetItem(filename))
            self.results_table.setItem(row, 2, QTableWidgetItem(f"{result.get('redshift', 0):.6f}"))
            
            transition_name = result.get('transition_name', 'Unknown')
            self.results_table.setItem(row, 3, QTableWidgetItem(transition_name))
            
            transition = result.get('transition', 0)
            self.results_table.setItem(row, 4, QTableWidgetItem(f"{transition:.2f}"))
            
            ew_value = result.get('W', 0)
            ew_error = result.get('W_e', 0)
            ew_text = f"{ew_value:.3f} ± {ew_error:.3f}"
            self.results_table.setItem(row, 5, QTableWidgetItem(ew_text))
            
            logn_value = result.get('logN', 0)
            logn_error = result.get('logN_e', 0)
            logn_text = f"{logn_value:.2f} ± {logn_error:.2f}"
            self.results_table.setItem(row, 6, QTableWidgetItem(logn_text))
            
            # Status with color coding
            status = result.get('status', 'unknown')
            status_item = QTableWidgetItem(status)
            
            if status == 'success':
                status_item.setBackground(QColor(200, 255, 200))  # Light green
            elif status == 'warning':
                status_item.setBackground(QColor(255, 255, 200))  # Light yellow
            elif status == 'error':
                status_item.setBackground(QColor(255, 200, 200))  # Light red
            
            self.results_table.setItem(row, 7, status_item)
    
    def on_table_item_clicked(self, item):
        """Handle clicks on table items."""
        row = item.row()
        
        if row < len(self.results):
            self.selected_item = self.results[row]
            self.update_preview()
            
            # Enable editing buttons
            self.edit_continuum_btn.setEnabled(True)
            self.edit_ew_btn.setEnabled(True)
    
    def update_preview(self):
        """Update the preview with the selected item's data."""
        if not self.selected_item:
            self.clear_preview()
            return
        
        # Update spectrum plot
        self.update_spectrum_plot()
        
        # Update details
        self.update_details()
    
    def update_spectrum_plot(self):
        """Update the spectrum plot with the selected item's data."""
        if not self.selected_item:
            return
        
        # Clear existing plot
        self.spectrum_canvas.axes.clear()
        
        # Get the rb_spec object from the selected item
        spec_object = self.selected_item.get('spec_object')
        
        if spec_object and hasattr(spec_object, 'velo') and hasattr(spec_object,     'fnorm'):
            # Plot the normalized flux from the rb_spec object
            self.spectrum_canvas.axes.step(
                spec_object.velo, 
                spec_object.fnorm, 
                'k-', 
                where='mid', 
                label='Normalized Flux'
            )
            
            # Plot the normalized error
            if hasattr(spec_object, 'enorm'):
                self.spectrum_canvas.axes.step(
                    spec_object.velo, 
                    spec_object.enorm, 
                    'r-', 
                    where='mid', 
                    alpha=0.5, 
                    label='Error'
                )
            
            # Reference line at y=1.0
            self.spectrum_canvas.axes.axhline(y=1.0, color='g', linestyle='--',     alpha=0.7)
            
            # Add EW measurement region
            vmin = spec_object.vmin if hasattr(spec_object, 'vmin') else     self.selected_item.get('vmin', -200)
            vmax = spec_object.vmax if hasattr(spec_object, 'vmax') else     self.selected_item.get('vmax', 200)
            
            self.spectrum_canvas.axes.axvspan(vmin, vmax, alpha=0.1, color='blue',     label='EW Region')
            self.spectrum_canvas.axes.axvline(x=vmin, color='b', linestyle=':',     alpha=0.7)
            self.spectrum_canvas.axes.axvline(x=vmax, color='b', linestyle=':',     alpha=0.7)
            
            # Set axis limits based on the data
            x_buffer = (max(spec_object.velo) - min(spec_object.velo)) * 0.1
            self.spectrum_canvas.axes.set_xlim(min(spec_object.velo) - x_buffer, max(    spec_object.velo) + x_buffer)
            
            # Try to set y limits nicely
            if hasattr(spec_object, 'fnorm'):
                y_data = spec_object.fnorm
                y_min = max(0, min(y_data) - 0.1)
                y_max = min(2.0, max(y_data) + 0.3)
                self.spectrum_canvas.axes.set_ylim(y_min, y_max)
            else:
                self.spectrum_canvas.axes.set_ylim(-0.1, 1.5)
        else:
            # Fallback to placeholder if rb_spec object is not available
            velo = range(-300, 301)
            flux = [0.9 + 0.2 * (1 - (i/100)**2) for i in velo]
            
            self.spectrum_canvas.axes.plot(velo, flux, 'k--', label='Preview (No     data)')
            self.spectrum_canvas.axes.axhline(y=1.0, color='r', linestyle='--',     alpha=0.7)
            
            # Add placeholder EW region
            vmin = self.selected_item.get('vmin', -200)
            vmax = self.selected_item.get('vmax', 200)
            self.spectrum_canvas.axes.axvspan(vmin, vmax, alpha=0.1, color='blue',     label='EW Region')
            
            # Set axis limits
            self.spectrum_canvas.axes.set_xlim(-300, 300)
            self.spectrum_canvas.axes.set_ylim(-0.1, 1.5)
        
        # Add labels and measurement info
        transition_name = self.selected_item.get('transition_name', 'Unknown')
        transition = self.selected_item.get('transition', 0)
        self.spectrum_canvas.axes.set_title(f"{transition_name} ({transition:.2f} Å)")
        self.spectrum_canvas.axes.set_xlabel('Velocity (km/s)')
        self.spectrum_canvas.axes.set_ylabel('Normalized Flux')
        
        # Add EW and column density info
        ew = self.selected_item.get('W', 0)
        ew_err = self.selected_item.get('W_e', 0)
        logn = self.selected_item.get('logN', 0)
        logn_err = self.selected_item.get('logN_e', 0)
        
        info_text = f"EW = {ew:.3f} ± {ew_err:.3f} Å\nlog N = {logn:.2f} ±     {logn_err:.2f}"
        
        # Add SNR if available
        if 'SNR' in self.selected_item and self.selected_item['SNR'] > 0:
            info_text += f"\nSNR = {self.selected_item['SNR']:.1f}"
        
        self.spectrum_canvas.axes.text(0.98, 0.05, info_text,     transform=self.spectrum_canvas.axes.transAxes, 
                                     horizontalalignment='right',     verticalalignment='bottom',
                                     bbox=dict(facecolor='white', alpha=0.7,     boxstyle='round'))
        
        # Add a legend
        #self.spectrum_canvas.axes.legend(loc='upper right')
        
        # Apply tight layout
        self.spectrum_canvas.fig.tight_layout()
        
        # Draw the updated plot
        self.spectrum_canvas.draw()    

    
    def update_details(self):
        """Update the details view with the selected item's data."""
        if not self.selected_item:
            self.details_label.setText("Select an item to view details.")
            return
        
        # Format details text
        details = f"<b>Filename:</b> {self.selected_item.get('filename', 'Unknown')}<br>"
        details += f"<b>Redshift:</b> {self.selected_item.get('redshift', 0):.6f}<br>"
        details += f"<b>Transition:</b> {self.selected_item.get('transition_name', 'Unknown')} ({self.selected_item.get('transition', 0):.2f} Å)<br>"
        details += f"<b>Velocity Range:</b> {self.selected_item.get('vmin', 0)} to {self.selected_item.get('vmax', 0)} km/s<br>"
        details += f"<b>Equivalent Width:</b> {self.selected_item.get('W', 0):.3f} ± {self.selected_item.get('W_e', 0):.3f} Å<br>"
        details += f"<b>Column Density (log):</b> {self.selected_item.get('logN', 0):.2f} ± {self.selected_item.get('logN_e', 0):.2f}<br>"
        details += f"<b>Status:</b> {self.selected_item.get('status', 'unknown')}<br>"
        
        # Add any additional details if available
        if 'SNR' in self.selected_item:
            details += f"<b>SNR:</b> {self.selected_item['SNR']:.1f}<br>"
        
        if 'vel_centroid' in self.selected_item:
            details += f"<b>Velocity Centroid:</b> {self.selected_item['vel_centroid']:.1f} km/s<br>"
        
        if 'vel_disp' in self.selected_item:
            details += f"<b>Velocity Dispersion:</b> {self.selected_item['vel_disp']:.1f} km/s<br>"
        
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
        """Select items with error or warning status."""
        for row in range(self.results_table.rowCount()):
            status_item = self.results_table.item(row, 7)
            checkbox_item = self.results_table.item(row, 0)
            
            if status_item and checkbox_item:
                status = status_item.text().lower()
                if status in ['error', 'warning']:
                    checkbox_item.setCheckState(Qt.Checked)
                else:
                    checkbox_item.setCheckState(Qt.Unchecked)
        
        self.emit_selected_indices()
    
    def emit_selected_indices(self):
        """Emit signal with currently selected indices."""
        selected = []
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                selected.append(row)
        
        self.items_selected.emit(selected)
    
    def process_selected_items(self):
        """Process the currently selected items."""
        selected = []
        for row in range(self.results_table.rowCount()):
            checkbox_item = self.results_table.item(row, 0)
            if checkbox_item and checkbox_item.checkState() == Qt.Checked:
                selected.append(row)
        
        if not selected:
            QMessageBox.warning(
                self, "No Selection", 
                "Please select items to process."
            )
            return
        
        # Emit signal with selected indices
        self.items_selected.emit(selected)
        
        # Switch to processing tab
        parent = self.parentWidget()
        if parent and isinstance(parent, QTabWidget):
            # Find the processing tab
            for i in range(parent.count()):
                if parent.tabText(i) == "Processing":
                    parent.setCurrentIndex(i)
                    break

    def edit_continuum(self):
        """Launch the interactive continuum editor for the selected item."""
        if not self.selected_item:
            return
        
        # Get the rb_spec object
        spec_object = self.selected_item.get('spec_object')
        
        if not spec_object or not hasattr(spec_object, 'velo'):
            QMessageBox.warning(
                self, "Cannot Edit", 
                "Cannot edit continuum: the spectrum data is not available."
            )
            return
        
        try:
            # Launch the interactive continuum fitter
            from rbcodes.GUIs.interactive_continuum_fit import launch_interactive_continuum_fit_dialog
            
            # Prepare input parameters
            input_params = {
                'wave': spec_object.wave_slice,
                'flux': spec_object.flux_slice,
                'error': spec_object.error_slice,
                'velocity': spec_object.velo,
                'existing_masks': getattr(spec_object, 'continuum_masks', []),
                'order': 3,  # Default order, could be improved
                'use_weights': False,
                'domain': [min(spec_object.velo), max(spec_object.velo)]
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
            
            # Update the rb_spec object with the new continuum
            spec_object.cont = result.get('continuum')
            spec_object.fnorm = spec_object.flux_slice / spec_object.cont
            spec_object.enorm = spec_object.error_slice / spec_object.cont
            
            # Store mask information if available
            if 'masks' in result:
                spec_object.continuum_masks = result.get('masks', [])
            
            # Re-compute EW with the new continuum
            spec_object.compute_EW(
                spec_object.transition,
                vmin=spec_object.vmin,
                vmax=spec_object.vmax,
                SNR=getattr(spec_object, 'SNR', 0) > 0,  # Use SNR if previously calculated
                _binsize=3  # Could be from settings
            )
            
            # Update the selected item with new EW values
            self.selected_item['W'] = spec_object.W
            self.selected_item['W_e'] = spec_object.W_e
            self.selected_item['logN'] = spec_object.logN
            self.selected_item['logN_e'] = spec_object.logN_e
            if hasattr(spec_object, 'SNR'):
                self.selected_item['SNR'] = spec_object.SNR
            
            # Update the plot and details
            self.update_preview()
            
            # Update the results table row
            row = -1
            for i in range(self.results_table.rowCount()):
                if self.results_table.item(i, 1).text() == os.path.basename(self.selected_item.get('filename', '')) and \
                   float(self.results_table.item(i, 2).text()) == self.selected_item.get('redshift', 0):
                    row = i
                    break
            
            if row >= 0:
                # Update EW cell
                ew_text = f"{self.selected_item['W']:.3f} ± {self.selected_item['W_e']:.3f}"
                self.results_table.setItem(row, 5, QTableWidgetItem(ew_text))
                
                # Update logN cell
                logn_text = f"{self.selected_item['logN']:.2f} ± {self.selected_item['logN_e']:.2f}"
                self.results_table.setItem(row, 6, QTableWidgetItem(logn_text))
            
            QMessageBox.information(
                self, "Continuum Updated", 
                "Continuum has been updated and measurements recalculated."
            )
            
        except Exception as e:
            QMessageBox.warning(
                self, "Error", 
                f"Error during continuum editing: {str(e)}"
            )
            import traceback
            traceback.print_exc()
    
    def edit_ew_range(self):
        """Launch the EW range editor for the selected item."""
        if not self.selected_item:
            return
        
        # Get the rb_spec object
        spec_object = self.selected_item.get('spec_object')
        
        if not spec_object or not hasattr(spec_object, 'velo'):
            QMessageBox.warning(
                self, "Cannot Edit", 
                "Cannot edit EW range: the spectrum data is not available."
            )
            return
        
        try:
            # Create a dialog to edit the velocity range
            dialog = QDialog(self)
            dialog.setWindowTitle("Edit EW Measurement Range")
            layout = QVBoxLayout(dialog)
            
            # Get current values
            vmin = spec_object.vmin if hasattr(spec_object, 'vmin') else     self.selected_item.get('vmin', -200)
            vmax = spec_object.vmax if hasattr(spec_object, 'vmax') else     self.selected_item.get('vmax', 200)
            
            # Form for velocity range
            form = QFormLayout()
            
            # Velocity min/max inputs
            vmin_spin = QSpinBox()
            vmin_spin.setRange(-5000, 0)
            vmin_spin.setValue(vmin)
            vmin_spin.setSingleStep(10)
            
            vmax_spin = QSpinBox()
            vmax_spin.setRange(0, 5000)
            vmax_spin.setValue(vmax)
            vmax_spin.setSingleStep(10)
            
            form.addRow("Velocity Min (km/s):", vmin_spin)
            form.addRow("Velocity Max (km/s):", vmax_spin)
            
            # SNR calculation option
            snr_check = QCheckBox("Calculate SNR")
            snr_check.setChecked(getattr(spec_object, 'SNR', 0) > 0)
            
            binsize_spin = QSpinBox()
            binsize_spin.setRange(1, 10)
            binsize_spin.setValue(3)
            binsize_spin.setEnabled(snr_check.isChecked())
            
            snr_check.toggled.connect(binsize_spin.setEnabled)
            
            snr_layout = QHBoxLayout()
            snr_layout.addWidget(snr_check)
            snr_layout.addWidget(QLabel("Bin Size:"))
            snr_layout.addWidget(binsize_spin)
            
            form.addRow("", snr_layout)
            
            layout.addLayout(form)
            
            # Add a small preview plot
            preview_label = QLabel("Current measurement range:")
            layout.addWidget(preview_label)
            
            canvas = MatplotlibCanvas(dialog, width=4, height=3, dpi=100)
            canvas.axes.step(spec_object.velo, spec_object.fnorm, 'k-', where='mid')
            canvas.axes.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
            canvas.axes.axvspan(vmin, vmax, alpha=0.2, color='blue')
            canvas.axes.set_xlabel('Velocity (km/s)')
            canvas.axes.set_ylabel('Normalized Flux')
            
            # Function to update the preview plot
            def update_preview():
                canvas.axes.clear()
                canvas.axes.step(spec_object.velo, spec_object.fnorm, 'k-',     where='mid')
                canvas.axes.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
                canvas.axes.axvspan(vmin_spin.value(), vmax_spin.value(), alpha=0.2,     color='blue')
                canvas.axes.set_xlabel('Velocity (km/s)')
                canvas.axes.set_ylabel('Normalized Flux')
                canvas.draw()
            
            # Connect value changes to update preview
            vmin_spin.valueChanged.connect(update_preview)
            vmax_spin.valueChanged.connect(update_preview)
            
            layout.addWidget(canvas)
            
            # Add buttons
            button_box = QDialogButtonBox(QDialogButtonBox.Ok |     QDialogButtonBox.Cancel)
            button_box.accepted.connect(dialog.accept)
            button_box.rejected.connect(dialog.reject)
            layout.addWidget(button_box)
            
            # Show dialog
            if dialog.exec_() == QDialog.Accepted:
                # Get new values
                new_vmin = vmin_spin.value()
                new_vmax = vmax_spin.value()
                calculate_snr = snr_check.isChecked()
                binsize = binsize_spin.value()
                
                # Re-compute EW with the new velocity range
                spec_object.compute_EW(
                    spec_object.transition,
                    vmin=new_vmin,
                    vmax=new_vmax,
                    SNR=calculate_snr,
                    _binsize=binsize
                )
                
                # Update the selected item with new values
                self.selected_item['vmin'] = new_vmin
                self.selected_item['vmax'] = new_vmax
                self.selected_item['W'] = spec_object.W
                self.selected_item['W_e'] = spec_object.W_e
                self.selected_item['logN'] = spec_object.logN
                self.selected_item['logN_e'] = spec_object.logN_e
                if hasattr(spec_object, 'SNR'):
                    self.selected_item['SNR'] = spec_object.SNR
                
                # Update the plot and details
                self.update_preview()
                
                # Update the results table row
                row = -1
                for i in range(self.results_table.rowCount()):
                    if self.results_table.item(i, 1).text() == os.path.basename(    self.selected_item.get('filename', '')) and \
                       float(self.results_table.item(i, 2).text()) ==     self.selected_item.get('redshift', 0):
                        row = i
                        break
                
                if row >= 0:
                    # Update EW cell
                    ew_text = f"{self.selected_item['W']:.3f} ±     {self.selected_item['W_e']:.3f}"
                    self.results_table.setItem(row, 5, QTableWidgetItem(ew_text))
                    
                    # Update logN cell
                    logn_text = f"{self.selected_item['logN']:.2f} ±     {self.selected_item['logN_e']:.2f}"
                    self.results_table.setItem(row, 6, QTableWidgetItem(logn_text))
                
                QMessageBox.information(
                    self, "EW Range Updated", 
                    "EW measurement range has been updated and measurements     recalculated."
                )
                
        except Exception as e:
            QMessageBox.warning(
                self, "Error", 
                f"Error during EW range editing: {str(e)}"
            )
            import traceback
            traceback.print_exc()    
        
        QMessageBox.information(
            self, "Not Implemented", 
            "EW range editing will be implemented in a future version."
        )