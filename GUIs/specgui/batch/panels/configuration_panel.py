# rbcodes/GUIs/specgui/batch/panels/configuration_panel.py
import os
import csv
import json
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QTableWidget, QTableWidgetItem,
                           QHeaderView, QComboBox, QGroupBox, QFormLayout,
                           QMessageBox, QFileDialog, QCheckBox, QSpinBox,
                           QDoubleSpinBox, QDialog, QDialogButtonBox, QLineEdit,
                           QAbstractItemView)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor
from datetime import datetime
import numpy as np
class ConfigurationPanel(QWidget):
    """Panel for configuring batch processing items using master table."""
    
    # Signals
    configuration_changed = pyqtSignal()  # Emitted when batch configuration changes
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.init_ui()
        
        # Connect to controller signals
        self.controller.table_changed.connect(self.refresh_table)
        
        # Connect table selection to update button states
        self.batch_table.itemSelectionChanged.connect(self.update_button_states)
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add instructions
        instructions = QLabel(
            "Configure batch processing by adding items to analyze. "
            "Each item can have different slice and EW measurement velocity ranges."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)
        
        # Input method selection
        input_group = QGroupBox("Add Items to Batch")
        input_layout = QHBoxLayout()
        
        # Single item button
        self.add_btn = QPushButton("Add Single Ion")
        self.add_btn.clicked.connect(self.add_batch_item)
        self.add_btn.setToolTip("Manually add one spectrum analysis item")
        
        # Bulk import button (new)
        self.import_json_btn = QPushButton("Import Multiple JSON Files")
        self.import_json_btn.clicked.connect(self.import_multiple_json_files)
        self.import_json_btn.setToolTip("Import multiple rb_spec JSON files to create batch items")
        
        # CSV import button
        self.import_btn = QPushButton("Import from CSV")
        self.import_btn.clicked.connect(self.import_from_csv)
        self.import_btn.setToolTip("Import batch items from a CSV template file")
        
        input_layout.addWidget(self.add_btn)
        input_layout.addWidget(self.import_json_btn)
        input_layout.addWidget(self.import_btn)
        
        input_group.setLayout(input_layout)
        main_layout.addWidget(input_group)
        
        # Batch items table
        self.table_group = QGroupBox("Batch Items")
        table_layout = QVBoxLayout()
        
        self.batch_table = QTableWidget()
        self.batch_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setup_table()
        
        # CREATE buttons first, THEN add them to layout
        self.export_template_btn = QPushButton("Export Template CSV")
        self.export_template_btn.clicked.connect(self.export_template_csv)
        
        # NEW: Edit button
        self.edit_btn = QPushButton("Edit Selected Item")
        self.edit_btn.clicked.connect(self.edit_batch_item)
        self.edit_btn.setEnabled(False)  # Initially disabled
        self.edit_btn.setToolTip("Edit the selected batch item parameters")
        
        self.remove_btn = QPushButton("Remove Selected")
        self.remove_btn.clicked.connect(self.remove_batch_item)
        self.remove_btn.setEnabled(False)  # Initially disabled
        
        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.clicked.connect(self.clear_batch_items)
        
        # Buttons for managing batch items
        table_buttons = QHBoxLayout()
        table_buttons.addWidget(self.export_template_btn)
        table_buttons.addWidget(self.edit_btn)  # Add edit button
        table_buttons.addWidget(self.remove_btn)
        table_buttons.addWidget(self.clear_btn)
        
        table_layout.addWidget(self.batch_table)
        table_layout.addLayout(table_buttons)
        
        self.table_group.setLayout(table_layout)
        main_layout.addWidget(self.table_group)
        
        # Configuration file management
        config_group = QGroupBox("Configuration File")
        config_layout = QHBoxLayout()
        
        self.save_config_btn = QPushButton("Save Configuration")
        self.save_config_btn.clicked.connect(self.save_configuration)
        
        self.load_config_btn = QPushButton("Load Configuration")
        self.load_config_btn.clicked.connect(self.load_configuration)
        
        config_layout.addWidget(self.save_config_btn)
        config_layout.addWidget(self.load_config_btn)
        
        config_group.setLayout(config_layout)
        main_layout.addWidget(config_group)

    def update_button_states(self):
        """Update button states based on table selection."""
        selected_rows = set(index.row() for index in self.batch_table.selectedIndexes())
        
        # Edit button: enabled only when exactly one row is selected
        self.edit_btn.setEnabled(len(selected_rows) == 1)
        
        # Remove button: enabled when at least one row is selected
        self.remove_btn.setEnabled(len(selected_rows) > 0)

    def edit_batch_item(self):
        """Edit the selected batch item."""
        selected_rows = set(index.row() for index in self.batch_table.selectedIndexes())
        
        if len(selected_rows) != 1:
            QMessageBox.warning(self, "Selection Error", "Please select exactly one item to edit.")
            return
        
        row_index = list(selected_rows)[0]
        
        # Get the item from master table
        item = self.controller.master_table.get_item(row_index)
        if not item:
            QMessageBox.warning(self, "Error", "Could not retrieve item data.")
            return
        
        # Launch edit dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Edit Ion Parameters")
        dialog.setMinimumSize(500, 450)
        layout = QVBoxLayout(dialog)
        
        # Main form - pre-populated with current values
        form_layout = QFormLayout()
        
        # Filename
        file_layout = QHBoxLayout()
        filename = QLineEdit()
        filename.setText(item.template.filename)  # Pre-populate
        filename.setReadOnly(True)
        
        browse_btn = QPushButton("Browse")
        
        def browse_file():
            options = QFileDialog.Options()
            filepath, _ = QFileDialog.getOpenFileName(
                dialog, "Select Spectrum File", "", 
                "All Files (*);;FITS Files (*.fits);;JSON Files (*.json)",
                options=options
            )
            if filepath:
                filename.setText(filepath)
        
        browse_btn.clicked.connect(browse_file)
        file_layout.addWidget(filename)
        file_layout.addWidget(browse_btn)
        form_layout.addRow("Spectrum File:", file_layout)
        
        # Redshift
        redshift = QDoubleSpinBox()
        redshift.setDecimals(6)
        redshift.setRange(0.0, 10.0)
        redshift.setValue(item.template.redshift)  # Pre-populate
        form_layout.addRow("Redshift:", redshift)
        
        # Transition
        transition_layout = QHBoxLayout()
        transition_combo = QComboBox()
        transitions = [
            "Custom",
            "Lyα (1215.67 Å)",
            "Lyβ (1025.72 Å)",
            "Lyγ (972.54 Å)",
            "CIV (1548.20 Å)",
            "CIV (1550.78 Å)",
            "MgII (2796.35 Å)",
            "MgII (2803.53 Å)",
            "SiIV (1393.76 Å)",
            "SiIV (1402.77 Å)"
        ]
        transition_combo.addItems(transitions)
        
        custom_wavelength = QDoubleSpinBox()
        custom_wavelength.setRange(1.0, 10000.0)
        custom_wavelength.setValue(item.template.transition)  # Pre-populate
        custom_wavelength.setDecimals(2)
        custom_wavelength.setEnabled(True)
        
        # Set transition combo to match current transition
        current_transition = item.template.transition
        found_match = False
        for i, trans_text in enumerate(transitions[1:], 1):  # Skip "Custom"
            try:
                wavelength_text = trans_text.split("(")[1].split(" Å")[0]
                wavelength = float(wavelength_text)
                if abs(wavelength - current_transition) < 0.01:
                    transition_combo.setCurrentIndex(i)
                    custom_wavelength.setEnabled(False)
                    found_match = True
                    break
            except (IndexError, ValueError):
                continue
        
        if not found_match:
            transition_combo.setCurrentIndex(0)  # Set to "Custom"
            custom_wavelength.setEnabled(True)
        
        transition_combo.currentTextChanged.connect(
            lambda text: custom_wavelength.setEnabled(text == "Custom")
        )
        
        transition_layout.addWidget(transition_combo)
        transition_layout.addWidget(custom_wavelength)
        form_layout.addRow("Transition:", transition_layout)
        
        # Line list selection
        linelist_combo = QComboBox()
        linelist_combo.addItems(["atom", "LLS", "LLS Small", "DLA", "LBG", "Gal"])
        # Set current linelist
        linelist_index = linelist_combo.findText(item.template.linelist)
        if linelist_index >= 0:
            linelist_combo.setCurrentIndex(linelist_index)
        form_layout.addRow("Line List:", linelist_combo)
        
        # Add separator
        separator = QLabel()
        separator.setStyleSheet("QLabel { border-bottom: 1px solid gray; }")
        form_layout.addRow(separator)
        
        # Velocity ranges section
        velocity_label = QLabel("<b>Velocity Ranges:</b>")
        form_layout.addRow(velocity_label)
        
        # Add explanatory text
        explanation = QLabel(
            "• Slice Range: Wide range for extracting spectrum (includes continuum regions)\n"
            "• EW Range: Narrow range for equivalent width integration (around absorption)"
        )
        explanation.setWordWrap(True)
        explanation.setStyleSheet("QLabel { color: gray; font-style: italic; margin-bottom: 10px; }")
        form_layout.addRow(explanation)
        
        # Slice velocity range - pre-populated
        slice_range_layout = QHBoxLayout()
        slice_vmin = QSpinBox()
        slice_vmin.setRange(-5000, 0)
        slice_vmin.setValue(item.template.slice_vmin)  # Pre-populate
        slice_vmin.setSuffix(" km/s")
        
        slice_vmax = QSpinBox()
        slice_vmax.setRange(0, 5000)
        slice_vmax.setValue(item.template.slice_vmax)  # Pre-populate
        slice_vmax.setSuffix(" km/s")
        
        slice_range_layout.addWidget(QLabel("Min:"))
        slice_range_layout.addWidget(slice_vmin)
        slice_range_layout.addWidget(QLabel("Max:"))
        slice_range_layout.addWidget(slice_vmax)
        
        form_layout.addRow("Slice Range:", slice_range_layout)
        
        # EW measurement range - pre-populated
        ew_range_layout = QHBoxLayout()
        ew_vmin = QSpinBox()
        ew_vmin.setRange(-2000, 0)
        ew_vmin.setValue(item.template.ew_vmin)  # Pre-populate
        ew_vmin.setSuffix(" km/s")
        
        ew_vmax = QSpinBox()
        ew_vmax.setRange(0, 2000)
        ew_vmax.setValue(item.template.ew_vmax)  # Pre-populate
        ew_vmax.setSuffix(" km/s")
        
        ew_range_layout.addWidget(QLabel("Min:"))
        ew_range_layout.addWidget(ew_vmin)
        ew_range_layout.addWidget(QLabel("Max:"))
        ew_range_layout.addWidget(ew_vmax)
        
        form_layout.addRow("EW Range:", ew_range_layout)
        
        # Validation function and label
        validation_label = QLabel("")
        form_layout.addRow(validation_label)
        
        def validate_ranges():
            """Validate that EW range is within slice range."""
            slice_min, slice_max = slice_vmin.value(), slice_vmax.value()
            ew_min, ew_max = ew_vmin.value(), ew_vmax.value()
            
            # Check if EW range is within slice range (with tolerance)
            tolerance = 50  # km/s
            if ew_min < (slice_min - tolerance) or ew_max > (slice_max + tolerance):
                validation_label.setText(f"⚠️ EW range should be within slice range ± {tolerance} km/s")
                validation_label.setStyleSheet("QLabel { color: red; }")
                return False
            else:
                validation_label.setText("✓ Velocity ranges are valid")
                validation_label.setStyleSheet("QLabel { color: green; }")
                return True
        
        # Connect validation to value changes
        slice_vmin.valueChanged.connect(validate_ranges)
        slice_vmax.valueChanged.connect(validate_ranges)
        ew_vmin.valueChanged.connect(validate_ranges)
        ew_vmax.valueChanged.connect(validate_ranges)
        
        # Initial validation
        validate_ranges()
        
        layout.addLayout(form_layout)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        
        def accept_with_validation():
            if not validate_ranges():
                QMessageBox.warning(dialog, "Validation Error", 
                                  "Please fix the velocity range validation errors before proceeding.")
                return
            dialog.accept()
        
        button_box.accepted.connect(accept_with_validation)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            # Get new values
            file_path = filename.text()
            if not file_path:
                QMessageBox.warning(self, "Input Error", "Please select a spectrum file.")
                return
                
            z = redshift.value()
            
            # Get transition details - UNIFIED rb_setline approach
            trans_text = transition_combo.currentText()
            linelist_text = linelist_combo.currentText()
            
            # Extract wavelength from either dropdown or custom input
            if trans_text == "Custom":
                trans_wavelength = custom_wavelength.value()
            else:
                # Extract wavelength from dropdown selection
                trans_wavelength = float(trans_text.split("(")[1].split(" Å")[0])
            
            # Always use rb_setline to get proper transition name from database
            try:
                from rbcodes.IGM import rb_setline as s
                line_info = s.rb_setline(trans_wavelength, 'closest', linelist=linelist_text)
                trans_name = line_info['name']
                actual_wavelength = line_info['wave']
                
                print(f"Input wavelength {trans_wavelength:.2f} Å matched to: {trans_name} ({actual_wavelength:.2f} Å)")
                
            except Exception as e:
                print(f"Warning: Could not find transition in database: {e}")
                # Fallback to current name or generic name if rb_setline fails
                if trans_text == "Custom":
                    trans_name = item.template.transition_name  # Keep current name
                else:
                    trans_name = trans_text.split(" (")[0]
            
            # Get velocity ranges
            slice_min, slice_max = slice_vmin.value(), slice_vmax.value()
            ew_min, ew_max = ew_vmin.value(), ew_vmax.value()
            
            # Update the item in master table
            success = self.controller.master_table.update_template(
                row_index,
                filename=file_path,
                redshift=z,
                transition=trans_wavelength,
                transition_name=trans_name,
                slice_vmin=slice_min,
                slice_vmax=slice_max,
                ew_vmin=ew_min,
                ew_vmax=ew_max,
                linelist=linelist_text
            )
            
            if success:
                # Reset status to ready since parameters changed
                self.controller.master_table.update_analysis(
                    row_index,
                    processing_status="ready",
                    error_message="",
                    last_modified=datetime.now().isoformat()
                )
                
                # Clear any existing results since they're no longer valid
                self.controller.master_table.update_results(
                    row_index,
                    W=0.0, W_e=0.0, N=0.0, N_e=0.0,
                    logN=0.0, logN_e=0.0, vel_centroid=0.0,
                    vel_disp=0.0, SNR=0.0,
                    calculation_timestamp=""
                )
                
                # Remove any cached rb_spec object since parameters changed
                if row_index in self.controller.master_table.rb_spec_objects:
                    del self.controller.master_table.rb_spec_objects[row_index]
                
                # Refresh display
                self.refresh_table()
                self.configuration_changed.emit()
                
                # Show success message
                self.controller.status_updated.emit(f"Item updated: {trans_name}. Status reset to 'ready' - reprocessing required.")
                
            else:
                QMessageBox.warning(self, "Update Error", "Failed to update item.")

    def import_multiple_json_files(self):
        """Import multiple rb_spec JSON files to create batch items."""
        # Select multiple JSON files
        options = QFileDialog.Options()
        filepaths, _ = QFileDialog.getOpenFileNames(
            self, "Select rb_spec JSON Files", "", 
            "JSON Files (*.json);;All Files (*)",
            options=options
        )
        
        if not filepaths:
            return
        
        imported_count = 0
        completed_count = 0  # Track how many were already analyzed
        skipped_count = 0
        
        for filepath in filepaths:
            try:
                # Load the complete rb_spec object
                from rbcodes.GUIs.rb_spec import load_rb_spec_object
                spec_object = load_rb_spec_object(filepath)
                
                # Extract basic parameters from rb_spec object
                if not all(hasattr(spec_object, attr) for attr in ['zabs', 'trans_wave']):
                    print(f"Skipping {os.path.basename(filepath)} - missing required rb_spec attributes")
                    skipped_count += 1
                    continue
                
                # Extract template parameters from the loaded rb_spec object
                redshift = spec_object.zabs
                transition = spec_object.trans_wave
                transition_name = getattr(spec_object, 'trans', f"λ {transition:.2f}")
                linelist = getattr(spec_object, 'linelist', 'atom')
                
                # Get slice ranges from velocity array if available
                if hasattr(spec_object, 'velo') and len(spec_object.velo) > 0:
                    slice_vmin = int(min(spec_object.velo))
                    slice_vmax = int(max(spec_object.velo))
                else:
                    # Default slice ranges
                    slice_vmin = -1500
                    slice_vmax = 1500
                
                # Get EW ranges if available, otherwise use defaults
                ew_vmin = getattr(spec_object, 'vmin', -200)
                ew_vmax = getattr(spec_object, 'vmax', 200)
                
                # Create template parameters
                template_params = {
                    'filename': filepath,
                    'redshift': redshift,
                    'transition': transition,
                    'transition_name': transition_name,
                    'slice_vmin': slice_vmin,
                    'slice_vmax': slice_vmax,
                    'ew_vmin': int(ew_vmin),
                    'ew_vmax': int(ew_vmax),
                    'linelist': linelist,
                    'method': 'closest'
                }
                
                # Add to master table
                row_index = self.controller.master_table.add_item(template_params)
                
                # Check if this rb_spec object has EW results (already analyzed)
                if (hasattr(spec_object, 'W') and hasattr(spec_object, 'W_e') and 
                    not np.isnan(spec_object.W) and not np.isnan(spec_object.W_e)):
                    
                    # Extract all results from the rb_spec object
                    results = {
                        'W': spec_object.W,
                        'W_e': spec_object.W_e,
                        'N': getattr(spec_object, 'N', 0.0),
                        'N_e': getattr(spec_object, 'N_e', 0.0),
                        'logN': getattr(spec_object, 'logN', 0.0),
                        'logN_e': getattr(spec_object, 'logN_e', 0.0),
                        'vel_centroid': getattr(spec_object, 'vel_centroid', 0.0),
                        'vel_disp': getattr(spec_object, 'vel_disp', 0.0),
                        'SNR': getattr(spec_object, 'SNR', 0.0)
                    }
                    
                    # Update results in master table
                    self.controller.master_table.update_results(row_index, **results)
                    
                    # Extract continuum parameters if available
                    continuum_method = "polynomial"  # Default
                    continuum_order = 3  # Default
                    continuum_masks = []  # Default
                    
                    # Try to extract continuum fitting parameters from rb_spec object
                    if hasattr(spec_object, 'fit_params') and spec_object.fit_params:
                        fit_params = spec_object.fit_params
                        if 'Legendre' in fit_params:
                            continuum_order = fit_params['Legendre']
                        if 'mask' in fit_params and fit_params['mask']:
                            # Convert flat mask list to tuples: [v1,v2,v3,v4] -> [(v1,v2), (v3,v4)]
                            mask_list = fit_params['mask']
                            for i in range(0, len(mask_list), 2):
                                if i + 1 < len(mask_list):
                                    continuum_masks.append((mask_list[i], mask_list[i+1]))
                    
                    # Check if flat continuum was used (continuum is all 1s)
                    if hasattr(spec_object, 'cont') and len(spec_object.cont) > 0:
                        if abs(spec_object.cont.mean() - 1.0) < 0.01 and spec_object.cont.std() < 0.01:
                            continuum_method = "flat"
                    
                    # Update analysis parameters
                    self.controller.master_table.update_analysis(
                        row_index,
                        continuum_method=continuum_method,
                        continuum_order=continuum_order,
                        continuum_masks=continuum_masks,
                        processing_status="complete",
                        calculate_snr=True,  # Assume SNR was calculated if results exist
                        binsize=3,  # Default
                        use_weights=False,  # Default
                        last_modified=datetime.now().isoformat()
                    )
                    
                    # Store the rb_spec object
                    self.controller.master_table.set_rb_spec_object(row_index, spec_object)
                    
                    completed_count += 1
                    print(f"Imported complete analysis: {os.path.basename(filepath)}")
                    
                else:
                    # No EW results - mark as ready for processing
                    self.controller.master_table.update_analysis(
                        row_index,
                        processing_status="ready"
                    )
                    print(f"Imported for processing: {os.path.basename(filepath)}")
                
                imported_count += 1
                
            except Exception as e:
                print(f"Error importing {os.path.basename(filepath)}: {e}")
                skipped_count += 1
        
        # Refresh UI
        self.refresh_table()
        self.configuration_changed.emit()
        
        # Show summary message
        if completed_count > 0:
            summary = f"Import completed: {imported_count} files imported ({completed_count} already analyzed), {skipped_count} files skipped"
            
            # Auto-navigate to Review panel if we have completed analyses
            main_window = self.window()  # Get the top-level window

            if hasattr(main_window, 'review_panel'):
                main_window.review_panel.importing_in_progress = True
            if hasattr(main_window, 'tabs') and completed_count > 0:
                # Enable and switch to Review tab (index 2)
                main_window.update_tab_states()  # This will enable review tab
                main_window.tabs.setCurrentIndex(2)  # Switch to review panel
                summary += " - Switched to Review panel"
        else:
            summary = f"Import completed: {imported_count} files imported, {skipped_count} files skipped"
        
        self.controller.status_updated.emit(summary)
        # At the very end of import_multiple_json_files(), after the status message:
        if hasattr(main_window, 'review_panel'):
            main_window.review_panel.importing_in_progress = False
    

    def setup_table(self):
        """Setup table columns."""
        self.batch_table.setColumnCount(7)
        self.batch_table.setHorizontalHeaderLabels([
            "Filename", "Redshift", "Transition", "Slice Range", "EW Range", "Line List", "Status"
        ])
        
        # Set column widths
        header = self.batch_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)  # Filename
        for i in range(1, 7):
            header.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        
        self.batch_table.setRowCount(0)
    


    def add_batch_item(self):
        """Add an item to the batch processing list with smart defaults."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Single Ion")
        dialog.setMinimumSize(500, 450)
        layout = QVBoxLayout(dialog)
        
        # Get smart defaults from existing table entries
        default_values = self._get_smart_defaults()
        
        # Main form
        form_layout = QFormLayout()
        
        # Filename
        file_layout = QHBoxLayout()
        filename = QLineEdit()
        filename.setReadOnly(True)
        
        # Pre-populate with smart default
        if default_values['filename']:
            filename.setText(default_values['filename'])
        
        browse_btn = QPushButton("Browse")
        
        def browse_file():
            options = QFileDialog.Options()
            filepath, _ = QFileDialog.getOpenFileName(
                dialog, "Select Spectrum File", "", 
                "All Files (*);;FITS Files (*.fits);;JSON Files (*.json)",
                options=options
            )
            if filepath:
                filename.setText(filepath)
        
        browse_btn.clicked.connect(browse_file)
        file_layout.addWidget(filename)
        file_layout.addWidget(browse_btn)
        form_layout.addRow("Spectrum File:", file_layout)
        
        # Redshift
        redshift = QDoubleSpinBox()
        redshift.setDecimals(6)
        redshift.setRange(0.0, 10.0)
        redshift.setValue(default_values['redshift'])  # Smart default
        form_layout.addRow("Redshift:", redshift)
        
        # Transition
        transition_layout = QHBoxLayout()
        transition_combo = QComboBox()
        transitions = [
            "Custom",
            "Lyα (1215.67 Å)",
            "Lyβ (1025.72 Å)",
            "Lyγ (972.54 Å)",
            "CIV (1548.20 Å)",
            "CIV (1550.78 Å)",
            "MgII (2796.35 Å)",
            "MgII (2803.53 Å)",
            "SiIV (1393.76 Å)",
            "SiIV (1402.77 Å)"
        ]
        transition_combo.addItems(transitions)
        
        custom_wavelength = QDoubleSpinBox()
        custom_wavelength.setRange(1.0, 10000.0)
        custom_wavelength.setValue(1215.67)
        custom_wavelength.setDecimals(2)
        custom_wavelength.setEnabled(True)  # ✅ Changed from False
        custom_wavelength.setFocus()        # ✅ Put cursor here
        custom_wavelength.selectAll()       # ✅ Select default value
        
        transition_combo.currentTextChanged.connect(
            lambda text: custom_wavelength.setEnabled(text == "Custom")
        )
        
        transition_layout.addWidget(transition_combo)
        transition_layout.addWidget(custom_wavelength)
        form_layout.addRow("Transition:", transition_layout)
        
        # Line list selection
        linelist_combo = QComboBox()
        linelist_combo.addItems(["atom", "LLS", "LLS Small", "DLA", "LBG", "Gal"])
        # Set smart default for linelist
        if default_values['linelist']:
            linelist_index = linelist_combo.findText(default_values['linelist'])
            if linelist_index >= 0:
                linelist_combo.setCurrentIndex(linelist_index)
        form_layout.addRow("Line List:", linelist_combo)
        
        # Add separator
        separator = QLabel()
        separator.setStyleSheet("QLabel { border-bottom: 1px solid gray; }")
        form_layout.addRow(separator)
        
        # Velocity ranges section
        velocity_label = QLabel("<b>Velocity Ranges:</b>")
        form_layout.addRow(velocity_label)
        
        # Add explanatory text
        explanation = QLabel(
            "• Slice Range: Wide range for extracting spectrum (includes continuum regions)\n"
            "• EW Range: Narrow range for equivalent width integration (around absorption)"
        )
        explanation.setWordWrap(True)
        explanation.setStyleSheet("QLabel { color: gray; font-style: italic; margin-bottom: 10px; }")
        form_layout.addRow(explanation)
        
        # Slice velocity range with smart defaults
        slice_range_layout = QHBoxLayout()
        slice_vmin = QSpinBox()
        slice_vmin.setRange(-5000, 0)
        slice_vmin.setValue(default_values['slice_vmin'])  # Smart default
        slice_vmin.setSuffix(" km/s")
        
        slice_vmax = QSpinBox()
        slice_vmax.setRange(0, 5000)
        slice_vmax.setValue(default_values['slice_vmax'])  # Smart default
        slice_vmax.setSuffix(" km/s")
        
        slice_range_layout.addWidget(QLabel("Min:"))
        slice_range_layout.addWidget(slice_vmin)
        slice_range_layout.addWidget(QLabel("Max:"))
        slice_range_layout.addWidget(slice_vmax)
        
        form_layout.addRow("Slice Range:", slice_range_layout)
        
        # EW measurement range with smart defaults
        ew_range_layout = QHBoxLayout()
        ew_vmin = QSpinBox()
        ew_vmin.setRange(-2000, 0)
        ew_vmin.setValue(default_values['ew_vmin'])  # Smart default
        ew_vmin.setSuffix(" km/s")
        
        ew_vmax = QSpinBox()
        ew_vmax.setRange(0, 2000)
        ew_vmax.setValue(default_values['ew_vmax'])  # Smart default
        ew_vmax.setSuffix(" km/s")
        
        ew_range_layout.addWidget(QLabel("Min:"))
        ew_range_layout.addWidget(ew_vmin)
        ew_range_layout.addWidget(QLabel("Max:"))
        ew_range_layout.addWidget(ew_vmax)
        
        form_layout.addRow("EW Range:", ew_range_layout)
        
        # Validation function and label
        validation_label = QLabel("")
        form_layout.addRow(validation_label)
        
        def validate_ranges():
            """Validate that EW range is within slice range."""
            slice_min, slice_max = slice_vmin.value(), slice_vmax.value()
            ew_min, ew_max = ew_vmin.value(), ew_vmax.value()
            
            # Check if EW range is within slice range (with tolerance)
            tolerance = 50  # km/s
            if ew_min < (slice_min - tolerance) or ew_max > (slice_max + tolerance):
                validation_label.setText(f"⚠️ EW range should be within slice range ± {tolerance} km/s")
                validation_label.setStyleSheet("QLabel { color: red; }")
                return False
            else:
                validation_label.setText("✓ Velocity ranges are valid")
                validation_label.setStyleSheet("QLabel { color: green; }")
                return True
        
        # Connect validation to value changes
        slice_vmin.valueChanged.connect(validate_ranges)
        slice_vmax.valueChanged.connect(validate_ranges)
        ew_vmin.valueChanged.connect(validate_ranges)
        ew_vmax.valueChanged.connect(validate_ranges)
        
        # Initial validation
        validate_ranges()
        
        layout.addLayout(form_layout)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        
        def accept_with_validation():
            if not validate_ranges():
                QMessageBox.warning(dialog, "Validation Error", 
                                  "Please fix the velocity range validation errors before proceeding.")
                return
            dialog.accept()
        
        button_box.accepted.connect(accept_with_validation)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        

        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            # Get values
            file_path = filename.text()
            if not file_path:
                QMessageBox.warning(self, "Input Error", "Please select a spectrum file.")
                return
                
            z = redshift.value()
            
            # Get transition details - UNIFIED rb_setline approach
            trans_text = transition_combo.currentText()
            linelist = linelist_combo.currentText()
            
            # Extract wavelength from either dropdown or custom input
            if trans_text == "Custom":
                trans_wavelength = custom_wavelength.value()
            else:
                # Extract wavelength from dropdown selection
                trans_wavelength = float(trans_text.split("(")[1].split(" Å")[0])
            
            # Always use rb_setline to get proper transition name from database
            try:
                from rbcodes.IGM import rb_setline as s
                line_info = s.rb_setline(trans_wavelength, 'closest', linelist=linelist)
                trans_name = line_info['name']
                actual_wavelength = line_info['wave']
                
                print(f"Input wavelength {trans_wavelength:.2f} Å matched to: {trans_name} ({actual_wavelength:.2f} Å)")
                
            except Exception as e:
                print(f"Warning: Could not find transition in database: {e}")
                # Fallback to generic name if rb_setline fails
                if trans_text == "Custom":
                    trans_name = f"Unknown ({trans_wavelength:.2f} Å)"
                else:
                    trans_name = trans_text.split(" (")[0]
            
            # Get velocity ranges
            slice_min, slice_max = slice_vmin.value(), slice_vmax.value()
            ew_min, ew_max = ew_vmin.value(), ew_vmax.value()
            
            # Create template parameters as dict (not object)
            template_params = {
                'filename': file_path,
                'redshift': z,
                'transition': trans_wavelength,
                'transition_name': trans_name,
                'slice_vmin': slice_min,
                'slice_vmax': slice_max,
                'ew_vmin': ew_min,
                'ew_vmax': ew_max,
                'linelist': linelist,
                'method': 'closest'
            }
            
            # Add to master table (returns row index)
            row_index = self.controller.master_table.add_item(template_params)
            
            # Refresh display
            self.refresh_table()
            self.configuration_changed.emit()

    def _get_smart_defaults(self):
        """Get smart default values from existing table entries."""
        # Default values (used if no existing entries)
        defaults = {
            'filename': '',
            'redshift': 0.0,
            'linelist': 'atom',
            'slice_vmin': -1500,
            'slice_vmax': 1500,
            'ew_vmin': -200,
            'ew_vmax': 200
        }
        
        # Get existing items from master table
        item_count = self.controller.master_table.get_item_count()
        
        if item_count > 0:
            # Use the last added item as template
            last_item = self.controller.master_table.get_item(item_count - 1)
            if last_item:
                defaults.update({
                    'filename': last_item.template.filename,
                    'redshift': last_item.template.redshift,
                    'linelist': last_item.template.linelist,
                    'slice_vmin': last_item.template.slice_vmin,
                    'slice_vmax': last_item.template.slice_vmax,
                    'ew_vmin': last_item.template.ew_vmin,
                    'ew_vmax': last_item.template.ew_vmax
                })
                
                print(f"Using smart defaults from: {last_item.template.transition_name}")
                print(f"  File: {os.path.basename(last_item.template.filename)}")
                print(f"  Redshift: {last_item.template.redshift}")
        
        return defaults  
    
    def refresh_table(self):
        """Refresh the table display from the master table."""
        # Clear existing rows
        self.batch_table.setRowCount(0)
        
        # Get all items from master table
        items = self.controller.master_table.get_all_items()
        
        # Populate table
        for item in items:
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            # Use basename for display
            basename = os.path.basename(item.template.filename)
            
            self.batch_table.setItem(row, 0, QTableWidgetItem(basename))
            self.batch_table.setItem(row, 1, QTableWidgetItem(f"{item.template.redshift:.6f}"))
            self.batch_table.setItem(row, 2, QTableWidgetItem(f"{item.template.transition_name} ({item.template.transition:.2f} Å)"))
            self.batch_table.setItem(row, 3, QTableWidgetItem(f"[{item.template.slice_vmin}, {item.template.slice_vmax}]"))
            self.batch_table.setItem(row, 4, QTableWidgetItem(f"[{item.template.ew_vmin}, {item.template.ew_vmax}]"))
            self.batch_table.setItem(row, 5, QTableWidgetItem(item.template.linelist))
            
            # Status with color coding
            status = item.analysis.processing_status
            status_item = QTableWidgetItem(status.title())
            
            if status == 'complete':
                status_item.setBackground(QColor(200, 255, 200))  # Light green
            elif status == 'processing':
                status_item.setBackground(QColor(255, 255, 200))  # Light yellow
            elif status == 'error':
                status_item.setBackground(QColor(255, 200, 200))  # Light red
            elif status in ['needs_processing', 'needs_ew_recalc']:
                status_item.setBackground(QColor(255, 240, 200))  # Light orange
            
            self.batch_table.setItem(row, 6, status_item)
    
    def remove_batch_item(self):
        """Remove the selected item from the batch."""
        selected_rows = set(index.row() for index in self.batch_table.selectedIndexes())
        
        if not selected_rows:
            QMessageBox.warning(self, "No Selection", "Please select items to remove.")
            return
        
        
        # Remove from bottom to top to avoid index shifting
        for row in sorted(selected_rows, reverse=True):
            if 0 <= row < self.controller.master_table.get_item_count():
                success = self.controller.master_table.remove_item(row)        
        self.refresh_table()
        self.configuration_changed.emit()
    
    def clear_batch_items(self):
        """Clear all batch items."""
        if self.controller.master_table.get_item_count() > 0:
            reply = QMessageBox.question(
                self, 'Clear All Items',
                'Are you sure you want to clear all batch items?',
                QMessageBox.Yes | QMessageBox.No
            )
            
            if reply == QMessageBox.Yes:
                self.controller.master_table.clear_all()
                self.refresh_table()
                self.configuration_changed.emit()
    
    def import_from_csv(self):
        """Import batch items from a CSV file."""
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Import Batch Items", "", 
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )
        
        if not filepath:
            return
        
        success, message = self.controller.master_table.import_template_csv(filepath)
        
        if success:
            self.refresh_table()
            self.configuration_changed.emit()
            QMessageBox.information(self, "Import Successful", message)
        else:
            QMessageBox.warning(self, "Import Error", f"Failed to import batch items:\n{message}")
    
    def export_template_csv(self):
        """Export current configuration as CSV template."""
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Export Template CSV", "batch_template.csv", 
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )
        
        if not filepath:
            return
        
        # Add .csv extension if not present
        if not filepath.lower().endswith('.csv'):
            filepath += '.csv'
        
        success = self.controller.master_table.export_template_csv(filepath)
        
        if success:
            QMessageBox.information(
                self, "Export Successful", 
                f"Template exported to {filepath}\n\nYou can edit this file and import it back using 'Import from CSV'."
            )
        else:
            QMessageBox.warning(
                self, "Export Error",
                "Failed to export template CSV."
            )
    
    def save_configuration(self):
        """Save the current batch configuration to a file."""
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save Batch Configuration", "", 
            "JSON Files (*.json);;All Files (*)",
            options=options
        )
        
        if not filepath:
            return
        
        # Add .json extension if not present
        if not filepath.lower().endswith('.json'):
            filepath += '.json'
        
        success = self.controller.save_batch_configuration(filepath)
        
        if success:
            QMessageBox.information(
                self, "Save Successful", 
                f"Batch configuration saved to {filepath}"
            )
        else:
            QMessageBox.warning(
                self, "Save Error",
                "Failed to save batch configuration."
            )
    
    def load_configuration(self):
        """Load a batch configuration from a file."""
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Load Batch Configuration", "", 
            "JSON Files (*.json);;All Files (*)",
            options=options
        )
        
        if not filepath:
            return
        
        success = self.controller.load_batch_configuration(filepath)
        
        if success:
            self.refresh_table()
            self.configuration_changed.emit()
            
            # FORCE REVIEW PANEL REFRESH - ADD THIS
            # Get the main window and refresh review panel
            main_window = self.window()  # Get the top-level window
            
            # Clear import flag
            if hasattr(main_window, 'review_panel'):
                main_window.review_panel.importing_in_progress = False
                # Force a final refresh now that import is complete
                main_window.review_panel.refresh_results_table()
            
            QMessageBox.information(
                self, "Load Successful", 
                f"Batch configuration loaded from {filepath}"
            )
        else:
            QMessageBox.warning(
                self, "Load Error",
                "Failed to load batch configuration."
            )
    
    def create_csv_template(self):
        """Create an empty CSV template file."""
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Create CSV Template", "batch_template.csv", 
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )
        
        if not filepath:
            return
        
        # Add .csv extension if not present
        if not filepath.lower().endswith('.csv'):
            filepath += '.csv'
        
        try:
            with open(filepath, 'w', newline='') as csvfile:
                fieldnames = ['filename', 'redshift', 'transition', 'transition_name', 
                             'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax', 'linelist']
                
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                # Add example rows
                examples = [
                    {
                        'filename': '/path/to/spectrum1.fits',
                        'redshift': 0.5,
                        'transition': 1215.67,
                        'transition_name': 'Lyα',
                        'slice_vmin': -1500,
                        'slice_vmax': 1500,
                        'ew_vmin': -200,
                        'ew_vmax': 200,
                        'linelist': 'atom'
                    },
                    {
                        'filename': '/path/to/spectrum1.fits',
                        'redshift': 0.5,
                        'transition': 1025.72,
                        'transition_name': 'Lyβ',
                        'slice_vmin': -1000,
                        'slice_vmax': 1000,
                        'ew_vmin': -150,
                        'ew_vmax': 150,
                        'linelist': 'atom'
                    }
                ]
                
                for example in examples:
                    writer.writerow(example)
            
            QMessageBox.information(
                self, "Template Created", 
                f"CSV template saved to {filepath}\n\n"
                f"Edit this file with your batch items and use 'Import from CSV' to load them.\n"
                f"The template includes example rows to show the format."
            )
            
        except Exception as e:
            QMessageBox.warning(
                self, "Template Error",
                f"Failed to create CSV template: {str(e)}"
            )
    
    def validate_configuration(self):
        """Validate the current batch configuration."""
        valid_ids, invalid_items = self.controller.master_table.validate_all()
        
        if not invalid_items:
            QMessageBox.information(
                self, "Validation Result",
                f"All {len(valid_ids)} batch items are valid and ready for processing."
            )
        else:
            error_details = "\n".join([f"Item {i+1}: {msg}" for i, (item_id, msg) in enumerate(invalid_items)])
            QMessageBox.warning(
                self, "Validation Errors",
                f"Found {len(invalid_items)} invalid items:\n\n{error_details}"
            )
    
    def get_item_count(self):
        """Get the current number of items."""
        return self.controller.master_table.get_item_count()
    
    def get_selected_item_ids(self):
        """Get IDs of currently selected items."""
        selected_rows = set(index.row() for index in self.batch_table.selectedIndexes())
        items = self.controller.master_table.get_all_items()
        
        selected_ids = []
        for row in selected_rows:
            if 0 <= row < len(items):
                selected_ids.append(items[row].id)
        
        return selected_ids