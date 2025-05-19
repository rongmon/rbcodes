# rbcodes/GUIs/specgui/batch/panels/configuration_panel.py
import os
import csv
import json
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QTableWidget, QTableWidgetItem,
                           QHeaderView, QComboBox, QGroupBox, QFormLayout,
                           QMessageBox, QFileDialog, QCheckBox, QSpinBox,
                           QDoubleSpinBox, QDialog, QDialogButtonBox, QLineEdit,QAbstractItemView)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor

class ConfigurationPanel(QWidget):
    """Panel for configuring batch processing items."""
    
    # Signals
    configuration_changed = pyqtSignal(list)  # Emitted when batch configuration changes
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.batch_items = []  # List to store batch items
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add instructions
        instructions = QLabel(
            "Configure batch processing by adding items to analyze. "
            "You can add multiple systems, transitions, or files to process in batch."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)
        
        # Batch mode selection
        mode_group = QGroupBox("Batch Mode")
        mode_layout = QHBoxLayout()
        
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Multiple Absorption Systems (Same File)", 
            "Multiple Transitions (Same Redshift)", 
            "Multiple Files (Single Transition Each)"
        ])
        self.mode_combo.currentIndexChanged.connect(self.update_ui_for_mode)
        
        mode_layout.addWidget(QLabel("Mode:"))
        mode_layout.addWidget(self.mode_combo)
        
        mode_group.setLayout(mode_layout)
        main_layout.addWidget(mode_group)
        
        # Batch items table
        self.table_group = QGroupBox("Batch Items")
        table_layout = QVBoxLayout()
        
        self.batch_table = QTableWidget()
        self.batch_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setup_table_for_mode(0)  # Initial setup for first mode
        
        # Buttons for managing batch items
        table_buttons = QHBoxLayout()
        
        self.add_btn = QPushButton("Add Item")
        self.add_btn.clicked.connect(self.add_batch_item)
        
        self.import_btn = QPushButton("Import from CSV")
        self.import_btn.clicked.connect(self.import_from_csv)
        
        self.remove_btn = QPushButton("Remove Selected")
        self.remove_btn.clicked.connect(self.remove_batch_item)
        
        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.clicked.connect(self.clear_batch_items)
        
        table_buttons.addWidget(self.add_btn)
        table_buttons.addWidget(self.import_btn)
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
    
    def setup_table_for_mode(self, mode_index):
        """Setup table columns based on the selected batch mode."""
        self.batch_table.clear()
        
        if mode_index == 0:  # Multiple systems
            self.batch_table.setColumnCount(5)
            self.batch_table.setHorizontalHeaderLabels([
                "Filename", "Redshift", "Transition", "Velocity Min", "Velocity Max"
            ])
        elif mode_index == 1:  # Multiple transitions
            self.batch_table.setColumnCount(5)
            self.batch_table.setHorizontalHeaderLabels([
                "Filename", "Redshift", "Transition", "Velocity Min", "Velocity Max"
            ])
        else:  # Multiple files
            self.batch_table.setColumnCount(5)
            self.batch_table.setHorizontalHeaderLabels([
                "Filename", "Redshift", "Transition", "Velocity Min", "Velocity Max"
            ])
        
        # Set column widths
        header = self.batch_table.horizontalHeader()
        for i in range(self.batch_table.columnCount()):
            header.setSectionResizeMode(i, QHeaderView.Stretch)
        
        self.batch_table.setRowCount(0)
    
    def update_ui_for_mode(self, mode_index):
        """Update UI elements based on the selected batch mode."""
        self.setup_table_for_mode(mode_index)
        self.clear_batch_items()  # Clear items when mode changes
    
    def add_batch_item(self):
        """Add an item to the batch processing list."""
        mode_index = self.mode_combo.currentIndex()
        
        if mode_index == 0:  # Multiple systems
            self.add_system_dialog()
        elif mode_index == 1:  # Multiple transitions
            self.add_transition_dialog()
        else:  # Multiple files
            self.add_file_dialog()
    
    def add_system_dialog(self):
        """Open dialog to add a multiple systems batch item."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Absorption System")
        layout = QVBoxLayout(dialog)
        
        # Form for system details
        form = QFormLayout()
        
        # Filename
        file_layout = QHBoxLayout()
        filename = QLineEdit()
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
        form.addRow("Spectrum File:", file_layout)
        
        # Redshift
        redshift = QDoubleSpinBox()
        redshift.setDecimals(6)
        redshift.setRange(0.0, 10.0)
        redshift.setValue(0.0)
        form.addRow("Redshift:", redshift)
        
        # Transition
        transition_layout = QHBoxLayout()
        transition_combo = QComboBox()
        transitions = [
            "Lyα (1215.67 Å)",
            "Lyβ (1025.72 Å)",
            "Lyγ (972.54 Å)",
            "CIV (1548.20 Å)",
            "CIV (1550.78 Å)",
            "MgII (2796.35 Å)",
            "MgII (2803.53 Å)",
            "SiIV (1393.76 Å)",
            "SiIV (1402.77 Å)",
            "Custom"
        ]
        transition_combo.addItems(transitions)
        
        custom_wavelength = QDoubleSpinBox()
        custom_wavelength.setRange(1.0, 10000.0)
        custom_wavelength.setValue(1215.67)
        custom_wavelength.setDecimals(2)
        custom_wavelength.setEnabled(False)
        
        transition_combo.currentTextChanged.connect(
            lambda text: custom_wavelength.setEnabled(text == "Custom")
        )
        
        transition_layout.addWidget(transition_combo)
        transition_layout.addWidget(custom_wavelength)
        form.addRow("Transition:", transition_layout)
        
        # Velocity range
        velocity_layout = QHBoxLayout()
        vmin = QSpinBox()
        vmin.setRange(-5000, 0)
        vmin.setValue(-200)
        
        vmax = QSpinBox()
        vmax.setRange(0, 5000)
        vmax.setValue(200)
        
        velocity_layout.addWidget(QLabel("Min:"))
        velocity_layout.addWidget(vmin)
        velocity_layout.addWidget(QLabel("Max:"))
        velocity_layout.addWidget(vmax)
        
        form.addRow("Velocity Range:", velocity_layout)
        
        layout.addLayout(form)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            # Get values
            file_path = filename.text()
            z = redshift.value()
            
            # Get transition details
            trans_text = transition_combo.currentText()
            if trans_text == "Custom":
                trans_wavelength = custom_wavelength.value()
                trans_name = f"Custom ({trans_wavelength:.2f} Å)"
            else:
                # Extract wavelength from the selected item
                trans_wavelength = float(trans_text.split("(")[1].split(" Å")[0])
                trans_name = trans_text.split(" (")[0]
            
            # Add to table
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            self.batch_table.setItem(row, 0, QTableWidgetItem(file_path))
            self.batch_table.setItem(row, 1, QTableWidgetItem(f"{z:.6f}"))
            self.batch_table.setItem(row, 2, QTableWidgetItem(f"{trans_name} ({trans_wavelength:.2f} Å)"))
            self.batch_table.setItem(row, 3, QTableWidgetItem(f"{vmin.value()}"))
            self.batch_table.setItem(row, 4, QTableWidgetItem(f"{vmax.value()}"))
            
            # Store item data
            self.batch_items.append({
                'type': 'multiple_systems',
                'filename': file_path,
                'redshift': z,
                'transition': trans_wavelength,
                'transition_name': trans_name,
                'vmin': vmin.value(),
                'vmax': vmax.value()
            })
            
            # Emit signal that configuration has changed
            self.configuration_changed.emit(self.batch_items)
    
    def add_transition_dialog(self):
        """Open dialog to add a transition to the batch."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Transition")
        layout = QVBoxLayout(dialog)
        
        # Form for transition details
        form = QFormLayout()
        
        # Filename
        file_layout = QHBoxLayout()
        filename = QLineEdit()
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
        form.addRow("Spectrum File:", file_layout)
        
        # Redshift
        redshift = QDoubleSpinBox()
        redshift.setDecimals(6)
        redshift.setRange(0.0, 10.0)
        redshift.setValue(0.0)
        form.addRow("Redshift:", redshift)
        
        # Transition
        transition_layout = QHBoxLayout()
        transition_combo = QComboBox()
        transitions = [
            "Lyα (1215.67 Å)",
            "Lyβ (1025.72 Å)",
            "Lyγ (972.54 Å)",
            "CIV (1548.20 Å)",
            "CIV (1550.78 Å)",
            "MgII (2796.35 Å)",
            "MgII (2803.53 Å)",
            "SiIV (1393.76 Å)",
            "SiIV (1402.77 Å)",
            "Custom"
        ]
        transition_combo.addItems(transitions)
        
        custom_wavelength = QDoubleSpinBox()
        custom_wavelength.setRange(1.0, 10000.0)
        custom_wavelength.setValue(1215.67)
        custom_wavelength.setDecimals(2)
        custom_wavelength.setEnabled(False)
        
        transition_combo.currentTextChanged.connect(
            lambda text: custom_wavelength.setEnabled(text == "Custom")
        )
        
        transition_layout.addWidget(transition_combo)
        transition_layout.addWidget(custom_wavelength)
        form.addRow("Transition:", transition_layout)
        
        # Line list selection
        linelist_combo = QComboBox()
        linelist_combo.addItems(["atom", "LLS", "LLS Small", "DLA", "LBG", "Gal"])
        form.addRow("Line List:", linelist_combo)
        
        # Velocity range
        velocity_layout = QHBoxLayout()
        vmin = QSpinBox()
        vmin.setRange(-5000, 0)
        vmin.setValue(-200)
        
        vmax = QSpinBox()
        vmax.setRange(0, 5000)
        vmax.setValue(200)
        
        velocity_layout.addWidget(QLabel("Min:"))
        velocity_layout.addWidget(vmin)
        velocity_layout.addWidget(QLabel("Max:"))
        velocity_layout.addWidget(vmax)
        
        form.addRow("Velocity Range:", velocity_layout)
        
        layout.addLayout(form)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            # Get values
            file_path = filename.text()
            z = redshift.value()
            
            # Get transition details
            trans_text = transition_combo.currentText()
            if trans_text == "Custom":
                trans_wavelength = custom_wavelength.value()
            else:
                # Extract wavelength from the selected item
                trans_wavelength = float(trans_text.split("(")[1].split(" Å")[0])
            
            # Use rb_setline to look up the transition
            try:
                # Import rb_setline
                try:
                    from rbcodes.IGM import rb_setline as s
                except ImportError:
                    try:
                        from IGM import rb_setline as s
                    except ImportError:
                        QMessageBox.warning(self, "Import Error", 
                                         "Could not import rb_setline. Using entered     wavelength without lookup.")
                        s = None
                
                # Get line list
                linelist = linelist_combo.currentText()
                
                # Get transition information if rb_setline is available
                if s:
                    try:
                        transition_info = s.rb_setline(trans_wavelength, 'closest',     linelist=linelist)
                        trans_wavelength = float(transition_info['wave'])
                        trans_name = transition_info['name']
                        fval = transition_info['fval']
                    except Exception as e:
                        print(f"Error in rb_setline: {str(e)}")
                        # Fallback if rb_setline lookup fails
                        if trans_text == "Custom":
                            trans_name = f"Custom ({trans_wavelength:.2f} Å)"
                        else:
                            trans_name = trans_text.split(" (")[0]
                        fval = None
                else:
                    # Fallback if rb_setline is not available
                    if trans_text == "Custom":
                        trans_name = f"Custom ({trans_wavelength:.2f} Å)"
                    else:
                        trans_name = trans_text.split(" (")[0]
                    fval = None
                    
            except Exception as e:
                QMessageBox.warning(self, "Transition Lookup Error", 
                                  f"Error looking up transition: {str(e)}\nUsing     entered wavelength.")
                if trans_text == "Custom":
                    trans_name = f"Custom ({trans_wavelength:.2f} Å)"
                else:
                    trans_name = trans_text.split(" (")[0]
                fval = None
            
            # Add to table
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            self.batch_table.setItem(row, 0, QTableWidgetItem(file_path))
            self.batch_table.setItem(row, 1, QTableWidgetItem(f"{z:.6f}"))
            self.batch_table.setItem(row, 2, QTableWidgetItem(f"{trans_name}     ({trans_wavelength:.2f} Å)"))
            self.batch_table.setItem(row, 3, QTableWidgetItem(f"{vmin.value()}"))
            self.batch_table.setItem(row, 4, QTableWidgetItem(f"{vmax.value()}"))
            
            # Store item data
            self.batch_items.append({
                'type': 'multiple_transitions',
                'filename': file_path,
                'redshift': z,
                'transition': trans_wavelength,
                'transition_name': trans_name,
                'fval': fval,  # Store the oscillator strength
                'linelist': linelist,  # Store the line list
                'vmin': vmin.value(),
                'vmax': vmax.value()
            })
            
            # Emit signal that configuration has changed
            self.configuration_changed.emit(self.batch_items)    
    
    def add_file_dialog(self):
        """Open dialog to add a multiple files batch item."""
        # First select the file
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select Spectrum File", "", 
            "All Files (*);;FITS Files (*.fits);;JSON Files (*.json)",
            options=options
        )
        
        if not filename:
            return
        
        # Now open dialog for additional parameters
        dialog = QDialog(self)
        dialog.setWindowTitle("File Parameters")
        layout = QVBoxLayout(dialog)
        
        # Form for file parameters
        form = QFormLayout()
        
        # Redshift
        redshift = QDoubleSpinBox()
        redshift.setDecimals(6)
        redshift.setRange(0.0, 10.0)
        redshift.setValue(0.0)
        form.addRow("Redshift:", redshift)
        
        # Transition
        transition_layout = QHBoxLayout()
        transition_combo = QComboBox()
        transitions = [
            "Lyα (1215.67 Å)",
            "Lyβ (1025.72 Å)",
            "Lyγ (972.54 Å)",
            "CIV (1548.20 Å)",
            "CIV (1550.78 Å)",
            "MgII (2796.35 Å)",
            "MgII (2803.53 Å)",
            "SiIV (1393.76 Å)",
            "SiIV (1402.77 Å)",
            "Custom"
        ]
        transition_combo.addItems(transitions)
        
        custom_wavelength = QDoubleSpinBox()
        custom_wavelength.setRange(1.0, 10000.0)
        custom_wavelength.setValue(1215.67)
        custom_wavelength.setDecimals(2)
        custom_wavelength.setEnabled(False)
        
        transition_combo.currentTextChanged.connect(
            lambda text: custom_wavelength.setEnabled(text == "Custom")
        )
        
        transition_layout.addWidget(transition_combo)
        transition_layout.addWidget(custom_wavelength)
        form.addRow("Transition:", transition_layout)
        
        # Velocity range
        velocity_layout = QHBoxLayout()
        vmin = QSpinBox()
        vmin.setRange(-5000, 0)
        vmin.setValue(-200)
        
        vmax = QSpinBox()
        vmax.setRange(0, 5000)
        vmax.setValue(200)
        
        velocity_layout.addWidget(QLabel("Min:"))
        velocity_layout.addWidget(vmin)
        velocity_layout.addWidget(QLabel("Max:"))
        velocity_layout.addWidget(vmax)
        
        form.addRow("Velocity Range:", velocity_layout)
        
        layout.addLayout(form)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            # Get values
            z = redshift.value()
            
            # Get transition details
            trans_text = transition_combo.currentText()
            if trans_text == "Custom":
                trans_wavelength = custom_wavelength.value()
                trans_name = f"Custom ({trans_wavelength:.2f} Å)"
            else:
                # Extract wavelength from the selected item
                trans_wavelength = float(trans_text.split("(")[1].split(" Å")[0])
                trans_name = trans_text.split(" (")[0]
            
            # Add to table
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            # Use basename for display but store full path
            basename = os.path.basename(filename)
            
            self.batch_table.setItem(row, 0, QTableWidgetItem(basename))
            self.batch_table.setItem(row, 1, QTableWidgetItem(f"{z:.6f}"))
            self.batch_table.setItem(row, 2, QTableWidgetItem(f"{trans_name} ({trans_wavelength:.2f} Å)"))
            self.batch_table.setItem(row, 3, QTableWidgetItem(f"{vmin.value()}"))
            self.batch_table.setItem(row, 4, QTableWidgetItem(f"{vmax.value()}"))
            
            # Store item data
            self.batch_items.append({
                'type': 'multiple_files',
                'filename': filename,  # Store full path
                'redshift': z,
                'transition': trans_wavelength,
                'transition_name': trans_name,
                'vmin': vmin.value(),
                'vmax': vmax.value()
            })
            
            # Emit signal that configuration has changed
            self.configuration_changed.emit(self.batch_items)
    
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
        
        try:
            with open(filepath, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                
                if not reader.fieldnames:
                    raise ValueError("CSV file has no headers")
                
                # Check for required fields based on mode
                mode_index = self.mode_combo.currentIndex()
                required_fields = ['filename', 'redshift', 'transition', 'vmin', 'vmax']
                
                # Ensure all required fields are in the CSV
                missing_fields = [field for field in required_fields if field not in reader.fieldnames]
                if missing_fields:
                    raise ValueError(f"Missing required fields: {', '.join(missing_fields)}")
                
                # Clear existing items if we're importing a new set
                self.clear_batch_items()
                
                # Add items from CSV
                for row in reader:
                    try:
                        # Parse row data
                        filename = row['filename']
                        redshift = float(row['redshift'])
                        transition = float(row['transition'])
                        vmin = int(row['vmin'])
                        vmax = int(row['vmax'])
                        
                        # Get transition name if available, otherwise use a generic name
                        if 'transition_name' in row and row['transition_name']:
                            transition_name = row['transition_name']
                        else:
                            transition_name = f"λ {transition:.2f}"
                        
                        # Add to table
                        table_row = self.batch_table.rowCount()
                        self.batch_table.insertRow(table_row)
                        
                        # Use basename for display but store full path
                        basename = os.path.basename(filename)
                        
                        self.batch_table.setItem(table_row, 0, QTableWidgetItem(basename))
                        self.batch_table.setItem(table_row, 1, QTableWidgetItem(f"{redshift:.6f}"))
                        self.batch_table.setItem(table_row, 2, QTableWidgetItem(f"{transition_name} ({transition:.2f} Å)"))
                        self.batch_table.setItem(table_row, 3, QTableWidgetItem(f"{vmin}"))
                        self.batch_table.setItem(table_row, 4, QTableWidgetItem(f"{vmax}"))
                        
                        # Determine item type based on mode
                        if mode_index == 0:
                            item_type = 'multiple_systems'
                        elif mode_index == 1:
                            item_type = 'multiple_transitions'
                        else:
                            item_type = 'multiple_files'
                        
                        # Store item data
                        self.batch_items.append({
                            'type': item_type,
                            'filename': filename,
                            'redshift': redshift,
                            'transition': transition,
                            'transition_name': transition_name,
                            'vmin': vmin,
                            'vmax': vmax
                        })
                    except (ValueError, KeyError) as e:
                        # Skip invalid rows but continue processing
                        print(f"Error parsing row: {e}")
                
                # Emit signal that configuration has changed
                self.configuration_changed.emit(self.batch_items)
                
                QMessageBox.information(
                    self, "Import Successful", 
                    f"Imported {len(self.batch_items)} batch items from {filepath}"
                )
        except Exception as e:
            QMessageBox.warning(
                self, "Import Error",
                f"Failed to import batch items: {str(e)}"
            )
    
    def remove_batch_item(self):
        """Remove the selected item from the batch."""
        selected_rows = set(index.row() for index in self.batch_table.selectedIndexes())
        
        # Remove from bottom to top to avoid index shifting
        for row in sorted(selected_rows, reverse=True):
            self.batch_table.removeRow(row)
            if 0 <= row < len(self.batch_items):
                self.batch_items.pop(row)
        
        # Emit signal that configuration has changed
        self.configuration_changed.emit(self.batch_items)
    
    def clear_batch_items(self):
        """Clear all batch items."""
        self.batch_table.setRowCount(0)
        self.batch_items = []
        
        # Emit signal that configuration has changed
        self.configuration_changed.emit(self.batch_items)
    
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
            # Update UI with loaded batch items
            self.update_ui_with_config(self.controller.batch_items)
            QMessageBox.information(
                self, "Load Successful", 
                f"Batch configuration loaded from {filepath}"
            )
        else:
            QMessageBox.warning(
                self, "Load Error",
                "Failed to load batch configuration."
            )
    
    def update_ui_with_config(self, batch_items):
        """Update the UI with loaded batch items."""
        # Clear existing items
        self.clear_batch_items()
        
        # Update items list
        self.batch_items = batch_items
        
        # Determine batch mode based on the first item
        if batch_items and 'type' in batch_items[0]:
            if batch_items[0]['type'] == 'multiple_systems':
                self.mode_combo.setCurrentIndex(0)
            elif batch_items[0]['type'] == 'multiple_transitions':
                self.mode_combo.setCurrentIndex(1)
            elif batch_items[0]['type'] == 'multiple_files':
                self.mode_combo.setCurrentIndex(2)
        
        # Update table
        for item in batch_items:
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            # Use basename for display but store full path
            basename = os.path.basename(item['filename'])
            
            self.batch_table.setItem(row, 0, QTableWidgetItem(basename))
            self.batch_table.setItem(row, 1, QTableWidgetItem(f"{item['redshift']:.6f}"))
            
            # Format transition display
            transition_display = f"{item.get('transition_name', 'Unknown')} ({item['transition']:.2f} Å)"
            self.batch_table.setItem(row, 2, QTableWidgetItem(transition_display))
            
            self.batch_table.setItem(row, 3, QTableWidgetItem(f"{item['vmin']}"))
            self.batch_table.setItem(row, 4, QTableWidgetItem(f"{item['vmax']}"))
        
        # Emit signal that configuration has changed
        self.configuration_changed.emit(self.batch_items)