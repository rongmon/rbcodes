# Modified batch_panel.py
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QTableWidget, QTableWidgetItem,
                           QHeaderView, QComboBox, QGroupBox, QFormLayout,
                           QMessageBox, QFileDialog, QCheckBox, QSpinBox,
                           QDoubleSpinBox, QAbstractItemView, QApplication,
                           QDialog, QDialogButtonBox, QLineEdit, QScrollArea,
                           QSplitter)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QBrush, QColor
import os
import numpy as np
import pandas as pd

class BatchPanel(QWidget):
    """Panel for batch processing multiple transitions or files."""
    
    # Signals
    batch_completed = pyqtSignal(list)  # Emitted when batch processing is complete
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.batch_items = []  # List to store batch items
        self.results = []      # List to store batch results
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add instructions
        instructions = QLabel(
            "This panel allows you to process multiple transitions or files in batch. "
            "You can add transitions to analyze, set parameters, and run the batch process."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)
        
        # Create a splitter for the main content
        main_splitter = QSplitter(Qt.Vertical)
        main_layout.addWidget(main_splitter, 1)  # Give it stretch factor
        
        # ==================== TOP SECTION: BATCH SETUP ====================
        top_widget = QWidget()
        top_layout = QVBoxLayout(top_widget)
        top_layout.setContentsMargins(0, 0, 0, 0)
        
        # Batch mode selection
        mode_group = QGroupBox("Batch Mode")
        mode_layout = QHBoxLayout()
        
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Multiple Transitions", "Multiple Files"])
        self.mode_combo.currentIndexChanged.connect(self.toggle_batch_mode)
        
        mode_layout.addWidget(QLabel("Mode:"))
        mode_layout.addWidget(self.mode_combo)
        
        mode_group.setLayout(mode_layout)
        top_layout.addWidget(mode_group)
        
        # Batch items section with scroll area
        table_scroll = QScrollArea()
        table_scroll.setWidgetResizable(True)
        table_scroll.setFrameShape(QScrollArea.NoFrame)
        
        self.table_group = QGroupBox("Batch Items")
        table_layout = QVBoxLayout()
        
        self.batch_table = QTableWidget()
        self.batch_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setup_transitions_table()  # Initial setup for transitions mode
        
        # Buttons for managing batch items
        table_buttons = QHBoxLayout()
        
        self.add_btn = QPushButton("Add Item")
        self.add_btn.clicked.connect(self.add_batch_item)
        
        self.load_file_btn = QPushButton("Load from File")
        self.load_file_btn.clicked.connect(self.load_from_file)
        
        self.export_template_btn = QPushButton("Export Template")
        self.export_template_btn.clicked.connect(self.export_template)
        
        self.remove_btn = QPushButton("Remove Selected")
        self.remove_btn.clicked.connect(self.remove_batch_item)
        
        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.clicked.connect(self.clear_batch_items)
        
        table_buttons.addWidget(self.add_btn)
        table_buttons.addWidget(self.load_file_btn)
        table_buttons.addWidget(self.export_template_btn)
        table_buttons.addWidget(self.remove_btn)
        table_buttons.addWidget(self.clear_btn)
        
        table_layout.addWidget(self.batch_table)
        table_layout.addLayout(table_buttons)
        
        self.table_group.setLayout(table_layout)
        table_scroll.setWidget(self.table_group)
        top_layout.addWidget(table_scroll)
        
        main_splitter.addWidget(top_widget)
        
        # ==================== MIDDLE SECTION: PARAMETERS ====================
        middle_widget = QWidget()
        middle_layout = QVBoxLayout(middle_widget)
        middle_layout.setContentsMargins(0, 0, 0, 0)
        
        # Parameters scroll area
        params_scroll = QScrollArea()
        params_scroll.setWidgetResizable(True)
        params_scroll.setFrameShape(QScrollArea.NoFrame)
        
        # Batch processing parameters
        params_group = QGroupBox("Processing Parameters")
        params_layout = QFormLayout()
        
        # Continuum fitting method
        self.cont_method = QComboBox()
        self.cont_method.addItems(["Polynomial", "Flat (Skip Fitting)"])
        params_layout.addRow("Continuum Method:", self.cont_method)
        
        # Polynomial order
        self.poly_order = QSpinBox()
        self.poly_order.setRange(1, 10)
        self.poly_order.setValue(3)
        params_layout.addRow("Polynomial Order:", self.poly_order)
        
        # Use weights
        self.use_weights = QCheckBox()
        params_layout.addRow("Use Weights:", self.use_weights)
        
        # Use default velocity range
        self.default_vrange = QCheckBox("Use Default")
        self.default_vrange.setChecked(True)
        self.default_vrange.toggled.connect(self.toggle_velocity_range)
        
        self.vmin = QSpinBox()
        self.vmin.setRange(-5000, 0)
        self.vmin.setValue(-200)
        self.vmin.setEnabled(False)
        
        self.vmax = QSpinBox()
        self.vmax.setRange(0, 5000)
        self.vmax.setValue(200)
        self.vmax.setEnabled(False)
        
        vrange_layout = QHBoxLayout()
        vrange_layout.addWidget(self.default_vrange)
        vrange_layout.addWidget(QLabel("Min:"))
        vrange_layout.addWidget(self.vmin)
        vrange_layout.addWidget(QLabel("Max:"))
        vrange_layout.addWidget(self.vmax)
        
        params_layout.addRow("Velocity Range:", vrange_layout)
        
        # Calculate SNR
        self.calc_snr = QCheckBox()
        self.calc_snr.setChecked(True)
        params_layout.addRow("Calculate SNR:", self.calc_snr)
        
        params_group.setLayout(params_layout)
        params_scroll.setWidget(params_group)
        middle_layout.addWidget(params_scroll)
        
        # Run batch and export buttons
        buttons_layout = QHBoxLayout()
        
        self.run_btn = QPushButton("Run Batch Process")
        self.run_btn.clicked.connect(self.run_batch)
        self.run_btn.setEnabled(False)  # Disabled until items are added
        
        self.export_btn = QPushButton("Export Results")
        self.export_btn.clicked.connect(self.export_results)
        self.export_btn.setEnabled(False)  # Disabled until batch is run
        
        self.advanced_btn = QPushButton("Launch Advanced Batch Mode")
        self.advanced_btn.clicked.connect(self.launch_advanced_mode)
        
        buttons_layout.addWidget(self.run_btn)
        buttons_layout.addWidget(self.export_btn)
        buttons_layout.addWidget(self.advanced_btn)
        
        middle_layout.addLayout(buttons_layout)
        main_splitter.addWidget(middle_widget)
        
        # ==================== BOTTOM SECTION: RESULTS ====================
        bottom_widget = QWidget()
        bottom_layout = QVBoxLayout(bottom_widget)
        bottom_layout.setContentsMargins(0, 0, 0, 0)
        
        # Results scroll area
        results_scroll = QScrollArea()
        results_scroll.setWidgetResizable(True)
        results_scroll.setFrameShape(QScrollArea.NoFrame)
        
        # Results group
        results_group = QGroupBox("Batch Results")
        results_layout = QVBoxLayout()
        
        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels(["Status", "Item", "Equivalent Width (Å)", "Error", "log(N)", "Error"])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.results_table.setRowCount(0)
        
        results_layout.addWidget(self.results_table)
        
        # Results summary
        self.results_label = QLabel("No batch processing results yet.")
        self.results_label.setWordWrap(True)
        results_layout.addWidget(self.results_label)
        
        results_group.setLayout(results_layout)
        results_scroll.setWidget(results_group)
        bottom_layout.addWidget(results_scroll)
        
        main_splitter.addWidget(bottom_widget)
        
        # Set initial splitter sizes (40% top, 30% middle, 30% bottom)
        main_splitter.setSizes([400, 300, 300])
        
        # Initial UI state
        self.toggle_batch_mode(0)  # Start with transitions mode
    
    def setup_transitions_table(self):
        """Setup table for transitions batch mode."""
        self.batch_table.clear()
        self.batch_table.setColumnCount(4)
        self.batch_table.setHorizontalHeaderLabels(["Transition", "Wavelength (Å)", "Description", "Status"])
        
        # Set column widths
        header = self.batch_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        
        self.batch_table.setRowCount(0)
        self.batch_items = []
    
    def setup_files_table(self):
        """Setup table for files batch mode."""
        self.batch_table.clear()
        self.batch_table.setColumnCount(4)
        self.batch_table.setHorizontalHeaderLabels(["Filename", "Transition", "Redshift", "Status"])
        
        # Set column widths
        header = self.batch_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        
        self.batch_table.setRowCount(0)
        self.batch_items = []
    
    def toggle_batch_mode(self, index):
        """Switch between transitions and files batch modes."""
        if index == 0:  # Transitions mode
            self.setup_transitions_table()
        else:  # Files mode
            self.setup_files_table()
        
        # Clear results
        self.results = []
        self.results_label.setText("No batch processing results yet.")
        self.results_table.setRowCount(0)
        self.export_btn.setEnabled(False)
    
    def toggle_velocity_range(self, use_default):
        """Enable/disable custom velocity range inputs."""
        self.vmin.setEnabled(not use_default)
        self.vmax.setEnabled(not use_default)
    
    def load_from_file(self):
        """Load batch items from a file."""
        options = QFileDialog.Options()
        file_filter = "All Files (*);;CSV Files (*.csv);;TSV Files (*.tsv);;Excel Files (*.xlsx);;JSON Files (*.json)"
        filename, selected_filter = QFileDialog.getOpenFileName(
            self, "Select Batch Input File", "", file_filter, options=options
        )
        
        if not filename:
            return
        
        try:
            # Determine file type based on extension
            file_ext = os.path.splitext(filename)[1].lower()
            
            if file_ext == '.json':
                self.load_from_json_file(filename)
            elif file_ext == '.xlsx':
                self.load_from_excel_file(filename)
            elif file_ext == '.csv':
                self.load_from_csv_file(filename, delimiter=',')
            elif file_ext == '.tsv':
                self.load_from_csv_file(filename, delimiter='\t')
            else:
                # Try to determine file type from content
                with open(filename, 'r') as f:
                    first_line = f.readline().strip()
                    
                if ',' in first_line:
                    self.load_from_csv_file(filename, delimiter=',')
                elif '\t' in first_line:
                    self.load_from_csv_file(filename, delimiter='\t')
                else:
                    QMessageBox.warning(self, "Unknown File Type", 
                                       f"Could not determine the format of {filename}.\n"
                                       "Please use CSV, TSV, Excel, or JSON format.")
                    return
            
            # Enable run button if items were loaded
            if len(self.batch_items) > 0:
                self.run_btn.setEnabled(True)
                
        except Exception as e:
            QMessageBox.critical(self, "Error Loading File", 
                                f"Failed to load batch items from {filename}:\n{str(e)}")
    
    def load_from_json_file(self, filename):
        """Load batch items from a JSON file containing rb_spec objects."""
        try:
            from rbcodes.GUIs.rb_spec import load_rb_spec_object
        except ImportError:
            try:
                from GUIs.rb_spec import load_rb_spec_object
            except ImportError:
                QMessageBox.warning(self, "Import Error", 
                                 "Could not import load_rb_spec_object. Make sure rbcodes is installed.")
                return
        
        try:
            # Load the rb_spec object
            spec = load_rb_spec_object(filename)
            
            # Check if we have the necessary attributes
            if not all(hasattr(spec, attr) for attr in ['trans', 'trans_wave', 'zabs']):
                QMessageBox.warning(self, "Incomplete JSON", 
                                   f"The JSON file {filename} does not contain all required attributes for batch processing.")
                return
            
            # Get the batch mode based on the file
            batch_mode = self.mode_combo.currentIndex()
            
            if batch_mode == 0:  # Transitions mode
                # Create a transition item
                item = {
                    'type': 'transition',
                    'name': spec.trans,
                    'wavelength': spec.trans_wave,
                    'description': f"Loaded from {os.path.basename(filename)}",
                    'linelist': getattr(spec, 'linelist', 'atom'),
                    'fvalue': getattr(spec, 'fval', None),
                    'vmin': getattr(spec, 'vmin', -200),
                    'vmax': getattr(spec, 'vmax', 200),
                    'source': filename
                }
                
                # Add to table
                row = self.batch_table.rowCount()
                self.batch_table.insertRow(row)
                
                self.batch_table.setItem(row, 0, QTableWidgetItem(item['name']))
                self.batch_table.setItem(row, 1, QTableWidgetItem(f"{item['wavelength']:.2f}"))
                self.batch_table.setItem(row, 2, QTableWidgetItem(item['description']))
                
                # Status column
                status_item = QTableWidgetItem("Loaded")
                status_item.setForeground(QBrush(QColor("blue")))
                self.batch_table.setItem(row, 3, status_item)
                
                # Add to batch items
                self.batch_items.append(item)
                
            else:  # Files mode
                # We need to consider this as a file to be processed
                item = {
                    'type': 'file',
                    'filename': filename,
                    'transition_name': spec.trans,
                    'wavelength': spec.trans_wave,
                    'redshift': spec.zabs,
                    'linelist': getattr(spec, 'linelist', 'atom'),
                    'fvalue': getattr(spec, 'fval', None),
                    'vmin': getattr(spec, 'vmin', -200),
                    'vmax': getattr(spec, 'vmax', 200),
                    'source': filename
                }
                
                # Add to table
                row = self.batch_table.rowCount()
                self.batch_table.insertRow(row)
                
                self.batch_table.setItem(row, 0, QTableWidgetItem(os.path.basename(filename)))
                self.batch_table.setItem(row, 1, QTableWidgetItem(item['transition_name']))
                self.batch_table.setItem(row, 2, QTableWidgetItem(f"{item['redshift']:.5f}"))
                
                # Status column
                status_item = QTableWidgetItem("Loaded")
                status_item.setForeground(QBrush(QColor("blue")))
                self.batch_table.setItem(row, 3, status_item)
                
                # Add to batch items
                self.batch_items.append(item)
            
            QMessageBox.information(self, "JSON Loaded", 
                                  f"Successfully loaded rb_spec object from {filename}")
            
        except Exception as e:
            QMessageBox.warning(self, "Error", 
                              f"Failed to load JSON file {filename}: {str(e)}")
    
    def load_from_excel_file(self, filename):
        """Load batch items from an Excel file."""
        try:
            # Import pandas for Excel reading
            import pandas as pd
            
            # Read the Excel file
            df = pd.read_excel(filename)
            
            # Process based on batch mode
            if self.mode_combo.currentIndex() == 0:  # Transitions mode
                self._process_transitions_dataframe(df)
            else:  # Files mode
                self._process_files_dataframe(df)
            
        except Exception as e:
            QMessageBox.warning(self, "Error", 
                              f"Failed to load Excel file: {str(e)}")
    
    def load_from_csv_file(self, filename, delimiter=','):
        """Load batch items from a CSV/TSV file."""
        try:
            # Import pandas for CSV reading
            import pandas as pd
            
            # Read the CSV file
            df = pd.read_csv(filename, delimiter=delimiter)
            
            # Process based on batch mode
            if self.mode_combo.currentIndex() == 0:  # Transitions mode
                self._process_transitions_dataframe(df)
            else:  # Files mode
                self._process_files_dataframe(df)
            
        except Exception as e:
            QMessageBox.warning(self, "Error", 
                              f"Failed to load CSV file: {str(e)}")
    
    def _process_transitions_dataframe(self, df):
        """Process a dataframe containing transition information."""
        # Check required columns
        required_cols = ['transition_name', 'wavelength']
        if not all(col in df.columns for col in required_cols):
            # Try alternative column names
            alt_cols = {
                'transition_name': ['transition', 'name', 'trans'],
                'wavelength': ['wave', 'lambda', 'wl']
            }
            
            for req_col, alternatives in alt_cols.items():
                if req_col not in df.columns:
                    for alt in alternatives:
                        if alt in df.columns:
                            df = df.rename(columns={alt: req_col})
                            break
        
        # Check again after trying alternatives
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            QMessageBox.warning(self, "Missing Columns", 
                               f"The file is missing required columns: {', '.join(missing)}")
            return
        
        # Process each row
        count = 0
        for _, row in df.iterrows():
            try:
                # Extract data with defaults for optional fields
                transition_name = str(row['transition_name'])
                wavelength = float(row['wavelength'])
                description = str(row.get('description', '')) if 'description' in row else ''
                linelist = str(row.get('linelist', 'atom')) if 'linelist' in row else 'atom'
                vmin = float(row.get('vmin', -200)) if 'vmin' in row else -200
                vmax = float(row.get('vmax', 200)) if 'vmax' in row else 200
                
                # Create transition item
                item = {
                    'type': 'transition',
                    'name': transition_name,
                    'wavelength': wavelength,
                    'description': description,
                    'linelist': linelist,
                    'fvalue': None,  # Will be looked up when processing
                    'vmin': vmin,
                    'vmax': vmax
                }
                
                # Add to table
                table_row = self.batch_table.rowCount()
                self.batch_table.insertRow(table_row)
                
                self.batch_table.setItem(table_row, 0, QTableWidgetItem(transition_name))
                self.batch_table.setItem(table_row, 1, QTableWidgetItem(f"{wavelength:.2f}"))
                self.batch_table.setItem(table_row, 2, QTableWidgetItem(description))
                
                # Status column
                status_item = QTableWidgetItem("Ready")
                status_item.setForeground(QBrush(QColor("blue")))
                self.batch_table.setItem(table_row, 3, status_item)
                
                # Add to batch items
                self.batch_items.append(item)
                count += 1
                
            except Exception as e:
                print(f"Error processing row: {str(e)}")
        
        QMessageBox.information(self, "Import Complete", 
                              f"Successfully imported {count} transitions.")
    
    def _process_files_dataframe(self, df):
        """Process a dataframe containing file information."""
        # Check required columns
        required_cols = ['filename', 'transition']
        if not all(col in df.columns for col in required_cols):
            # Try alternative column names
            alt_cols = {
                'filename': ['file', 'path', 'file_path'],
                'transition': ['transition_name', 'trans', 'line']
            }
            
            for req_col, alternatives in alt_cols.items():
                if req_col not in df.columns:
                    for alt in alternatives:
                        if alt in df.columns:
                            df = df.rename(columns={alt: req_col})
                            break
        
        # Check again after trying alternatives
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            QMessageBox.warning(self, "Missing Columns", 
                               f"The file is missing required columns: {', '.join(missing)}")
            return
        
        # Process each row
        count = 0
        for _, row in df.iterrows():
            try:
                # Extract data with defaults for optional fields
                filename = str(row['filename'])
                transition = str(row['transition'])
                redshift = float(row.get('redshift', 0.0)) if 'redshift' in row else 0.0
                wavelength = float(row.get('wavelength', 0.0)) if 'wavelength' in row else 0.0
                linelist = str(row.get('linelist', 'atom')) if 'linelist' in row else 'atom'
                
                # If wavelength is not provided, try to extract it from transition name
                if wavelength == 0.0 and '(' in transition and ')' in transition:
                    try:
                        wl_text = transition.split('(')[1].split(')')[0].strip()
                        if 'Å' in wl_text:
                            wl_text = wl_text.replace('Å', '').strip()
                        wavelength = float(wl_text)
                    except:
                        # Use a default and look it up later
                        pass
                
                # Check if file exists
                if not os.path.exists(filename):
                    # Try relative to the input file
                    base_dir = os.path.dirname(df.file_path if hasattr(df, 'file_path') else '')
                    alt_path = os.path.join(base_dir, filename)
                    if os.path.exists(alt_path):
                        filename = alt_path
                    else:
                        print(f"Warning: File not found: {filename}")
                
                # Create file item
                item = {
                    'type': 'file',
                    'filename': filename,
                    'transition_name': transition,
                    'wavelength': wavelength,
                    'redshift': redshift,
                    'linelist': linelist,
                    'fvalue': None  # Will be looked up when processing
                }
                
                # Add to table
                table_row = self.batch_table.rowCount()
                self.batch_table.insertRow(table_row)
                
                self.batch_table.setItem(table_row, 0, QTableWidgetItem(os.path.basename(filename)))
                self.batch_table.setItem(table_row, 1, QTableWidgetItem(transition))
                self.batch_table.setItem(table_row, 2, QTableWidgetItem(f"{redshift:.5f}"))
                
                # Status column
                status_item = QTableWidgetItem("Ready")
                status_item.setForeground(QBrush(QColor("blue")))
                self.batch_table.setItem(table_row, 3, status_item)
                
                # Add to batch items
                self.batch_items.append(item)
                count += 1
                
            except Exception as e:
                print(f"Error processing row: {str(e)}")
        
        QMessageBox.information(self, "Import Complete", 
                              f"Successfully imported {count} file entries.")
    
    def export_template(self):
        """Export a template file for batch input."""
        options = QFileDialog.Options()
        file_filter = "CSV Files (*.csv);;TSV Files (*.tsv);;Excel Files (*.xlsx)"
        filename, selected_filter = QFileDialog.getSaveFileName(
            self, "Save Template File", "batch_template", file_filter, options=options
        )
        
        if not filename:
            return
        
        try:
            # Determine file type based on selection
            if "CSV" in selected_filter:
                ext = ".csv"
                delim = ","
            elif "TSV" in selected_filter:
                ext = ".tsv"
                delim = "\t"
            elif "Excel" in selected_filter:
                ext = ".xlsx"
            else:
                ext = ".csv"
                delim = ","
            
            # Ensure correct extension
            if not filename.lower().endswith(ext):
                filename += ext
            
            # Create template data
            if self.mode_combo.currentIndex() == 0:  # Transitions mode
                # Create transitions template
                headers = ["transition_name", "wavelength", "description", "linelist", "vmin", "vmax"]
                
                # Sample data rows
                sample_data = [
                    ["Lyα", 1215.67, "Lyman Alpha Example", "atom", -200, 200],
                    ["CIV", 1548.20, "Carbon IV Example", "atom", -200, 200],
                    ["MgII", 2796.35, "Magnesium II Example", "atom", -200, 200],
                    ["SiIV", 1393.76, "Silicon IV Example", "atom", -200, 200]
                ]
                
                # Create dataframe
                df = pd.DataFrame(sample_data, columns=headers)
                
            else:  # Files mode
                # Create files template
                headers = ["filename", "transition", "wavelength", "redshift", "linelist"]
                
                # Sample data rows
                sample_data = [
                    ["spectrum1.fits", "Lyα", 1215.67, 0.5, "atom"],
                    ["spectrum2.fits", "CIV", 1548.20, 1.0, "atom"],
                    ["spectrum3.fits", "MgII", 2796.35, 0.75, "atom"],
                    ["spectrum4.fits", "SiIV", 1393.76, 1.2, "atom"]
                ]
                
                # Create dataframe
                df = pd.DataFrame(sample_data, columns=headers)
            
            # Save to file
            if ext == ".xlsx":
                df.to_excel(filename, index=False)
            else:
                df.to_csv(filename, index=False, sep=delim)
            
            QMessageBox.information(self, "Template Exported", 
                                  f"Template file saved to {filename}")
            
        except Exception as e:
            QMessageBox.warning(self, "Export Error", 
                              f"Failed to export template: {str(e)}")
    
    def add_batch_item(self):
        """Add an item to the batch processing list."""
        if self.mode_combo.currentIndex() == 0:
            # Add transition
            self.add_transition_dialog()
        else:
            # Add file
            self.add_file_dialog()
    
    def add_transition_dialog(self):
        """Open dialog to add a transition to the batch."""
        from PyQt5.QtWidgets import QDialog, QVBoxLayout, QDialogButtonBox
        
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Transition")
        layout = QVBoxLayout(dialog)
        
        # Form for transition details
        form = QFormLayout()
        
        # Common transitions dropdown
        transitions_combo = QComboBox()
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
        transitions_combo.addItems(transitions)
        transitions_combo.currentTextChanged.connect(
            lambda text: custom_wavelength.setEnabled(text == "Custom"))
        
        form.addRow("Transition:", transitions_combo)
        
        # Custom wavelength input
        custom_wavelength = QDoubleSpinBox()
        custom_wavelength.setRange(1.0, 10000.0)
        custom_wavelength.setValue(1215.67)
        custom_wavelength.setDecimals(2)
        custom_wavelength.setEnabled(False)
        
        form.addRow("Custom Wavelength (Å):", custom_wavelength)
        
        # Line list selection
        linelist_combo = QComboBox()
        linelist_combo.addItems(["atom", "LLS", "LLS Small", "DLA", "LBG", "Gal"])
        form.addRow("Line List:", linelist_combo)
        
        # Description
        description = QLineEdit()
        description.setPlaceholderText("Optional description")
        
        form.addRow("Description:", description)
        
        layout.addLayout(form)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            # Get the line list
            linelist = linelist_combo.currentText()
            
            # Process the selected transition
            transition_text = transitions_combo.currentText()
            
            if transition_text == "Custom":
                # Get wavelength and look up transition using rb_setline
                wavelength = custom_wavelength.value()
                
                try:
                    # Import rb_setline
                    from rbcodes.IGM import rb_setline as s
                except ImportError:
                    # Try alternate import path
                    try:
                        from IGM import rb_setline as s
                    except ImportError:
                        QMessageBox.warning(self, "Import Error", 
                                         "Could not import rb_setline. Using custom wavelength without lookup.")
                        transition_name = f"Custom ({wavelength:.2f} Å)"
                        transition_info = None
                
                # Look up the transition
                try:
                    transition_info = s.rb_setline(wavelength, 'closest', linelist=linelist)
                    wavelength = float(transition_info['wave'])
                    transition_name = transition_info['name']
                except Exception as e:
                    print(f"Error looking up transition: {str(e)}")
                    transition_name = f"Custom ({wavelength:.2f} Å)"
                    transition_info = None
            else:
                # Extract wavelength from the selected item
                wavelength = float(transition_text.split("(")[1].split(" Å")[0])
                transition_name = transition_text.split(" (")[0]
                
                # Look up the transition for consistency
                try:
                    from rbcodes.IGM import rb_setline as s
                except ImportError:
                    try:
                        from IGM import rb_setline as s
                    except ImportError:
                        QMessageBox.warning(self, "Import Error", 
                                         "Could not import rb_setline. Using predefined wavelength.")
                        transition_info = None
                
                try:
                    transition_info = s.rb_setline(wavelength, 'closest', linelist=linelist)
                    wavelength = float(transition_info['wave'])
                    transition_name = transition_info['name']
                except Exception as e:
                    print(f"Error looking up transition: {str(e)}")
                    transition_info = None
            
            # Add to table
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            self.batch_table.setItem(row, 0, QTableWidgetItem(transition_name))
            self.batch_table.setItem(row, 1, QTableWidgetItem(f"{wavelength:.2f}"))
            self.batch_table.setItem(row, 2, QTableWidgetItem(description.text()))
            
            # Status column
            status_item = QTableWidgetItem("Ready")
            status_item.setForeground(QBrush(QColor("blue")))
            self.batch_table.setItem(row, 3, status_item)
            
            # Store item data
            self.batch_items.append({
                'type': 'transition',
                'name': transition_name,
                'wavelength': wavelength,
                'description': description.text(),
                'linelist': linelist,
                'fvalue': transition_info['fval'] if transition_info else None
            })
            
            # Enable run button
            self.run_btn.setEnabled(True)
        
    def add_file_dialog(self):
        """Open dialog to add a file to the batch."""
        # First select the file
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select Spectrum File", "", 
            "All Files (*);;FITS Files (*.fits);;JSON Files (*.json)",
            options=options
        )
        
        if not filename:
            return
        
        # Now open dialog for transition and redshift
        from PyQt5.QtWidgets import QDialog, QVBoxLayout, QDialogButtonBox, QLineEdit
        
        dialog = QDialog(self)
        dialog.setWindowTitle("File Parameters")
        layout = QVBoxLayout(dialog)
        
        # Form for file parameters
        form = QFormLayout()
        
        # Transition
        transitions_combo = QComboBox()
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
        transitions_combo.addItems(transitions)
        transitions_combo.currentTextChanged.connect(
            lambda text: custom_wavelength.setEnabled(text == "Custom"))
        
        form.addRow("Transition:", transitions_combo)
        
        # Custom wavelength input
        custom_wavelength = QDoubleSpinBox()
        custom_wavelength.setRange(1.0, 10000.0)
        custom_wavelength.setValue(1215.67)
        custom_wavelength.setDecimals(2)
        custom_wavelength.setEnabled(False)
        
        form.addRow("Custom Wavelength (Å):", custom_wavelength)
        
        # Line list selection
        linelist_combo = QComboBox()
        linelist_combo.addItems(["atom", "LLS", "LLS Small", "DLA", "LBG", "Gal"])
        form.addRow("Line List:", linelist_combo)
        
        # Redshift
        redshift = QDoubleSpinBox()
        redshift.setRange(0.0, 10.0)
        redshift.setValue(0.0)
        redshift.setDecimals(5)
        
        form.addRow("Redshift:", redshift)
        
        layout.addLayout(form)
        
        # Add buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        layout.addWidget(button_box)
        
        # Show dialog
        if dialog.exec_() == QDialog.Accepted:
            # Get the line list
            linelist = linelist_combo.currentText()
            
            # Process the selected transition
            transition_text = transitions_combo.currentText()
            
            if transition_text == "Custom":
                # Get wavelength and look up transition using rb_setline
                wavelength = custom_wavelength.value()
                
                try:
                    # Import rb_setline
                    from rbcodes.IGM import rb_setline as s
                except ImportError:
                    # Try alternate import path
                    try:
                        from IGM import rb_setline as s
                    except ImportError:
                        QMessageBox.warning(self, "Import Error", 
                                         "Could not import rb_setline. Using custom wavelength without lookup.")
                        transition_name = f"Custom ({wavelength:.2f} Å)"
                        transition_info = None
                
                # Look up the transition
                try:
                    transition_info = s.rb_setline(wavelength, 'closest', linelist=linelist)
                    wavelength = float(transition_info['wave'])
                    transition_name = transition_info['name']
                except Exception as e:
                    print(f"Error looking up transition: {str(e)}")
                    transition_name = f"Custom ({wavelength:.2f} Å)"
                    transition_info = None
            else:
                # Extract wavelength from the selected item
                wavelength = float(transition_text.split("(")[1].split(" Å")[0])
                transition_name = transition_text.split(" (")[0]
                
                # Look up the transition for consistency
                try:
                    from rbcodes.IGM import rb_setline as s
                except ImportError:
                    try:
                        from IGM import rb_setline as s
                    except ImportError:
                        QMessageBox.warning(self, "Import Error", 
                                         "Could not import rb_setline. Using predefined wavelength.")
                        transition_info = None
                
                try:
                    transition_info = s.rb_setline(wavelength, 'closest', linelist=linelist)
                    wavelength = float(transition_info['wave'])
                    transition_name = transition_info['name']
                except Exception as e:
                    print(f"Error looking up transition: {str(e)}")
                    transition_info = None
            
            # Add to table
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            self.batch_table.setItem(row, 0, QTableWidgetItem(os.path.basename(filename)))
            self.batch_table.setItem(row, 1, QTableWidgetItem(transition_name))
            self.batch_table.setItem(row, 2, QTableWidgetItem(f"{redshift.value():.5f}"))
            
            # Status column
            status_item = QTableWidgetItem("Ready")
            status_item.setForeground(QBrush(QColor("blue")))
            self.batch_table.setItem(row, 3, status_item)
            
            # Store item data
            self.batch_items.append({
                'type': 'file',
                'filename': filename,
                'transition_name': transition_name,
                'wavelength': wavelength,
                'redshift': redshift.value(),
                'linelist': linelist,
                'fvalue': transition_info['fval'] if transition_info else None
            })
            
            # Enable run button
            self.run_btn.setEnabled(True)
        
    def remove_batch_item(self):
        """Remove the selected item from the batch."""
        selected_rows = set(index.row() for index in self.batch_table.selectedIndexes())
        
        # Remove from bottom to top to avoid index shifting
        for row in sorted(selected_rows, reverse=True):
            self.batch_table.removeRow(row)
            if 0 <= row < len(self.batch_items):
                self.batch_items.pop(row)
        
        # Disable run button if no items left
        if len(self.batch_items) == 0:
            self.run_btn.setEnabled(False)
    
    def clear_batch_items(self):
        """Clear all batch items."""
        self.batch_table.setRowCount(0)
        self.batch_items = []
        self.run_btn.setEnabled(False)
    
    def run_batch(self):
        """Run the batch processing tasks."""
        if not self.batch_items:
            QMessageBox.warning(self, "No Items", "Add items to the batch first.")
            return
        
        # Get processing parameters
        params = {
            'cont_method': self.cont_method.currentText(),
            'poly_order': self.poly_order.value(),
            'use_weights': self.use_weights.isChecked(),
            'use_default_vrange': self.default_vrange.isChecked(),
            'vmin': self.vmin.value(),
            'vmax': self.vmax.value(),
            'calc_snr': self.calc_snr.isChecked()
        }
        
        # Clear previous results
        self.results = []
        self.results_table.setRowCount(0)
        
        # Process each item
        for i, item in enumerate(self.batch_items):
            try:
                # Update status
                self.results_label.setText(f"Processing item {i+1} of {len(self.batch_items)}...")
                QApplication.processEvents()  # Update UI
                
                # Update status in batch table
                status_item = QTableWidgetItem("Processing...")
                status_item.setForeground(QBrush(QColor("orange")))
                self.batch_table.setItem(i, self.batch_table.columnCount()-1, status_item)
                QApplication.processEvents()  # Update UI
                
                if item['type'] == 'transition':
                    # Process transition
                    result = self.process_transition(item, params)
                else:
                    # Process file
                    result = self.process_file(item, params)
                
                # Add to results
                if result:
                    self.results.append(result)
                    
                    # Update table with result info in batch table
                    status_item = QTableWidgetItem("✓ Completed")
                    status_item.setForeground(QBrush(QColor("green")))
                    self.batch_table.setItem(i, self.batch_table.columnCount()-1, status_item)
                    
                    # Add to results table
                    self.add_result_to_table(result)
                else:
                    # Update table with error info
                    status_item = QTableWidgetItem("❌ Error")
                    status_item.setForeground(QBrush(QColor("red")))
                    self.batch_table.setItem(i, self.batch_table.columnCount()-1, status_item)
                
                QApplication.processEvents()  # Update UI
                
            except Exception as e:
                print(f"Error processing batch item {i+1}: {str(e)}")
                
                # Update table with error info
                status_item = QTableWidgetItem("❌ Error")
                status_item.setForeground(QBrush(QColor("red")))
                self.batch_table.setItem(i, self.batch_table.columnCount()-1, status_item)
                QApplication.processEvents()  # Update UI
        
        # Update results summary
        if self.results:
            summary = f"Processed {len(self.results)} of {len(self.batch_items)} items successfully."
            self.results_label.setText(summary)
            self.export_btn.setEnabled(True)
            
            # Signal completion
            self.batch_completed.emit(self.results)
        else:
            self.results_label.setText("Batch processing completed with no successful results.")
            self.export_btn.setEnabled(False)
    
    def add_result_to_table(self, result):
        """Add a result to the results table."""
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        
        # Status column with icon
        status_item = QTableWidgetItem("✓")
        status_item.setForeground(QBrush(QColor("green")))
        self.results_table.setItem(row, 0, status_item)
        
        # Item name/identifier
        if result['type'] == 'transition':
            item_text = f"{result['name']} ({result['wavelength']:.2f} Å)"
        else:
            item_text = f"{os.path.basename(result['filename'])}: {result['name']}"
        
        self.results_table.setItem(row, 1, QTableWidgetItem(item_text))
        
        # Equivalent width
        if 'W' in result and 'W_e' in result:
            self.results_table.setItem(row, 2, QTableWidgetItem(f"{result['W']:.3f}"))
            self.results_table.setItem(row, 3, QTableWidgetItem(f"{result['W_e']:.3f}"))
        else:
            self.results_table.setItem(row, 2, QTableWidgetItem("N/A"))
            self.results_table.setItem(row, 3, QTableWidgetItem("N/A"))
        
        # Column density
        if 'logN' in result and 'logN_e' in result:
            self.results_table.setItem(row, 4, QTableWidgetItem(f"{result['logN']:.2f}"))
            self.results_table.setItem(row, 5, QTableWidgetItem(f"{result['logN_e']:.2f}"))
        else:
            self.results_table.setItem(row, 4, QTableWidgetItem("N/A"))
            self.results_table.setItem(row, 5, QTableWidgetItem("N/A"))
        
    def process_transition(self, item, params):
        """Process a transition batch item."""
        # This assumes we're using the current loaded spectrum
        if not self.controller.has_spectrum():
            raise ValueError("No spectrum loaded.")
        
        # Save current state
        original_state = self.controller.save_state()
        
        try:
            # Get parameters
            wavelength = item['wavelength']
            linelist = item.get('linelist', 'atom')
            
            # Get velocity ranges
            vmin = item.get('vmin', params['vmin'])
            vmax = item.get('vmax', params['vmax'])
            
            # Slice the spectrum
            success = self.controller.slice_spectrum(
                wavelength, vmin, vmax, 
                use_vel=True, 
                linelist=linelist
            )
            
            if not success:
                raise ValueError(f"Failed to slice spectrum for {item['name']}")
            
            # Fit continuum
            if params['cont_method'] == 'Polynomial':
                # Use polynomial fitting
                success = self.controller.fit_continuum(
                    Legendre=params['poly_order'],
                    use_weights=params['use_weights']
                )
            else:
                # Use flat continuum
                success = self.controller.apply_flat_continuum()
            
            if not success:
                raise ValueError(f"Failed to fit continuum for {item['name']}")
            
            # Compute EW
            success, ew_results = self.controller.compute_equivalent_width(
                vmin=vmin, 
                vmax=vmax, 
                snr=params['calc_snr']
            )
            
            if not success or not ew_results:
                raise ValueError(f"Failed to compute EW for {item['name']}")
            
            # Create result that includes original item info and EW results
            result = item.copy()  # Start with the original item
            result.update(ew_results)  # Add EW results
            
            return result
            
        except Exception as e:
            print(f"Error processing transition {item['name']}: {str(e)}")
            return None
        finally:
            # Restore original state
            self.controller.restore_state(original_state)
    
    def process_file(self, item, params):
        """Process a file batch item."""
        # Save current state
        original_state = self.controller.save_state()
        
        try:
            # Get parameters
            filename = item['filename']
            wavelength = item['wavelength']
            redshift = item['redshift']
            linelist = item.get('linelist', 'atom')
            
            # Get velocity ranges
            vmin = item.get('vmin', params['vmin'])
            vmax = item.get('vmax', params['vmax'])
            
            # Check if JSON file
            is_json = filename.lower().endswith('.json')
            
            # Load the file
            if is_json:
                success, _, _ = self.controller.load_from_json(filename)
            else:
                success, _, _ = self.controller.load_from_file(filename)
                
            if not success:
                raise ValueError(f"Failed to load file {filename}")
            
            # Apply redshift (skip if loaded from JSON and already has redshift)
            if not (is_json and hasattr(self.controller.spec, 'zabs')):
                success = self.controller.apply_redshift(redshift)
                if not success:
                    raise ValueError(f"Failed to apply redshift {redshift} to {filename}")
            
            # Slice the spectrum
            success = self.controller.slice_spectrum(
                wavelength, vmin, vmax, 
                use_vel=True, 
                linelist=linelist
            )
            
            if not success:
                raise ValueError(f"Failed to slice spectrum for {filename}")
            
            # Fit continuum
            if params['cont_method'] == 'Polynomial':
                # Use polynomial fitting
                success = self.controller.fit_continuum(
                    Legendre=params['poly_order'],
                    use_weights=params['use_weights']
                )
            else:
                # Use flat continuum
                success = self.controller.apply_flat_continuum()
            
            if not success:
                raise ValueError(f"Failed to fit continuum for {filename}")
            
            # Compute EW
            success, ew_results = self.controller.compute_equivalent_width(
                vmin=vmin, 
                vmax=vmax, 
                snr=params['calc_snr']
            )
            
            if not success or not ew_results:
                raise ValueError(f"Failed to compute EW for {filename}")
            
            # Create result that includes original item info and EW results
            result = item.copy()  # Start with the original item
            result.update(ew_results)  # Add EW results
            
            return result
            
        except Exception as e:
            print(f"Error processing file {item['filename']}: {str(e)}")
            return None
        finally:
            # Restore original state
            self.controller.restore_state(original_state)
    
    def export_results(self):
        """Export the batch processing results."""
        if not self.results:
            QMessageBox.warning(self, "No Results", "Run batch processing first.")
            return
        
        # Open file dialog to select save location
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Batch Results", "batch_results.csv", 
            "CSV Files (*.csv);;Excel Files (*.xlsx);;All Files (*)",
            options=options
        )
        
        if not filename:
            return
        
        try:
            # Determine file type
            if filename.lower().endswith('.xlsx'):
                self.export_to_excel(filename)
            else:
                # Default to CSV
                if not filename.lower().endswith('.csv'):
                    filename += '.csv'
                self.export_to_csv(filename)
            
            QMessageBox.information(self, "Export Complete", f"Results saved to {filename}")
        except Exception as e:
            print(f"Error exporting results: {str(e)}")
            QMessageBox.warning(self, "Export Error", f"Failed to export results: {str(e)}")
    
    def export_to_csv(self, filename):
        """Export results to a CSV file."""
        with open(filename, 'w') as f:
            # Write header
            if self.mode_combo.currentIndex() == 0:
                # Transitions mode
                f.write("Transition,Wavelength,Description,EW,EW_err,logN,logN_err,Velocity_Centroid,Velocity_Dispersion,SNR\n")
            else:
                # Files mode
                f.write("Filename,Transition,Redshift,EW,EW_err,logN,logN_err,Velocity_Centroid,Velocity_Dispersion,SNR\n")
            
            # Write data
            for result in self.results:
                if result['type'] == 'transition':
                    # Transitions mode
                    f.write(f"{result['name']},{result['wavelength']:.2f},{result.get('description', '')},")
                else:
                    # Files mode
                    f.write(f"{os.path.basename(result['filename'])},{result['transition_name']},{result['redshift']:.5f},")
                
                # Common fields
                f.write(f"{result.get('W', 0.0):.3f},{result.get('W_e', 0.0):.3f},")
                f.write(f"{result.get('logN', 0.0):.2f},{result.get('logN_e', 0.0):.2f},")
                f.write(f"{result.get('vel_centroid', 0.0):.1f},{result.get('vel_disp', 0.0):.1f},")
                f.write(f"{result.get('SNR', -99):.1f}\n")
    
    def export_to_excel(self, filename):
        """Export results to an Excel file."""
        # Create a pandas DataFrame
        data = []
        
        for result in self.results:
            row = {}
            
            if result['type'] == 'transition':
                # Transitions mode
                row['Item_Type'] = 'Transition'
                row['Transition'] = result['name']
                row['Wavelength'] = result['wavelength']
                row['Description'] = result.get('description', '')
            else:
                # Files mode
                row['Item_Type'] = 'File'
                row['Filename'] = os.path.basename(result['filename'])
                row['Transition'] = result['transition_name']
                row['Redshift'] = result['redshift']
            
            # Common fields
            row['EW'] = result.get('W', 0.0)
            row['EW_err'] = result.get('W_e', 0.0)
            row['logN'] = result.get('logN', 0.0)
            row['logN_e'] = result.get('logN_e', 0.0)
            row['Velocity_Centroid'] = result.get('vel_centroid', 0.0)
            row['Velocity_Dispersion'] = result.get('vel_disp', 0.0)
            row['SNR'] = result.get('SNR', -99)
            
            data.append(row)
        
        # Create DataFrame
        df = pd.DataFrame(data)
        
        # Export to Excel
        df.to_excel(filename, index=False)
    
    def launch_advanced_mode(self):
        """Launch the advanced batch processing dialog."""
        try:
            # Import the advanced batch dialog
            from advanced_batch_dialog import AdvancedBatchDialog
            
            # Create and show the dialog
            dialog = AdvancedBatchDialog(self.controller, self.batch_items, self.results)
            dialog.exec_()
            
            # Update our batch items and results if the dialog modified them
            if dialog.batch_items:
                self.batch_items = dialog.batch_items
                self.update_batch_table()
            
            if dialog.results:
                self.results = dialog.results
                self.update_results_table()
                
        except ImportError:
            QMessageBox.information(self, "Not Available", 
                                  "Advanced batch mode is not yet implemented.\n"
                                  "This feature will be available in a future update.")
    
    def update_batch_table(self):
        """Update the batch table with current batch items."""
        # Clear the table
        self.batch_table.setRowCount(0)
        
        # Re-populate with current items
        for item in self.batch_items:
            row = self.batch_table.rowCount()
            self.batch_table.insertRow(row)
            
            if item['type'] == 'transition':
                self.batch_table.setItem(row, 0, QTableWidgetItem(item['name']))
                self.batch_table.setItem(row, 1, QTableWidgetItem(f"{item['wavelength']:.2f}"))
                self.batch_table.setItem(row, 2, QTableWidgetItem(item.get('description', '')))
            else:
                self.batch_table.setItem(row, 0, QTableWidgetItem(os.path.basename(item['filename'])))
                self.batch_table.setItem(row, 1, QTableWidgetItem(item['transition_name']))
                self.batch_table.setItem(row, 2, QTableWidgetItem(f"{item['redshift']:.5f}"))
            
            # Status column
            status_item = QTableWidgetItem("Ready")
            status_item.setForeground(QBrush(QColor("blue")))
            self.batch_table.setItem(row, 3, status_item)
    
    def update_results_table(self):
        """Update the results table with current results."""
        # Clear the table
        self.results_table.setRowCount(0)
        
        # Re-populate with current results
        for result in self.results:
            self.add_result_to_table(result)
        
        # Update summary
        if self.results:
            summary = f"Processed {len(self.results)} of {len(self.batch_items)} items successfully."
            self.results_label.setText(summary)
            self.export_btn.setEnabled(True)
        else:
            self.results_label.setText("No batch processing results yet.")
            self.export_btn.setEnabled(False)

