# panels/input_panel.py
import os
import numpy as np
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
                           QFileDialog, QComboBox, QGroupBox, QFormLayout, QLineEdit)
from PyQt5.QtCore import pyqtSignal

class InputPanel(QWidget):
    """Panel for loading spectrum data from file or arrays."""
    
    # Signals
    spectrum_loaded = pyqtSignal(bool, bool, bool)  # success, has_redshift, has_transition

    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # File loading section
        file_group = QGroupBox("Load from File")
        file_layout = QVBoxLayout()
        
        file_form = QFormLayout()
        self.file_path = QLineEdit()
        self.file_path.setReadOnly(True)
        self.browse_btn = QPushButton("Browse")
        self.browse_btn.clicked.connect(self.browse_file)
        
        file_path_layout = QHBoxLayout()
        file_path_layout.addWidget(self.file_path)
        file_path_layout.addWidget(self.browse_btn)
        
        file_form.addRow("File:", file_path_layout)
        
        self.file_type = QComboBox()
        self.file_type.addItems(["Auto", "FITS", "ASCII", "JSON", "linetools"])
        file_form.addRow("File Type:", self.file_type)
        
        self.load_file_btn = QPushButton("Load File")
        self.load_file_btn.clicked.connect(self.load_from_file)
        
        file_layout.addLayout(file_form)
        file_layout.addWidget(self.load_file_btn)
        
        file_group.setLayout(file_layout)
        main_layout.addWidget(file_group)
        
        # Array input section
        array_group = QGroupBox("Load from Data Arrays")
        array_layout = QVBoxLayout()
        
        # In a real implementation, we might have text areas or CSV import
        # For now, just a placeholder button for demonstration
        self.load_data_btn = QPushButton("Import Data Arrays")
        self.load_data_btn.clicked.connect(self.import_arrays)
        array_layout.addWidget(self.load_data_btn)
        
        array_group.setLayout(array_layout)
        main_layout.addWidget(array_group)
        
        # Add stretch to push everything to the top
        main_layout.addStretch(1)
    
    def browse_file(self):
        """Open file dialog to browse for spectrum file."""
        file_filter = "All Files (*);;FITS Files (*.fits);;Text Files (*.txt);;JSON Files (*.json)"
        filename, _ = QFileDialog.getOpenFileName(self, "Open Spectrum File", "", file_filter)
        
        if filename:
            self.file_path.setText(filename)
    

    def load_from_file(self):
        """Load spectrum from the selected file."""
        filename = self.file_path.text()
        if not filename:
            QMessageBox.warning(self, "Input Error", "Please select a file first.")
            return
            
        # Get file type (None for auto)
        filetype = None if self.file_type.currentText() == "Auto" else self.file_type.currentText().lower()
        
        # Check if JSON format
        if filetype == "json" or (filetype is None and filename.lower().endswith('.json')):
            success, has_redshift, has_transition = self.controller.load_from_json(filename)
            self.spectrum_loaded.emit(success, has_redshift, has_transition)
        else:
            success, _, _ = self.controller.load_from_file(filename, filetype)
            self.spectrum_loaded.emit(success, False, False)  # Regular file - no redshift or transition info


    def import_arrays(self):
        """
        Import data arrays - this is a placeholder.
        In a real implementation, this would open a dialog to import CSV, paste data, etc.
        For now, we'll just create some test data.
        """
        try:
            # Create some test data
            wave = np.linspace(4000, 5000, 1000)
            flux = np.ones_like(wave) + 0.1 * np.sin(wave/100)
            error = 0.05 * np.ones_like(wave)
            
            success = self.controller.load_from_data(wave, flux, error)
            self.spectrum_loaded.emit(success)
        except Exception as e:
            print(f"Error creating test data: {str(e)}")
            self.spectrum_loaded.emit(False)

    def create_test_data(self):
        """Create and load test data."""
        try:
            # Create some test data
            wave = np.linspace(4000, 5000, 1000)
            flux = np.ones_like(wave) + 0.1 * np.sin(wave/100)
            error = 0.05 * np.ones_like(wave)
            
            success = self.controller.load_from_data(wave, flux, error)
            self.spectrum_loaded.emit(success, False, False)  # No redshift or transition info for test data
            
            if success:
                QMessageBox.information(self, "Success", "Test data created and loaded successfully.")
        except Exception as e:
            print(f"Error creating test data: {str(e)}")
            self.spectrum_loaded.emit(False, False, False)
            QMessageBox.warning(self, "Error", f"Failed to create test data: {str(e)}")    