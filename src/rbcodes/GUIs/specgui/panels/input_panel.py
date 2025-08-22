# panels/input_panel.py
import os
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
                           QFileDialog, QComboBox, QGroupBox, QFormLayout, QLineEdit,
                           QMessageBox)
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
        self.browse_btn.setToolTip("Select a spectrum file from your computer")
        
        file_path_layout = QHBoxLayout()
        file_path_layout.addWidget(self.file_path)
        file_path_layout.addWidget(self.browse_btn)
        
        file_form.addRow("File:", file_path_layout)
        
        self.file_type = QComboBox()
        self.file_type.addItems(["Auto","linetools", "FITS", "JSON", "ASCII"])
        file_form.addRow("File Type:", self.file_type)
        
        self.load_file_btn = QPushButton("Load File")
        self.load_file_btn.clicked.connect(self.load_from_file)
        self.load_file_btn.setToolTip("Load the selected file into the application")
        self.load_file_btn.setMaximumWidth(120)
        
        file_layout.addLayout(file_form)
        button_layout = QHBoxLayout()
        button_layout.addStretch()  # Left stretch
        button_layout.addWidget(self.load_file_btn)
        button_layout.addStretch()  # Right stretch
        
        file_layout.addLayout(button_layout)
        
        file_group.setLayout(file_layout)
        main_layout.addWidget(file_group)
        
 
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
        elif filetype is None and filename.lower().endswith('.fits'):
            # Use 'linetools' as the default filetype for FITS files
            filetype = 'linetools'
            success, _, _ = self.controller.load_from_file(filename, filetype)
            self.spectrum_loaded.emit(success, False, False)
        else:
            success, _, _ = self.controller.load_from_file(filename, filetype)
            self.spectrum_loaded.emit(success, False, False)