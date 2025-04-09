"""
AbsTools Launcher - A wrapper GUI for the Absorption Line Analysis Toolbox.

This GUI provides a user-friendly interface for setting up the absorption line analysis
tool, allowing users to select input files, enter parameters, and launch the main
analysis tool. Now with JSON support for loading saved analyses.
"""

import sys
import os
import numpy as np
import pickle
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QLabel, QPushButton, QLineEdit,
    QFileDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox,
    QCheckBox, QMessageBox, QScrollArea, QFrame, QSplitter, QSpacerItem,
    QSizePolicy, QComboBox
)
from PyQt5.QtGui import QFont, QColor, QPalette
from PyQt5.QtCore import Qt, QSettings

# Import JSON utilities if available
try:
    from json_utils import load_from_json
    JSON_SUPPORT = True
except ImportError:
    JSON_SUPPORT = False
    print("JSON utilities not found. JSON loading support will be disabled.")

class AbsToolsLauncher(QMainWindow):
    """
    Main wrapper GUI window for launching the absorption line analysis tool.
    """
    
    def __init__(self):
        super().__init__()
        
        # Initialize settings
        self.settings = QSettings("AbsTools", "Launcher")
        
        # Set up the main window
        self.setWindowTitle("AbsTools Launcher")
        self.setMinimumSize(900, 700)  # Increased minimum size
        
        # Restore window size and position
        if self.settings.contains("geometry"):
            self.restoreGeometry(self.settings.value("geometry"))
        else:
            self.resize(1100, 800)  # Increased default size if no saved geometry
        
        # Create the central widget and main layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)
        
        # Create UI components
        self.create_ui_components()
        
        # Load settings
        self.load_settings()
        
        # Apply dark mode if enabled
        if self.dark_mode.isChecked():
            self.apply_dark_mode()
        
        # Show the window
        self.show()
    
    def create_ui_components(self):
        """Create all the UI components for the launcher."""
        # Create header with title and description
        self.create_header()
        
        # Create input file selection area
        self.create_file_selection()
        
        # Create parameter input area
        self.create_parameter_inputs()
        
        # Create optional settings area
        self.create_optional_settings()
        
        # Create saved file loading area (NEW)
        self.create_saved_file_loading()
        
        # Create action buttons (Launch, Save Settings, etc.)
        self.create_action_buttons()
        
        # Create status bar
        self.statusBar().showMessage("Ready")
    
    def create_header(self):
        """Create the header section with title and description."""
        header_layout = QVBoxLayout()
        
        # Title label
        title_label = QLabel("AbsTools Launcher")
        title_font = QFont()
        title_font.setPointSize(16)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)
        
        # Description label
        desc_label = QLabel(
            "This tool helps set up and launch the Absorption Line Analysis Toolbox. "
            "Select input files and parameters below to get started."
        )
        desc_label.setWordWrap(True)
        desc_label.setAlignment(Qt.AlignCenter)
        
        # Add to layout
        header_layout.addWidget(title_label)
        header_layout.addWidget(desc_label)
        header_layout.addSpacing(10)
        
        # Add a horizontal line
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        header_layout.addWidget(line)
        
        self.main_layout.addLayout(header_layout)
    
    def create_file_selection(self):
        """Create the file selection area."""
        # Create group box
        file_group = QGroupBox("Input Files")
        file_layout = QGridLayout()
        
        # Spectrum file selection
        file_layout.addWidget(QLabel("Spectrum File:"), 0, 0)
        self.spectrum_path = QLineEdit()
        file_layout.addWidget(self.spectrum_path, 0, 1)
        spectrum_browse = QPushButton("Browse...")
        spectrum_browse.clicked.connect(self.browse_spectrum)
        file_layout.addWidget(spectrum_browse, 0, 2)
        
        # File format selection
        file_layout.addWidget(QLabel("File Format:"), 1, 0)
        self.file_format = QComboBox()
        self.file_format.addItems(["FITS", "ASCII", "HDF5", "Custom"])
        file_layout.addWidget(self.file_format, 1, 1)
        
        # Optional: Intervening line list
        file_layout.addWidget(QLabel("Intervening Line List (Optional):"), 2, 0)
        self.intervening_path = QLineEdit()
        file_layout.addWidget(self.intervening_path, 2, 1)
        intervening_browse = QPushButton("Browse...")
        intervening_browse.clicked.connect(self.browse_intervening)
        file_layout.addWidget(intervening_browse, 2, 2)
        
        file_group.setLayout(file_layout)
        self.main_layout.addWidget(file_group)
    
    def create_parameter_inputs(self):
        """Create the parameter input area."""
        # Create group box
        param_group = QGroupBox("Analysis Parameters")
        param_layout = QGridLayout()
        
        # Redshift
        param_layout.addWidget(QLabel("Redshift:"), 0, 0)
        self.redshift = QLineEdit()
        self.redshift.setPlaceholderText("e.g., 0.348")
        param_layout.addWidget(self.redshift, 0, 1)
        
        # Rest-frame wavelengths
        param_layout.addWidget(QLabel("Rest-frame Wavelengths:"), 1, 0)
        self.wavelengths = QLineEdit()
        self.wavelengths.setPlaceholderText("e.g., 977.6, 1025.5, 1031.2, 1036.2")
        param_layout.addWidget(self.wavelengths, 1, 1)
        
        # Common wavelength presets
        param_layout.addWidget(QLabel("Wavelength Presets:"), 2, 0)
        self.wavelength_presets = QComboBox()
        self.wavelength_presets.addItems([
            "Select a preset...",
            "Lyman Series (HI)",
            "CIV Doublet",
            "SiIV Doublet",
            "OVI Doublet",
            "MgII Doublet"
        ])
        self.wavelength_presets.currentIndexChanged.connect(self.apply_wavelength_preset)
        param_layout.addWidget(self.wavelength_presets, 2, 1)
        
        # Window limits
        param_layout.addWidget(QLabel("Velocity Window Limits (km/s):"), 3, 0)
        window_layout = QHBoxLayout()
        self.vmin = QLineEdit("-2000")
        self.vmax = QLineEdit("2000")
        window_layout.addWidget(self.vmin)
        window_layout.addWidget(QLabel(" to "))
        window_layout.addWidget(self.vmax)
        temp_widget = QWidget()
        temp_widget.setLayout(window_layout)
        param_layout.addWidget(temp_widget, 3, 1)
        
        param_group.setLayout(param_layout)
        self.main_layout.addWidget(param_group)
    
    def create_optional_settings(self):
        """Create the optional settings area."""
        # Create group box
        opt_group = QGroupBox("Optional Settings")
        opt_layout = QGridLayout()
        
        # Initial polynomial order
        opt_layout.addWidget(QLabel("Initial Polynomial Order:"), 0, 0)
        self.poly_order = QComboBox()
        self.poly_order.addItems(["2", "3", "4", "5", "6"])
        self.poly_order.setCurrentIndex(2)  # Default to 4
        opt_layout.addWidget(self.poly_order, 0, 1)
        
        # Initial mask range
        opt_layout.addWidget(QLabel("Initial Mask Range (km/s):"), 1, 0)
        mask_layout = QHBoxLayout()
        self.mask_min = QLineEdit("-200")
        self.mask_max = QLineEdit("200")
        mask_layout.addWidget(self.mask_min)
        mask_layout.addWidget(QLabel(" to "))
        mask_layout.addWidget(self.mask_max)
        temp_widget = QWidget()
        temp_widget.setLayout(mask_layout)
        opt_layout.addWidget(temp_widget, 1, 1)
        
        # Dark mode toggle
        opt_layout.addWidget(QLabel("Dark Mode:"), 2, 0)
        self.dark_mode = QCheckBox("Enable Dark Mode")
        self.dark_mode.setChecked(True)
        self.dark_mode.stateChanged.connect(self.toggle_dark_mode)
        opt_layout.addWidget(self.dark_mode, 2, 1)
        
        opt_group.setLayout(opt_layout)
        self.main_layout.addWidget(opt_group)
    
    def create_saved_file_loading(self):
        """Create the saved file loading area (NEW)."""
        # Create group box
        load_group = QGroupBox("Load Saved Analysis")
        load_layout = QGridLayout()
        
        # Saved file selection
        load_layout.addWidget(QLabel("Saved Analysis File:"), 0, 0)
        self.saved_file_path = QLineEdit()
        self.saved_file_path.setPlaceholderText("Path to previously saved .p (pickle) or .json file")
        load_layout.addWidget(self.saved_file_path, 0, 1)
        saved_browse = QPushButton("Browse...")
        saved_browse.clicked.connect(self.browse_saved_file)
        load_layout.addWidget(saved_browse, 0, 2)
        
        # Load button
        load_btn = QPushButton("Load Saved Analysis")
        load_btn.setMinimumHeight(35)
        font = load_btn.font()
        font.setBold(True)
        load_btn.setFont(font)
        load_btn.clicked.connect(self.load_saved_file)
        load_layout.addWidget(load_btn, 1, 1)
        
        # Note about file types
        note_label = QLabel("Supported formats: .p (Pickle) and .json (JSON)")
        note_label.setStyleSheet("font-style: italic; color: gray;")
        load_layout.addWidget(note_label, 2, 1)
        
        # Add separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        load_layout.addWidget(separator, 3, 0, 1, 3)  # span across all columns
        
        # Add OR label
        or_label = QLabel("OR")
        or_label.setAlignment(Qt.AlignCenter)
        font = or_label.font()
        font.setBold(True)
        font.setPointSize(12)
        or_label.setFont(font)
        load_layout.addWidget(or_label, 4, 0, 1, 3)  # span across all columns
        
        load_group.setLayout(load_layout)
        self.main_layout.addWidget(load_group)
        
    def create_action_buttons(self):
        """Create action buttons for launching the tool and managing settings."""
        button_layout = QHBoxLayout()
        
        # Save settings button
        self.save_settings_btn = QPushButton("Save Settings")
        self.save_settings_btn.clicked.connect(self.save_settings)
        button_layout.addWidget(self.save_settings_btn)
        
        # Spacer
        button_layout.addItem(QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))
        
        # Clear button
        self.clear_btn = QPushButton("Clear Fields")
        self.clear_btn.clicked.connect(self.clear_fields)
        button_layout.addWidget(self.clear_btn)
        
        # Launch button
        self.launch_btn = QPushButton("Launch Analysis Tool")
        self.launch_btn.setMinimumHeight(40)
        font = self.launch_btn.font()
        font.setBold(True)
        self.launch_btn.setFont(font)
        self.launch_btn.clicked.connect(self.launch_tool)
        button_layout.addWidget(self.launch_btn)
        
        self.main_layout.addLayout(button_layout)
    
    def browse_spectrum(self):
        """Open file dialog to browse for spectrum file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Spectrum File", "",
            "FITS Files (*.fits);;ASCII Files (*.txt *.dat);;All Files (*)"
        )
        if file_path:
            self.spectrum_path.setText(file_path)
            self.update_file_format(file_path)
    
    def browse_intervening(self):
        """Open file dialog to browse for intervening line list file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Intervening Line List", "",
            "Text Files (*.txt);;CSV Files (*.csv);;All Files (*)"
        )
        if file_path:
            self.intervening_path.setText(file_path)
    
    def browse_saved_file(self):
        """Open file dialog to browse for saved analysis files (pickle or JSON)."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Saved Analysis File", "",
            "Analysis Files (*.p *.json);;Pickle Files (*.p);;JSON Files (*.json);;All Files (*)"
        )
        if file_path:
            self.saved_file_path.setText(file_path)
    
    def update_file_format(self, file_path):
        """Update file format combobox based on file extension."""
        _, ext = os.path.splitext(file_path)
        ext = ext.lower()
        
        if ext == '.fits':
            self.file_format.setCurrentText("FITS")
        elif ext in ['.txt', '.dat']:
            self.file_format.setCurrentText("ASCII")
        elif ext == '.hdf5':
            self.file_format.setCurrentText("HDF5")
        else:
            self.file_format.setCurrentText("Custom")
    
    def apply_wavelength_preset(self, index):
        """Apply the selected wavelength preset."""
        if index == 0:
            return  # "Select a preset..." option
        
        presets = {
            1: "1215.67, 1025.72, 972.54, 949.74, 937.80",  # Lyman Series
            2: "1548.20, 1550.78",  # CIV
            3: "1393.76, 1402.77",  # SiIV
            4: "1031.93, 1037.62",  # OVI
            5: "2796.35, 2803.53"   # MgII
        }
        
        if index in presets:
            self.wavelengths.setText(presets[index])
            
        # Reset the combobox to "Select a preset..." after applying
        self.wavelength_presets.setCurrentIndex(0)
    
    def apply_dark_mode(self):
        """Apply dark mode styling to the application."""
        app = QApplication.instance()
        if app:
            # Create dark palette
            palette = QPalette()
            palette.setColor(QPalette.Window, QColor(53, 53, 53))
            palette.setColor(QPalette.WindowText, Qt.white)        
            palette.setColor(QPalette.Base, QColor(25, 25, 25))
            palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
            palette.setColor(QPalette.ToolTipBase, Qt.white)
            palette.setColor(QPalette.ToolTipText, Qt.white)
            palette.setColor(QPalette.Text, Qt.white)
            palette.setColor(QPalette.Button, QColor(53, 53, 53))
            palette.setColor(QPalette.ButtonText, Qt.white)
            palette.setColor(QPalette.BrightText, Qt.red)
            palette.setColor(QPalette.Link, QColor(42, 130, 218))
            palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
            palette.setColor(QPalette.HighlightedText, Qt.black)
            
            # Apply the palette
            app.setPalette(palette)
    
    def apply_light_mode(self):
        """Apply light mode styling to the application."""
        app = QApplication.instance()
        if app:
            app.setStyle("Fusion")  # Reset to default style
            app.setPalette(app.style().standardPalette())  # Reset to default palette
    
    def toggle_dark_mode(self, state):
        """Toggle between dark mode and light mode."""
        if state == Qt.Checked:
            self.apply_dark_mode()
        else:
            self.apply_light_mode()
    
    def save_settings(self):
        """Save current settings."""
        try:
            self.settings.setValue("geometry", self.saveGeometry())
            self.settings.setValue("spectrum_path", self.spectrum_path.text())
            self.settings.setValue("file_format", self.file_format.currentText())
            self.settings.setValue("intervening_path", self.intervening_path.text())
            self.settings.setValue("redshift", self.redshift.text())
            self.settings.setValue("wavelengths", self.wavelengths.text())
            self.settings.setValue("vmin", self.vmin.text())
            self.settings.setValue("vmax", self.vmax.text())
            self.settings.setValue("poly_order", self.poly_order.currentText())
            self.settings.setValue("mask_min", self.mask_min.text())
            self.settings.setValue("mask_max", self.mask_max.text())
            self.settings.setValue("dark_mode", self.dark_mode.isChecked())
            self.settings.setValue("saved_file_path", self.saved_file_path.text())
            
            QMessageBox.information(self, "Settings Saved", "Your settings have been saved successfully.")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to save settings: {str(e)}")
    
    def load_settings(self):
        """Load previously saved settings."""
        try:
            if self.settings.contains("spectrum_path"):
                self.spectrum_path.setText(self.settings.value("spectrum_path"))
            
            if self.settings.contains("file_format"):
                self.file_format.setCurrentText(self.settings.value("file_format"))
            
            if self.settings.contains("intervening_path"):
                self.intervening_path.setText(self.settings.value("intervening_path"))
            
            if self.settings.contains("redshift"):
                self.redshift.setText(self.settings.value("redshift"))
            
            if self.settings.contains("wavelengths"):
                self.wavelengths.setText(self.settings.value("wavelengths"))
            
            if self.settings.contains("vmin"):
                self.vmin.setText(self.settings.value("vmin"))
                
            if self.settings.contains("vmax"):
                self.vmax.setText(self.settings.value("vmax"))
            
            if self.settings.contains("poly_order"):
                self.poly_order.setCurrentText(self.settings.value("poly_order"))
            
            if self.settings.contains("mask_min"):
                self.mask_min.setText(self.settings.value("mask_min"))
                
            if self.settings.contains("mask_max"):
                self.mask_max.setText(self.settings.value("mask_max"))
            
            if self.settings.contains("dark_mode"):
                self.dark_mode.setChecked(self.settings.value("dark_mode", type=bool))
                
            if self.settings.contains("saved_file_path"):
                self.saved_file_path.setText(self.settings.value("saved_file_path"))
                
        except Exception as e:
            print(f"Error loading settings: {e}")
    
    def clear_fields(self):
        """Clear all input fields."""
        self.spectrum_path.clear()
        self.intervening_path.clear()
        self.redshift.clear()
        self.wavelengths.clear()
        self.vmin.setText("-2000")
        self.vmax.setText("2000")
        self.poly_order.setCurrentIndex(2)  # Reset to 4
        self.mask_min.setText("-200")
        self.mask_max.setText("200")
        self.saved_file_path.clear()
    
    def load_saved_file(self):
        """Load a saved analysis file (pickle or JSON) and launch the analysis tool."""
        file_path = self.saved_file_path.text().strip()
        
        if not file_path:
            QMessageBox.warning(self, "No File Selected", "Please select a saved analysis file.")
            return
            
        if not os.path.exists(file_path):
            QMessageBox.warning(self, "File Not Found", f"File not found: {file_path}")
            return
            
        try:
            ions = None
            
            # Determine file type based on extension
            if file_path.lower().endswith('.json'):
                if not JSON_SUPPORT:
                    QMessageBox.warning(
                        self, "JSON Support Not Available", 
                        "JSON utilities not found. Please install the json_utils.py module."
                    )
                    return
                    
                # Load from JSON
                ions = load_from_json(file_path)
                
            elif file_path.lower().endswith('.p'):
                # Load from pickle
                with open(file_path, 'rb') as f:
                    ions = pickle.load(f)
            else:
                QMessageBox.warning(
                    self, "Unknown File Type", 
                    "Unrecognized file extension. Please use .p for pickle files or .json for JSON files."
                )
                return
                
            if ions is None:
                QMessageBox.critical(self, "Loading Error", "Failed to load the analysis data.")
                return
                
            # Launch the analysis tool with the loaded data
            self.statusBar().showMessage("Launching analysis tool with loaded data...")
            
            # Hide the launcher window before starting the analysis tool
            self.hide()
            
            try:
                from GUIs.abstools import Metal_Plot as M
                M.Transitions(ions)
                
                # Show the launcher again after the analysis tool is closed
                self.show()
                self.statusBar().showMessage("Analysis tool closed")
                
            except Exception as e:
                # Make sure to show the launcher again if there's an error
                self.show()
                QMessageBox.critical(self, "Launch Error", f"Failed to launch analysis tool: {str(e)}")
                return
                
        except Exception as e:
            QMessageBox.critical(self, "Loading Error", f"Error loading file: {str(e)}")
    
    def launch_tool(self):
        """Validate inputs and launch the absorption line analysis tool."""
        # Validate required inputs
        if not self.spectrum_path.text().strip():
            QMessageBox.warning(self, "Missing Input", "Please select a spectrum file.")
            return
        
        if not self.redshift.text().strip():
            QMessageBox.warning(self, "Missing Input", "Please enter a redshift value.")
            return
        
        if not self.wavelengths.text().strip():
            QMessageBox.warning(self, "Missing Input", "Please enter rest-frame wavelengths.")
            return
        
        try:
            # Parse inputs
            z = float(self.redshift.text())
            
            # Parse wavelengths (handle various input formats)
            wavelengths_text = self.wavelengths.text().replace('[', '').replace(']', '')
            lines = [float(x.strip()) for x in wavelengths_text.split(',') if x.strip()]
            
            window_lim = [float(self.vmin.text()), float(self.vmax.text())]
            mask_init = [float(self.mask_min.text()), float(self.mask_max.text())]
            
            # Check if files exist
            spectrum_file = self.spectrum_path.text()
            if not os.path.exists(spectrum_file):
                QMessageBox.warning(self, "File Not Found", f"Spectrum file not found: {spectrum_file}")
                return
            
            intervening_file = None
            if self.intervening_path.text().strip():
                intervening_file = self.intervening_path.text()
                if not os.path.exists(intervening_file):
                    QMessageBox.warning(self, "File Not Found", f"Intervening line list file not found: {intervening_file}")
                    return
            
            # Generate and show the command to run
            if not self.generate_command(spectrum_file, z, lines, window_lim, mask_init, intervening_file):
                return
            
            # Execute the command
            self.execute_analysis(spectrum_file, z, lines, window_lim, mask_init, intervening_file)
            
        except ValueError as e:
            QMessageBox.warning(self, "Invalid Input", f"Please check your numerical inputs: {str(e)}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
    
    def generate_command(self, spectrum_file, z, lines, window_lim, mask_init, intervening_file):
        """Generate and display the command that will be run."""
        lines_str = str(lines).replace(' ', '')
        cmd = f"""
from linetools.spectra.xspectrum1d import XSpectrum1D
from GUIs.abstools import Absorber as A
from GUIs.abstools import Metal_Plot as M

# Read spectrum file
spectrum_file = '{spectrum_file}'
sp = XSpectrum1D.from_file(spectrum_file)
wave = sp.wavelength.value
flux = sp.flux.value
error = sp.sig.value

# Analysis parameters
z = {z}
lines = {lines_str}
window_lim = {window_lim}
mask_init = {mask_init}

# Create absorber and prepare for analysis
absys = A.Absorber(z, wave, flux, error, lines=lines, window_lim=window_lim, mask_init=mask_init)
Abs = absys.ions

# Launch analysis tool
"""
        
        if intervening_file:
            cmd += f"intervening_file = '{intervening_file}'\n"
            cmd += "M.Transitions(Abs, intervening=intervening_file)"
        else:
            cmd += "M.Transitions(Abs)"
        
        # Display the command in a message box for reference
        cmd_box = QMessageBox(self)
        cmd_box.setWindowTitle("Command to Execute")
        cmd_box.setText("The following command will be executed:")
        cmd_box.setDetailedText(cmd)
        cmd_box.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        
        if cmd_box.exec_() == QMessageBox.Cancel:
            return False
            
        return True
    
    def execute_analysis(self, spectrum_file, z, lines, window_lim, mask_init, intervening_file):
        """Execute the analysis using the provided parameters."""
        try:
            # Import required modules
            self.statusBar().showMessage("Importing required modules...")
            
            try:
                from linetools.spectra.xspectrum1d import XSpectrum1D
                from GUIs.abstools import Absorber as A
                from GUIs.abstools import Metal_Plot as M
                    
            except ImportError as e:
                QMessageBox.critical(self, "Import Error", 
                    f"Failed to import required modules: {str(e)}\n\n"
                    "Please ensure that linetools and the abstools packages are correctly installed."
                )
                return
            
            self.statusBar().showMessage("Reading spectrum file...")
            
            # Read spectrum file based on format
            file_format = self.file_format.currentText()
            try:
                if file_format == "FITS":
                    sp = XSpectrum1D.from_file(spectrum_file)
                    wave = sp.wavelength.value
                    flux = sp.flux.value
                    error = sp.sig.value
                elif file_format == "ASCII":
                    # Simple ASCII file reading
                    data = np.loadtxt(spectrum_file)
                    if data.shape[1] < 3:
                        QMessageBox.warning(self, "Format Error", 
                            "ASCII file must have at least 3 columns (wavelength, flux, error)")
                        return
                    wave = data[:, 0]
                    flux = data[:, 1]
                    error = data[:, 2]
                else:
                    # For other formats, prompt user for custom loading
                    QMessageBox.warning(self, "Custom Format",
                        "Custom format loading is not yet implemented. Please use FITS or ASCII formats.")
                    return
            except Exception as e:
                QMessageBox.critical(self, "File Reading Error", 
                    f"Failed to read spectrum file: {str(e)}")
                return
            
            self.statusBar().showMessage("Creating absorber...")
            
            # Create absorber
            try:
                absys = A.Absorber(z, wave, flux, error, lines=lines, 
                                   window_lim=window_lim, mask_init=mask_init)
                Abs = absys.ions
            except Exception as e:
                QMessageBox.critical(self, "Absorber Error", 
                    f"Failed to create absorber: {str(e)}")
                return
            
            self.statusBar().showMessage("Launching analysis tool...")
            
            # Hide the launcher window before starting the analysis tool
            self.hide()

            # Launch the analysis tool
            try:
                if intervening_file:
                    M.Transitions(Abs, intervening=intervening_file)
                else:
                    M.Transitions(Abs)
                
                # Show the launcher again after the analysis tool is closed
                self.show()
            except Exception as e:
                # Make sure to show the launcher again if there's an error
                self.show()
                QMessageBox.critical(self, "Launch Error", 
                    f"Failed to launch analysis tool: {str(e)}")
                return
            
            self.statusBar().showMessage("Analysis tool launched successfully")
            
        except Exception as e:
            QMessageBox.critical(self, "Execution Error", f"An error occurred: {str(e)}")
    
        
    def closeEvent(self, event):
        """Handle window close event to save settings."""
        self.save_settings()
        event.accept()