#!/usr/bin/env python
"""
AbsTools Launcher GUI

A PyQt5 GUI for setting up and launching the AbsTools absorption line analysis toolbox.
This launcher simplifies the process of loading spectra, setting redshift and lines to analyze,
and launching the Metal_Plot visualization tool.

Usage:
    python abstools_launcher.py
"""

import os
import sys
import pickle
import numpy as np
import time
import subprocess
import tempfile
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QLabel, QLineEdit, QPushButton,
    QVBoxLayout, QHBoxLayout, QGridLayout, QFileDialog, QMessageBox,
    QCheckBox, QGroupBox, QComboBox, QTextEdit, QSplitter, QFrame,
    QTableWidget, QTableWidgetItem, QHeaderView
)
from PyQt5.QtCore import Qt, QSettings, QTimer

# Global constants
DEFAULT_LINES = "1031.93, 1037.62, 1215.67, 1548.20, 1550.78"
DEFAULT_REDSHIFT = "0.348"
DEFAULT_WINDOW_LIMITS = "-2000, 2000"

class AbsToolsLauncher(QMainWindow):
    """
    Main window for the AbsTools launcher GUI.
    """
    
    def __init__(self):
        super().__init__()
        
        # Set up window properties
        self.setWindowTitle("AbsTools Launcher")
        self.setMinimumSize(600, 780)
        
        # Initialize settings
        self.settings = QSettings("AbsTools", "Launcher")
        
        # Set up the main UI
        self.setup_ui()
        
        # Load saved settings
        self.load_settings()
        
        # Show the window
        self.show()
    
    def setup_ui(self):
        """Set up the user interface components."""
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QVBoxLayout(central_widget)
        
        # Create a splitter for flexible layout
        splitter = QSplitter(Qt.Vertical)
        main_layout.addWidget(splitter)
        
        # TOP SECTION: Input parameters
        top_widget = QWidget()
        top_layout = QVBoxLayout(top_widget)
        
        # Add title
        title_label = QLabel("AbsTools Launcher")
        title_label.setStyleSheet("font-size: 18pt; font-weight: bold; margin-bottom: 10px;")
        title_label.setAlignment(Qt.AlignCenter)
        top_layout.addWidget(title_label)
        
        # Add description
        desc_label = QLabel(
            "Configure and launch the AbsTools absorption line analysis toolbox.\n"
            "Load a spectrum file, set the redshift, specify absorption lines to analyze, and start the analysis."
        )
        desc_label.setStyleSheet("font-size: 10pt; margin-bottom: 20px;")
        desc_label.setAlignment(Qt.AlignCenter)
        top_layout.addWidget(desc_label)
        
        # Create groups for different input methods
        tabs = QtWidgets.QTabWidget()
        top_layout.addWidget(tabs)
        
        # === Tab 1: Load New Spectrum ===
        new_spectrum_tab = QWidget()
        tabs.addTab(new_spectrum_tab, "Load New Spectrum")
        
        new_spectrum_layout = QGridLayout(new_spectrum_tab)
        
        # Spectrum file selection
        new_spectrum_layout.addWidget(QLabel("Spectrum File:"), 0, 0)
        self.spectrum_path_edit = QLineEdit()
        new_spectrum_layout.addWidget(self.spectrum_path_edit, 0, 1)
        
        browse_button = QPushButton("Browse...")
        browse_button.clicked.connect(self.browse_spectrum)
        new_spectrum_layout.addWidget(browse_button, 0, 2)
        
        # Format selection
        new_spectrum_layout.addWidget(QLabel("Format:"), 1, 0)
        self.format_combo = QComboBox()
        self.format_combo.addItems(["Auto Detect", "FITS", "ASCII", "HDF5"])
        new_spectrum_layout.addWidget(self.format_combo, 1, 1)
        
        # Redshift
        new_spectrum_layout.addWidget(QLabel("Redshift:"), 2, 0)
        self.redshift_edit = QLineEdit(DEFAULT_REDSHIFT)
        new_spectrum_layout.addWidget(self.redshift_edit, 2, 1)
        
        # Absorption lines
        new_spectrum_layout.addWidget(QLabel("Absorption Lines:"), 3, 0)
        self.lines_edit = QLineEdit(DEFAULT_LINES)
        new_spectrum_layout.addWidget(self.lines_edit, 3, 1)
        
        # Window limits
        new_spectrum_layout.addWidget(QLabel("Window Limits (km/s):"), 4, 0)
        self.window_limits_edit = QLineEdit(DEFAULT_WINDOW_LIMITS)
        new_spectrum_layout.addWidget(self.window_limits_edit, 4, 1)
        
        # Intervening absorbers file
        new_spectrum_layout.addWidget(QLabel("Intervening Absorbers File:"), 5, 0)
        self.intervening_edit = QLineEdit()
        new_spectrum_layout.addWidget(self.intervening_edit, 5, 1)
        
        intervening_browse = QPushButton("Browse...")
        intervening_browse.clicked.connect(self.browse_intervening)
        new_spectrum_layout.addWidget(intervening_browse, 5, 2)
        
        # Add line presets with improved interface
        preset_group = QGroupBox("Absorption Line Selection")
        preset_layout = QVBoxLayout(preset_group)
        
        # Current selection area with text box
        selection_layout = QHBoxLayout()
        selection_layout.addWidget(QLabel("Current Lines:"))
        
        # This text box replaces the original lines_edit
        self.lines_display = QLineEdit(self.lines_edit.text())
        self.lines_display.setReadOnly(False)  # Allow direct editing
        self.lines_display.textChanged.connect(self._sync_line_fields)
        selection_layout.addWidget(self.lines_display)
        
        preset_layout.addLayout(selection_layout)
        
        # Preset dropdown
        preset_layout.addWidget(QLabel("Add from presets:"))
        
        dropdown_button_layout = QHBoxLayout()
        
        self.preset_combo = QComboBox()
        self.preset_combo.setMinimumWidth(300)
        self.preset_combo.addItem("Select a preset...", "")
        self.preset_combo.addItem("OVI Doublet", "1031.93, 1037.62")
        self.preset_combo.addItem("CIV Doublet", "1548.20, 1550.78")
        self.preset_combo.addItem("Lyman-α", "1215.67")
        self.preset_combo.addItem("MgII Doublet", "2796.35, 2803.53")
        self.preset_combo.addItem("SiIV Doublet", "1393.76, 1402.77")
        self.preset_combo.addItem("NV Doublet", "1238.82, 1242.80")
        self.preset_combo.addItem("FeII Lines", "2344.21, 2374.46, 2382.77, 2586.65, 2600.17")
        self.preset_combo.addItem("Common IGM Lines", "1215.67, 1031.93, 1037.62, 1548.20, 1550.78, 1393.76, 1402.77")
        
        dropdown_button_layout.addWidget(self.preset_combo)
        
        add_preset_button = QPushButton("Add")
        add_preset_button.clicked.connect(lambda: self.apply_preset(append=True))
        dropdown_button_layout.addWidget(add_preset_button)
        
        replace_preset_button = QPushButton("Replace")
        replace_preset_button.clicked.connect(self.apply_preset)
        dropdown_button_layout.addWidget(replace_preset_button)
        
        preset_layout.addLayout(dropdown_button_layout)
        
        # Add a help text
        help_text = QLabel("You can directly edit the lines above or select from presets.")
        help_text.setStyleSheet("font-style: italic; color: gray;")
        preset_layout.addWidget(help_text)
        
        new_spectrum_layout.addWidget(preset_group, 6, 0, 1, 3)
        
        # Hide the original lines edit (we'll keep it synced with the new display)
        self.lines_edit.setVisible(False)
        
        # === Tab 2: Load Saved Analysis ===
        saved_analysis_tab = QWidget()
        tabs.addTab(saved_analysis_tab, "Load Saved Analysis")
        
        saved_analysis_layout = QGridLayout(saved_analysis_tab)
        
        # Analysis file selection
        saved_analysis_layout.addWidget(QLabel("Analysis File:"), 0, 0)
        self.analysis_path_edit = QLineEdit()
        saved_analysis_layout.addWidget(self.analysis_path_edit, 0, 1)
        
        analysis_browse_button = QPushButton("Browse...")
        analysis_browse_button.clicked.connect(self.browse_analysis)
        saved_analysis_layout.addWidget(analysis_browse_button, 0, 2)
        
        # File format info
        saved_analysis_layout.addWidget(QLabel("File Format:"), 1, 0)
        self.format_label = QLabel("None selected")
        saved_analysis_layout.addWidget(self.format_label, 1, 1)
        
        # Intervening absorbers file
        saved_analysis_layout.addWidget(QLabel("Intervening Absorbers File:"), 2, 0)
        self.saved_intervening_edit = QLineEdit()
        saved_analysis_layout.addWidget(self.saved_intervening_edit, 2, 1)
        
        saved_intervening_browse = QPushButton("Browse...")
        saved_intervening_browse.clicked.connect(self.browse_saved_intervening)
        saved_analysis_layout.addWidget(saved_intervening_browse, 2, 2)
        
        # Add info about file formats
        format_info = QTextEdit()
        format_info.setReadOnly(True)
        format_info.setMaximumHeight(100)
        format_info.setText(
            "Supported analysis file formats:\n"
            "- Pickle (.p): Binary format saved from AbsTools\n"
            "- JSON (.json): JSON format saved from AbsTools\n\n"
            "These files contain previous analysis sessions that can be resumed."
        )
        saved_analysis_layout.addWidget(format_info, 3, 0, 1, 3)
        
        # Add the top widget to the splitter
        splitter.addWidget(top_widget)
        
        # BOTTOM SECTION: Controls and status
        bottom_widget = QWidget()
        bottom_layout = QVBoxLayout(bottom_widget)
        
        # Add a separator
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        bottom_layout.addWidget(line)
        
        # Status area
        status_group = QGroupBox("Status")
        status_layout = QVBoxLayout(status_group)
        
        self.status_text = QTextEdit()
        self.status_text.setReadOnly(True)
        self.status_text.setMaximumHeight(100)
        self.status_text.setText("Ready to start AbsTools. Configure options above and click Launch.")
        status_layout.addWidget(self.status_text)
        
        bottom_layout.addWidget(status_group)
        
        # Launch button area
        button_layout = QHBoxLayout()
        
        # Save settings button
        save_settings_button = QPushButton("Save Settings")
        save_settings_button.clicked.connect(self.save_settings)
        button_layout.addWidget(save_settings_button)
        
        # Spacer
        button_layout.addStretch()
        
        # Launch button
        self.launch_button = QPushButton("Launch AbsTools")
        self.launch_button.setStyleSheet("font-size: 12pt; font-weight: bold; padding: 8px;")
        self.launch_button.clicked.connect(self.launch_abstools)
        button_layout.addWidget(self.launch_button)
        
        bottom_layout.addLayout(button_layout)
        
        # Add the bottom widget to the splitter
        splitter.addWidget(bottom_widget)
        
        # Set the initial sizes for the splitter
        splitter.setSizes([400, 200])
    
    def browse_spectrum(self):
        """Open a file dialog to select a spectrum file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Spectrum File", "",
            "Spectrum Files (*.fits *.fit *.txt *.dat *.hdf5);;All Files (*)"
        )
        
        if file_path:
            self.spectrum_path_edit.setText(file_path)
            self.update_status(f"Selected spectrum file: {file_path}")
    
    def browse_intervening(self):
        """Open a file dialog to select an intervening absorbers file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Intervening Absorbers File", "",
            "Text Files (*.txt *.dat);;All Files (*)"
        )
        
        if file_path:
            self.intervening_edit.setText(file_path)
            self.update_status(f"Selected intervening absorbers file: {file_path}")
    
    def browse_saved_intervening(self):
        """Open a file dialog to select an intervening absorbers file for saved analysis."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Intervening Absorbers File", "",
            "Text Files (*.txt *.dat);;All Files (*)"
        )
        
        if file_path:
            self.saved_intervening_edit.setText(file_path)
            self.update_status(f"Selected intervening absorbers file: {file_path}")
    
    def browse_analysis(self):
        """Open a file dialog to select a saved analysis file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Analysis File", "",
            "Analysis Files (*.p *.json);;Pickle Files (*.p);;JSON Files (*.json);;All Files (*)"
        )
        
        if file_path:
            self.analysis_path_edit.setText(file_path)
            
            # Determine and display the format
            if file_path.lower().endswith('.p'):
                self.format_label.setText("Pickle (.p)")
            elif file_path.lower().endswith('.json'):
                self.format_label.setText("JSON (.json)")
            else:
                self.format_label.setText("Unknown format")
            
            self.update_status(f"Selected analysis file: {file_path}")
    
    def _sync_line_fields(self):
        """Keep the hidden lines_edit and visible lines_display in sync."""
        # Update the original field from the display field
        self.lines_edit.setText(self.lines_display.text())
    
    def apply_preset(self, append=False):
        """Apply the selected preset to the absorption lines field."""
        preset = self.preset_combo.currentData()
        if not preset:
            return
            
        if append:
            # Get current lines
            current_lines = self.lines_display.text().strip()
            if current_lines:
                # Append the new preset
                combined = current_lines + ", " + preset
                self.lines_display.setText(combined)
                self.update_status(f"Added {self.preset_combo.currentText()} lines to current selection")
            else:
                # If no current lines, just set the preset
                self.lines_display.setText(preset)
                self.update_status(f"Set absorption lines to: {preset}")
        else:
            # Replace current lines with the preset
            self.lines_display.setText(preset)
            self.update_status(f"Set absorption lines to: {preset}")
            
        # Reset the dropdown to the first item
        self.preset_combo.setCurrentIndex(0)
    
    def update_status(self, message):
        """Update the status text area with a new message."""
        current_text = self.status_text.toPlainText()
        self.status_text.setText(f"{message}\n{current_text}")
    
    def save_settings(self):
        """Save current settings to QSettings."""
        self.settings.setValue("spectrum_path", self.spectrum_path_edit.text())
        self.settings.setValue("format", self.format_combo.currentText())
        self.settings.setValue("redshift", self.redshift_edit.text())
        self.settings.setValue("lines", self.lines_display.text())  # Use the visible field
        self.settings.setValue("window_limits", self.window_limits_edit.text())
        self.settings.setValue("intervening", self.intervening_edit.text())
        self.settings.setValue("analysis_path", self.analysis_path_edit.text())
        self.settings.setValue("saved_intervening", self.saved_intervening_edit.text())
        
        self.update_status("Settings saved")
    
    def load_settings(self):
        """Load saved settings from QSettings."""
        if self.settings.contains("spectrum_path"):
            self.spectrum_path_edit.setText(self.settings.value("spectrum_path"))
        
        if self.settings.contains("format"):
            index = self.format_combo.findText(self.settings.value("format"))
            if index >= 0:
                self.format_combo.setCurrentIndex(index)
        
        if self.settings.contains("redshift"):
            self.redshift_edit.setText(self.settings.value("redshift"))
        
        if self.settings.contains("lines"):
            lines_value = self.settings.value("lines")
            self.lines_edit.setText(lines_value)
            self.lines_display.setText(lines_value)  # Update both fields
        
        if self.settings.contains("window_limits"):
            self.window_limits_edit.setText(self.settings.value("window_limits"))
        
        if self.settings.contains("intervening"):
            self.intervening_edit.setText(self.settings.value("intervening"))
        
        if self.settings.contains("analysis_path"):
            self.analysis_path_edit.setText(self.settings.value("analysis_path"))
            
            # Update format label based on the file extension
            path = self.settings.value("analysis_path")
            if path.lower().endswith('.p'):
                self.format_label.setText("Pickle (.p)")
            elif path.lower().endswith('.json'):
                self.format_label.setText("JSON (.json)")
            elif path:
                self.format_label.setText("Unknown format")
        
        if self.settings.contains("saved_intervening"):
            self.saved_intervening_edit.setText(self.settings.value("saved_intervening"))
    
    def clean_up_ui(self):
        """Disconnect signals and prepare for clean exit."""
        # Save settings before exit
        self.save_settings()
        
        # Disconnect all signals
        # This helps prevent issues when closing
        for child in self.findChildren(QtCore.QObject):
            try:
                # Try to disconnect all signals for each widget
                if hasattr(child, 'disconnect'):
                    child.disconnect()
            except:
                pass
    
    def launch_abstools(self):
        """Launch AbsTools with the configured settings."""
        # Get the current tab index to determine the launch mode
        current_tab = self.centralWidget().findChild(QtWidgets.QTabWidget).currentIndex()
        
        try:
            # Launch mode: New spectrum
            if current_tab == 0:
                spectrum_path = self.spectrum_path_edit.text()
                if not spectrum_path:
                    QMessageBox.warning(self, "Input Error", "Please select a spectrum file.")
                    return
                
                if not os.path.exists(spectrum_path):
                    QMessageBox.warning(self, "File Error", f"Spectrum file not found: {spectrum_path}")
                    return
                
                try:
                    redshift = float(self.redshift_edit.text())
                except ValueError:
                    QMessageBox.warning(self, "Input Error", "Redshift must be a valid number.")
                    return
                
                # Parse lines
                try:
                    lines_text = self.lines_edit.text().strip()
                    if not lines_text:
                        QMessageBox.warning(self, "Input Error", "Please specify at least one absorption line.")
                        return
                    
                    lines = [float(line.strip()) for line in lines_text.split(',')]
                    if not lines:
                        QMessageBox.warning(self, "Input Error", "Please specify at least one absorption line.")
                        return
                except ValueError:
                    QMessageBox.warning(self, "Input Error", "Invalid absorption line format. Use comma-separated wavelengths (e.g., 1031.93, 1037.62).")
                    return
                
                # Parse window limits
                try:
                    window_limits_text = self.window_limits_edit.text().strip()
                    window_limits = [int(limit.strip()) for limit in window_limits_text.split(',')]
                    if len(window_limits) != 2:
                        QMessageBox.warning(self, "Input Error", "Window limits must be two comma-separated values (e.g., -2000, 2000).")
                        return
                except ValueError:
                    QMessageBox.warning(self, "Input Error", "Invalid window limit format. Use comma-separated integers (e.g., -2000, 2000).")
                    return
                
                # Check intervening absorbers file
                intervening_path = self.intervening_edit.text()
                if intervening_path and not os.path.exists(intervening_path):
                    QMessageBox.warning(self, "File Error", f"Intervening absorbers file not found: {intervening_path}")
                    return
                
                # Update status
                self.update_status("Loading spectrum and initializing AbsTools...")
                
                # Clean up the UI before launching
                self.clean_up_ui()
                
                # Launch with delay to allow clean shutdown of this window
                self.update_status("Preparing to launch Metal_Plot...")
                QTimer.singleShot(100, lambda: self.run_abstools_new_spectrum(
                    spectrum_path,
                    self.format_combo.currentText(),
                    redshift,
                    lines,
                    window_limits,
                    intervening_path
                ))
            
            # Launch mode: Load saved analysis
            elif current_tab == 1:
                analysis_path = self.analysis_path_edit.text()
                if not analysis_path:
                    QMessageBox.warning(self, "Input Error", "Please select an analysis file.")
                    return
                
                if not os.path.exists(analysis_path):
                    QMessageBox.warning(self, "File Error", f"Analysis file not found: {analysis_path}")
                    return
                
                # Check intervening absorbers file
                intervening_path = self.saved_intervening_edit.text()
                if intervening_path and not os.path.exists(intervening_path):
                    QMessageBox.warning(self, "File Error", f"Intervening absorbers file not found: {intervening_path}")
                    return
                
                # Update status
                self.update_status("Loading saved analysis and initializing AbsTools...")
                
                # Clean up the UI before launching
                self.clean_up_ui()
                
                # Launch with delay to allow clean shutdown of this window
                self.update_status("Preparing to launch Metal_Plot...")
                QTimer.singleShot(100, lambda: self.run_abstools_load_analysis(analysis_path, intervening_path))
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while launching AbsTools: {str(e)}")
            self.update_status(f"Error: {str(e)}")
    
    def launch_metal_plot_separately(self, data_file, is_json=False, intervening_path=None):
        """
        Launch Metal_Plot in a separate process to avoid event loop conflicts.
        
        Parameters:
        -----------
        data_file : str
            Path to the data file (either a temporary pickle file or a JSON file)
        is_json : bool
            Whether the file is a JSON file
        intervening_path : str, optional
            Path to intervening absorbers file
        """
        # Create a temporary Python script file instead of passing a long command
        with tempfile.NamedTemporaryFile(suffix='.py', delete=False, mode='w') as script_file:
            script_path = script_file.name
            
            # Write the launcher script
            if is_json:
                script_file.write(f"""
import sys
import os
import json
print('Starting Metal_Plot with JSON file: {data_file}')

# Always use package imports
try:
    from rbcodes.GUIs.abstools.json_utils import load_from_json
    from rbcodes.GUIs.abstools import Metal_Plot as M
    print('Successfully imported from rbcodes package')
    
    try:
        ions = load_from_json('{data_file}')
        print('JSON data loaded successfully')
        print('Launching Metal_Plot...')
        
        # Launch Metal_Plot
        M.Transitions(ions, intervening='{intervening_path}' if '{intervening_path}' != 'None' else False)
    except Exception as e:
        print(f'Error loading or running with JSON: {{e}}')
        import traceback
        traceback.print_exc()
        sys.exit(1)
        
except ImportError as e1:
    print(f'Package import failed: {{e1}}')
    sys.exit(1)
""")
            else:
                script_file.write(f"""
import sys
import os
import pickle
print('Starting Metal_Plot with pickle file: {data_file}')

# Always use package imports
try:
    from rbcodes.GUIs.abstools import Metal_Plot as M
    print('Successfully imported from rbcodes package')
    
    try:
        with open('{data_file}', 'rb') as f:
            ions = pickle.load(f)
        print('Pickle data loaded successfully')
        print('Launching Metal_Plot...')
        
        # Launch Metal_Plot
        M.Transitions(ions, intervening='{intervening_path}' if '{intervening_path}' != 'None' else False)
    except Exception as e:
        print(f'Error loading or running with Pickle: {{e}}')
        import traceback
        traceback.print_exc()
        sys.exit(1)
        
except ImportError as e1:
    print(f'Package import failed: {{e1}}')
    sys.exit(1)
""")
        
        # Print the script path for debugging
        print(f"Created launcher script at: {script_path}")
        
        # Hide the launcher
        self.hide()
        
        # Allow time for the window to hide before executing script
        QTimer.singleShot(100, lambda: self._execute_script(script_path, data_file))
    
    def _execute_script(self, script_path, data_file):
        """Execute a Python script file in a separate process."""
        try:
            print(f"Executing Python script: {script_path}")
            
            # Use a more reliable method to launch the subprocess
            # Run in the foreground (wait=True) to catch any immediate errors
            result = subprocess.run([sys.executable, script_path], 
                                     capture_output=True, text=True)
            
            # Check for errors
            if result.returncode != 0:
                print(f"Error launching Metal_Plot: {result.stderr}")
                print(f"Output: {result.stdout}")
                # Don't exit immediately, show the error
                QMessageBox.critical(None, "Error", 
                                   f"Failed to launch Metal_Plot:\n{result.stderr}\n\nOutput:\n{result.stdout}")
                self.show()  # Show the launcher again
            else:
                print(f"Successfully launched Metal_Plot")
                print(f"Output: {result.stdout}")
                # If successful, now exit this process
                sys.exit(0)
                
        except Exception as e:
            print(f"Error launching Metal_Plot: {e}")
            import traceback
            traceback.print_exc()
            # Don't exit immediately, show the error
            QMessageBox.critical(None, "Error", 
                               f"Failed to launch Metal_Plot:\n{str(e)}")
            self.show()  # Show the launcher again
    
    def run_abstools_new_spectrum(self, spectrum_path, format_name, redshift, lines, window_limits, intervening_path=None):
        """
        Run AbsTools with a new spectrum by creating an Absorber object and saving to a temporary file.
        Then launch Metal_Plot in a separate process.
        
        Parameters:
        -----------
        spectrum_path : str
            Path to the spectrum file
        format_name : str
            Format of the spectrum file
        redshift : float
            Redshift at which to perform analysis
        lines : list
            List of absorption line wavelengths to analyze
        window_limits : list
            List of [min, max] velocity window limits in km/s
        intervening_path : str, optional
            Path to intervening absorbers file
        """
        try:
            # Create a temporary file to store the data
            with tempfile.NamedTemporaryFile(suffix='.p', delete=False) as tmp_file:
                temp_path = tmp_file.name
            
            # Print debug info
            print(f"Creating temporary file: {temp_path}")
            
            # Create the Absorber object
            try:
                # Always use package imports
                from linetools.spectra.xspectrum1d import XSpectrum1D
                from rbcodes.GUIs.abstools import Absorber as A
                
                # Read spectrum data using XSpectrum1D
                print(f"Reading spectrum file: {spectrum_path}")
                sp = XSpectrum1D.from_file(spectrum_path)
                wave = sp.wavelength.value
                flux = sp.flux.value
                error = sp.sig.value
                
                # Create absorber object and save to file
                print(f"Creating Absorber object with redshift={redshift}, lines={lines}, window_limits={window_limits}")
                absys = A.Absorber(redshift, wave, flux, error, lines=lines, window_lim=window_limits)
                
                print(f"Saving ions to {temp_path}")
                with open(temp_path, 'wb') as f:
                    pickle.dump(absys.ions, f)
                print("Successful pickle save")
                
                # Launch Metal_Plot in a separate process
                print("Launching Metal_Plot with the saved data")
                self.launch_metal_plot_separately(temp_path, is_json=False, intervening_path=intervening_path)
                
            except Exception as e:
                self.show()  # Show the window again on error
                self.update_status(f"Error creating Absorber object: {str(e)}")
                import traceback
                traceback.print_exc()
                raise e
            
        except Exception as e:
            self.show()  # Ensure window is visible on error
            self.update_status(f"Error: {str(e)}")
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
    
    def run_abstools_load_analysis(self, analysis_path, intervening_path=None):
        """
        Run AbsTools with a saved analysis file.
        Launch Metal_Plot in a separate process to avoid event loop conflicts.
        
        Parameters:
        -----------
        analysis_path : str
            Path to the saved analysis file
        intervening_path : str, optional
            Path to intervening absorbers file
        """
        try:
            # For JSON files, we can use the file directly
            if analysis_path.lower().endswith('.json'):
                print(f"Loading JSON file: {analysis_path}")
                self.launch_metal_plot_separately(analysis_path, is_json=True, intervening_path=intervening_path)
            
            # For pickle files, we can use the file directly
            elif analysis_path.lower().endswith('.p'):
                print(f"Loading pickle file: {analysis_path}")
                self.launch_metal_plot_separately(analysis_path, is_json=False, intervening_path=intervening_path)
            
            else:
                self.show()  # Show the window again on error
                self.update_status(f"Unsupported file format: {analysis_path}")
                QMessageBox.warning(self, "Format Error", f"Unsupported file format: {analysis_path}")
            
        except Exception as e:
            self.show()  # Ensure window is visible on error
            self.update_status(f"Error: {str(e)}")
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")

def run_launcher():
    """Run the AbsTools launcher application."""
    # Ensure only one QApplication instance exists
    app = QApplication.instance() or QApplication(sys.argv)
    launcher = AbsToolsLauncher()
    sys.exit(app.exec_())

if __name__ == "__main__":
    run_launcher()