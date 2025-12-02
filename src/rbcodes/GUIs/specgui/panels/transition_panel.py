import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QDoubleSpinBox, QFormLayout, QGroupBox,
                           QMessageBox, QComboBox, QSpinBox, QCheckBox)
from PyQt5.QtCore import pyqtSignal

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots."""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)


class TransitionPanel(QWidget):
    """Panel for selecting a transition and slicing the spectrum."""
    
    # Signals
    spectrum_sliced = pyqtSignal(float, float, float)  # Emitted when spectrum is sliced (transition, vmin, vmax)
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        
        # Default values
        self.transition = 1215.67  # Lyman-alpha by default
        self.vmin = -1500
        self.vmax = 1500
        
        self.init_ui()
        
        # Connect to controller signals
        self.controller.spectrum_changed.connect(self.update_plot)
        self.controller.spectrum_changed.connect(self.reset)

    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Transition selection
        transition_group = QGroupBox("Transition Selection")
        transition_layout = QFormLayout()
        
        # Common transitions dropdown
        self.transition_combo = QComboBox()
        self.populate_transition_combo()
        self.transition_combo.currentTextChanged.connect(self.on_transition_selected)
        
        # Custom transition input
        self.transition_spinbox = QDoubleSpinBox()
        self.transition_spinbox.setRange(1.0, 10000.0)
        self.transition_spinbox.setDecimals(2)
        self.transition_spinbox.setValue(self.transition)
        self.transition_spinbox.setSingleStep(1.0)
        self.transition_spinbox.valueChanged.connect(self.on_transition_changed)
        
        transition_input_layout = QHBoxLayout()
        transition_input_layout.addWidget(self.transition_combo)
        transition_input_layout.addWidget(self.transition_spinbox)
        
        transition_layout.addRow("Select Transition:", transition_input_layout)
        
        # Velocity range
        velocity_layout = QHBoxLayout()
        
        self.vmin_spinbox = QSpinBox()
        self.vmin_spinbox.setRange(-10000, 0)
        self.vmin_spinbox.setValue(self.vmin)
        self.vmin_spinbox.setSingleStep(100)
        self.vmin_spinbox.valueChanged.connect(self.on_vmin_changed)
        
        self.vmax_spinbox = QSpinBox()
        self.vmax_spinbox.setRange(0, 10000)
        self.vmax_spinbox.setValue(self.vmax)
        self.vmax_spinbox.setSingleStep(100)
        self.vmax_spinbox.valueChanged.connect(self.on_vmax_changed)
        
        velocity_layout.addWidget(QLabel("Min:"))
        velocity_layout.addWidget(self.vmin_spinbox)
        velocity_layout.addWidget(QLabel("Max:"))
        velocity_layout.addWidget(self.vmax_spinbox)
        
        transition_layout.addRow("Velocity Range (km/s):", velocity_layout)
        
        # Use velocity checkbox
        self.use_vel_checkbox = QCheckBox("Use Velocity")
        self.use_vel_checkbox.setChecked(True)
        transition_layout.addRow("", self.use_vel_checkbox)
        
        # Linelist selection
        self.linelist_combo = QComboBox()
        self.linelist_combo.addItems(["atom", "LLS", "LLS Small", "DLA", "LBG", "Gal"])
        transition_layout.addRow("Line List:", self.linelist_combo)
        
        # Slice button
        self.slice_btn = QPushButton("Slice Spectrum")
        self.slice_btn.clicked.connect(self.slice_spectrum)
        self.slice_btn.setToolTip("Extract a portion of the spectrum around the selected transition")
    
        transition_layout.addRow("", self.slice_btn)
        
        # Status label
        self.status_label = QLabel("No slicing applied")
        transition_layout.addRow("Status:", self.status_label)
   
        # Add New Transition button
        new_transition_layout = QHBoxLayout()
        self.new_transition_btn = QPushButton("New Transition")
        self.new_transition_btn.setToolTip("Reset current transition selection and analyze a different transition at the same redshift")
        self.new_transition_btn.clicked.connect(self.reset_transition)
        
        # Set a pale red color for the new transition button
        self.new_transition_btn.setStyleSheet("""
            QPushButton {
                background-color: #ffcccc; 
                color: #990000; 
                border: 1px solid #cc0000;
                padding: 5px;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #ffb3b3;
            }
            QPushButton:pressed {
                background-color: #ff9999;
            }
        """)
        
        # Add the button to your layout, positioned to the right
        transition_button_layout = QHBoxLayout()
        transition_button_layout.addWidget(self.slice_btn)
        transition_button_layout.addWidget(self.new_transition_btn)


        new_transition_layout.addWidget(self.new_transition_btn)
        
        transition_layout.addRow("", new_transition_layout)

        # Tooltips for dropdown and spinboxes
        self.transition_combo.setToolTip("Select a common transition or choose 'Custom'")
        self.transition_spinbox.setToolTip("Wavelength of the transition in Angstroms")
        self.vmin_spinbox.setToolTip("Minimum velocity for slicing (km/s)")
        self.vmax_spinbox.setToolTip("Maximum velocity for slicing (km/s)")
        self.use_vel_checkbox.setToolTip("Use velocity space for slicing instead of wavelength")
        self.linelist_combo.setToolTip("Select the atomic line list to use")
        
        transition_group.setLayout(transition_layout)
        main_layout.addWidget(transition_group)
        
        # Spectrum plot
        plot_group = QGroupBox("Spectrum Preview")
        plot_layout = QVBoxLayout()
        
        self.canvas = MatplotlibCanvas(self, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        
        plot_group.setLayout(plot_layout)
        main_layout.addWidget(plot_group)

    def reset(self):
        """Reset panel state when a new file is loaded."""
        self.status_label.setText("No slicing applied")
        # Only reset transition selector if we loaded from a JSON file with transition info
        if hasattr(self.controller.spec, 'trans_wave') and hasattr(self.controller.spec, 'trans'):
            # JSON file with transition info - keep what was loaded
            pass
        # For regular files, don't reset the dropdown - keep user's selection
        self.update_plot()

    def reset_transition(self):
        """Reset the transition selection and clear downstream analysis."""
        # Ask for confirmation
        reply = QMessageBox.question(self, 'Confirm Reset', 
                                   'Reset transition? This will clear all downstream analysis steps.',
                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        
        if reply == QMessageBox.No:
            return
        
        # Reset UI elements to default values
        self.transition_combo.setCurrentIndex(0)  # Reset to first transition
        self.vmin_spinbox.setValue(self.vmin)     # Reset to default velocity range
        self.vmax_spinbox.setValue(self.vmax)
        self.status_label.setText("Transition reset")
        
        # Signal to the main application that transition was reset
        self.spectrum_sliced.emit(-99999, 0, 0)  # Special values to indicate reset
    

    
    def populate_transition_combo(self):
        """Populate the transition combo box with common transitions."""
        transitions = [
            "Custom",
            "Lyα (1215.67 Å)",
            "Lyβ (1025.72 Å)",
            "Lyγ (972.54 Å)",
            "CII (1036.3367 Å)",
            "CII (1334.532 Å)",
            "CIV (1548.20 Å)",
            "CIV (1550.78 Å)",
            "MgII (2796.35 Å)",
            "MgII (2803.53 Å)",
            "SiIV (1393.76 Å)",
            "SiIV (1402.77 Å)",
            "OVI (1031.926 Å)",
            "OVI (1037.617 Å)"
        ]
        self.transition_combo.addItems(transitions)
        
        # Set default to Lyα
        self.transition_combo.setCurrentIndex(0)
    
    def on_transition_selected(self, text):
        """Handle selection from common transitions dropdown."""
        if text == "Custom":
            return  # Keep the current custom value
        
        # Extract wavelength from text (format: "Name (XXXX.XX Å)")
        try:
            wavelength_text = text.split("(")[1].split(" Å")[0]
            wavelength = float(wavelength_text)
            self.transition_spinbox.setValue(wavelength)
            self.transition = wavelength
        except (IndexError, ValueError) as e:
            print(f"Error parsing transition: {str(e)}")
    
    def on_transition_changed(self, value):
        """Handle direct transition wavelength changes."""
        self.transition = value
        
        # Set combo box to "Custom" if the value doesn't match a predefined transition
        found = False
        for i in range(self.transition_combo.count()):
            text = self.transition_combo.itemText(i)
            if "Custom" in text:
                continue
                
            try:
                wavelength_text = text.split("(")[1].split(" Å")[0]
                wavelength = float(wavelength_text)
                if abs(wavelength - value) < 0.01:  # Small tolerance for floating point comparison
                    self.transition_combo.setCurrentIndex(i)
                    found = True
                    break
            except (IndexError, ValueError):
                pass
        
        if not found:
            # Find and select the "Custom" item
            for i in range(self.transition_combo.count()):
                if self.transition_combo.itemText(i) == "Custom":
                    self.transition_combo.setCurrentIndex(i)
                    break
    
    def on_vmin_changed(self, value):
        """Handle velocity minimum changes."""
        self.vmin = value
        # Ensure vmin < vmax
        if value >= self.vmax:
            self.vmax_spinbox.setValue(value + 100)
    
    def on_vmax_changed(self, value):
        """Handle velocity maximum changes."""
        self.vmax = value
        # Ensure vmax > vmin
        if value <= self.vmin:
            self.vmin_spinbox.setValue(value - 100)
    
    def slice_spectrum(self):
        """Slice the spectrum around the selected transition."""
        if not self.controller.has_spectrum():
            QMessageBox.warning(self, "Error", "No spectrum loaded or no redshift applied.")
            return
        
        try:
            # Get parameters
            transition = self.transition
            vmin = self.vmin
            vmax = self.vmax
            use_vel = self.use_vel_checkbox.isChecked()
            linelist = self.linelist_combo.currentText()
            
            # Call slice_spectrum method
            success = self.controller.slice_spectrum(
                transition, vmin, vmax, 
                use_vel=use_vel, 
                linelist=linelist
            )
            
            if success:
                self.status_label.setText(f"Spectrum sliced around {transition:.2f} Å")
                self.update_plot()
                self.spectrum_sliced.emit(transition, vmin, vmax)
            else:
                QMessageBox.warning(self, "Error", "Failed to slice spectrum.")
        except Exception as e:
            print(f"Error slicing spectrum: {str(e)}")
            QMessageBox.warning(self, "Error", f"Failed to slice spectrum: {str(e)}")


    def set_transition(self, wavelength, name=None):
        """Set the transition wavelength and update UI (e.g., when loading from JSON)."""
        self.transition = wavelength
        self.transition_spinbox.setValue(wavelength)
        
        # Update status label with the proper name if provided
        if name:
            self.status_label.setText(f"Transition: {name} ({wavelength:.2f} Å)")
    
    def set_velocity_limits(self, vmin, vmax):
        """Set the velocity limits (e.g., when loading from JSON)."""
        self.vmin = vmin
        self.vmax = vmax
        self.vmin_spinbox.setValue(vmin)
        self.vmax_spinbox.setValue(vmax)
        
    def update_plot(self):
        """Update the spectrum plot."""
        if not self.controller.has_spectrum():
            return
        
        try:
            # Get sliced data from controller if available
            velocity, flux, error = self.controller.get_sliced_data()
            
            if velocity is None:
                # No sliced data yet, get full spectrum
                wave, flux, error = self.controller.get_spectrum_data()
                x_data = wave
                x_label = 'Wavelength (Å)'
                title = 'Full Spectrum'
            else:
                # Plot sliced data
                x_data = velocity
                x_label = 'Velocity (km/s)' if self.use_vel_checkbox.isChecked() else 'Wavelength (Å)'
                
                # Get transition info for title
                wavelength, name = self.controller.get_transition_info()
                if name:
                    title = f'Transition: {name} ({wavelength:.2f} Å)'
                else:
                    title = f'Transition: {self.transition:.2f} Å'
            
            # Clear the plot
            self.canvas.axes.clear()
            
            # Plot the data
            self.canvas.axes.step(x_data, flux, 'k-', where='mid', label='Flux')
            self.canvas.axes.step(x_data, error, 'r-', where='mid', alpha=0.5, label='Error')
            
            # Set labels and title
            self.canvas.axes.set_xlabel(x_label)
            self.canvas.axes.set_ylabel('Relative Flux')
            self.canvas.axes.set_title(title)
            #self.canvas.axes.legend()
            
            # Redraw
            self.canvas.fig.tight_layout()  # Add tight_layout to ensure proper spacing
            self.canvas.draw()
        except Exception as e:
            print(f"Error updating plot: {str(e)}")