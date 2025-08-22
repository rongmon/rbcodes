import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QDoubleSpinBox, QFormLayout, QGroupBox,
                           QMessageBox, QComboBox)
from PyQt5.QtCore import pyqtSignal

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots."""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)

class RedshiftPanel(QWidget):
    """Panel for setting and applying the absorber redshift."""
    
    # Signals
    redshift_applied = pyqtSignal(float)  # Emitted when redshift is applied
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.redshift = 0.0  # Initial redshift
        self.init_ui()
        
        # Connect to controller signals
        self.controller.spectrum_changed.connect(self.update_plot)
        self.controller.spectrum_changed.connect(self.reset)


    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Redshift controls
        redshift_group = QGroupBox("Redshift Settings")
        redshift_layout = QFormLayout()
        
        # Redshift input
        redshift_input_layout = QHBoxLayout()
        self.redshift_spinbox = QDoubleSpinBox()
        self.redshift_spinbox.setDecimals(6)  # 6 decimal places for precision
        self.redshift_spinbox.setRange(0.0, 10.0)  # Reasonable range for redshifts
        self.redshift_spinbox.setSingleStep(0.01)
        self.redshift_spinbox.setValue(self.redshift)
        self.redshift_spinbox.valueChanged.connect(self.on_redshift_changed)
        
        # Common redshift values (optional)
        self.common_redshifts = QComboBox()
        self.common_redshifts.addItem("Custom")
        self.common_redshifts.addItems(["0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0"])
        self.common_redshifts.currentTextChanged.connect(self.on_common_redshift_selected)
        
        # Apply button
        self.apply_btn = QPushButton("Apply Redshift")
        self.apply_btn.clicked.connect(self.apply_redshift)
        self.apply_btn.setToolTip("Apply the current redshift value to the spectrum")

        
        redshift_input_layout.addWidget(self.redshift_spinbox)
        redshift_input_layout.addWidget(self.common_redshifts)
        redshift_input_layout.addWidget(self.apply_btn)
        
        redshift_layout.addRow("Absorber Redshift (z):", redshift_input_layout)
        
        # Add Reset button
        reset_layout = QHBoxLayout()
        self.reset_btn = QPushButton("Reset Redshift")
        self.reset_btn.setToolTip("Clear current redshift and reset analysis. Use this to analyze the spectrum with a different redshift.")
        self.reset_btn.clicked.connect(self.reset_redshift)
        

        # Set a pale red color for the reset button
        self.reset_btn.setStyleSheet("""
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

        # Add the button to your layout
        redshift_button_layout = QHBoxLayout()
        redshift_button_layout.addWidget(self.apply_btn)
        redshift_button_layout.addWidget(self.reset_btn)

        reset_layout.addWidget(self.reset_btn)
        
        redshift_layout.addRow("", reset_layout)
        
        self.common_redshifts.setToolTip("Select from common redshift values or choose 'Custom'")

        # Status label
        self.status_label = QLabel("No redshift applied")
        redshift_layout.addRow("Status:", self.status_label)
        
        redshift_group.setLayout(redshift_layout)
        main_layout.addWidget(redshift_group)
        
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
        self.status_label.setText("No redshift applied")
        # Only reset if we don't have redshift info from a JSON file
        if not hasattr(self.controller.spec, 'zabs'):
            self.redshift_spinbox.setValue(0.0)  # Reset to default redshift
        self.update_plot()

    
    def on_redshift_changed(self, value):
        """Handle redshift spinbox value changes."""
        self.redshift = value
        
        # Update the common redshifts dropdown if it's a custom value
        found = False
        for i in range(self.common_redshifts.count()):
            if self.common_redshifts.itemText(i) == str(value):
                self.common_redshifts.setCurrentIndex(i)
                found = True
                break
        
        if not found:
            self.common_redshifts.setCurrentIndex(0)  # Set to "Custom"
    
    def reset_redshift(self):
        """Reset the redshift and clear downstream analysis."""
        # Ask for confirmation
        reply = QMessageBox.question(self, 'Confirm Reset', 
                                   'Reset redshift? This will clear all downstream analysis steps.',
                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        
        if reply == QMessageBox.No:
            return
        
        # Reset the spinbox to default value
        self.redshift_spinbox.setValue(0.0)
        self.redshift = 0.0
        self.status_label.setText("Redshift reset")
        
        # Signal to the main application that redshift was reset
        self.redshift_applied.emit(-99999)  # Use a special value to indicate reset   
    

    def on_common_redshift_selected(self, text):
        """Handle selection from common redshifts dropdown."""
        if text == "Custom":
            return  # Keep the current custom value
        
        try:
            z = float(text)
            self.redshift_spinbox.setValue(z)
        except ValueError:
            pass


    def set_redshift(self, zabs):
        """Set the redshift value (e.g., when loading from JSON)."""
        self.redshift = zabs
        self.redshift_spinbox.setValue(zabs)
        
    def apply_redshift(self):
        """Apply the redshift to the spectrum."""
        if not self.controller.has_spectrum():
            QMessageBox.warning(self, "Error", "No spectrum loaded.")
            return
        
        try:
            # Apply redshift using rb_spec
            success = self.controller.apply_redshift(self.redshift)
            
            if success:
                self.status_label.setText(f"Redshift z = {self.redshift} applied")
                self.update_plot()
                self.redshift_applied.emit(self.redshift)
            else:
                QMessageBox.warning(self, "Error", "Failed to apply redshift.")
        except Exception as e:
            print(f"Error applying redshift: {str(e)}")
            QMessageBox.warning(self, "Error", f"Failed to apply redshift: {str(e)}")
    
    def update_plot(self):
        """Update the spectrum plot."""
        if not self.controller.has_spectrum():
            return
        
        try:
            # Get data from controller
            wave, flux, error = self.controller.get_spectrum_data()
            
            # Clear the plot
            self.canvas.axes.clear()
            
            # Plot the data
            self.canvas.axes.step(wave, flux, 'k-', where='mid', label='Flux')
            self.canvas.axes.step(wave, error, 'r-', where='mid', alpha=0.5, label='Error')
            
            # Set labels and legend
            self.canvas.axes.set_xlabel('Wavelength (Å)')
            self.canvas.axes.set_ylabel('Flux')
            #self.canvas.axes.legend()
            
            # Redraw
            self.canvas.draw()
        except Exception as e:
            print(f"Error updating plot: {str(e)}")