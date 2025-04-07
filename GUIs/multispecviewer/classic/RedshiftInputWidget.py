import sys
from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QComboBox, QPushButton, 
                             QVBoxLayout, QHBoxLayout, QMessageBox, QApplication)
from PyQt5.QtCore import pyqtSignal

class RedshiftInputWidget(QWidget):
    """
    A standalone PyQt5 widget that provides UI elements for entering a redshift value
    and selecting a line list (LLS or DLA).
    
    Signals:
        submitted(float, str): Emitted when valid data is submitted, containing the redshift value and selected line list
        linelist_changed(str): Emitted when linelist selection changes
    """
    
    # Define custom signals
    submitted = pyqtSignal(float, str)
    linelist_changed = pyqtSignal(float, str)  # New signal for linelist changes
    
    def __init__(self, parent=None):
        super().__init__(parent)
        # Set default values
        self.default_redshift = 0.0
        self.default_linelist = "LLS"
        self.initUI()
        
    def initUI(self):
        # Create widgets
        redshift_label = QLabel("Redshift Guess:")
        self.redshift_input = QLineEdit()
        # Set default redshift value
        self.redshift_input.setText(f"{self.default_redshift:.6f}")
        
        linelist_label = QLabel("Line List:")
        self.linelist_combo = QComboBox()
        self.linelist_combo.addItems(["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong"])  # Updated linelist options
        
        # Set default line list
        default_index = self.linelist_combo.findText(self.default_linelist)
        if default_index >= 0:
            self.linelist_combo.setCurrentIndex(default_index)
        
        # Connect linelist combo box change event
        self.linelist_combo.currentIndexChanged.connect(self.on_linelist_changed)
        
        submit_button = QPushButton("Submit")
        submit_button.clicked.connect(self.validate_and_submit)
        
        # Create layouts
        redshift_layout = QHBoxLayout()
        redshift_layout.addWidget(redshift_label)
        redshift_layout.addWidget(self.redshift_input)
        
        linelist_layout = QHBoxLayout()
        linelist_layout.addWidget(linelist_label)
        linelist_layout.addWidget(self.linelist_combo)
        
        # Main layout
        main_layout = QVBoxLayout()
        main_layout.addLayout(redshift_layout)
        main_layout.addLayout(linelist_layout)
        main_layout.addWidget(submit_button)
        
        self.setLayout(main_layout)
        
        # Emit the initial values to ensure everything is synchronized at startup
        try:
            self.linelist_changed.emit(self.default_redshift, self.default_linelist)
        except Exception as e:
            print(f"Error emitting initial values: {str(e)}")
    
    def on_linelist_changed(self, index):
        """Handle linelist selection change"""
        try:
            # Get the current redshift value
            redshift_text = self.redshift_input.text().strip()
            linelist = self.linelist_combo.currentText()
            
            try:
                redshift = float(redshift_text)
                # Emit the linelist_changed signal with current redshift and new linelist
                self.linelist_changed.emit(redshift, linelist)
            except ValueError:
                # If the redshift is invalid, don't emit the signal
                # Optionally, could show a warning here, but that might interrupt the workflow
                pass
        except Exception as e:
            print(f"Error in linelist change handler: {str(e)}")
    
    def validate_and_submit(self):
        """Validate inputs and emit the submitted signal if valid"""
        redshift_text = self.redshift_input.text().strip()
        linelist = self.linelist_combo.currentText()
        
        # Validate redshift
        try:
            redshift = float(redshift_text)
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid redshift value (number).")
            return
            
        # Validate line list selection
        if not linelist:
            QMessageBox.warning(self, "Missing Selection", "Please select a line list.")
            return
        
        # All validations passed, emit the signal with the data
        self.submitted.emit(redshift, linelist)
        
    def get_values(self):
        """Get the current values from the inputs (can be used as an alternative to the signal)"""
        try:
            redshift = float(self.redshift_input.text().strip())
            linelist = self.linelist_combo.currentText()
            return redshift, linelist
        except ValueError:
            return None, self.linelist_combo.currentText()
    
    def set_values(self, redshift=None, linelist=None):
        """Set the values of the input fields"""
        if redshift is not None:
            self.redshift_input.setText(str(redshift))
        
        if linelist is not None and linelist in ["LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong"]:
            index = self.linelist_combo.findText(linelist)
            if index >= 0:
                self.linelist_combo.setCurrentIndex(index)
                
    # Add the set_redshift method to update the widget from external sources
    def set_redshift(self, redshift):
        """
        Sets the redshift value in the input field programmatically
        
        :param redshift: The redshift value to set
        """
        try:
            # Format the redshift to 6 decimal places and set in the text field
            formatted_redshift = f"{float(redshift):.6f}"
            self.redshift_input.setText(formatted_redshift)
            
            # Note: This only updates the display - it doesn't trigger the submitted signal
            # If you want to automatically trigger calculation, you'd need to add:
            # self.submitted.emit(float(formatted_redshift), self.linelist_combo.currentText())
        except Exception as e:
            print(f"Error setting redshift value: {str(e)}")
            # If there's an error, set a fallback value
            self.redshift_input.setText("0.000000")