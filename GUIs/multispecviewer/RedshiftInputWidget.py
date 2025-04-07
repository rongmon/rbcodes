import sys
from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QComboBox, QPushButton, 
                             QVBoxLayout, QHBoxLayout, QMessageBox, QApplication,)
from PyQt5.QtCore import pyqtSignal

# Import color utility from rb_utility
from utils import rb_utility as rt

class RedshiftInputWidget(QWidget):
    """
    A standalone PyQt5 widget that provides UI elements for entering a redshift value,
    selecting a line list, choosing a color, and a catalog button.
    
    Signals:
        submitted(float, str, str): Emitted when valid data is submitted, containing the redshift value, 
                                    selected line list, and color
        linelist_changed(float, str, str): Emitted when linelist selection changes
        catalog_clicked(float, str, str): Emitted when catalog button is clicked
    """
    
    # Define custom signals
    submitted = pyqtSignal(float, str, str)
    linelist_changed = pyqtSignal(float, str, str)
    catalog_clicked = pyqtSignal(float, str, str)  # New signal for catalog button
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Get color list from rb_utility
        self.clr = rt.rb_set_color()
        self.color_list = list(self.clr.keys())[1:]  # Skip the first color (usually background)
        
        # Set default values
        self.default_redshift = 0.0
        self.default_linelist = "LLS"
        self.default_color = "white"
        
        self.initUI()
        
    def initUI(self):
        # Dark theme colors
        dark_background = '#353535'      # Main dark background
        darker_background = '#252525'    # Slightly darker background for inputs
        white_text = '#FFFFFF'           # White text
        highlight_color = '#252525'  #'#2A82DA'      # Highlight/accent color
        
        # Apply dark theme stylesheet
        dark_stylesheet = f"""
        QWidget {{
            background-color: {dark_background};
            color: {white_text};
        }}
        
        QLineEdit, QComboBox {{
            background-color: {darker_background};
            color: {white_text};
            border: 1px solid {highlight_color};
            padding: 5px;
            min-width: 200px;  
            min-height: 30px;  
        }}
        
        QPushButton {{
            background-color: {dark_background};
            color: {white_text};
            border: 1px solid {highlight_color};
            padding: 5px;
        }}
        
        QPushButton:hover {{
            background-color: {highlight_color};
        }}
        """
        
        # Create widgets
        redshift_label = QLabel("Redshift Guess:")
        self.redshift_input = QLineEdit()
        # Set default redshift value
        self.redshift_input.setText(f"{self.default_redshift:.6f}")
        
        linelist_label = QLabel("Line List:")
        self.linelist_combo = QComboBox()
        self.linelist_combo.addItems(["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong"])
        
        # Add color selection dropdown using the color list from rb_utility
        color_label = QLabel("Color:")
        self.color_combo = QComboBox()
        self.color_combo.addItems(self.color_list)
        
        # Set default line list and color
        default_index = self.linelist_combo.findText(self.default_linelist)
        if default_index >= 0:
            self.linelist_combo.setCurrentIndex(default_index)
        
        default_color_index = self.color_combo.findText(self.default_color)
        if default_color_index >= 0:
            self.color_combo.setCurrentIndex(default_color_index)
        
        # Connect linelist combo box change event
        self.linelist_combo.currentIndexChanged.connect(self.on_linelist_changed)
        
        # Create buttons
        submit_button = QPushButton("Submit")
        submit_button.clicked.connect(self.validate_and_submit)
        
        # Add catalog button
        catalog_button = QPushButton("Catalog")
        catalog_button.clicked.connect(self.on_catalog_clicked)
        
        # Create layouts
        redshift_layout = QHBoxLayout()
        redshift_layout.addWidget(redshift_label)
        redshift_layout.addWidget(self.redshift_input)
        
        linelist_layout = QHBoxLayout()
        linelist_layout.addWidget(linelist_label)
        linelist_layout.addWidget(self.linelist_combo)
        
        color_layout = QHBoxLayout()
        color_layout.addWidget(color_label)
        color_layout.addWidget(self.color_combo)
        
        button_layout = QHBoxLayout()
        button_layout.addWidget(submit_button)
        button_layout.addWidget(catalog_button)
        
        # Main layout
        main_layout = QVBoxLayout()
        main_layout.addLayout(redshift_layout)
        main_layout.addLayout(linelist_layout)
        main_layout.addLayout(color_layout)
        main_layout.addLayout(button_layout)
        
        # Apply dark theme stylesheet
        self.setStyleSheet(dark_stylesheet)
        
        self.setLayout(main_layout)
        
        # Emit the initial values to ensure everything is synchronized at startup
        try:
            self.linelist_changed.emit(self.default_redshift, self.default_linelist, self.default_color)
        except Exception as e:
            print(f"Error emitting initial values: {str(e)}")
    
    
    def on_color_changed(self, index):
        """Handle color selection change"""
        try:
            # Get the current redshift value and linelist
            redshift_text = self.redshift_input.text().strip()
            linelist = self.linelist_combo.currentText()
            color = self.color_combo.currentText()
            
            try:
                redshift = float(redshift_text)
                # Emit the linelist_changed signal with current redshift, linelist, and new color
                self.linelist_changed.emit(redshift, linelist, color)
            except ValueError:
                # If the redshift is invalid, don't emit the signal
                pass
        except Exception as e:
            print(f"Error in color change handler: {str(e)}")
    
    def on_linelist_changed(self, index):
        """Handle linelist selection change"""
        try:
            # Get the current redshift value and color
            redshift_text = self.redshift_input.text().strip()
            linelist = self.linelist_combo.currentText()
            color = self.color_combo.currentText()
            
            try:
                redshift = float(redshift_text)
                # Emit the linelist_changed signal with current redshift, new linelist, and color
                self.linelist_changed.emit(redshift, linelist, color)
            except ValueError:
                # If the redshift is invalid, don't emit the signal
                pass
        except Exception as e:
            print(f"Error in linelist change handler: {str(e)}")
    
    def on_catalog_clicked(self):
        """Handle catalog button click"""
        try:
            redshift_text = self.redshift_input.text().strip()
            linelist = self.linelist_combo.currentText()
            color = self.color_combo.currentText()
            
            try:
                redshift = float(redshift_text)
                # Emit the catalog_clicked signal
                self.catalog_clicked.emit(redshift, linelist, color)
            except ValueError:
                QMessageBox.warning(self, "Invalid Input", "Please enter a valid redshift value (number).")
        except Exception as e:
            print(f"Error in catalog button handler: {str(e)}")
    
    def validate_and_submit(self):
        """Validate inputs and emit the submitted signal if valid"""
        redshift_text = self.redshift_input.text().strip()
        linelist = self.linelist_combo.currentText()
        color = self.color_combo.currentText()
        
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
        self.submitted.emit(redshift, linelist, color)
        
    def get_values(self):
        """Get the current values from the inputs (can be used as an alternative to the signal)"""
        try:
            redshift = float(self.redshift_input.text().strip())
            linelist = self.linelist_combo.currentText()
            color = self.color_combo.currentText()
            return redshift, linelist, color
        except ValueError:
            return None, self.linelist_combo.currentText(), self.color_combo.currentText()
    
    def set_values(self, redshift=None, linelist=None, color=None):
        """Set the values of the input fields"""
        if redshift is not None:
            self.redshift_input.setText(str(redshift))
        
        if linelist is not None and linelist in ["LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong"]:
            index = self.linelist_combo.findText(linelist)
            if index >= 0:
                self.linelist_combo.setCurrentIndex(index)
                
        if color is not None:
            index = self.color_combo.findText(color)
            if index >= 0:
                self.color_combo.setCurrentIndex(index)
                
    def set_redshift(self, redshift):
        """
        Sets the redshift value in the input field programmatically
        
        :param redshift: The redshift value to set
        """
        try:
            # Format the redshift to 6 decimal places and set in the text field
            formatted_redshift = f"{float(redshift):.6f}"
            self.redshift_input.setText(formatted_redshift)
        except Exception as e:
            print(f"Error setting redshift value: {str(e)}")
            # If there's an error, set a fallback value
            self.redshift_input.setText("0.000000")