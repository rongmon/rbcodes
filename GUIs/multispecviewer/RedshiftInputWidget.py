import sys
from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QComboBox, QPushButton, 
                             QVBoxLayout, QHBoxLayout, QMessageBox, QApplication,
                             QGridLayout, QSpacerItem, QSizePolicy)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor

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
    catalog_clicked = pyqtSignal(float, str, str)
    
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
        highlight_color = '#2A82DA'      # Highlight/accent color
        
        # Set up size policy
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.setMinimumWidth(300)
        self.setMaximumHeight(150)
        
        # Apply stylesheet with improved button and input field styling
        self.setStyleSheet(f"""
            QWidget {{
                background-color: {dark_background};
                color: {white_text};
            }}
            
            QLabel {{
                color: {white_text};
                font-size: 12px;
                min-height: 20px;
            }}
            
            QLineEdit {{
                background-color: {darker_background};
                color: {white_text};
                border: 1px solid {highlight_color};
                border-radius: 3px;
                padding: 5px;
                min-height: 25px;
                font-size: 12px;
            }}
            
            QComboBox {{
                background-color: {darker_background};
                color: {white_text};
                border: 1px solid {highlight_color};
                border-radius: 3px;
                padding: 5px;
                min-height: 25px;
                font-size: 12px;
            }}
            
            QComboBox::drop-down {{
                background-color: {dark_background};
                width: 20px;
                border-left: 1px solid {highlight_color};
            }}
            
            QComboBox QAbstractItemView {{
                background-color: {darker_background};
                color: {white_text};
                selection-background-color: {highlight_color};
            }}
            
            QPushButton {{
                background-color: {dark_background};
                color: {white_text};
                border: 1px solid {highlight_color};
                border-radius: 3px;
                padding: 5px;
                min-height: 30px;
                min-width: 80px;
                font-size: 12px;
            }}
            
            QPushButton:hover {{
                background-color: {highlight_color};
                color: {white_text};
            }}
            
            QPushButton:pressed {{
                background-color: #1C5A99;
            }}
        """)
        
        # Use a grid layout for better alignment of labels and inputs
        main_layout = QGridLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setVerticalSpacing(8)
        main_layout.setHorizontalSpacing(10)
        
        # Create labels with bold font
        redshift_label = QLabel("Redshift:")
        redshift_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        linelist_label = QLabel("Line List:")
        linelist_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        color_label = QLabel("Color:")
        color_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        # Create widgets with appropriate sizing
        self.redshift_input = QLineEdit()
        self.redshift_input.setText(f"{self.default_redshift:.6f}")
        self.redshift_input.setMinimumWidth(100)
        
        self.linelist_combo = QComboBox()
        self.linelist_combo.addItems(["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong"])
        self.linelist_combo.setMinimumWidth(120)
        
        self.color_combo = QComboBox()
        self.color_combo.addItems(self.color_list)
        self.color_combo.setMinimumWidth(120)
        
        # Set default selections
        default_index = self.linelist_combo.findText(self.default_linelist)
        if default_index >= 0:
            self.linelist_combo.setCurrentIndex(default_index)
        
        default_color_index = self.color_combo.findText(self.default_color)
        if default_color_index >= 0:
            self.color_combo.setCurrentIndex(default_color_index)
        
        # Connect linelist combo box change event
        self.linelist_combo.currentIndexChanged.connect(self.on_linelist_changed)
        
        # Create buttons with appropriate sizing
        submit_button = QPushButton("Submit")
        submit_button.clicked.connect(self.validate_and_submit)
        submit_button.setMinimumWidth(80)
        
        catalog_button = QPushButton("Catalog")
        catalog_button.clicked.connect(self.on_catalog_clicked)
        catalog_button.setMinimumWidth(80)
        
        # Add widgets to the grid layout
        main_layout.addWidget(redshift_label, 0, 0)
        main_layout.addWidget(self.redshift_input, 0, 1)
        
        main_layout.addWidget(linelist_label, 1, 0)
        main_layout.addWidget(self.linelist_combo, 1, 1)
        
        main_layout.addWidget(color_label, 2, 0)
        main_layout.addWidget(self.color_combo, 2, 1)
        
        # Create a horizontal layout for buttons
        button_layout = QHBoxLayout()
        button_layout.addWidget(submit_button)
        button_layout.addWidget(catalog_button)
        button_layout.addStretch(1)  # Add stretch to align buttons left
        
        # Add button layout to the grid
        main_layout.addLayout(button_layout, 3, 0, 1, 2)
        
        # Set stretch factors to ensure fields expand properly
        main_layout.setColumnStretch(1, 1)
        
        # Set column widths the grid layout
        main_layout.setColumnMinimumWidth(0, 80)  # Width for labels
        main_layout.setColumnMinimumWidth(1, 180) # Width for input fields
        
        # Emit the initial values to ensure everything is synchronized at startup
        try:
            self.linelist_changed.emit(self.default_redshift, self.default_linelist, self.default_color)
        except Exception as e:
            print(f"Error emitting initial values: {str(e)}")
    
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

# For testing the widget in isolation
if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Apply dark theme to entire application
    app.setStyle("Fusion")
    dark_palette = QPalette()
    dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.WindowText, Qt.white)
    dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
    dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
    dark_palette.setColor(QPalette.ToolTipText, Qt.white)
    dark_palette.setColor(QPalette.Text, Qt.white)
    dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
    dark_palette.setColor(QPalette.ButtonText, Qt.white)
    dark_palette.setColor(QPalette.BrightText, Qt.red)
    dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(dark_palette)
    
    # Create widget and show it
    widget = RedshiftInputWidget()
    widget.show()
    
    sys.exit(app.exec_())