import sys
from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit, QComboBox, QPushButton, 
                             QVBoxLayout, QHBoxLayout, QMessageBox, QApplication,
                             QGridLayout, QSpacerItem, QSizePolicy)
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QColor

# Import color utility from rb_utility
try:
    from rbcodes.utils import rb_utility as rt
except:
    from utils import rb_utility as rt

class RedshiftInputWidget(QWidget):
    """
    A standalone PyQt5 widget that provides UI elements for entering a redshift value
    and selecting a line list and color.
    
    Signals:
        submitted(float, str, str): Emitted when valid data is submitted
        linelist_changed(float, str, str): Emitted when linelist or color changes
        catalog_clicked(float, str, str): Emitted when catalog button is clicked
    """
    
    # Define custom signals
    submitted = pyqtSignal(float, str, str)
    linelist_changed = pyqtSignal(float, str, str)
    catalog_clicked = pyqtSignal(float, str, str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Set default values
        self.default_redshift = 0.0
        self.default_linelist = "LLS"
        self.default_color = "white"
        
        # Get color options from utility
        clr = rt.rb_set_color()
        self.color_options = list(clr.keys())[1:]  # Skip the first color (usually background)
        
        self.initUI()
        
    def initUI(self):
        # Create main layout
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(5, 5, 5, 5)
        main_layout.setSpacing(5)
        
        # Redshift input layout
        redshift_layout = QHBoxLayout()
        redshift_label = QLabel("Redshift:")
        redshift_label.setStyleSheet("color: white;")  

        self.redshift_input = QLineEdit()
        self.redshift_input.setText(f"{self.default_redshift:.6f}")
        self.redshift_input.setStyleSheet("""
            QLineEdit {
                background-color: #3C3C3C;
                color: white;
                border: 1px solid #555555;
                border-radius: 1px;
                padding: 1px;
                selection-background-color: #2A82DA;
                selection-color: white;
            }
        """)
        redshift_layout.addWidget(redshift_label)
        redshift_layout.addWidget(self.redshift_input)
        main_layout.addLayout(redshift_layout)
        
        # Line list layout
        linelist_layout = QHBoxLayout()
        linelist_label = QLabel("Line List:")
        linelist_label.setStyleSheet("color: white;")  
        self.linelist_combo = QComboBox()
        self.linelist_combo.addItems(["None", "LLS", "LLS Small", "DLA", "LBG", "Gal", "Eiger_Strong","AGN"])
        # AGN, DLA, Eiger_Strong, Gal, Gal_Abs, Gal_Em, Gal_long, HI_recomb, HI_recomb_light, LBG, LLS, LLS Small, atom
        self.linelist_combo.setCurrentText(self.default_linelist)

        self.linelist_combo.setStyleSheet("""
            QComboBox { 
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 1px;
                padding: 1px;
                font-size: 12px;
            }
            QComboBox::drop-down {
                background-color: #48484A;
                border-left: 1px solid #636366;
                width: 30px;
            }
            QComboBox QAbstractItemView {
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                selection-background-color: #0A84FF;
                selection-color: white;
            }                
        """)

        linelist_layout.addWidget(linelist_label)
        linelist_layout.addWidget(self.linelist_combo)
        main_layout.addLayout(linelist_layout)
        
        # Color layout
        color_layout = QHBoxLayout()
        color_label = QLabel("Color:")
        color_label.setStyleSheet("color: white;")  
        self.color_combo = QComboBox()
        self.color_combo.addItems(self.color_options)
        self.color_combo.setCurrentText(self.default_color)
        self.color_combo.setStyleSheet("""
            QComboBox { 
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 1px;
                padding: 1px;
                font-size: 14px;
            }
            QComboBox::drop-down {
                background-color: #48484A;
                border-left: 1px solid #636366;
                width: 30px;
            }
            QComboBox QAbstractItemView {
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                selection-background-color: #0A84FF;
                selection-color: white;
            }
        """)
        color_layout.addWidget(color_label)
        color_layout.addWidget(self.color_combo)
        main_layout.addLayout(color_layout)
        
        # Buttons layout
        buttons_layout = QHBoxLayout()
        
        # Submit button
        self.submit_button = QPushButton("Submit")
        self.submit_button.setStyleSheet("""
            QPushButton {
                background-color: #474747;  
                color: #F2F2F2;             
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050; 
            }
            QPushButton:pressed {
                background-color: #2A2A2A;  
            }
        """)
        self.submit_button.clicked.connect(self.validate_and_submit)
        buttons_layout.addWidget(self.submit_button)
        
        # Catalog button
        self.catalog_button = QPushButton("Catalog")
        self.catalog_button.setStyleSheet("""
            QPushButton {
                background-color: #474747;  
                color: #F2F2F2;             
                border: none;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #505050; 
            }
            QPushButton:pressed {
                background-color: #2A2A2A;  
            }
        """)

        self.catalog_button.clicked.connect(self.on_catalog_clicked)
        buttons_layout.addWidget(self.catalog_button)
        
        main_layout.addLayout(buttons_layout)
        
        # Connect signals for changes
        #This plots changes on the fly
        #self.linelist_combo.currentIndexChanged.connect(self.on_linelist_changed)
        #self.color_combo.currentIndexChanged.connect(self.on_color_changed)
        
    def on_catalog_clicked(self):
        """Handle catalog button click"""
        try:
            redshift_text = self.redshift_input.text().strip()
            redshift = float(redshift_text)
            linelist = self.linelist_combo.currentText()
            color = self.color_combo.currentText()
            
            self.catalog_clicked.emit(redshift, linelist, color)
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid redshift value (number).")
    
    def on_linelist_changed(self, index):
        """Handle linelist selection change"""
        try:
            redshift_text = self.redshift_input.text().strip()
            redshift = float(redshift_text)
            linelist = self.linelist_combo.currentText()
            color = self.color_combo.currentText()
            
            self.linelist_changed.emit(redshift, linelist, color)
        except ValueError:
            pass  # Silently handle invalid redshift during changes
    
    def on_color_changed(self, index):
        """Handle color selection change"""
        try:
            redshift_text = self.redshift_input.text().strip()
            redshift = float(redshift_text)
            linelist = self.linelist_combo.currentText()
            color = self.color_combo.currentText()
            
            self.linelist_changed.emit(redshift, linelist, color)
        except ValueError:
            pass  # Silently handle invalid redshift during changes
    
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
        if not linelist or linelist == "None":
            QMessageBox.warning(self, "Missing Selection", "Please select a line list.")
            return
        
        # All validations passed, emit the signal with the data
        self.submitted.emit(redshift, linelist, color)
    
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