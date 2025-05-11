from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QListWidget, QListWidgetItem,
                             QLabel, QPushButton, QHBoxLayout, QWidget)
from PyQt5.QtCore import Qt, pyqtSignal

class LineSelectionDialog(QDialog):
    """
    A non-modal dialog for selecting spectral lines from a list.
    Uses signals to communicate back with the parent window instead of blocking execution.
    """
    # Define a signal for line selection
    lineSelected = pyqtSignal(str, float, float)  # transition_name, rest_wavelength, redshift
    
    def __init__(self, parent=None, observed_wavelength=None, line_list=None):
        super().__init__(parent)
        self.setWindowTitle("Line Identification")
        self.setWindowFlags(Qt.Dialog | Qt.WindowStaysOnTopHint)
        self.setMinimumWidth(400)
        
        self.observed_wavelength = observed_wavelength
        self.line_list = line_list
        
        # Setup UI
        self.initUI()
        
        # Position the dialog near the cursor but not directly under it
        # and make sure it's visible
        if parent:
            self.setGeometry(
                parent.x() + parent.width() // 4,
                parent.y() + parent.height() // 4,
                400, 900
            )
            
        # Apply dark theme styling
        self.apply_theme()
    
    def initUI(self):
        """Initialize the user interface"""
        # Main layout
        layout = QVBoxLayout(self)
        
        # Display the observed wavelength
        if self.observed_wavelength is not None:
            info_label = QLabel(f"Observed Wavelength: {self.observed_wavelength:.2f} Å")
            layout.addWidget(info_label)
        
        instruction_label = QLabel("Select a spectral line to compute redshift:")
        layout.addWidget(instruction_label)
        
        # Create list widget
        self.list_widget = QListWidget()
        self.list_widget.setAlternatingRowColors(False) # Change this from True to False
        
        # Populate the list if we have data
        if self.line_list and self.observed_wavelength is not None:
            for ix in range(len(self.line_list)):
                rest_wavelength = self.line_list[ix]['wrest']
                transition_name = self.line_list[ix]['ion']
                # Calculate implied redshift if this line is selected
                implied_redshift = (self.observed_wavelength / rest_wavelength) - 1.0
                list_item_text = f"{transition_name} ({rest_wavelength:.2f} Å) → z = {implied_redshift:.4f}"
                
                item = QListWidgetItem(list_item_text)
                # Store the data with the item
                item.setData(Qt.UserRole, {
                    'transition_name': transition_name,
                    'rest_wavelength': rest_wavelength,
                    'implied_redshift': implied_redshift
                })
                self.list_widget.addItem(item)
        
        layout.addWidget(self.list_widget)
        
        # Buttons layout
        button_layout = QHBoxLayout()
        
        self.select_button = QPushButton("Select")
        self.select_button.clicked.connect(self.on_select_clicked)
        
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.close)
        
        button_layout.addWidget(self.select_button)
        button_layout.addWidget(self.cancel_button)
        
        layout.addLayout(button_layout)
        
        # Connect double-click on list item to selection
        self.list_widget.itemDoubleClicked.connect(self.on_item_double_clicked)
    
    def apply_theme(self):
        """Apply the dark theme styling consistent with other widgets"""
        # Dialog background
        self.setStyleSheet("""
            QDialog {
                background-color: #353535;
                color: #F2F2F7;
            }
            
            QLabel {
                color: #F2F2F7;
                font-size: 14px;
            }
            
            QListWidget {
                background-color: #252525;
                alternate-background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 6px;
                padding: 4px;
                font-size: 14px;
            }
            
            QListWidget::item {
                padding: 6px;
            }
            
            QListWidget::item:selected {
                background-color: #0A84FF;
                color: white;
            }
            
            QPushButton {
                background-color: #3A3A3C;
                color: #F2F2F7;
                border: 1px solid #636366;
                border-radius: 6px;
                padding: 6px 12px;
                font-size: 14px;
            }
            
            QPushButton:hover {
                background-color: #48484A;
            }
            
            /* Style the Select button specifically */
            QPushButton#select_button {
                background-color: #0A84FF;
                color: white;
                border: none;
            }
            
            QPushButton#select_button:hover {
                background-color: #409CFF;
            }
        """)
        
        # Set object names to allow specific styling
        self.select_button.setObjectName("select_button")
    
    def on_select_clicked(self):
        """Handle select button click"""
        current_item = self.list_widget.currentItem()
        if current_item:
            self.process_selection(current_item)
    
    def on_item_double_clicked(self, item):
        """Handle double-click on an item"""
        self.process_selection(item)
    
    def process_selection(self, item):
        """Process the selected item"""
        if not item:
            return
            
        # Get the stored data
        data = item.data(Qt.UserRole)
        if data:
            # Emit the signal with the selection data
            self.lineSelected.emit(
                data['transition_name'],
                data['rest_wavelength'],
                data['implied_redshift']
            )
        
        # Close the dialog
        self.close()