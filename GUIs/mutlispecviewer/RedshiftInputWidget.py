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
    """
    
    # Define a custom signal that emits the redshift (float) and line list (str) when submitted
    submitted = pyqtSignal(float, str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        
    def initUI(self):
        # Create widgets
        redshift_label = QLabel("Redshift Guess:")
        self.redshift_input = QLineEdit()
        
        linelist_label = QLabel("Line List:")
        self.linelist_combo = QComboBox()
        self.linelist_combo.addItems(["None", "LLS","LLS_Small" ,"DLA"])  # Empty string as default "Select" option
        
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
        
        if linelist is not None and linelist in ["LLS", "DLA"]:
            index = self.linelist_combo.findText(linelist)
            if index >= 0:
                self.linelist_combo.setCurrentIndex(index)


# Example of how to use this widget in an existing application
if __name__ == "__main__":
    # This is just for demonstration - not needed when embedding in your application
    app = QApplication(sys.argv)
    
    # Create a main window
    main_window = QWidget()
    main_window.setWindowTitle("Redshift Input Example")
    main_layout = QVBoxLayout()
    
    # Create and add our widget
    redshift_widget = RedshiftInputWidget()
    
    # Connect to the signal
    def on_submitted(redshift, linelist):
        print(f"Submitted: Redshift = {redshift}, Line List = {linelist}")
        QMessageBox.information(
            main_window, 
            "Submission Received", 
            f"Redshift: {redshift}\nLine List: {linelist}"
        )
    
    redshift_widget.submitted.connect(on_submitted)
    
    # Add to layout
    main_layout.addWidget(redshift_widget)
    main_window.setLayout(main_layout)
    
    # Show the window
    main_window.show()
    sys.exit(app.exec_())