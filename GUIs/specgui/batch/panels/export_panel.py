# rbcodes/GUIs/specgui/batch/panels/export_panel.py
import os
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QFileDialog, QGroupBox, QFormLayout,
                           QCheckBox, QLineEdit, QMessageBox, QComboBox)
from PyQt5.QtCore import pyqtSignal, Qt

class ExportPanel(QWidget):
    """Panel for exporting batch processing results."""
    
    # Signals
    export_completed = pyqtSignal(bool, str)  # success, path
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add instructions
        instructions = QLabel(
            "Export batch processing results as a CSV file or individual JSON files. "
            "You can also export the error log for troubleshooting."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)
        
        # CSV export section
        csv_group = QGroupBox("Export Summary CSV")
        csv_layout = QFormLayout()
        
        # CSV file path
        csv_path_layout = QHBoxLayout()
        self.csv_path = QLineEdit()
        self.csv_path.setPlaceholderText("Enter path for CSV file...")
        
        self.csv_browse_btn = QPushButton("Browse")
        self.csv_browse_btn.clicked.connect(self.browse_csv_path)
        
        csv_path_layout.addWidget(self.csv_path)
        csv_path_layout.addWidget(self.csv_browse_btn)
        
        csv_layout.addRow("CSV File:", csv_path_layout)
        
        # CSV export button
        self.export_csv_btn = QPushButton("Export CSV")
        self.export_csv_btn.clicked.connect(self.export_csv)
        
        csv_layout.addRow("", self.export_csv_btn)
        
        csv_group.setLayout(csv_layout)
        main_layout.addWidget(csv_group)
        
        # JSON export section
        json_group = QGroupBox("Export Individual JSON Files")
        json_layout = QFormLayout()
        
        # Output directory
        json_dir_layout = QHBoxLayout()
        self.json_dir = QLineEdit()
        self.json_dir.setPlaceholderText("Enter directory for JSON files...")
        
        self.json_browse_btn = QPushButton("Browse")
        self.json_browse_btn.clicked.connect(self.browse_json_dir)
        
        json_dir_layout.addWidget(self.json_dir)
        json_dir_layout.addWidget(self.json_browse_btn)
        
        json_layout.addRow("Output Directory:", json_dir_layout)
        
        # File name template
        self.json_template = QLineEdit("{transition}_{redshift:.3f}.json")
        json_layout.addRow("Filename Template:", self.json_template)
        
        # Export all option
        self.export_all_json = QCheckBox("Export All Items")
        self.export_all_json.setChecked(True)
        json_layout.addRow("", self.export_all_json)
        
        # JSON export button
        self.export_json_btn = QPushButton("Export JSON Files")
        self.export_json_btn.clicked.connect(self.export_json)
        
        json_layout.addRow("", self.export_json_btn)
        
        json_group.setLayout(json_layout)
        main_layout.addWidget(json_group)
        
        # Error log export section
        error_group = QGroupBox("Export Error Log")
        error_layout = QFormLayout()
        
        # Error log path
        error_path_layout = QHBoxLayout()
        self.error_path = QLineEdit()
        self.error_path.setPlaceholderText("Enter path for error log file...")
        
        self.error_browse_btn = QPushButton("Browse")
        self.error_browse_btn.clicked.connect(self.browse_error_path)
        
        error_path_layout.addWidget(self.error_path)
        error_path_layout.addWidget(self.error_browse_btn)
        
        error_layout.addRow("Error Log File:", error_path_layout)
        
        # Error log export button
        self.export_error_btn = QPushButton("Export Error Log")
        self.export_error_btn.clicked.connect(self.export_error_log)
        
        error_layout.addRow("", self.export_error_btn)
        
        error_group.setLayout(error_layout)
        main_layout.addWidget(error_group)
    
    def browse_csv_path(self):
        """Browse for CSV file path."""
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save CSV File", "", 
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )
        
        if filepath:
            if not filepath.lower().endswith('.csv'):
                filepath += '.csv'
            self.csv_path.setText(filepath)
    
    def browse_json_dir(self):
        """Browse for JSON output directory."""
        directory = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", "", 
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks
        )
        
        if directory:
            self.json_dir.setText(directory)
    
    def browse_error_path(self):
        """Browse for error log file path."""
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save Error Log", "", 
            "Text Files (*.txt);;All Files (*)",
            options=options
        )
        
        if filepath:
            if not filepath.lower().endswith('.txt'):
                filepath += '.txt'
            self.error_path.setText(filepath)
    
    def export_csv(self):
        """Export results to a CSV file."""
        filepath = self.csv_path.text()
        
        if not filepath:
            QMessageBox.warning(
                self, "No File Path", 
                "Please specify a file path for the CSV file."
            )
            return
        
        # Create directory if it doesn't exist
        directory = os.path.dirname(filepath)
        if directory and not os.path.isdir(directory):
            try:
                os.makedirs(directory)
            except Exception as e:
                QMessageBox.warning(
                    self, "Directory Error", 
                    f"Failed to create directory: {str(e)}"
                )
                return
        
        # Export the CSV
        success = self.controller.export_results_csv(filepath)
        
        if success:
            QMessageBox.information(
                self, "Export Successful", 
                f"Results exported to {filepath}"
            )
            self.export_completed.emit(True, filepath)
        else:
            QMessageBox.warning(
                self, "Export Failed", 
                "Failed to export results to CSV."
            )
            self.export_completed.emit(False, filepath)
    
    def export_json(self):
        """Export results to individual JSON files."""
        directory = self.json_dir.text()
        
        if not directory:
            QMessageBox.warning(
                self, "No Directory", 
                "Please specify an output directory for the JSON files."
            )
            return
        
        # Create directory if it doesn't exist
        if not os.path.isdir(directory):
            try:
                os.makedirs(directory)
            except Exception as e:
                QMessageBox.warning(
                    self, "Directory Error", 
                    f"Failed to create directory: {str(e)}"
                )
                return
        
        # Update controller settings
        template = self.json_template.text()
        export_all = self.export_all_json.isChecked()
        
        # This would typically call a method on the controller to export JSON files
        # For now, just show a message
        QMessageBox.information(
            self, "Not Fully Implemented", 
            f"JSON export functionality will be fully implemented in a future version.\n\n"
            f"Would export to: {directory}\n"
            f"Using template: {template}\n"
            f"Export all: {export_all}"
        )
        
        self.export_completed.emit(True, directory)
    
    def export_error_log(self):
        """Export the error log to a file."""
        filepath = self.error_path.text()
        
        if not filepath:
            QMessageBox.warning(
                self, "No File Path", 
                "Please specify a file path for the error log."
            )
            return
        
        # Create directory if it doesn't exist
        directory = os.path.dirname(filepath)
        if directory and not os.path.isdir(directory):
            try:
                os.makedirs(directory)
            except Exception as e:
                QMessageBox.warning(
                    self, "Directory Error", 
                    f"Failed to create directory: {str(e)}"
                )
                return
        
        # Export the error log
        success = self.controller.export_error_log(filepath)
        
        if success:
            QMessageBox.information(
                self, "Export Successful", 
                f"Error log exported to {filepath}"
            )
            self.export_completed.emit(True, filepath)
        else:
            QMessageBox.warning(
                self, "Export Failed", 
                "Failed to export error log."
            )
            self.export_completed.emit(False, filepath)