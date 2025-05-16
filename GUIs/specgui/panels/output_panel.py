import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QFormLayout, QGroupBox, QMessageBox,
                           QLineEdit, QFileDialog, QComboBox, QCheckBox,
                           QDialog, QDialogButtonBox)
from PyQt5.QtCore import pyqtSignal
class OutputPanel(QWidget):
    """Panel for saving analysis results and generating output."""
    
    # Signals
    analysis_saved = pyqtSignal(str)  # Emitted when analysis is saved (filepath)
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add instructions
        instructions = QLabel(
            "Save your analysis results as a JSON file that can be loaded later. "
            "You can also export plots of the spectrum and measurements."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)
        
        # Save JSON section
        save_group = QGroupBox("Save Analysis")
        save_layout = QFormLayout()
        
        # File path
        file_path_layout = QHBoxLayout()
        self.file_path = QLineEdit()
        # Create a good initial guess for the filename
        default_filename = "analysis_results.json"
        # If we have transition information, use it to create a more specific filename
        if hasattr(self.controller, 'spec') and self.controller.has_spectrum():
            if hasattr(self.controller.spec, 'trans'):
                trans_name = self.controller.spec.trans.replace(' ', '_')
                transition_wave = getattr(self.controller.spec, 'trans_wave', 0)
                zabs = getattr(self.controller.spec, 'zabs', 0)
                default_filename = f"{trans_name}_{transition_wave:.1f}_z{zabs:.3f}.json"
        
        self.file_path.setText(default_filename)
        self.file_path.setPlaceholderText("Enter a file path or use Browse...")
        
        self.browse_btn = QPushButton("Browse")
        self.browse_btn.clicked.connect(self.browse_save_path)
        
        file_path_layout.addWidget(self.file_path)
        file_path_layout.addWidget(self.browse_btn)
        
        save_layout.addRow("Save Path:", file_path_layout)
        
        # Format selector
        self.format_combo = QComboBox()
        self.format_combo.addItems(["JSON", "Pickle"])
        save_layout.addRow("Format:", self.format_combo)
        
        # Save button
        self.save_btn = QPushButton("Save Analysis")
        self.save_btn.clicked.connect(self.save_analysis)
        save_layout.addRow("", self.save_btn)
        
        save_group.setLayout(save_layout)
        main_layout.addWidget(save_group)
        
        # Export plots section
        export_group = QGroupBox("Export Plots")
        export_layout = QVBoxLayout()
        
        # Plot options
        self.export_continuum = QCheckBox("Export Continuum Fit Plot")
        self.export_continuum.setChecked(True)
        
        self.export_measurement = QCheckBox("Export Measurement Plot")
        self.export_measurement.setChecked(True)
        
        export_layout.addWidget(self.export_continuum)
        export_layout.addWidget(self.export_measurement)
        
        # Format selector for plots
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("Format:"))
        
        self.plot_format = QComboBox()
        self.plot_format.addItems(["PNG", "PDF", "SVG"])
        format_layout.addWidget(self.plot_format)
        
        export_layout.addLayout(format_layout)
        
        # Export button
        self.export_btn = QPushButton("Export Plots")
        self.export_btn.clicked.connect(self.export_plots)
        export_layout.addWidget(self.export_btn)
        
        export_group.setLayout(export_layout)
        main_layout.addWidget(export_group)
        
        # Summary display
        summary_group = QGroupBox("Analysis Summary")
        summary_layout = QVBoxLayout()
        
        self.summary_label = QLabel("No analysis to summarize yet.")
        self.summary_label.setWordWrap(True)
        summary_layout.addWidget(self.summary_label)
        
        summary_group.setLayout(summary_layout)
        main_layout.addWidget(summary_group)
        
        # Update summary when this panel is shown
        self.update_summary()
    
    def browse_save_path(self):
        """Open file dialog to select save path."""
        # Get format
        fmt = self.format_combo.currentText().lower()
        extension = ".json" if fmt == "json" else ".p"
        
        # Get current filename as default
        current_filename = self.file_path.text()
        if not current_filename:
            current_filename = "analysis_results" + extension
        
        # Ensure it has the right extension
        if not current_filename.lower().endswith(extension):
            current_filename = os.path.splitext(current_filename)[0] + extension
        
        file_filter = f"{fmt.upper()} Files (*{extension});;All Files (*)"
        
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Analysis", current_filename, file_filter
        )
        
        if filename:
            # Add extension if not present
            if not filename.lower().endswith(extension):
                filename += extension
                
            self.file_path.setText(filename)    
    def save_analysis(self):
        """Save the analysis to a file."""
        filepath = self.file_path.text()
        
        if not filepath:
            QMessageBox.warning(self, "Input Required", "Please specify a file path.")
            return
        
        # Get format
        fmt = self.format_combo.currentText().lower()
        
        try:
            success = self.controller.save_analysis(filepath, format=fmt)
            
            if success:
                QMessageBox.information(
                    self, "Success", f"Analysis saved to {filepath}"
                )
                self.update_summary()
                self.analysis_saved.emit(filepath)
            else:
                QMessageBox.warning(
                    self, "Error", "Failed to save analysis."
                )
        except Exception as e:
            print(f"Error saving analysis: {str(e)}")
            QMessageBox.warning(
                self, "Error", f"Failed to save analysis: {str(e)}"
            )
    
    
    def update_summary(self):
        """Update the analysis summary display."""
        if not self.controller.has_spectrum():
            self.summary_label.setText("No analysis to summarize yet.")
            return
        
        try:
            # Get information from controller
            zabs, transition, name, ew, col_dens = self.controller.get_analysis_summary()
            
            if transition is None:
                self.summary_label.setText("Analysis incomplete. Some steps have not been performed.")
                return
            
            # Format the summary
            summary = (
                f"<b>Absorber Redshift:</b> {zabs:.6f}<br>"
                f"<b>Transition:</b> {name} ({transition:.2f} Å)<br>"
            )
            
            if ew is not None:
                summary += f"<b>Equivalent Width:</b> {ew['value']:.3f} ± {ew['error']:.3f} Å<br>"
            
            if col_dens is not None:
                summary += f"<b>Column Density (log):</b> {col_dens['value']:.2f} ± {col_dens['error']:.2f}<br>"
            
            self.summary_label.setText(summary)
        except Exception as e:
            print(f"Error updating summary: {str(e)}")
            self.summary_label.setText(f"Error updating summary: {str(e)}")


    def export_plots(self):
        """Export plots of the analysis using the same path as the JSON file."""
        # Get format
        fmt = self.plot_format.currentText().lower()
        
        # Check if any plots are selected
        if not self.export_continuum.isChecked() and not self.export_measurement.isChecked():
            QMessageBox.warning(
                self, "No Plots Selected", 
                "Please select at least one plot to export."
            )
            return
        
        # Check if we have a save path
        filepath = self.file_path.text()
        
        if not filepath:
            # If no path set, prompt user to select one
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getSaveFileName(
                self, "Save Results", "", 
                f"JSON Files (*.json);;All Files (*)", 
                options=options
            )
            
            if not filename:
                return  # User cancelled
            
            # Add .json extension if not present
            if not filename.lower().endswith('.json'):
                filename += '.json'
                
            self.file_path.setText(filename)
            filepath = filename
        
        # Extract directory path from filepath
        directory = os.path.dirname(filepath)
        if not directory:
            directory = '.'  # Use current directory if not specified
        
        # Get base name without extension
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        
        try:
            # Export continuum plot if selected
            if self.export_continuum.isChecked():
                continuum_path = os.path.join(directory, f"{base_name}_continuum.{fmt}")
                success = self.controller.export_continuum_plot(continuum_path)
                if not success:
                    QMessageBox.warning(
                        self, "Error", "Failed to export continuum plot."
                    )
                else:
                    print(f"Continuum plot saved to: {continuum_path}")
            
            # Export measurement plot if selected
            if self.export_measurement.isChecked():
                ew_path = os.path.join(directory, f"{base_name}_ew.{fmt}")
                success = self.controller.export_measurement_plot(ew_path)
                if not success:
                    QMessageBox.warning(
                        self, "Error", "Failed to export measurement plot."
                    )
                else:
                    print(f"EW plot saved to: {ew_path}")
            
            QMessageBox.information(
                self, "Success", f"Plots saved to {directory}"
            )
        except Exception as e:
            print(f"Error exporting plots: {str(e)}")
            QMessageBox.warning(
                self, "Error", f"Failed to export plots: {str(e)}"
            )    