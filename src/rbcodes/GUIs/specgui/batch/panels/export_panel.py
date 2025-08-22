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
            "You can also export figures and error logs for review and troubleshooting."
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
        
        # Export options for JSON
        self.export_all_json = QCheckBox("Export All Items")
        self.export_all_json.setChecked(True)
        self.export_complete_only = QCheckBox("Export Only Completed Items")
        self.export_complete_only.setChecked(True)
        
        json_layout.addRow("", self.export_all_json)
        json_layout.addRow("", self.export_complete_only)
        
        # JSON export button
        self.export_json_btn = QPushButton("Export JSON Files")
        self.export_json_btn.clicked.connect(self.export_json)
        
        json_layout.addRow("", self.export_json_btn)
        
        json_group.setLayout(json_layout)
        main_layout.addWidget(json_group)
        
        # Figure export section (placeholder for now)
        figure_group = QGroupBox("Export Figures")
        figure_layout = QFormLayout()
        
        # Figure directory
        figure_dir_layout = QHBoxLayout()
        self.figure_dir = QLineEdit()
        self.figure_dir.setPlaceholderText("Enter directory for figure files...")
        
        self.figure_browse_btn = QPushButton("Browse")
        self.figure_browse_btn.clicked.connect(self.browse_figure_dir)
        
        figure_dir_layout.addWidget(self.figure_dir)
        figure_dir_layout.addWidget(self.figure_browse_btn)
        
        figure_layout.addRow("Output Directory:", figure_dir_layout)
        
        # Figure format
        self.figure_format = QComboBox()
        self.figure_format.addItems(["PDF", "PNG"])
        figure_layout.addRow("Format:", self.figure_format)
        
        # Figure type (enable/disable based on format)
        self.figure_type = QComboBox()
        self.figure_type.addItems(["Individual Files", "Multi-page PDF"])
        
        # Connect format change to update available types
        self.figure_format.currentTextChanged.connect(self._update_figure_type_options)
        
        figure_layout.addRow("Output Type:", self.figure_type)
        
        # Figure export button
        self.export_figure_btn = QPushButton("Export Figures")
        self.export_figure_btn.clicked.connect(self.export_figures)
        # Enable the button now that functionality is implemented
        self.export_figure_btn.setEnabled(True)
        
        figure_layout.addRow("", self.export_figure_btn)
        
        figure_group.setLayout(figure_layout)
        main_layout.addWidget(figure_group)
        
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
        
        # Set default directories
        self._set_default_directories()
        
        # Connect to controller signals to update defaults when batch changes
        self.controller.table_changed.connect(self._update_defaults_on_change)
    
    def _update_defaults_on_change(self):
        """Update default directories when batch configuration changes - simplified."""
        # Always update all defaults to be consistent
        self._set_default_directories()
    
    def _is_csv_path_default(self):
        """Check if CSV path looks like our auto-generated name."""
        current_path = self.csv_path.text()
        if not current_path:
            return True
        
        filename = os.path.basename(current_path)
        # Check if it matches our naming pattern
        return filename.startswith('batch_results_') and filename.endswith('.csv')
    
    def _are_directories_default(self):
        """Check if directories haven't been manually changed by user."""
        # Simple heuristic: if all three directories are the same, they're probably defaults
        json_dir = self.json_dir.text()
        figure_dir = self.figure_dir.text()
        
        if not json_dir or not figure_dir:
            return True
        
        return json_dir == figure_dir  # Same directory suggests they're defaults
    
    def _update_csv_filename(self):
        """Update just the CSV filename, keeping the consistent directory."""
        # Get the consistent default directory
        default_dir = self._get_default_directory()
        csv_filename = self._generate_smart_csv_filename()
        self.csv_path.setText(os.path.join(default_dir, csv_filename))
    
    def _update_default_directories(self):
        """Update the directory fields only."""
        default_dir = self._get_default_directory()
        
        # Update directories
        self.json_dir.setText(default_dir)
        self.figure_dir.setText(default_dir)
        
        # Update error log path
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y-%m-%d")
        error_filename = f"batch_errors_{timestamp}.txt"
        self.error_path.setText(os.path.join(default_dir, error_filename))
    
    def _set_default_directories(self):
        """Set default output directories and filenames - all use same directory."""
        try:
            # Get the consistent default directory
            default_dir = self._get_default_directory()
            
            # Set smart CSV filename
            csv_filename = self._generate_smart_csv_filename()
            csv_path = os.path.join(default_dir, csv_filename)
            
            # Set smart error log filename
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y-%m-%d")
            error_filename = f"batch_errors_{timestamp}.txt"
            error_path = os.path.join(default_dir, error_filename)
            
            # Update ALL UI fields with the SAME directory
            self.csv_path.setText(csv_path)
            self.json_dir.setText(default_dir)
            self.figure_dir.setText(default_dir)
            self.error_path.setText(error_path)
                        
        except Exception as e:
            print(f"Could not set default directories: {e}")

    
    def _get_default_directory(self):
        """Get the best default directory - simplified for consistency."""
        # Priority 1: Configuration loading directory (where batch config was loaded from)
        if hasattr(self.controller, 'config_load_directory') and self.controller.config_load_directory:
            if os.path.exists(self.controller.config_load_directory):
                return self.controller.config_load_directory
        
        # Priority 2: Current working directory (fallback)
        return os.getcwd()
    
    def _generate_smart_csv_filename(self):
        """Generate smart CSV filename based on batch content."""
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y-%m-%d")
        
        items = self.controller.master_table.get_all_items()
        if not items:
            return f"batch_results_{timestamp}.csv"
        
        # Analyze batch content
        transitions = set()
        redshifts = set()
        
        for item in items:
            transitions.add(item.template.transition_name)
            # Round redshift to 3 decimal places for grouping
            redshifts.add(round(item.template.redshift, 3))
        
        # Build filename based on content
        filename_parts = ["batch_results"]
        
        # Add transition info
        if len(transitions) == 1:
            # Single transition type
            transition_name = list(transitions)[0].replace(' ', '').replace('IV', '4')  # Clean up name
            filename_parts.append(transition_name)
        elif len(transitions) <= 3:
            # Few transitions - list them
            sorted_transitions = sorted(list(transitions))
            transition_str = "_".join([t.replace(' ', '').replace('IV', '4') for t in sorted_transitions])
            filename_parts.append(transition_str)
        else:
            # Many transitions
            filename_parts.append("mixed")
        
        # Add redshift info
        if len(redshifts) == 1:
            # Single redshift
            z_value = list(redshifts)[0]
            filename_parts.append(f"z{z_value:.3f}")
        elif len(redshifts) <= 3:
            # Few redshifts - add range info
            min_z = min(redshifts)
            max_z = max(redshifts)
            if min_z == max_z:
                filename_parts.append(f"z{min_z:.3f}")
            else:
                filename_parts.append(f"z{min_z:.3f}-{max_z:.3f}")
        # For many redshifts, don't add redshift info to keep filename reasonable
        
        # Add timestamp and extension
        filename_parts.append(timestamp)
        
        return "_".join(filename_parts) + ".csv"
    
    def _update_figure_type_options(self, format_text):
        """Update available figure type options based on selected format."""
        current_selection = self.figure_type.currentText()
        
        self.figure_type.clear()
        self.figure_type.addItem("Individual Files")
        
        # Multi-page is only available for PDF
        if format_text.lower() == 'pdf':
            self.figure_type.addItem("Multi-page PDF")
            
            # Restore selection if it was multi-page
            if current_selection == "Multi-page PDF":
                self.figure_type.setCurrentText("Multi-page PDF")
    
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
    
    def browse_figure_dir(self):
        """Browse for figure output directory."""
        directory = QFileDialog.getExistingDirectory(
            self, "Select Output Directory", "", 
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks
        )
        
        if directory:
            self.figure_dir.setText(directory)
    
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
        """Export results to a CSV file with enhanced fields."""
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
        
        # Export the CSV with enhanced fields
        success = self._export_enhanced_csv(filepath)
        
        if success:
            # Send status message instead of popup
            self.controller.status_updated.emit(f"CSV results exported successfully to {os.path.basename(filepath)}")
            self.export_completed.emit(True, filepath)
        else:
            # Send error message to status bar
            self.controller.status_updated.emit("Failed to export CSV results")
            self.export_completed.emit(False, filepath)
    
    def update_defaults_after_processing(self):
        """Update defaults after batch processing completes (use processing output directory)."""
        # This can be called by the main window when processing completes
        self._set_default_directories()
    
    def _export_enhanced_csv(self, filepath):
        """Export enhanced CSV with velocity centroid and dispersion."""
        try:
            import csv
            
            items = self.controller.master_table.get_all_items()
            if not items:
                self.controller.status_updated.emit("No results to export.")
                return False
            
            with open(filepath, 'w', newline='') as csvfile:
                fieldnames = ['filename', 'redshift', 'transition', 'transition_name', 
                             'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax',
                             'W', 'W_e', 'logN', 'logN_e', 'SNR', 
                             'vel_centroid', 'vel_disp', 'status']
                
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for item in items:
                    row = {
                        'filename': os.path.basename(item.template.filename),
                        'redshift': item.template.redshift,
                        'transition': item.template.transition,
                        'transition_name': item.template.transition_name,
                        'slice_vmin': item.template.slice_vmin,
                        'slice_vmax': item.template.slice_vmax,
                        'ew_vmin': item.template.ew_vmin,
                        'ew_vmax': item.template.ew_vmax,
                        'W': item.results.W,
                        'W_e': item.results.W_e,
                        'logN': item.results.logN,
                        'logN_e': item.results.logN_e,
                        'SNR': item.results.SNR,
                        'vel_centroid': item.results.vel_centroid,
                        'vel_disp': item.results.vel_disp,
                        'status': item.analysis.processing_status
                    }
                    writer.writerow(row)
            
            self.controller.status_updated.emit(f"Enhanced CSV exported to {filepath}")
            return True
        except Exception as e:
            self.controller.status_updated.emit(f"Error exporting enhanced CSV: {str(e)}")
            return False
    
    def export_json(self):
        """Export individual JSON files for selected items."""
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
        
        # Get items to export
        items = self.controller.master_table.get_all_items()
        if not items:
            QMessageBox.warning(
                self, "No Items", 
                "No batch items found to export."
            )
            return
        
        # Filter items based on options
        export_items = []
        
        if self.export_complete_only.isChecked():
            # Only export completed items
            for item in items:
                if item.analysis.processing_status == 'complete':
                    export_items.append(item)
        else:
            # Export all items if requested
            if self.export_all_json.isChecked():
                export_items = items
        
        if not export_items:
            QMessageBox.warning(
                self, "No Items to Export", 
                "No items match the export criteria."
            )
            return
        
        # Export individual JSON files
        success_count, error_count = self._export_individual_jsons(export_items, directory)
        
        # Show results
        if error_count == 0:
            # Send status message instead of popup
            self.controller.status_updated.emit(f"Successfully exported {success_count} JSON files to {os.path.basename(directory)}")
            self.export_completed.emit(True, directory)
        else:
            # Send warning message to status bar
            self.controller.status_updated.emit(f"JSON export completed: {success_count} successful, {error_count} failed")
            self.export_completed.emit(False, directory)
    
    def _export_individual_jsons(self, items, directory):
        """Export individual JSON files for the given items."""
        success_count = 0
        error_count = 0
        
        for row_index, item in enumerate(items):
            try:
                # Get rb_spec object from master table using row index
                spec_object = self.controller.master_table.get_rb_spec_object(row_index)
                
                if spec_object is None:
                    print(f"Cannot export {item.template.transition_name}: No rb_spec object available")
                    error_count += 1
                    continue
                
                # Simple logic: if JSON input, keep same name; otherwise create new name
                if item.template.filename.lower().endswith('.json'):
                    # For JSON files, use the exact same filename
                    output_filename = os.path.basename(item.template.filename)
                else:
                    # For other files, create the standard naming convention
                    basename = os.path.splitext(os.path.basename(item.template.filename))[0]
                    transition_id = f"{item.template.transition_name}_{item.template.transition:.0f}"
                    output_filename = f"{basename}_{transition_id}_z{item.template.redshift:.3f}.json"
                output_path = os.path.join(directory, output_filename)
                
                # Check for filename conflicts and resolve
                output_path = self._resolve_filename_conflict(output_path)
                
                # Save the rb_spec object
                spec_object.save_slice(output_path, file_format='json', verbose=False)
                
                print(f"Exported: {output_filename}")
                success_count += 1
                
            except Exception as e:
                print(f"Error exporting {item.template.transition_name}: {str(e)}")
                error_count += 1
        
        self.controller.status_updated.emit(
            f"JSON export completed: {success_count} successful, {error_count} failed"
        )
        
        return success_count, error_count    

    def _resolve_filename_conflict(self, filepath):
        """Resolve filename conflicts by appending numbers."""
        if not os.path.exists(filepath):
            return filepath
        
        base, ext = os.path.splitext(filepath)
        counter = 1
        
        while os.path.exists(f"{base}_{counter}{ext}"):
            counter += 1
        
        return f"{base}_{counter}{ext}"
    
    def export_figures(self):
        """Export figures using the batch figure generator."""
        directory = self.figure_dir.text()
        
        if not directory:
            QMessageBox.warning(
                self, "No Directory", 
                "Please specify an output directory for the figures."
            )
            return
        
        # Get items to export
        items = self.controller.master_table.get_all_items()
        if not items:
            QMessageBox.warning(
                self, "No Items", 
                "No batch items found to export."
            )
            return
        
        # Filter to only completed items
        completed_items = [item for item in items if item.analysis.processing_status == 'complete']
        
        if not completed_items:
            QMessageBox.warning(
                self, "No Completed Items", 
                "No completed items found to export figures for."
            )
            return
        
        # Get export options
        file_format = self.figure_format.currentText().lower()
        figure_type_text = self.figure_type.currentText()
        
        if figure_type_text == "Individual Files":
            figure_type = 'individual'
        elif figure_type_text == "Multi-page PDF":
            if file_format != 'pdf':
                QMessageBox.warning(
                    self, "Format Error", 
                    "Multi-page option is only available for PDF format."
                )
                return
            figure_type = 'multipage'
        else:
            figure_type = 'individual'
        
        # Import and use the figure generator
        try:
            from rbcodes.GUIs.specgui.batch.panels.batch_figure_generator import export_batch_figures
            
            success, message = export_batch_figures(
                    completed_items, 
                    directory, 
                    file_format, 
                    figure_type,
                    master_table=self.controller.master_table  # Add this
                )
            
            if success:
                # Send status message instead of popup
                self.controller.status_updated.emit(f"Figures exported successfully to {os.path.basename(directory)}")
                self.export_completed.emit(True, directory)
            else:
                # Send error message to status bar
                self.controller.status_updated.emit(f"Figure export failed: {message}")
                self.export_completed.emit(False, directory)
                
        except ImportError as e:
            QMessageBox.warning(
                self, "Import Error", 
                f"Could not import figure generator module:\n{str(e)}"
            )
        except Exception as e:
            QMessageBox.warning(
                self, "Export Error", 
                f"Error during figure export:\n{str(e)}"
            )
            self.export_completed.emit(False, directory)
    
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
            # Send status message instead of popup
            self.controller.status_updated.emit(f"Error log exported to {os.path.basename(filepath)}")
            self.export_completed.emit(True, filepath)
        else:
            # Send error message to status bar  
            self.controller.status_updated.emit("Failed to export error log")
            self.export_completed.emit(False, filepath)