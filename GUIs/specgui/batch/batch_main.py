# rbcodes/GUIs/specgui/batch/batch_main.py
import sys
import os
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTabWidget, QVBoxLayout, 
                           QWidget, QStatusBar, QMenuBar, QAction, QMessageBox,
                           QProgressBar, QLabel, QHBoxLayout)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QKeySequence, QIcon

from rbcodes.GUIs.specgui.batch.batch_controller import BatchController
from rbcodes.GUIs.specgui.batch.panels.configuration_panel import ConfigurationPanel
from rbcodes.GUIs.specgui.batch.panels.processing_panel import ProcessingPanel
from rbcodes.GUIs.specgui.batch.panels.review_panel import ReviewPanel
from rbcodes.GUIs.specgui.batch.panels.export_panel import ExportPanel


class BatchSpecGUI(QMainWindow):
    """Main window for batch processing of absorption line spectra."""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("rb_spec Batch Processing")
        self.setMinimumSize(1000, 700)
        
        # Application state
        self.has_unsaved_changes = False
        self.current_config_file = None
        self.batch_results = []
        
        # Create the batch controller
        self.controller = BatchController()
        
        # Initialize UI
        self.init_ui()
        self.create_menus()
        self.create_status_bar()
        
        # Set initial tab states
        self.update_tab_states()
        
    def init_ui(self):
        """Initialize the user interface."""
        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Create tabs with custom styling
        self.tabs = QTabWidget()
        
        # Set custom stylesheet for better tab visibility
        self.tabs.setStyleSheet("""
            QTabWidget::pane {
                border: 1px solid #c0c0c0;
                background-color: white;
            }
            QTabWidget::tab-bar {
                alignment: left;
            }
            QTabBar::tab {
                background-color: #f0f0f0;
                border: 1px solid #c0c0c0;
                border-bottom-color: #c0c0c0;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
                min-width: 150px;
                padding: 8px 12px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background-color: white;
                border-bottom-color: white;
                color: #000000;
                font-weight: bold;
            }
            QTabBar::tab:hover {
                background-color: #e0e0e0;
            }
            QTabBar::tab:!selected {
                margin-top: 2px;
            }
        """)
        
        main_layout.addWidget(self.tabs)
        
        # Add panels to tabs
        self.config_panel = ConfigurationPanel(self.controller)
        self.tabs.addTab(self.config_panel, "1. Configuration")
        
        self.processing_panel = ProcessingPanel(self.controller)
        self.tabs.addTab(self.processing_panel, "2. Processing")
        
        self.review_panel = ReviewPanel(self.controller)
        self.tabs.addTab(self.review_panel, "3. Review Results")
        
        self.export_panel = ExportPanel(self.controller)
        self.tabs.addTab(self.export_panel, "4. Export")
        
        # Connect signals between panels and controller
        self.connect_signals()
    
    def create_menus(self):
        """Create the menu bar."""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('&File')
        
        # New configuration
        new_action = QAction('&New Configuration', self)
        new_action.setShortcut(QKeySequence.New)
        new_action.setStatusTip('Create a new batch configuration')
        new_action.triggered.connect(self.new_configuration)
        file_menu.addAction(new_action)
        
        # Open configuration
        open_action = QAction('&Open Configuration...', self)
        open_action.setShortcut(QKeySequence.Open)
        open_action.setStatusTip('Open an existing batch configuration')
        open_action.triggered.connect(self.open_configuration)
        file_menu.addAction(open_action)
        
        # Save configuration
        save_action = QAction('&Save Configuration', self)
        save_action.setShortcut(QKeySequence.Save)
        save_action.setStatusTip('Save the current batch configuration')
        save_action.triggered.connect(self.save_configuration)
        file_menu.addAction(save_action)
        
        # Save configuration as
        save_as_action = QAction('Save Configuration &As...', self)
        save_as_action.setShortcut(QKeySequence.SaveAs)
        save_as_action.setStatusTip('Save the batch configuration with a new name')
        save_as_action.triggered.connect(self.save_configuration_as)
        file_menu.addAction(save_as_action)
        
        file_menu.addSeparator()
        
        # Create CSV template
        template_action = QAction('Create CSV &Template...', self)
        template_action.setStatusTip('Create a CSV template for batch import')
        template_action.triggered.connect(self.create_csv_template)
        file_menu.addAction(template_action)
        
        file_menu.addSeparator()
        
        # Exit
        exit_action = QAction('E&xit', self)
        exit_action.setShortcut(QKeySequence.Quit)
        exit_action.setStatusTip('Exit the application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Edit menu
        edit_menu = menubar.addMenu('&Edit')
        
        # Clear all items
        clear_action = QAction('&Clear All Items', self)
        clear_action.setStatusTip('Clear all batch items')
        clear_action.triggered.connect(self.clear_all_items)
        edit_menu.addAction(clear_action)
        
        # Validate configuration
        validate_action = QAction('&Validate Configuration', self)
        validate_action.setStatusTip('Validate the current batch configuration')
        validate_action.triggered.connect(self.validate_configuration)
        edit_menu.addAction(validate_action)
        
        # Processing menu
        processing_menu = menubar.addMenu('&Processing')
        
        # Run batch
        run_action = QAction('&Run Batch Process', self)
        run_action.setShortcut('Ctrl+R')
        run_action.setStatusTip('Start batch processing')
        run_action.triggered.connect(self.run_batch_processing)
        processing_menu.addAction(run_action)
        
        # Stop processing
        stop_action = QAction('&Stop Processing', self)
        stop_action.setShortcut('Ctrl+T')
        stop_action.setStatusTip('Stop batch processing')
        stop_action.triggered.connect(self.stop_batch_processing)
        stop_action.setEnabled(False)
        processing_menu.addAction(stop_action)
        
        self.stop_action = stop_action  # Keep reference for enabling/disabling
        
        # Help menu
        help_menu = menubar.addMenu('&Help')
        
        # About
        about_action = QAction('&About', self)
        about_action.setStatusTip('About rb_spec Batch Processing')
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
        
        # User guide
        guide_action = QAction('&User Guide', self)
        guide_action.setStatusTip('Open the user guide')
        guide_action.triggered.connect(self.show_user_guide)
        help_menu.addAction(guide_action)
    
    def create_status_bar(self):
        """Create an enhanced status bar."""
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        
        # Main status label
        self.status_label = QLabel("Ready. Configure batch items to begin.")
        self.statusBar.addWidget(self.status_label)
        
        # Progress bar (initially hidden)
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximumWidth(200)
        self.statusBar.addPermanentWidget(self.progress_bar)
        
        # Items count label 
        self.items_count_label = QLabel("Items: 0")
        self.statusBar.addPermanentWidget(self.items_count_label)
        
        # Errors count label
        self.errors_count_label = QLabel("")
        self.errors_count_label.setStyleSheet("QLabel { color: red; }")
        self.statusBar.addPermanentWidget(self.errors_count_label)
    
    def connect_signals(self):
        """Connect signals between panels and controller."""
        # Configuration panel signals - updated for master table architecture
        self.config_panel.configuration_changed.connect(self.on_configuration_changed)

        #Connect config changes to review panel refresh
        self.config_panel.configuration_changed.connect(self.review_panel.refresh_results_table)
    
        
        # Controller signals
        self.controller.status_updated.connect(self.update_status)
        self.controller.table_changed.connect(self.on_configuration_changed)  # Updated signal name
        self.controller.batch_progress.connect(self.update_progress)
        self.controller.batch_item_completed.connect(self.on_batch_item_completed)
        self.controller.batch_item_failed.connect(self.on_batch_item_failed)
        
        # Processing panel signals
        self.processing_panel.batch_started.connect(self.on_batch_started)
        self.processing_panel.batch_completed.connect(self.on_batch_completed)
        
        # Review panel signals
        self.review_panel.items_selected.connect(self.processing_panel.set_selected_items)
        
        # Export panel signals
        self.export_panel.export_completed.connect(self.on_export_completed)
        
        # Tab change signal
        self.tabs.currentChanged.connect(self.on_tab_changed)
    
    def update_tab_states(self):
        """Update which tabs are enabled based on current state."""
        has_items = self.controller.get_item_count() > 0
        has_results = len(self.batch_results) > 0
        
        # Configuration tab is always enabled
        self.tabs.setTabEnabled(0, True)
        
        # Processing tab enabled if we have items
        self.tabs.setTabEnabled(1, has_items)
        
        # Review tab enabled if we have results
        self.tabs.setTabEnabled(2, has_results)
        
        # Export tab enabled if we have results
        self.tabs.setTabEnabled(3, has_results)
        
        # Update tab tooltips
        if not has_items:
            self.tabs.setTabToolTip(1, "Add batch items in Configuration tab first")
        else:
            self.tabs.setTabToolTip(1, "Configure and run batch processing")
            
        if not has_results:
            self.tabs.setTabToolTip(2, "Run batch processing first")
            self.tabs.setTabToolTip(3, "Run batch processing first")
        else:
            self.tabs.setTabToolTip(2, "Review and edit batch results")
            self.tabs.setTabToolTip(3, "Export results to files")
    
    def on_configuration_changed(self):
        """Handle configuration changes - updated for master table architecture."""
        self.has_unsaved_changes = True
        self.update_window_title()
        self.update_tab_states()
        
        # Update items count using new controller methods
        item_count = self.controller.get_item_count()
        self.items_count_label.setText(f"Items: {item_count}")
        
        # Clear error count when configuration changes
        self.errors_count_label.setText("")
        
        # Update status based on item count
        if item_count == 0:
            self.update_status("No batch items configured. Add items to begin.")
        else:
            valid_count = self.controller.get_valid_item_count()
            if valid_count == item_count:
                self.update_status(f"Configuration ready with {item_count} valid items.")
            else:
                invalid_count = item_count - valid_count
                self.update_status(f"Configuration has {invalid_count} invalid items. Check configuration.")
    
    def update_status(self, message):
        """Update the status bar with a message."""
        # Convert message to string if it's not already
        if isinstance(message, list):
            message = ", ".join(str(item) for item in message)
        elif not isinstance(message, str):
            message = str(message)
        
        self.status_label.setText(message)
    
    def update_progress(self, current, total):
        """Update the progress bar."""
        if total > 0:
            self.progress_bar.setVisible(True)
            self.progress_bar.setMaximum(total)
            self.progress_bar.setValue(current)
        else:
            self.progress_bar.setVisible(False)
    
    def on_batch_item_completed(self, index, result):
        """Handle batch item completion."""
        # Could add per-item status updates here
        pass
    
    def on_batch_item_failed(self, index, error):
        """Handle batch item failure."""
        # Update error count
        current_errors = len(self.controller.error_log)
        self.errors_count_label.setText(f"Errors: {current_errors}")
    
    def on_batch_started(self):
        """Handle batch processing started signal."""
        self.update_status("Batch processing started...")
        # Disable configuration tab during processing
        self.tabs.setTabEnabled(0, False)
        # Enable stop action
        if hasattr(self, 'stop_action'):
            self.stop_action.setEnabled(True)
    
    def on_batch_completed(self, success, results):
        """Handle batch processing completed signal."""
        self.batch_results = results
        
        if success:
            self.update_status(f"Batch processing completed. {len(results)} items processed.")
            # Update review panel with results
            self.review_panel.set_results(results)
            # Switch to review tab
            self.tabs.setCurrentIndex(2)
        else:
            self.update_status("Batch processing failed or was cancelled.")
        
        # Re-enable configuration tab
        self.tabs.setTabEnabled(0, True)
        # Update tab states
        self.update_tab_states()
        # Hide progress bar
        self.progress_bar.setVisible(False)
        # Disable stop action
        if hasattr(self, 'stop_action'):
            self.stop_action.setEnabled(False)
    
    def on_export_completed(self, success, path):
        """Handle export completed signal."""
        if success:
            self.update_status(f"Export completed successfully to {path}")
        else:
            self.update_status("Export failed.")
    
    def on_tab_changed(self, index):
        """Handle tab changes to provide guidance."""
        tab_names = ["Configuration", "Processing", "Review Results", "Export"]
        if index < len(tab_names):
            tab_name = tab_names[index]
            
            # Provide contextual guidance
            if index == 0:  # Configuration
                item_count = self.controller.get_item_count()
                if item_count == 0:
                    self.update_status("Add batch items using the buttons below the table.")
                else:
                    self.update_status(f"Configuration loaded with {item_count} items. Ready for processing.")
            elif index == 1:  # Processing
                item_count = self.controller.get_item_count()
                if item_count == 0:
                    self.update_status("Configure batch items first before processing.")
                else:
                    self.update_status("Configure processing settings and run batch processing.")
            elif index == 2:  # Review
                if len(self.batch_results) == 0:
                    self.update_status("No results to review. Run batch processing first.")
                else:
                    self.update_status(f"Reviewing {len(self.batch_results)} batch results. Select items to edit or reprocess.")
            elif index == 3:  # Export
                if len(self.batch_results) == 0:
                    self.update_status("No results to export. Run batch processing first.")
                else:
                    self.update_status("Export batch results to CSV or individual JSON files.")
    
    def update_window_title(self):
        """Update the window title based on current state."""
        title = "rb_spec Batch Processing"
        
        if self.current_config_file:
            title += f" - {os.path.basename(self.current_config_file)}"
        
        if self.has_unsaved_changes:
            title += " *"
        
        self.setWindowTitle(title)
    
    # Menu action implementations
    def new_configuration(self):
        """Create a new batch configuration."""
        if self.has_unsaved_changes:
            reply = QMessageBox.question(
                self, 'Unsaved Changes',
                'You have unsaved changes. Do you want to save before creating a new configuration?',
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            )
            
            if reply == QMessageBox.Yes:
                if not self.save_configuration():
                    return  # Save was cancelled or failed
            elif reply == QMessageBox.Cancel:
                return
        
        # Clear current configuration
        self.config_panel.clear_batch_items()
        self.batch_results = []
        self.current_config_file = None
        self.has_unsaved_changes = False
        self.update_window_title()
        self.update_tab_states()
        
        # Switch to configuration tab
        self.tabs.setCurrentIndex(0)
        self.update_status("New configuration created. Add batch items to begin.")
    
    def open_configuration(self):
        """Open an existing batch configuration."""
        if self.has_unsaved_changes:
            reply = QMessageBox.question(
                self, 'Unsaved Changes',
                'You have unsaved changes. Do you want to save before opening a new configuration?',
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            )
            
            if reply == QMessageBox.Yes:
                if not self.save_configuration():
                    return  # Save was cancelled or failed
            elif reply == QMessageBox.Cancel:
                return
        
        # Use the config panel's load functionality
        self.config_panel.load_configuration()
        
        # If successful, update state
        if len(self.controller.batch_items) > 0:
            self.has_unsaved_changes = False
            self.update_window_title()
            self.update_tab_states()
    
    def save_configuration(self):
        """Save the current batch configuration."""
        if self.current_config_file:
            success = self.controller.save_batch_configuration(self.current_config_file)
            if success:
                self.has_unsaved_changes = False
                self.update_window_title()
                return True
            return False
        else:
            return self.save_configuration_as()
    
    def save_configuration_as(self):
        """Save the batch configuration with a new name."""
        from PyQt5.QtWidgets import QFileDialog
        
        options = QFileDialog.Options()
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save Batch Configuration", "", 
            "JSON Files (*.json);;All Files (*)",
            options=options
        )
        
        if not filepath:
            return False
        
        # Add .json extension if not present
        if not filepath.lower().endswith('.json'):
            filepath += '.json'
        
        success = self.controller.save_batch_configuration(filepath)
        
        if success:
            self.current_config_file = filepath
            self.has_unsaved_changes = False
            self.update_window_title()
            QMessageBox.information(
                self, "Save Successful", 
                f"Batch configuration saved to {filepath}"
            )
            return True
        else:
            QMessageBox.warning(
                self, "Save Error",
                "Failed to save batch configuration."
            )
            return False
    
    def create_csv_template(self):
        """Create a CSV template for batch import."""
        if hasattr(self.config_panel, 'create_csv_template'):
            self.config_panel.create_csv_template()
        else:
            QMessageBox.information(
                self, "CSV Template",
                "CSV template functionality will be available in the next update."
            )
    
    def clear_all_items(self):
        """Clear all batch items."""
        item_count = self.controller.get_item_count()
        if item_count > 0:
            reply = QMessageBox.question(
                self, 'Clear All Items',
                'Are you sure you want to clear all batch items?',
                QMessageBox.Yes | QMessageBox.No
            )
            
            if reply == QMessageBox.Yes:
                self.controller.clear_all_items()
                self.batch_results = []
                self.has_unsaved_changes = True
                self.update_window_title()
                self.update_tab_states()
                self.update_status("All batch items cleared.")
    
    def validate_configuration(self):
        """Validate the current batch configuration."""
        item_count = self.controller.get_item_count()
        if item_count == 0:
            QMessageBox.information(
                self, "Validation Result",
                "No batch items to validate. Add items first."
            )
            return
        
        is_valid, message = self.controller.validate_all_items()
        
        if is_valid:
            QMessageBox.information(
                self, "Validation Result",
                message
            )
        else:
            QMessageBox.warning(
                self, "Validation Errors",
                message
            )
    
    def run_batch_processing(self):
        """Start batch processing."""
        item_count = self.controller.get_item_count()
        if item_count == 0:
            QMessageBox.warning(
                self, "No Items",
                "No batch items configured. Add items in the Configuration tab first."
            )
            return
        
        # Switch to processing tab and start processing
        self.tabs.setCurrentIndex(1)
        
        # Use QTimer to allow UI to update before starting processing
        QTimer.singleShot(100, self.processing_panel.run_batch)
    
    def stop_batch_processing(self):
        """Stop batch processing."""
        self.processing_panel.cancel_processing()
    
    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self, "About rb_spec Batch Processing",
            "<h3>rb_spec Batch Processing</h3>"
            "<p>A tool for batch analysis of absorption line spectra.</p>"
            "<p>This application allows you to:</p>"
            "<ul>"
            "<li>Configure multiple absorption systems or transitions for analysis</li>"
            "<li>Automatically fit continua and measure equivalent widths</li>"
            "<li>Review and edit results interactively</li>"
            "<li>Export results to CSV or individual JSON files</li>"
            "</ul>"
            "<p>Version 1.0<br>"
            "Part of the rbcodes package</p>"
        )
    
    def show_user_guide(self):
        """Show user guide."""
        guide_text = """
        <h3>rb_spec Batch Processing User Guide</h3>
        
        <h4>Workflow:</h4>
        <ol>
        <li><b>Configuration:</b> Add batch items specifying files, redshifts, transitions, and velocity ranges</li>
        <li><b>Processing:</b> Configure processing settings and run batch analysis</li>
        <li><b>Review:</b> Examine results, edit problematic fits interactively</li>
        <li><b>Export:</b> Save results to CSV or individual JSON files</li>
        </ol>
        
        <h4>Batch Modes:</h4>
        <ul>
        <li><b>Multiple Systems:</b> Same file, same transition, different redshifts</li>
        <li><b>Multiple Transitions:</b> Same file, same redshift, different transitions</li>
        <li><b>Multiple Files:</b> Different files, each with their own parameters</li>
        </ul>
        
        <h4>Velocity Ranges:</h4>
        <ul>
        <li><b>Slice Range:</b> Wide range for extracting spectrum (e.g., ±1500 km/s)</li>
        <li><b>EW Range:</b> Narrow range for measuring equivalent width (e.g., ±200 km/s)</li>
        </ul>
        
        <h4>Keyboard Shortcuts:</h4>
        <ul>
        <li>Ctrl+N: New configuration</li>
        <li>Ctrl+O: Open configuration</li>
        <li>Ctrl+S: Save configuration</li>
        <li>Ctrl+R: Run batch processing</li>
        <li>Ctrl+T: Stop processing</li>
        </ul>
        """
        
        QMessageBox.about(self, "User Guide", guide_text)
    
    def closeEvent(self, event):
        """Handle window close event."""
        if self.has_unsaved_changes:
            reply = QMessageBox.question(
                self, 'Unsaved Changes',
                'You have unsaved changes. Do you want to save before exiting?',
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            )
            
            if reply == QMessageBox.Yes:
                if not self.save_configuration():
                    event.ignore()
                    return
            elif reply == QMessageBox.Cancel:
                event.ignore()
                return
        
        # Stop any ongoing processing
        if hasattr(self.processing_panel, 'cancel_processing'):
            self.processing_panel.cancel_processing()
        
        event.accept()


def launch_batch_gui():
    """Launch the batch processing GUI."""
    app = QApplication(sys.argv)
    
    # Set application properties
    app.setApplicationName("rb_spec Batch Processing")
    app.setApplicationVersion("1.0")
    app.setOrganizationName("rbcodes")
    
    window = BatchSpecGUI()
    window.show()
    
    return app.exec_()


def main():
    """Main entry point for the batch processing application."""
    import argparse
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Launch rb_spec Batch Processing GUI')
    parser.add_argument('config', nargs='?', help='Batch configuration file to load on startup')
    parser.add_argument('--version', action='version', version='rb_spec Batch Processing 1.0')
    
    args = parser.parse_args()
    
    # Create application
    app = QApplication(sys.argv)
    app.setApplicationName("rb_spec Batch Processing")
    app.setApplicationVersion("1.0")
    app.setOrganizationName("rbcodes")
    
    # Create main window
    window = BatchSpecGUI()
    window.show()
    
    # Load configuration file if specified
    if args.config and os.path.exists(args.config):
        def load_config():
            success = window.controller.load_batch_configuration(args.config)
            if success:
                window.config_panel.update_ui_with_config(window.controller.batch_items)
                window.current_config_file = args.config
                window.has_unsaved_changes = False
                window.update_window_title()
                window.update_tab_states()
                window.update_status(f"Loaded configuration from {args.config}")
            else:
                window.update_status(f"Failed to load configuration from {args.config}")
        
        # Schedule config loading after GUI is fully initialized
        QTimer.singleShot(100, load_config)
    
    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())