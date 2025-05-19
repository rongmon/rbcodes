# rbcodes/GUIs/specgui/batch/batch_main.py
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget, QStatusBar
from PyQt5.QtCore import Qt

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
        
        # Create the batch controller
        self.controller = BatchController()
        
        # Initialize UI
        self.init_ui()
        
    def init_ui(self):
        """Initialize the user interface."""
        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Create tabs
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)
        
        # Add panels to tabs
        self.config_panel = ConfigurationPanel(self.controller)
        self.tabs.addTab(self.config_panel, "Configuration")
        
        self.processing_panel = ProcessingPanel(self.controller)
        self.tabs.addTab(self.processing_panel, "Processing")
        
        self.review_panel = ReviewPanel(self.controller)
        self.tabs.addTab(self.review_panel, "Review Results")
        
        self.export_panel = ExportPanel(self.controller)
        self.tabs.addTab(self.export_panel, "Export")
        
        # Connect signals between panels and controller
        self.connect_signals()
        
        # Create status bar
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.showMessage("Ready. Configure batch items to begin.")
    
    def connect_signals(self):
        """Connect signals between panels and controller."""
        # Configuration panel signals
        self.config_panel.configuration_changed.connect(self.controller.update_batch_items)
        self.config_panel.configuration_changed.connect(self.update_status)
        
        # Controller signals
        self.controller.status_updated.connect(self.update_status)
        
        # Processing panel signals
        self.processing_panel.batch_started.connect(self.on_batch_started)
        self.processing_panel.batch_completed.connect(self.on_batch_completed)
        
        # Review panel signals
        self.review_panel.items_selected.connect(self.processing_panel.set_selected_items)
        
        # Export panel signals
        self.export_panel.export_completed.connect(self.on_export_completed)
    
    # batch_main.py - update_status method
    def update_status(self, message):
        """Update the status bar with a message."""
        # Convert message to string if it's not already
        if isinstance(message, list):
            message = ", ".join(str(item) for item in message)
        elif not isinstance(message, str):
            message = str(message)
        
        self.statusBar.showMessage(message)
        
    def on_batch_started(self):
        """Handle batch processing started signal."""
        self.update_status("Batch processing started...")
        # Enable/disable tabs as needed
        self.tabs.setTabEnabled(0, False)  # Disable configuration during processing
    
    def on_batch_completed(self, success, results):
        """Handle batch processing completed signal."""
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
    
    def on_export_completed(self, success, path):
        """Handle export completed signal."""
        if success:
            self.update_status(f"Export completed successfully to {path}")
        else:
            self.update_status("Export failed.")


def launch_batch_gui():
    """Launch the batch processing GUI."""
    app = QApplication(sys.argv)
    window = BatchSpecGUI()
    window.show()
    return app.exec_()


if __name__ == "__main__":
    sys.exit(launch_batch_gui())