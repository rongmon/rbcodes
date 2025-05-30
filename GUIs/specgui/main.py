# main.py
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget, QMessageBox
from PyQt5.QtCore import Qt

from rbcodes.GUIs.specgui.spectrum_controller import SpectrumController
from rbcodes.GUIs.specgui.panels.input_panel import InputPanel
from rbcodes.GUIs.specgui.panels.redshift_panel import RedshiftPanel
from rbcodes.GUIs.specgui.panels.transition_panel import TransitionPanel
from rbcodes.GUIs.specgui.panels.continuum_panel import ContinuumPanel
from rbcodes.GUIs.specgui.panels.measurement_panel import MeasurementPanel
from rbcodes.GUIs.specgui.panels.output_panel import OutputPanel

class RbSpecGUI(QMainWindow):
    """Main application window for the rb_spec GUI."""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("rb_spec GUI")
        self.setMinimumSize(800, 600)
        
        # Create the spectrum controller
        self.controller = SpectrumController()
        
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
                min-width: 120px;
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
        self.input_panel = InputPanel(self.controller)
        self.tabs.addTab(self.input_panel, "Input Spectrum")
        
        self.redshift_panel = RedshiftPanel(self.controller)
        self.tabs.addTab(self.redshift_panel, "Set Redshift")
        
        self.transition_panel = TransitionPanel(self.controller)
        self.tabs.addTab(self.transition_panel, "Select Transition")

        self.continuum_panel = ContinuumPanel(self.controller)
        self.tabs.addTab(self.continuum_panel, "Fit Continuum")

        self.measurement_panel = MeasurementPanel(self.controller)
        self.tabs.addTab(self.measurement_panel, "Measurements")

        self.output_panel = OutputPanel(self.controller)
        self.tabs.addTab(self.output_panel, "Save Results")


        
        # Initially disable redshift and transition tabs
        self.tabs.setTabEnabled(1, False)
        self.tabs.setTabEnabled(2, False)
        self.tabs.setTabEnabled(3, False)  # Continuum tab
        self.tabs.setTabEnabled(4, False)  # Measurement tab
        self.tabs.setTabEnabled(5, False)  # Output tab

        # Connect signals
        self.input_panel.spectrum_loaded.connect(self.on_spectrum_loaded)
        self.redshift_panel.redshift_applied.connect(self.on_redshift_applied)
        self.transition_panel.spectrum_sliced.connect(self.on_spectrum_sliced)
        self.continuum_panel.continuum_fitted.connect(self.on_continuum_fitted)
        self.measurement_panel.measurement_completed.connect(self.on_measurement_completed)



        
        # Show status message
        self.statusBar().showMessage("Ready. Load a spectrum to begin.")
    

    def on_spectrum_loaded(self, success, has_redshift=False, has_transition=False):
        """Handle spectrum loaded signal."""
        if success:
            self.statusBar().showMessage("Spectrum loaded successfully.")
            
            # Reset all tab states when loading a new file
            # Disable all tabs except input and redshift initially
            self.tabs.setTabEnabled(1, True)  # Redshift tab
            self.tabs.setTabEnabled(2, False)  # Transition tab
            self.tabs.setTabEnabled(3, False)  # Continuum tab
            self.tabs.setTabEnabled(4, False)  # Measurement tab
            self.tabs.setTabEnabled(5, False)  # Output tab
            
            if has_redshift:
                # Pre-populate redshift panel with value from JSON
                zabs, trans_wave, trans_name, vel_limits = self.controller.get_json_info()
                if zabs is not None:
                    self.redshift_panel.set_redshift(zabs)
                    self.redshift_panel.apply_redshift()  # Automatically apply the loaded redshift
                
                # Enable transition tab if we also have transition info
                if has_transition:
                    self.tabs.setTabEnabled(2, True)
                    
                    # Pre-populate transition panel
                    if trans_wave is not None:
                        self.transition_panel.set_transition(trans_wave, trans_name)
                        
                        # Automatically slice with loaded values
                        self.transition_panel.slice_spectrum()
                        
                    # Set tab to the furthest we can go
                    self.tabs.setCurrentIndex(2)
                else:
                    # Just go to redshift tab
                    self.tabs.setCurrentIndex(1)
            else:
                # Just go to redshift tab
                self.tabs.setCurrentIndex(1)
        else:
            self.statusBar().showMessage("Failed to load spectrum.")
            QMessageBox.warning(self, "Load Error", "Failed to load spectrum. Check console for details.")    
    
    def on_redshift_applied(self, redshift):
        """Handle redshift applied signal."""
        if redshift == -99999:  # Special value indicating reset
            # Reset downstream tabs
            self.statusBar().showMessage("Redshift has been reset. Select a new redshift.")
            self.tabs.setTabEnabled(2, False)  # Disable transition tab
            self.tabs.setTabEnabled(3, False)  # Disable continuum tab
            self.tabs.setTabEnabled(4, False)  # Disable measurement tab
            self.tabs.setTabEnabled(5, False)  # Disable output tab
            
            # Reset controller state
            self.controller.reset_after_redshift()
        else:
            # Original code for normal redshift application
            self.statusBar().showMessage(f"Redshift z = {redshift} applied successfully.")
            # Enable the transition tab
            self.tabs.setTabEnabled(2, True)
            # Switch to the transition tab
            self.tabs.setCurrentIndex(2)
    


    def on_spectrum_sliced(self, transition, vmin, vmax):
        """Handle spectrum sliced signal."""
        if transition == -99999:  # Special value indicating reset
            # Reset downstream tabs
            self.statusBar().showMessage("Transition has been reset. Select a new transition.")
            self.tabs.setTabEnabled(3, False)  # Disable continuum tab
            self.tabs.setTabEnabled(4, False)  # Disable measurement tab
            self.tabs.setTabEnabled(5, False)  # Disable output tab
            
            # Reset controller state
            self.controller.reset_after_transition()
        else:
            # Original code for normal transition slicing
            self.statusBar().showMessage(f"Spectrum sliced around {transition:.2f} Å")
            # Enable the continuum tab
            self.tabs.setTabEnabled(3, True)
            # Switch to the continuum tab
            self.tabs.setCurrentIndex(3)
            
    def on_continuum_fitted(self, success):
        """Handle continuum fitted signal."""
        if success:
            self.statusBar().showMessage("Continuum fitted successfully")
            # Enable the measurement tab
            self.tabs.setTabEnabled(4, True)
            # Switch to the measurement tab
            self.tabs.setCurrentIndex(4)
        else:
            self.statusBar().showMessage("Failed to fit continuum")

    def on_measurement_completed(self, results):
        """Handle measurement completed signal."""
        if results:
            self.statusBar().showMessage("Measurements computed successfully")
            # Enable the output tab
            self.tabs.setTabEnabled(5, True)
            # Switch to the output tab
            self.tabs.setCurrentIndex(5)


