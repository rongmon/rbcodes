# panels/continuum_panel.py
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QGroupBox, QMessageBox,QCheckBox)
from PyQt5.QtCore import pyqtSignal
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Simple matplotlib canvas for the spectrum preview."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)
   
class ContinuumPanel(QWidget):
    """Panel that launches the existing interactive continuum fitter."""
    
    # Signals
    continuum_fitted = pyqtSignal(bool)  # Emitted when continuum is fitted (success/failure)
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.init_ui()
        
        # Connect to controller signals
        self.controller.spectrum_changed.connect(self.update_plot)
        self.controller.spectrum_changed.connect(self.reset)

    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Add instructions
        instructions = QLabel(
            "This panel allows you to fit a continuum to the spectrum. "
            "You can use the interactive fitter or set a flat continuum."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)
        
        # Launch button group
        launch_group = QGroupBox("Continuum Options")
        launch_layout = QVBoxLayout()

        # Add option to skip continuum normalization
        self.skip_continuum = QCheckBox("Skip continuum fitting")
        self.skip_continuum.setChecked(False)
        self.skip_continuum.toggled.connect(self.toggle_skip_continuum)
        self.skip_continuum.setToolTip("Skip continuum fitting and use a flat continuum instead")

        launch_layout.addWidget(self.skip_continuum)

        # Launch interactive fitter button
        self.launch_btn = QPushButton("Launch Interactive Continuum Fitter")
        self.launch_btn.clicked.connect(self.launch_fitter)
        self.launch_btn.setToolTip("Open the interactive tool to fit a continuum to the spectrum")
        launch_layout.addWidget(self.launch_btn)

        # Apply flat continuum button
        self.flat_btn = QPushButton("Apply Flat Continuum")
        self.flat_btn.clicked.connect(self.apply_flat_continuum)
        self.flat_btn.setEnabled(False)  # Initially disabled
        self.flat_btn.setToolTip("Use a flat line at value 1.0 as the continuum (skips fitting)")
        launch_layout.addWidget(self.flat_btn)

        # Apply current continuum button
        self.current_btn = QPushButton("Apply Current Continuum")
        self.current_btn.clicked.connect(self.apply_current_continuum)
        self.current_btn.setEnabled(False)  # Initially disabled
        self.current_btn.setToolTip("Apply the existing fitted continuum (if available)")
        launch_layout.addWidget(self.current_btn)

        #launch_layout.addWidget(self.flat_btn)
        
        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(self.launch_btn)
        buttons_layout.addSpacing(10)  # Add some space between buttons
        buttons_layout.addWidget(self.flat_btn)
        buttons_layout.addStretch()  # Push buttons to the left, or remove for center
        buttons_layout.addSpacing(10)  # Add some space between buttons
        buttons_layout.addWidget(self.current_btn)
        buttons_layout.addStretch()  # Push buttons to the left, or remove for center

        launch_layout.addLayout(buttons_layout)
        self.launch_btn.setMaximumWidth(320)
        self.flat_btn.setMaximumWidth(320)
        self.current_btn.setMaximumWidth(320)

        # Status label
        self.status_label = QLabel("No continuum fitted yet")
        launch_layout.addWidget(self.status_label)
        
        launch_group.setLayout(launch_layout)
        main_layout.addWidget(launch_group)
        
        # Plot area
        plot_group = QGroupBox("Preview")
        plot_layout = QVBoxLayout()
        
        self.canvas = MatplotlibCanvas(self)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        
        plot_group.setLayout(plot_layout)
        main_layout.addWidget(plot_group)    


    def reset(self):
        """Reset panel state when a new file is loaded."""
        self.status_label.setText("No continuum fitted yet")
        # Reset any continuum-related UI elements
        self.skip_continuum.setChecked(False)
        self.update_plot()



    def launch_fitter(self):
        """Launch the interactive continuum fitter."""
        if not self.controller.has_spectrum():
            QMessageBox.warning(self, "Error", "No spectrum loaded or not sliced.")
            return
        
        try:
            # Call the interactive fitter through rb_spec
            success = self.controller.fit_continuum()
            
            if success:
                self.status_label.setText("Continuum fitted successfully")
                self.update_plot()
                self.continuum_fitted.emit(True)
            else:
                self.status_label.setText("Fitting was cancelled or failed")
                self.continuum_fitted.emit(False)
        except Exception as e:
            print(f"Error launching fitter: {str(e)}")
            QMessageBox.warning(self, "Error", f"Failed to launch fitter: {str(e)}")
            self.continuum_fitted.emit(False)
    
    def update_plot(self):
        """Update the preview plot with current spectrum and continuum."""
        if not self.controller.has_spectrum():
            return
        
        try:
            # Get data from controller
            velocity, flux, error = self.controller.get_sliced_data()
            
            if velocity is None:
                return
            
            # Get continuum if available
            continuum = self.controller.get_continuum()
            
            # Clear existing plot
            self.canvas.fig.clear()
            
            if continuum is not None:
                # Create two subplots for original and normalized
                ax1 = self.canvas.fig.add_subplot(211)
                ax2 = self.canvas.fig.add_subplot(212, sharex=ax1)
                
                # Plot original data with continuum
                ax1.step(velocity, flux, 'k-', where='mid', label='Flux')
                ax1.plot(velocity, continuum, 'r-', linewidth=2, label='Continuum')
                #ax1.legend()
                ax1.set_ylabel('Flux')
                
                # Plot normalized flux
                norm_flux = flux / continuum
                ax2.step(velocity, norm_flux, 'k-', where='mid')
                ax2.axhline(y=1.0, color='r', linestyle='--')
                ax2.set_ylabel('Normalized Flux')
                ax2.set_xlabel('Velocity (km/s)')
            else:
                # Just plot the original spectrum
                ax = self.canvas.fig.add_subplot(111)
                ax.step(velocity, flux, 'k-', where='mid', label='Flux')
                ax.set_xlabel('Velocity (km/s)')
                ax.set_ylabel('Flux')
            
            # Use tight_layout to avoid overlap
            self.canvas.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            print(f"Error updating plot: {str(e)}")

    def toggle_skip_continuum(self, checked):
        """Toggle between interactive fitter and flat continuum."""
        self.launch_btn.setEnabled(not checked)
        self.flat_btn.setEnabled(checked)
        self.current_btn.setEnabled(checked)

    def apply_flat_continuum(self):
        """Apply a flat continuum (equal to 1 everywhere)."""
        if not self.controller.has_spectrum():
            QMessageBox.warning(self, "Error", "No spectrum loaded or not sliced.")
            return
        
        try:
            # Apply flat continuum
            success = self.controller.apply_flat_continuum()
            
            if success:
                self.status_label.setText("Flat continuum applied (value = 1.0)")
                self.update_plot()
                self.continuum_fitted.emit(True)
            else:
                self.status_label.setText("Failed to apply flat continuum")
                self.continuum_fitted.emit(False)
        except Exception as e:
            print(f"Error applying flat continuum: {str(e)}")
            QMessageBox.warning(self, "Error", f"Failed to apply flat continuum: {str(e)}")
            self.continuum_fitted.emit(False)    

    def apply_current_continuum(self):
        """Apply the currently stored continuum (no re-fitting)."""
        if not self.controller.has_spectrum():
            QMessageBox.warning(self, "Error", "No spectrum loaded or not sliced.")
            return

        try:
            success = self.controller.apply_current_continuum()

            if success:
                self.status_label.setText("Applied existing continuum to normalize spectrum")
                self.update_plot()
                self.continuum_fitted.emit(True)
            else:
                self.status_label.setText("Failed to apply existing continuum")
                self.continuum_fitted.emit(False)

        except Exception as e:
            print(f"Error applying current continuum: {str(e)}")
            QMessageBox.warning(self, "Error", f"Failed to apply existing continuum: {str(e)}")
            self.continuum_fitted.emit(False)