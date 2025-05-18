import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                           QPushButton, QDoubleSpinBox, QFormLayout, QGroupBox,
                           QMessageBox, QSpinBox, QCheckBox, QTableWidget,
                           QTableWidgetItem, QHeaderView)
from PyQt5.QtCore import pyqtSignal

class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib plots."""
    
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)

class MeasurementPanel(QWidget):
    """Panel for computing equivalent width and column density."""
    
    # Signals
    measurement_completed = pyqtSignal(dict)  # Emitted when measurement is complete
    
    def __init__(self, controller):
        super().__init__()
        self.controller = controller
        self.measurement_results = None  # Store the most recent results

    
        # Initialize with default values
        self.vmin_default = -200
        self.vmax_default = 200
        
        # Check if controller already has vmin/vmax values (e.g., from loaded JSON)
        if controller.has_spectrum() and hasattr(controller.spec, 'vmin') and hasattr(controller.spec, 'vmax'):
            # Use existing values if available
            self.vmin_default = controller.spec.vmin
            self.vmax_default = controller.spec.vmax
            print(f"Using existing velocity range from loaded data: [{self.vmin_default}, {self.vmax_default}]")
    
        self.init_ui()
        
        # Connect to controller signals
        self.controller.spectrum_changed.connect(self.update_plot)
        self.controller.spectrum_changed.connect(self.update_velocity_range)
        
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Top section with two columns (40% of total height)
        top_widget = QWidget()
        top_layout = QHBoxLayout(top_widget)
        
        # Left side: controls (smaller width)
        controls_widget = QWidget()
        controls_layout = QVBoxLayout(controls_widget)
        
        # Measurement controls
        measurement_group = QGroupBox("Equivalent Width Measurement")
        measurement_layout = QFormLayout()
        
        # Velocity range
        velocity_layout = QHBoxLayout()
        
        self.vmin_spinbox = QSpinBox()
        self.vmin_spinbox.setRange(-10000, 10000)
        self.vmin_spinbox.setValue(-200)
        self.vmin_spinbox.setSingleStep(10)
        
        self.vmax_spinbox = QSpinBox()
        self.vmax_spinbox.setRange(-10000, 10000)
        self.vmax_spinbox.setValue(200)
        self.vmax_spinbox.setSingleStep(10)
        
        velocity_layout.addWidget(QLabel("Min:"))
        velocity_layout.addWidget(self.vmin_spinbox)
        velocity_layout.addWidget(QLabel("Max:"))
        velocity_layout.addWidget(self.vmax_spinbox)
        
        measurement_layout.addRow("Velocity Range (km/s):", velocity_layout)
        
        # SNR calculation option
        self.snr_checkbox = QCheckBox("Calculate Signal-to-Noise Ratio")
        self.snr_checkbox.setChecked(True)
        
        # Binsize for SNR
        self.binsize_spinbox = QSpinBox()
        self.binsize_spinbox.setRange(1, 20)
        self.binsize_spinbox.setValue(3)
        self.binsize_spinbox.setEnabled(True)
        
        # Connect SNR checkbox to binsize enabled state
        self.snr_checkbox.toggled.connect(self.binsize_spinbox.setEnabled)
        
        snr_layout = QHBoxLayout()
        snr_layout.addWidget(self.snr_checkbox)
        snr_layout.addWidget(QLabel("Bin Size:"))
        snr_layout.addWidget(self.binsize_spinbox)
        
        measurement_layout.addRow("", snr_layout)
        
        # Show plot option
        self.plot_checkbox = QCheckBox("Show Diagnostic Plot")
        self.plot_checkbox.setChecked(True)
        measurement_layout.addRow("", self.plot_checkbox)
        
        # Compute button
        self.compute_btn = QPushButton("Compute Equivalent Width")
        self.compute_btn.clicked.connect(self.compute_equivalent_width)
        self.compute_btn.setToolTip("Calculate the equivalent width and column density within the selected velocity range")

        measurement_layout.addRow("", self.compute_btn)
    
        # Interactive selection option
        self.interactive_selection = QCheckBox("Enable Interactive Selection")
        self.interactive_selection.setChecked(False)
        self.interactive_selection.toggled.connect(self.toggle_interactive_selection)
        self.interactive_selection.setToolTip("Enable clicking on the plot to set velocity limits")

        measurement_layout.addRow("", self.interactive_selection) 

        # Other input elements
        self.vmin_spinbox.setToolTip("Minimum velocity for EW calculation (km/s)")
        self.vmax_spinbox.setToolTip("Maximum velocity for EW calculation (km/s)")
        self.snr_checkbox.setToolTip("Calculate the signal-to-noise ratio")
        self.binsize_spinbox.setToolTip("Bin size for SNR calculation")
        self.plot_checkbox.setToolTip("Show diagnostic plot after calculation")
    
        # Tooltip help:
        help_text = QLabel("When enabled, click on plot to set integration limits. Press 'r' to reset.")
        help_text.setWordWrap(True)
        help_text.setStyleSheet("color: gray; font-style: italic;")
        measurement_layout.addRow("", help_text)
        
        measurement_group.setLayout(measurement_layout)
        controls_layout.addWidget(measurement_group)
        controls_layout.addStretch(1)  # Add stretch to push everything to the top
        
        # Right side: results table (larger width)
        results_group = QGroupBox("Measurement Results")
        results_layout = QVBoxLayout()
        
        # Table for results
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(2)
        self.results_table.setHorizontalHeaderLabels(["Measurement", "Value"])
        self.results_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.results_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.results_table.setRowCount(0)  # No rows initially
        
        results_layout.addWidget(self.results_table)
        results_group.setLayout(results_layout)
        
        # Add widgets to top layout with desired proportions (40:60)
        top_layout.addWidget(controls_widget, 3)  # 40%
        top_layout.addWidget(results_group, 7)    # 60%
        
        # Bottom section: plot (60% of total height)
        plot_group = QGroupBox("Measurement Preview")
        plot_layout = QVBoxLayout()
        
        self.canvas = MatplotlibCanvas(self, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)


        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)


        
        plot_group.setLayout(plot_layout)
        
        # Add sections to main layout with desired proportions
        main_layout.addWidget(top_widget, 3)    # 35% for top section
        main_layout.addWidget(plot_group, 7)    # 65% for bottom section

    def update_velocity_range(self):
        """Update velocity range inputs if spectrum has changed and contains vmin/vmax values."""
        if self.controller.has_spectrum() and hasattr(self.controller.spec, 'vmin') and hasattr(self.controller.spec, 'vmax'):
            # Only update if the values differ from current spinbox values
            new_vmin = self.controller.spec.vmin
            new_vmax = self.controller.spec.vmax
            
            if new_vmin != self.vmin_spinbox.value() or new_vmax != self.vmax_spinbox.value():
                print(f"Updating velocity range to: [{new_vmin}, {new_vmax}]")
                self.vmin_spinbox.setValue(new_vmin)
                self.vmax_spinbox.setValue(new_vmax)
    

    def toggle_interactive_selection(self, enabled):
        """Enable or disable interactive velocity range selection."""
        if enabled:
            # Connect click and key events
            self.canvas.mpl_connect('button_press_event', self.on_plot_click)
            self.canvas.mpl_connect('key_press_event', self.on_plot_key_press)

            # Give focus to the canvas
            self.canvas.setFocus()
            
            # Update status with instructions
            self.setToolTip("Click on plot to set vmin/vmax. Press 'r' to reset.")
            
            # Update plot to show current selection points more prominently
            self.update_plot()
        else:
            # Disconnect events (not strictly necessary as the connections
            # would be automatically garbage collected)
            # Disconnect events if they exist
            if hasattr(self, 'click_cid'):
                self.canvas.mpl_disconnect(self.click_cid)
                delattr(self, 'click_cid')
            if hasattr(self, 'key_cid'):
                self.canvas.mpl_disconnect(self.key_cid)
                delattr(self, 'key_cid')            
            self.setToolTip("")
    
    def on_plot_click(self, event):
        """Handle mouse clicks on the plot for setting vmin/vmax."""
        if not event.inaxes or not self.interactive_selection.isChecked():
            return
            
        # Get the x-coordinate (velocity)
        velocity = event.xdata
        
        # Determine whether to set vmin or vmax based on which is closer
        current_vmin = self.vmin_spinbox.value()
        current_vmax = self.vmax_spinbox.value()
        
        # If velocity is outside current range, set the closest endpoint
        if velocity < current_vmin:
            self.vmin_spinbox.setValue(int(velocity))
        elif velocity > current_vmax:
            self.vmax_spinbox.setValue(int(velocity))
        else:
            # If inside current range, set whichever endpoint is closer
            if abs(velocity - current_vmin) < abs(velocity - current_vmax):
                self.vmin_spinbox.setValue(int(velocity))
            else:
                self.vmax_spinbox.setValue(int(velocity))
        
        # Update the plot
        self.update_plot()

    def on_plot_key_press(self, event):
        """Handle key presses on the plot."""
        if not self.interactive_selection.isChecked():
            return
            
        if event.key == 'r':
            # Reset velocity limits to default
            self.vmin_spinbox.setValue(-200)
            self.vmax_spinbox.setValue(200)
            self.update_plot()
        
    def compute_equivalent_width(self):
        """Compute equivalent width using rb_spec's compute_EW method."""
        if not self.controller.has_spectrum():
            QMessageBox.warning(self, "Error", "No spectrum loaded or not prepared.")
            return
            
        if not self.controller.has_continuum():
            QMessageBox.warning(self, "Error", "Continuum must be fitted before measuring equivalent width.")
            return
        
        try:
            # Get parameters
            vmin = self.vmin_spinbox.value()
            vmax = self.vmax_spinbox.value()
            snr = self.snr_checkbox.isChecked()
            binsize = self.binsize_spinbox.value() if snr else 1
            plot = self.plot_checkbox.isChecked()
            
            # Compute equivalent width
            success, results = self.controller.compute_equivalent_width(
                vmin=vmin,
                vmax=vmax,
                snr=snr,
                binsize=binsize,
                plot=plot
            )
            
            if success:
                self.measurement_results = results
                
                # Update results table
                self.update_results_table(results)
                
                # Update plot
                self.update_plot()
                
                # Emit signal
                self.measurement_completed.emit(results)
            else:
                QMessageBox.warning(self, "Error", "Failed to compute equivalent width.")
        except Exception as e:
            print(f"Error computing equivalent width: {str(e)}")
            QMessageBox.warning(self, "Error", f"Failed to compute equivalent width: {str(e)}")
        
    def update_results_table(self, results):
        """Update the results table with the computation results."""
        # Clear existing rows
        self.results_table.setRowCount(0)
        
        # Define the measurements to display
        measurements = [
            ("Transition", results.get('transition_name', 'Unknown')),
            ("Wavelength", f"{results.get('transition_wave', 0):.2f} Å"),
            ("Velocity Range", f"{results.get('vmin', 0):.1f} to {results.get('vmax', 0):.1f} km/s"),
            ("Equivalent Width", f"{results.get('W', 0):.3f} ± {results.get('W_e', 0):.3f} Å"),
            ("Column Density (log)", f"{results.get('logN', 0):.2f} ± {results.get('logN_e', 0):.2f}"),
            ("Velocity Centroid", f"{results.get('vel_centroid', 0):.1f} km/s"),
            ("Velocity Dispersion", f"{results.get('vel_disp', 0):.1f} km/s"),
        ]
        
        # Add SNR if calculated
        if 'SNR' in results and results['SNR'] > 0:
            measurements.append(("Signal-to-Noise Ratio", f"{results.get('SNR', 0):.1f}"))
        
        # Add saturation information if available
        if 'line_saturation' in results:
            if results['line_saturation']:
                saturation = f"Yes ({results.get('saturation_fraction', 0)*100:.1f}%)"
            else:
                saturation = "No"
            measurements.append(("Line Saturation", saturation))
        
        # Populate the table
        self.results_table.setRowCount(len(measurements))
        
        for i, (name, value) in enumerate(measurements):
            self.results_table.setItem(i, 0, QTableWidgetItem(name))
            self.results_table.setItem(i, 1, QTableWidgetItem(str(value)))
    
    def update_plot(self):
        """Update the preview plot with normalized spectrum and measurement results."""
        if not self.controller.has_spectrum():
            return
        
        try:
            # Get normalized data from controller
            velocity, flux, error, continuum = self.controller.get_normalized_data()
            
            if velocity is None or continuum is None:
                return
            
            # Calculate normalized flux and error
            norm_flux = flux / continuum
            norm_error = error / continuum
            
            # Clear the existing figure
            self.canvas.axes.clear()
            
            # Plot normalized flux
            self.canvas.axes.step(velocity, norm_flux, 'k-', where='mid', label='Normalized Flux')
            self.canvas.axes.step(velocity, norm_error, 'r-', where='mid', alpha=0.5, label='Error')
            self.canvas.axes.axhline(y=1.0, color='g', linestyle='--', alpha=0.7)
            self.canvas.axes.set_ylabel('Normalized Flux')
            self.canvas.axes.set_xlabel('Velocity (km/s)')
            self.canvas.axes.grid(True, linestyle='--', alpha=0.5)
            
            # Get current velocity range from spinboxes
            vmin = self.vmin_spinbox.value()
            vmax = self.vmax_spinbox.value()
            
            # Highlight the EW measurement region
            if vmin < vmax:
                self.canvas.axes.axvspan(vmin, vmax, alpha=0.1, color='blue', label='EW Region')
                
                # Add vertical lines at vmin and vmax
                self.canvas.axes.axvline(x=vmin, color='b', linestyle=':', alpha=0.7)
                self.canvas.axes.axvline(x=vmax, color='b', linestyle=':', alpha=0.7)
            
            # If measurement results exist, add them to the plot
            if self.measurement_results:
                # Construct measurement text
                measurement_text = ""
                if 'W' in self.measurement_results and 'W_e' in self.measurement_results:
                    measurement_text += f"EW = {self.measurement_results['W']:.3f} ± {self.measurement_results['W_e']:.3f} Å\n"
                if 'logN' in self.measurement_results and 'logN_e' in self.measurement_results:
                    measurement_text += f"log N = {self.measurement_results['logN']:.2f} ± {self.measurement_results['logN_e']:.2f}"
                
                if measurement_text:
                    # Add title with transition name
                    transition_name = self.measurement_results.get('transition_name', 'Unknown')
                    title = f"{transition_name}: Equivalent Width Measurement"
                    
                    # Add saturation warning if line is saturated
                    if self.measurement_results.get('line_saturation', False):
                        title += " [SATURATED]"
                    
                    self.canvas.axes.set_title(title)
                    
                    # Add detailed text box
                    self.canvas.axes.text(0.98, 0.05, measurement_text,
                                         transform=self.canvas.axes.transAxes,
                                         horizontalalignment='right',
                                         verticalalignment='bottom',
                                         bbox=dict(facecolor='white', alpha=0.7, boxstyle='round'),
                                         fontsize=10)
            
            # Add markers at vmin and vmax if interactive selection is enabled
            if hasattr(self, 'interactive_selection') and self.interactive_selection.isChecked():
                # Add circle markers at vmin and vmax with labels
                self.canvas.axes.plot([vmin], [1.0], 'bo', ms=6, zorder=10)
                self.canvas.axes.plot([vmax], [1.0], 'bo', ms=6, zorder=10)
                
                # Add labels
                self.canvas.axes.text(vmin, 1.05, f'vmin: {vmin}', 
                                     ha='center', va='bottom', fontsize=8, color='blue')
                self.canvas.axes.text(vmax, 1.05, f'vmax: {vmax}', 
                                     ha='center', va='bottom', fontsize=8, color='blue')
            
            x_limits = self.canvas.axes.get_xlim()
            y_limits = self.canvas.axes.get_ylim()
            
            #self.canvas.fig.subplots_adjust(left=0.10, right=0.95, top=0.92, bottom=0.12)

            
            # Restore the original axis limits
            self.canvas.axes.set_xlim(x_limits)
            self.canvas.axes.set_ylim(y_limits)
            
            self.canvas.fig.tight_layout() 
            # Draw the canvas
            self.canvas.draw()
            
            # Connect interactive events and set focus when interactive selection is enabled
            if hasattr(self, 'interactive_selection') and self.interactive_selection.isChecked():
                # First disconnect any existing connections to avoid duplicates
                if hasattr(self, 'click_cid') and self.click_cid:
                    self.canvas.mpl_disconnect(self.click_cid)
                if hasattr(self, 'key_cid') and self.key_cid:
                    self.canvas.mpl_disconnect(self.key_cid)
                    
                # Connect the events
                self.click_cid = self.canvas.mpl_connect('button_press_event', self.on_plot_click)
                self.key_cid = self.canvas.mpl_connect('key_press_event', self.on_plot_key_press)
                
                # Set focus to the canvas so it can receive keyboard events
                self.canvas.setFocus()
            
        except Exception as e:
            print(f"Error updating measurement plot: {str(e)}")    