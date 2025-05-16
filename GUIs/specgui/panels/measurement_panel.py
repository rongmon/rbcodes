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
        self.init_ui()
        
        # Connect to controller signals
        self.controller.spectrum_changed.connect(self.update_plot)
    
    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Measurement controls
        measurement_group = QGroupBox("Equivalent Width Measurement")
        measurement_layout = QFormLayout()
        
        # Velocity range
        velocity_layout = QHBoxLayout()
        
        self.vmin_spinbox = QSpinBox()
        self.vmin_spinbox.setRange(-10000, 0)
        self.vmin_spinbox.setValue(-200)
        self.vmin_spinbox.setSingleStep(10)
        
        self.vmax_spinbox = QSpinBox()
        self.vmax_spinbox.setRange(0, 10000)
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
        measurement_layout.addRow("", self.compute_btn)
        
        measurement_group.setLayout(measurement_layout)
        main_layout.addWidget(measurement_group)
        
        # Results display
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
        
        # Save button will be added in the next panel (Output Panel)
        
        results_group.setLayout(results_layout)
        main_layout.addWidget(results_group)
        
        # Plot area
        plot_group = QGroupBox("Measurement Preview")
        plot_layout = QVBoxLayout()
        
        self.canvas = MatplotlibCanvas(self, width=8, height=9, dpi=100)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        
        plot_group.setLayout(plot_layout)
        main_layout.addWidget(plot_group)
    
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
        """Update the preview plot with current data and measurement results."""
        if not self.controller.has_spectrum():
            return
        
        try:
            # Use rb_spec's built-in plot if a measurement has been made and plot option is checked
            if self.measurement_results is not None and self.plot_checkbox.isChecked():
                # Clear the existing figure
                self.canvas.fig.clear()
                
                # Get rb_spec's figure for the measurement
                rb_spec_fig = self.controller.get_measurement_figure()
                
                if rb_spec_fig:
                    # Copy the rb_spec figure to our canvas
                    for i, ax in enumerate(rb_spec_fig.axes):
                        # Create a new axis in our figure
                        if i == 0:
                            new_ax = self.canvas.fig.add_subplot(2, 1, 1)
                        else:
                            new_ax = self.canvas.fig.add_subplot(2, 1, 2, sharex=self.canvas.fig.axes[0])
                        
                        # Copy the content from rb_spec's figure
                        for line in ax.lines:
                            new_ax.plot(line.get_xdata(), line.get_ydata(), 
                                      color=line.get_color(), linestyle=line.get_linestyle(),
                                      linewidth=0.5,marker=line.get_marker())
                        
                        # Copy axis labels and limits
                        new_ax.set_xlabel(ax.get_xlabel())
                        new_ax.set_ylabel(ax.get_ylabel())
                        new_ax.set_xlim(ax.get_xlim())
                        new_ax.set_ylim(ax.get_ylim())
                        
                        # Add a grid
                        new_ax.grid(True, linestyle='--', alpha=0.7)
                    
                    # Add title with key measurement results
                    if 'W' in self.measurement_results and 'W_e' in self.measurement_results:
                        ew = self.measurement_results['W']
                        ew_err = self.measurement_results['W_e']
                        transition = self.measurement_results.get('transition_name', 'Unknown')
                        
                        title = f"{transition}: W = {ew:.3f} ± {ew_err:.3f} Å"
                        
                        if 'line_saturation' in self.measurement_results and self.measurement_results['line_saturation']:
                            title += " [SATURATED]"
                            
                        self.canvas.fig.suptitle(title)
                
                # Apply tight layout
                self.canvas.fig.tight_layout()
                if self.measurement_results:
                    self.canvas.fig.subplots_adjust(top=0.9)  # Make room for the title
                
                # Redraw
                self.canvas.draw()
            else:
                # If no measurement, or plot option unchecked, just show the normalized spectrum
                velocity = None
                norm_flux = None
                norm_error = None
                
                if self.controller.has_continuum():
                    velocity, flux, error, continuum = self.controller.get_normalized_data()
                    if velocity is not None:
                        norm_flux = flux / continuum
                        norm_error = error / continuum
                
                if velocity is not None and norm_flux is not None:
                    # Clear the plot
                    self.canvas.axes.clear()
                    
                    # Plot normalized data
                    self.canvas.axes.step(velocity, norm_flux, 'k-', where='mid', label='Normalized Flux')
                    self.canvas.axes.step(velocity, norm_error, 'r-', where='mid', alpha=0.5, label='Error')
                    
                    # Add a horizontal line at y=1
                    self.canvas.axes.axhline(y=1.0, color='g', linestyle='--', alpha=0.7)
                    
                    # Plot velocity range if it's been set
                    vmin = self.vmin_spinbox.value()
                    vmax = self.vmax_spinbox.value()
                    
                    if vmin < vmax:
                        self.canvas.axes.axvline(x=vmin, color='b', linestyle=':', alpha=0.7)
                        self.canvas.axes.axvline(x=vmax, color='b', linestyle=':', alpha=0.7)
                        self.canvas.axes.axvspan(vmin, vmax, alpha=0.1, color='blue')
                    
                    # Labels and title
                    self.canvas.axes.set_xlabel('Velocity (km/s)')
                    self.canvas.axes.set_ylabel('Normalized Flux')
                    self.canvas.axes.set_title('Normalized Spectrum')
                    #self.canvas.axes.legend()
                    
                    # Redraw
                    self.canvas.fig.tight_layout()
                    self.canvas.draw()
        except Exception as e:
            print(f"Error updating measurement plot: {str(e)}")