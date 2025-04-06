import sys
import os
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, 
                             QVBoxLayout, QWidget, QFileDialog, QSplitter,
                             QLabel, QHBoxLayout, QSizePolicy, QLineEdit, QComboBox)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from linetools.spectra.io import readspec
from your_module import RedshiftInputWidget

class SpectralPlot(FigureCanvas):
    """
    A Matplotlib canvas for displaying spectral plots.
    Supports interactive key commands for adjusting axis limits and quitting the application.
    """
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(SpectralPlot, self).__init__(self.fig)
        self.setParent(parent)
        self.spectra = []  # List to store spectra data
        self.axes = None  # Stores the Matplotlib subplot axes
        self.parent_window = parent  # Reference to the parent window for quitting
        
        # Configure size policy
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        # Enable focus to receive key events
        self.setFocusPolicy(Qt.StrongFocus)
        
        # Connect Matplotlib's key press event
        self.mpl_connect('key_press_event', self.on_key_press)
        
    def plot_spectra(self, spectra):
        """
        Plots the given list of spectra as subplots.
        
        :param spectra: List of spectra objects from linetools.XSpectrum1D
        """
        self.spectra = spectra
        self.fig.clear()
        
        num_spectra = len(spectra)
        self.axes = self.fig.subplots(num_spectra, 1, sharex=True)
        
        if num_spectra == 1:
            self.axes = [self.axes]  # Ensure axes is always a list
        
        for i, spec in enumerate(spectra):
            wave = spec.wavelength.value
            flux = spec.flux.value
            error = spec.sig.value if hasattr(spec, 'sig') else None
            
            self.axes[i].plot(wave, flux, 'k-', lw=1)
            if error is not None:
                self.axes[i].plot(wave, error, 'r-', lw=0.5, alpha=0.5)
            
            if i == num_spectra - 1:
                self.axes[i].set_xlabel('Wavelength')
            self.axes[i].set_ylabel('Flux')
            
            filename = os.path.basename(spec.filename)
            self.axes[i].set_title(filename, fontsize=10)
        
        self.fig.tight_layout()
        self.draw()
        
    def on_key_press(self, event):
        """
        Handles key press events for interactive adjustments of the plots.
        """
        if not self.spectra or self.axes is None or len(self.axes) == 0:
            if event.key == 'q':
                print("Exiting application...")
                self.parent_window.close()
            return
        
        x, y = event.xdata, event.ydata
        
        if event.key == 'q':
            print("Exiting application...")
            self.parent_window.close()
            return
        
        if x is None or y is None:
            print("No valid coordinates at cursor position")
            return
        
        ax_index = 0
        for i, ax in enumerate(self.axes):
            if event.inaxes == ax:
                ax_index = i
                break
        
        if event.key == 'x':
            current_xlim = self.axes[0].get_xlim()
            self.axes[0].set_xlim(x, current_xlim[1])
            self.draw()
        elif event.key == 'X':
            current_xlim = self.axes[0].get_xlim()
            self.axes[0].set_xlim(current_xlim[0], x)
            self.draw()
        elif event.key == 't':
            current_ylim = self.axes[ax_index].get_ylim()
            self.axes[ax_index].set_ylim(current_ylim[0], y)
            self.draw()
        elif event.key == 'b':
            current_ylim = self.axes[ax_index].get_ylim()
            self.axes[ax_index].set_ylim(y, current_ylim[1])
            self.draw()


class MainWindow(QMainWindow):
    """
    Main application window for displaying spectral plots.
    Includes a file selection button and a Matplotlib canvas.
    """
    def __init__(self):
        super(MainWindow, self).__init__()
        
        self.setWindowTitle("FITS File Viewer")
        self.setGeometry(100, 100, 1000, 800)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Create input box and dropdown menu
        input_layout = QHBoxLayout()
        self.zabs_input = QLineEdit()
        self.zabs_input.setPlaceholderText("Enter zabs")
        self.linelist_dropdown = QComboBox()
        self.linelist_dropdown.addItems(["None", "LLS", "LLS Small", "DLA"])
        input_layout.addWidget(QLabel("zabs:"))
        input_layout.addWidget(self.zabs_input)
        input_layout.addWidget(QLabel("Linelist:"))
        input_layout.addWidget(self.linelist_dropdown)
        main_layout.addLayout(input_layout)
        
        button_layout = QHBoxLayout()
        self.select_button = QPushButton("Select FITS Files")
        self.select_button.clicked.connect(self.select_fits_files)
        self.file_label = QLabel("No files selected")
        button_layout.addWidget(self.select_button)
        button_layout.addWidget(self.file_label)
        main_layout.addLayout(button_layout)
        
        self.canvas = SpectralPlot(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.canvas)

        #----
        # Create our redshift input widget
        self.redshift_widget = RedshiftInputWidget()
    
        # Connect to the submitted signal
        self.redshift_widget.submitted.connect(self.handle_redshift_submission)
    
        # Add it to your layout
        main_layout.addWidget(self.redshift_widget)
        #-----
        
        self.spectra = []  # Stores loaded spectra
        
        status_text = "Ready - Click on plot then use keys: x = set min x-limit, X = set max x-limit, "
        status_text += "t = set max y-limit, b = set min y-limit, q = quit application"
        self.statusBar().showMessage(status_text)

    def handle_redshift_submission(self, redshift, linelist):
        # Process the submitted values
        print(f"Processing: z={redshift}, list={linelist}")
    
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
