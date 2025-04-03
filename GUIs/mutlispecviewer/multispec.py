import sys
import os
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, 
                             QVBoxLayout, QWidget, QFileDialog, QSplitter,
                             QLabel, QHBoxLayout, QSizePolicy,QInputDialog)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from linetools.spectra.io import readspec
from RedshiftInputWidget import RedshiftInputWidget
from MessageBox import MessageBox
from astropy.convolution import convolve, Box1DKernel, Gaussian2DKernel
from rbcodes.IGM import rb_setline as rb_setline
class SpectralPlot(FigureCanvas):
    """
    A Matplotlib canvas for displaying spectral plots.
    Supports interactive key commands for adjusting axis limits and quitting the application.
    """
    # Add this line to the __init__ method of SpectralPlot:
    def __init__(self, parent=None, width=12, height=10, dpi=100, message_box=None):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(SpectralPlot, self).__init__(self.fig)
        self.setParent(parent)
        self.spectra = []  # List to store spectra data
        self.axes = None  # Stores the Matplotlib subplot axes
        self.parent_window = parent  # Reference to the parent window for quitting
        self.scale = 1. # 1D spec convolution kernel size
        self.redshift_lines = []  # Add this line to store references to plotted lines
    
        # Store reference to the message box
        self.message_box = message_box
        
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
        self.num_spectra=num_spectra
        self.axes = self.fig.subplots(num_spectra, 1, sharex=True)
        
        if num_spectra == 1:
            self.axes = [self.axes]  # Ensure axes is always a list
        
        for i, spec in enumerate(spectra):
            wave = spec.wavelength.value
            flux = spec.flux.value
            error = spec.sig.value if hasattr(spec, 'sig') else None
            filename = os.path.basename(spec.filename)

            self.plot_one_spec(wave,flux,error,i,filename)
            
            #self.axes[i].plot(wave, flux, 'k-', lw=1)
            #if error is not None:
            #    self.axes[i].plot(wave, error, 'r-', lw=0.5, alpha=0.5)
            
            #if i == num_spectra - 1:
            #    self.axes[i].set_xlabel('Wavelength')
            #self.axes[i].set_ylabel('Flux')            
            #self.axes[i].set_title(filename, fontsize=10)
        
        self.fig.tight_layout()
        self.draw()
        
        # Send message if message box is available
        if self.message_box:
            self.message_box.on_sent_message(f"Plotted {num_spectra} spectra", "#008000")
        
    def on_key_press(self, event):
        """
        Handles key press events for interactive adjustments of the plots.
        """
        if not self.spectra or self.axes is None or len(self.axes) == 0:
            if event.key == 'q':
                print("Exiting application...")
                if self.message_box:
                    self.message_box.on_sent_message("Exiting application...", "#FF0000")
                self.parent_window.close()
            return
        
        x, y = event.xdata, event.ydata
        
        if event.key == 'q':
            print("Exiting application...")
            if self.message_box:
                self.message_box.on_sent_message("Exiting application...", "#FF0000")
            self.parent_window.close()
            return
        
        if x is None or y is None:
            print("No valid coordinates at cursor position")
            if self.message_box:
                self.message_box.on_sent_message("No valid coordinates at cursor position", "#FF0000")
            return
        
        ax_index = 0
        for i, ax in enumerate(self.axes):
            if event.inaxes == ax:
                ax_index = i
                break
        
        if event.key == 'x':
            current_xlim = self.axes[0].get_xlim()
            self.axes[0].set_xlim(x, current_xlim[1])
            if self.message_box:
                self.message_box.on_sent_message(f"Set minimum x-limit to {x:.2f}", "#0000FF")
            self.draw()
        elif event.key == 'X':
            current_xlim = self.axes[0].get_xlim()
            self.axes[0].set_xlim(current_xlim[0], x)
            if self.message_box:
                self.message_box.on_sent_message(f"Set maximum x-limit to {x:.2f}", "#0000FF")
            self.draw()
        elif event.key == 't':
            current_ylim = self.axes[ax_index].get_ylim()
            self.axes[ax_index].set_ylim(current_ylim[0], y)
            if self.message_box:
                self.message_box.on_sent_message(f"Set maximum y-limit to {y:.2f} on panel {ax_index+1}", "#0000FF")
            self.draw()
        elif event.key == 'b':
            current_ylim = self.axes[ax_index].get_ylim()
            self.axes[ax_index].set_ylim(y, current_ylim[1])
            if self.message_box:
                self.message_box.on_sent_message(f"Set minimum y-limit to {y:.2f} on panel {ax_index+1}", "#0000FF")
            self.draw()

        elif event.key=='S':
            self.scale += 2
            # First grab the corresponding spectra 
            spec = self.spectra[ax_index]
            filename = os.path.basename(spec.filename)
            wave = spec.wavelength.value
            flux = spec.flux.value
            error = spec.sig.value if hasattr(spec, 'sig') else None
            new_flux = convolve(flux, Box1DKernel(self.scale))
            new_err = convolve(error, Box1DKernel(self.scale))
            self.replot(wave, new_flux,new_err,ax_index,filename)
            self.draw()

            self.message_box.on_sent_message(f'Convolutional kernel size = {int(self.scale)}.')

        elif event.key=='U':
            self.scale -= 2
            # First grab the corresponding spectra 
            spec = self.spectra[ax_index]
            filename = os.path.basename(spec.filename)
            wave = spec.wavelength.value
            flux = spec.flux.value
            error = spec.sig.value if hasattr(spec, 'sig') else None
            new_flux = convolve(flux, Box1DKernel(self.scale))
            new_err = convolve(error, Box1DKernel(self.scale))
            self.replot(wave, new_flux,new_err,ax_index,filename)
            self.draw()

            self.message_box.on_sent_message(f'Convolutional kernel size = {int(self.scale)}.')



        elif event.key=='[':
            xlim = self.axes[ax_index].get_xlim()
            delx = (xlim[-1] - xlim[0])
            self.axes[ax_index].set_xlim([xlim[0] - delx, xlim[0]])
            if self.message_box:
                self.message_box.on_sent_message(f"Shifted view left on panel {ax_index+1}", "#0000FF")
            self.draw()
        elif event.key==']':
            xlim = self.axes[ax_index].get_xlim()
            delx = (xlim[-1] - xlim[0])
            self.axes[ax_index].set_xlim([xlim[1], xlim[1] + delx])
            if self.message_box:
                self.message_box.on_sent_message(f"Shifted view right on panel {ax_index+1}", "#0000FF")
            self.draw()
        elif event.key == 'Y':
            Windowname='Manual y-Limits'
            instruction='Input range (e.g. 0.,2.)'
            ylim, ok = QInputDialog.getText(self,Windowname,instruction)
            if ok:
                ylimit = ylim.split(',')
                ylimit = np.array(ylimit).astype('float32')
                self.axes[ax_index].set_ylim(ylimit)
                self.draw()

        #zoom out of xrange
        elif (event.key == 'o'):
            xlim=self.axes[ax_index].get_xlim()
            ylim=self.axes[ax_index].get_ylim()
            xcen = (xlim[0]+xlim[1])/2.0
            delx   = xlim[1] - xcen
            self.axes[ax_index].set_xlim([xcen - 1.5*delx,xcen + 1.5*delx])
            self.draw() 



        



    def replot(self, wave, new_spec, new_err,ax_index,filename):
        '''Re-plot smoothed/unsmoothed spectrum
        '''


        self.cur_xlims = self.axes[ax_index].get_xlim()
        self.cur_ylims = self.axes[ax_index].get_ylim()

        self.axes[ax_index].clear() 

        self.plot_one_spec(wave,new_spec,new_err,ax_index,filename)
        
        #self.axes[ax_index].plot(wave, new_err, color='red')# label='Error')
        #self.axes[ax_index].plot(wave, new_spec, color='black')#, label='Flux')
        # for a better y range
        ytmp = np.nan_to_num(new_spec, nan=0., posinf=0., neginf=0.)
        self.axes[ax_index].set_ylim(self.cur_ylims)
        self.axes[ax_index].set_xlim(self.cur_xlims)

        #del self.axes[ax_index].lines[2:]
        self.draw()

    def plot_one_spec(self,wave,flux,error,index,filename):
        
        self.axes[index].plot(wave, flux, 'k-', lw=1)
        if error is not None:
            self.axes[index].plot(wave, error, 'r-', lw=0.5, alpha=0.5)
            
        if index == self.num_spectra - 1:
            self.axes[index].set_xlabel('Wavelength')
        self.axes[index].set_ylabel('Flux')
            
        self.axes[index].set_title(filename, fontsize=10)

    def set_redshift_data(self, redshift, linelist):
        """
        Receives redshift and linelist data from the main window.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        """
        self.redshift = redshift
        self.linelist = linelist
        
        # If spectra are already loaded, update the plot with the new redshift lines
        if self.spectra and len(self.spectra) > 0:
            # Add debug message
            if self.message_box:
                self.message_box.on_sent_message(f"Setting redshift z={redshift} with {linelist} line list. Plotting lines...", "#0000FF")
            self.plot_redshift_lines()
            if self.message_box:
                self.message_box.on_sent_message(f"Applied redshift z={redshift} with {linelist} line list", "#008000")
        else:
            if self.message_box:
                self.message_box.on_sent_message(f"No spectra loaded. Load spectra before applying redshift.", "#FF0000")

    def plot_redshift_lines(self):
        """
        Plots spectral lines based on the redshift and linelist.
        This method draws vertical lines at the expected wavelengths of 
        common spectral features adjusted by the redshift.
        """
        if not hasattr(self, 'redshift') or not hasattr(self, 'linelist'):
            if self.message_box:
                self.message_box.on_sent_message("Redshift or linelist not set", "#FF0000")
            return
        
        # Clear any existing redshift lines first
        self.clear_redshift_lines()
        
        # Read linelist - use a different variable name (not 'line')
        line_list = rb_setline.read_line_list(self.linelist)
        
        self.redshift_lines = []  # Store references to lines for later removal
        
          
        # Draw vertical lines for each spectral feature
        lines_plotted = 0
        
        for ix in range(0, len(line_list)):
            # Loop through each transition
            rest_wavelength = line_list[ix]['wrest']
            transition_name = line_list[ix]['ion']
            # Calculate observed wavelength using the redshift
            observed_wavelength = rest_wavelength * (1. + self.redshift)
            
            # Draw on all axes
            for i, ax in enumerate(self.axes):
                # Get current axis limits
                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                
                # Always draw the line (even if outside current view)
                line_obj = ax.axvline(x=observed_wavelength, color='blue', linestyle='--', alpha=0.7)
                self.redshift_lines.append(line_obj)
                
                # Only add label if it's in view
                if xlim[0] <= observed_wavelength <= xlim[1]:
                    y_pos = ylim[0] + 0.85 * (ylim[1] - ylim[0])  # 85% of the way up
                    text = ax.text(observed_wavelength, y_pos, transition_name, rotation=90, 
                          horizontalalignment='right', verticalalignment='top',
                          fontsize=8, color='blue')
                    self.redshift_lines.append(text)
                    lines_plotted += 1
        
        self.fig.canvas.draw_idle()  # Ensure the figure updates
        
        if self.message_box:
            if lines_plotted > 0:
                self.message_box.on_sent_message(f"Successfully plotted {lines_plotted} spectral lines", "#008000")
            else:
                self.message_box.on_sent_message("No lines visible in current view. Try adjusting x-axis limits.", "#FFA500")

    def clear_redshift_lines(self):
        """
        Removes any previously plotted redshift lines.
        """
        if hasattr(self, 'redshift_lines') and self.redshift_lines:
            for line_obj in self.redshift_lines:
                try:
                    line_obj.remove()
                except:
                    pass
            self.redshift_lines = []

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
        
        # Layout of a file selection Button
        button_layout = QHBoxLayout()
        self.select_button = QPushButton("Select FITS Files")
        self.select_button.clicked.connect(self.select_fits_files)
        self.file_label = QLabel("No files selected")
        button_layout.addWidget(self.select_button)
        button_layout.addWidget(self.file_label)
        main_layout.addLayout(button_layout)
        
        # Create message box first so we can pass it to the SpectralPlot
        self.message_box = MessageBox()
        
        # Create SpectralPlot with reference to the message box
        self.canvas = SpectralPlot(self, message_box=self.message_box)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.canvas)

        # Create bottom widget container with redshift inputs and message box
        bottom_widget = QWidget()
        bottom_layout = QHBoxLayout(bottom_widget)
        
        # Create redshift input widget
        self.redshift_widget = RedshiftInputWidget()
        # Connect to the submitted signal
        self.redshift_widget.submitted.connect(self.handle_redshift_submission)
        
        # Add both widgets to the bottom layout
        bottom_layout.addWidget(self.redshift_widget)
        bottom_layout.addWidget(self.message_box)
        
        # Add the bottom widget to the main layout
        main_layout.addWidget(bottom_widget)
        
        self.spectra = []  # Stores loaded spectra
        
        status_text = "Ready - Click on plot then use keys: x = set min x-limit, X = set max x-limit, "
        status_text += "t = set max y-limit, b = set min y-limit, q = quit application"
        self.statusBar().showMessage(status_text)
        self.message_box.on_sent_message("Application ready. Select FITS files to begin.", "#000000")

    def handle_redshift_submission(self, redshift, linelist):
        """
        Processes the redshift and linelist data submitted by the RedshiftInputWidget
        and sends it to the SpectralPlot instance.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        """
        try:
            # Convert redshift to float (in case it's received as string)
            redshift = float(redshift)
            
            print(f"Processing: z={redshift}, list={linelist}")
            
            # Send a message to the message box
            self.message_box.on_sent_message(f"Processing redshift z={redshift} with {linelist} line list", "#008000")
            
            # Pass the data to the SpectralPlot instance
            self.canvas.set_redshift_data(redshift, linelist)
        except ValueError:
            error_msg = f"Invalid redshift value: {redshift}. Please enter a number."
            print(error_msg)
            self.message_box.on_sent_message(error_msg, "#FF0000")
        except Exception as e:
            error_msg = f"Error processing redshift: {str(e)}"
            print(error_msg)
            self.message_box.on_sent_message(error_msg, "#FF0000")



    def select_fits_files(self):
        """
        Opens a file dialog for selecting FITS files and loads them for plotting.
        """
        options = QFileDialog.Options()
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select FITS Files", "", 
                                                     "FITS Files (*.fits *.fit);;All Files (*)", 
                                                     options=options)
        
        if file_paths:
            self.statusBar().showMessage(f"Loading {len(file_paths)} files...")
            self.message_box.on_sent_message(f"Loading {len(file_paths)} files...", "#0000FF")
            self.load_fits_files(file_paths)
            self.file_label.setText(f"{len(file_paths)} files selected")
            self.statusBar().showMessage(f"Loaded {len(file_paths)} files - Use keyboard shortcuts for adjustments")
            self.message_box.on_sent_message(f"Successfully loaded {len(file_paths)} files", "#008000")
            self.canvas.setFocus()
    
    def load_fits_files(self, file_paths):
        """
        Loads and plots the selected FITS files.
        """
        try:
            self.spectra = []
            for file_path in file_paths:
                spec = readspec(file_path)
                self.spectra.append(spec)
            self.canvas.plot_spectra(self.spectra)
        except Exception as e:
            error_message = f"Error: {str(e)}"
            self.statusBar().showMessage(error_message)
            self.message_box.on_sent_message(error_message, "#FF0000")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())