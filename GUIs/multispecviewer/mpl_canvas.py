from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt, pyqtSignal
import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np

mpl.use('Qt5Agg')

class MplCanvas(FigureCanvasQTAgg):
    send_message = pyqtSignal(str)
    send_z_est = pyqtSignal(list)
    send_scale_limits = pyqtSignal(list)
    send_gauss_num = pyqtSignal(int)
    send_extract1d = pyqtSignal(dict)

    def __init__(self, parent=None, width=5, height=3, dpi=100):
        self.figsize = [width, height]
        self.fig = Figure(figsize=(self.figsize[0], self.figsize[-1]), dpi=dpi)
        pad = 0.05
        self.fig.subplots_adjust(left=pad+0.02, bottom=pad*2, right=0.995, top=0.95)

        self.wave, self.flux, self.error = [], [], []
        self.init_xlims, self.init_ylims = [], []
        self.cur_xlims, self.cur_ylims = [], []

        # Connect functions to events
        super().__init__(self.fig)

        # Set focus policy on the canvas (not on the underlying figure)
        self.setFocusPolicy(Qt.ClickFocus)
        self.setFocus()  # Set focus on the canvas for keyboard events

        # Connect key press event
        self.cid_k = self.fig.canvas.mpl_connect('key_press_event', self.ontype)

        # Connect mouse click event (fix the issue with missing 'onclick' method)
        self.cid_m = self.fig.canvas.mpl_connect('button_press_event', self.onclick)

    def plot_specs(self, spectra_dict):
        """Plot multiple spectra in subplots with shared x-axis."""
        num_files = len(spectra_dict)
        self.fig.clf()

        # Create subplots
        self.axes = []
        for i in range(num_files):
            ax = self.fig.add_subplot(num_files, 1, i + 1, sharex=self.fig.axes[0] if i > 0 else None)
            self.axes.append(ax)
            ax.set_xlabel('Wavelength (Angstrom)')
            ax.set_ylabel('Flux')

        for i, (filename, spec_data) in enumerate(spectra_dict.items()):
            wave, flux, error = spec_data['wave'], spec_data['flux'], spec_data['error']
            self.axes[i].plot(wave, flux, color='black', label='Flux')
            self.axes[i].plot(wave, error, color='red', label='Error')
            self.axes[i].set_title(filename)

        self.draw()

    def replot(self, wave, new_spec, new_err):
        axes = self.fig.gca()
        self.cur_xlims = axes.get_xlim()
        self.cur_ylims = axes.get_ylim()

        self.axes.lines[0] = self.axes.plot(wave, new_err, color='red')
        self.axes.lines[1] = self.axes.plot(wave, new_spec, color='black')

        ytmp = np.nan_to_num(new_spec, nan=0., posinf=0., neginf=0.)
        self.axes.set_ylim(self.cur_ylims)
        self.axes.set_xlim(self.cur_xlims)

        del self.axes.lines[2:]
        self.draw()

    def ontype(self, event):
        if event.key == 'x':  # Set minimum xlim based on the click position
            if event.inaxes:
                min_val = event.xdata
                self.cur_xlims[0] = min_val
                self.axes.set_xlim(self.cur_xlims)
                self.draw()

        elif event.key == 'X':  # Set maximum xlim based on the click position
            if event.inaxes:
                max_val = event.xdata
                self.cur_xlims[1] = max_val
                self.axes.set_xlim(self.cur_xlims)
                self.draw()

        elif event.key == 'y':  # Set minimum ylim based on the click position
            if event.inaxes:
                min_val = event.ydata
                self.cur_ylims[0] = min_val
                self.axes.set_ylim(self.cur_ylims)
                self.draw()

        elif event.key == 'Y':  # Set maximum ylim based on the click position
            if event.inaxes:
                max_val = event.ydata
                self.cur_ylims[1] = max_val
                self.axes.set_ylim(self.cur_ylims)
                self.draw()

        elif event.key == 'r':  # Reset axes to original limits
            self.axes.set_xlim(self.init_xlims)
            self.axes.set_ylim(self.init_ylims)
            self.cur_xlims = self.init_xlims
            self.cur_ylims = self.init_ylims
            self.draw()

        elif event.key == 'h':  # Display help message
            help_msg = (
                "Available keys:\n"
                "'x'  : Set min xlim\n"
                "'X'  : Set max xlim\n"
                "'y'  : Set min ylim\n"
                "'Y'  : Set max ylim\n"
                "'r'  : Reset x and y limits\n"
                "'h'  : Display this help message"
            )
            self.send_message.emit(help_msg)

    def onclick(self, event):
        """
        Mouse click event handler to capture click coordinates.
        """
        if event.inaxes:
            x, y = event.xdata, event.ydata
            print(f"Mouse clicked at x={x}, y={y}")
            # You can perform additional actions here based on the click location.
