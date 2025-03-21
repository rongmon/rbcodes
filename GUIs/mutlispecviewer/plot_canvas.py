from PyQt5.QtWidgets import QWidget, QVBoxLayout
from PyQt5.QtCore import Qt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class PlotCanvas(QWidget):
    def __init__(self, spectra_dict, parent=None):
        super().__init__(parent)
        self.spectra_dict = spectra_dict
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        self.setFocusPolicy(Qt.StrongFocus)  # Allow key events
        self.setFocus()  # Set initial focus on the widget
        self.plot_spectra()

    
    def keyPressEvent(self, event):
        """Ensure keypress events are detected at the widget level."""
        print(f"Qt Key Pressed: {event.text()}")  # Debug print
        if hasattr(self, "key_handler"):
            self.key_handler.on_key_press(event)
        super().keyPressEvent(event)  # Allow Qt's default key press handling if needed
    
    
    

    def plot_spectra(self):
        """Plots all spectra as subplots with a shared x-axis."""
        num_spectra = len(self.spectra_dict)
        self.figure.clear()

        if num_spectra == 0:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "No Data Loaded", ha='center', va='center', fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            axes = self.figure.subplots(num_spectra, 1, sharex=True)
            if num_spectra == 1:
                axes = [axes]

            for ax, (filename, data) in zip(axes, self.spectra_dict.items()):
                ax.plot(data['wave'], data['flux'], label=filename)
                ax.set_ylabel("Flux")
                ax.legend()

            axes[-1].set_xlabel("Wavelength")

        self.canvas.draw()
