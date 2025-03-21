import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QFileDialog
from io_module import read_fits_files
from mpl_canvas import MplCanvas

class MultiSpecViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Multi Spectra Viewer')

        # Canvas for displaying spectra
        self.canvas = MplCanvas(self)
        self.setCentralWidget(self.canvas)

        # Button to load files
        self.load_data_button = QPushButton("Load FITS Files")
        self.load_data_button.clicked.connect(self.load_data)

        # Layout for central widget
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.load_data_button)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def load_data(self):
        # Open file dialog to choose multiple FITS files
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Open FITS Files", "", "FITS Files (*.fits)")
        
        if file_paths:
            # Read the selected FITS files and store spectra in dictionary
            self.spectra_dict = read_fits_files(file_paths)
            
            # Plot each spectrum in the dictionary
            for filename, data in self.spectra_dict.items():
                self.canvas.plot_spec(data['wave'], data['flux'], data['error'], filename)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MultiSpecViewer()
    viewer.show()
    sys.exit(app.exec_())
