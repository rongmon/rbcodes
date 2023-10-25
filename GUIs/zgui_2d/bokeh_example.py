import sys
import os
from bokeh.plotting import figure, output_file, save
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl

from astropy.io import fits


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Create a QWebEngineView widget to display the Plotly plot
        view = QWebEngineView()
        view.setGeometry(100, 100, 800, 600)


        # create a simple Plotly image plot
        fitsfile = fits.open('./eiger-mock-data/spec2d_coadd_QSO_J0100_sID010242.fits')
        image = fitsfile['SCI'].data

        # get the absolute path of the current directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # define the full path to the HTML file
        html_file_path = os.path.join(script_dir, 'bokeh_plot.html')
        output_file(html_file_path, mode='cdn')
        
        # create a simple Bokeh image plot
        
        '''
        p = figure(title='Bokeh Image Plot: spec2d_coadd_QSO_J0100_sID010242.fits',
                   width=800, height=400)
        img = p.image(image=[image], # flip the image
                      x=0, y=0, # set the origin at the bottom left corner
                      dw=800, dh=400, # set the width and height of the image
                      palette="Spectral11", # set the color palette
                      )
        '''
        p = figure(title="Bokeh Plot", width=800, height=400)
        r = p.line([1, 2, 3], [4, 1, 2], line_width=2)
        save(p, html_file_path)

        # Load the HTML file in the QWebEngineView
        view.setUrl(QUrl.fromLocalFile(html_file_path))
        
        self.setCentralWidget(view)


# function to delete the HTML file when the application is closed
def delete_html_file():
    if os.path.exists('bokeh_plot.html'):
        os.remove('bokeh_plot.html')

app = QApplication(sys.argv)
window = MainWindow()
window.show()

# delete the HTML file when the application is closed
app.aboutToQuit.connect(delete_html_file)
sys.exit(app.exec_())