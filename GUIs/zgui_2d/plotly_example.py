import sys
import os
import plotly.graph_objs as go
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl

from astropy.io import fits


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Create a QWebEngineView widget to display the Plotly plot
        view = QWebEngineView()
        view.setGeometry(100, 100, 2000, 600)

        # Create a simple Plotly scatter plot
        #data = [go.Scatter(x=[1, 2, 3], y=[4, 1, 2], mode='markers', text=["Point 1", "Point 2", "Point 3"])]
        #go_layout = go.Layout(title='Plotly Scatter Plot')

        # create a simple Plotly image plot
        fitsfile = fits.open('./eiger-mock-data/spec2d_coadd_QSO_J0100_sID010242.fits')
        image = fitsfile['SCI'].data
        #print(image)
        data = [go.Heatmap(z=image)]
        go_layout = go.Layout(title='Plotly Image Plot: spec2d_coadd_QSO_J0100_sID010242.fits')

        fig = go.Figure(data=data, layout=go_layout)

        # Set hoverinfo to display point values
        fig.update_traces(hoverinfo='text+x+y')

        # get the absolute path of the current directory
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # define the full path to the HTML file
        html_file_path = os.path.join(script_dir, 'plotly_plot.html')

        # Export the plot as an HTML file
        fig.write_html(html_file_path)

        # Load the HTML file in the QWebEngineView
        view.setUrl(QUrl.fromLocalFile(html_file_path))
        



        self.setCentralWidget(view)


# function to delete the HTML file when the application is closed
def delete_html_file():
    if os.path.exists('plotly_plot.html'):
        os.remove('plotly_plot.html')

app = QApplication(sys.argv)
window = MainWindow()
window.show()

# delete the HTML file when the application is closed
app.aboutToQuit.connect(delete_html_file)
sys.exit(app.exec_())