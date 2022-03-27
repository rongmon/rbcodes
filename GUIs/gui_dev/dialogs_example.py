import sys

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QDialog
from spec_fit_gauss2d import Gaussfit_2d
from linetools.spectra.xspectrum1d import XSpectrum1D  
import numpy as np

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('My simple dialog app')

        button = QPushButton('Press me for a dialog!')
        button.clicked.connect(self.button_clicked)

        self.setCentralWidget(button)


        #Add your flux extraction code in here
        
        sp=XSpectrum1D.from_file('../../example-data/test.fits')

        wave=sp.wavelength.value
        flux=sp.flux.value
        error=sp.sig.value
        q=np.where((wave>1330) & (wave <1340))
        self.wave=wave[q]
        self.flux=flux[q]
        self.error=error[q]
        

    def button_clicked(self, s):
        dlg = Gaussfit_2d(self.wave, self.flux, self.error)
        dlg.setWindowTitle('HELLO!!')
        dlg.exec_()
        #print('click', s)

app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec_()