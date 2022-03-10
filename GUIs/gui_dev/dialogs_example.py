import sys

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QDialog
from spec_fit_gauss2d import Gaussfit_2d

class MainWindow(QMainWindow):
	def __init__(self):
		super().__init__()

		self.setWindowTitle('My simple dialog app')

		button = QPushButton('Press me for a dialog!')
		button.clicked.connect(self.button_clicked)

		self.setCentralWidget(button)

		#Add your flux extraction code in here
		self.wave = [1,2,3] 
		self.flux = [20, 40, 20]
		self.error = [5,4,2]


	def button_clicked(self, s):
		dlg = Gaussfit_2d(self.wave, self.flux, self.error)
		dlg.setWindowTitle('HELLO!!')
		dlg.exec_()
		#print('click', s)

app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec_()