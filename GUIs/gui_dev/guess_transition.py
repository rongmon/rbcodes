import sys
import numpy as np
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QComboBox, QLineEdit, QListWidget
from PyQt5.QtCore import pyqtSignal

class GuessTransition(QWidget):
	send_z_cal = pyqtSignal(float)

	def __init__(self, linelist, wavelength):
		super().__init__()
		self.linelist = linelist
		self.wave = wavelength
		#self.z_cal = 0.
		self.resize(200, 1000)
		self.layout = QVBoxLayout()

		self.transitions = QListWidget()

		self.transitions.addItems(self.linelist['name'].to_numpy())
		#print(self.linelist['name'].to_numpy())
		#print(type(self.linelist['name']))
		self.layout.addWidget(self.transitions)
		self.setLayout(self.layout)

		self.transitions.itemClicked.connect(self._select_ion)

	def _select_ion(self):
		ion_num = self.transitions.currentRow()
		ion_wrest = self.linelist['wave'][ion_num]
		z_cal = np.round((self.wave - ion_wrest)/ion_wrest, 4)
		self.send_z_cal.emit(z_cal)
		#self.z_cal = z_cal
