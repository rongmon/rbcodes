import sys
import numpy as np
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QComboBox, QLineEdit, QListWidget
from PyQt5.QtCore import pyqtSignal

class GuessTransition(QWidget):
	'''Pop-up Widget to select ions'''

	# Exporting signals
	send_z_cal = pyqtSignal(list)

	# Initialization - Widget main layout
	def __init__(self, linelist, wavelength, wave_std):
		super().__init__()

		# internal values
		self.linelist = linelist
		self.wave = wavelength
		self.std = wave_std
		#self.z_cal = 0.

		# widget layout
		self.resize(200, 1000)
		self.layout = QVBoxLayout()
		self.transitions = QListWidget()
		self.transitions.addItems(self.linelist['name'].to_numpy())
		#print(self.linelist['name'].to_numpy())
		#print(type(self.linelist['name']))
		self.layout.addWidget(self.transitions)
		self.setLayout(self.layout)

		# link clicking movement to more actions
		self.transitions.itemClicked.connect(self._select_ion)

	# define more actions after selecting an ion from widget
	def _select_ion(self):
		ion_num = self.transitions.currentRow()
		ion_wrest = self.linelist['wave'][ion_num]
		z_cal = (self.wave - ion_wrest)/ion_wrest
		if self.std is np.nan:
			z_std = np.nan
		else:
			z_std = self.std / self.wave * z_cal
		self.send_z_cal.emit([z_cal, z_std])
