from PyQt5.QtWidgets import QWidget, QDialog, QLabel, QVBoxLayout, QHBoxLayout, QGridLayout
from PyQt5.QtCore import QSize, Qt
# This contains the complicated user manual
class UserManualDialog(QWidget):
	def __init__(self, method = 0):
		super().__init__()
		self.setWindowTitle('GUI User Manual')

		if method == 0:
			content = QLabel('''MAIN GUI HELP:
				--------------Keyboard Events----------        
				'r':    Resets/clears the axes and replots the spectra
				'b':    Sets bottom y limit (ymin)
				't':    Sets top y limit (ymax)
				'x':    Sets left x limit (xmin)
				'X':    Sets right x limit (xmax)
				'Y':    Set y limits with precise user-input values
				'W':    Set x limits with precise user-input values
				'S':    Smoothes the spectra
				'U':    Unsmooth spectra
				']':    Shifts canvas to the right
				'[':    Shifts canvas to the left

				'G':    Fit a Gaussian profile with 3 user-selected points
				'D':    Delete previous selected points''')
			self.layout = QVBoxLayout()
			self.layout.addWidget(content)

		if method == 1:
			self.layout = QVBoxLayout()

			title = QLabel('MAIN GUI USER MANUAL')
			self.layout.addWidget(title, Qt.AlignCenter)
			self.layout.addWidget(QLabel('Key Presses'))
			self.layout.addWidget(QLabel('For 1D spectrum plot,'))
			self._create_event_doc('r', 'Resets/clears the axes and replots the spectra')
			self._create_event_doc('t', 'Sets bottom y limit (ymin)')
			self._create_event_doc('b', 'Sets top y limit (ymax)')
			self._create_event_doc('x', 'Sets left x limit (xmin)')
			self._create_event_doc('X', 'Sets right x limit (xmax)')
			self._create_event_doc('Y', 'Set y limits with precise user-input values')
			self._create_event_doc('W', 'Set x limits with precise user-input values')
			self._create_event_doc('S', 'Smoothes the spectra')
			self._create_event_doc('U', 'Unsmooth spectra')
			self._create_event_doc(']', 'Shifts canvas to the right')
			self._create_event_doc('[', 'Shifts canvas to the left')
			self._create_event_doc('G', 'Fit a Gaussian profile with 3 user-selected points')
			self._create_event_doc('D', 'Delete previous selected points')
			self.layout.addWidget(QLabel('For 2D spectrum plot,'))
			self._create_event_doc('C', 'Plots the cumulative data along wavelength axis')
			self._create_event_doc('P', 'Plots the histogram of entire 2D spectrum range')
			self.layout.addWidget(QLabel('-'*100), Qt.AlignCenter)
			self.layout.addWidget(QLabel('Mouse Clicks'))
			self._create_event_doc('Right', 'Bring up the ion list')
			self.layout.addWidget(QLabel('-'*100), Qt.AlignCenter)
			self.layout.addWidget(QLabel('Steps:'))
			self._create_step_doc('1.', 'Load a FITS-format spectrum')
			self._create_step_doc('2.', 'Select a line-list to work on')
			self._create_step_doc('3.', 'Zoom to region of interest and fit a Gaussian')
			self._create_step_doc('4.', 'Right click to select a potential line')
			self._create_step_doc('5.', 'Estimate error and Write your confidence level and additional flags')
			self._create_step_doc('6.', 'Save to the table below')



		self.setLayout(self.layout)

	def _create_event_doc(self, e_name='', e_doc=''):
		name = QLabel(' Key ' + e_name)
		name.setFixedWidth(80)
		doc = QLabel(e_doc)

		l = QHBoxLayout()
		l.addWidget(name, Qt.AlignLeft)
		l.addWidget(doc, Qt.AlignLeft)
		self.layout.addLayout(l)

	def _create_step_doc(self, e_name='', e_doc=''):
		name = QLabel(' ' + e_name)
		name.setFixedWidth(80)
		doc = QLabel(e_doc)
		l = QHBoxLayout()
		l.addWidget(name, Qt.AlignLeft)
		l.addWidget(doc, Qt.AlignLeft)
		self.layout.addLayout(l)