from PyQt5.QtWidgets import QDialog, QLabel, QVBoxLayout, QGridLayout
from PyQt5.QtCore import QSize
# This contains the complicated user manual
class UserManualDialog(QDialog):
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
			self.layout = QGridLayout()
			title = QLabel('MAIN GUI USER MANUAL')
			self.layout.addWidget(title, 0,0)
			self.layout.addWidget(QLabel('Key Presses'), 1,0)
			self._create_event_doc('r', 'Resets/clears the axes and replots the spectra', 2)
			self._create_event_doc('t', 'Sets bottom y limit (ymin)', 3)
			self._create_event_doc('b', 'Sets top y limit (ymax)', 4)
			self._create_event_doc('x', 'Sets left x limit (xmin)', 5)
			self._create_event_doc('X', 'Sets right x limit (xmax)', 6)
			self._create_event_doc('Y', 'Set y limits with precise user-input values', 7)
			self._create_event_doc('W', 'Set x limits with precise user-input values', 8)
			self._create_event_doc('S', 'Smoothes the spectra', 9)
			self._create_event_doc('U', 'Unsmooth spectra', 10)
			self._create_event_doc(']', 'Shifts canvas to the right', 11)
			self._create_event_doc('[', 'Shifts canvas to the left', 12)
			self._create_event_doc('G', 'Fit a Gaussian profile with 3 user-selected points', 13)
			self._create_event_doc('D', 'Delete previous selected points', 14)
			self.layout.addWidget(QLabel('------------'), 15,0)
			self.layout.addWidget(QLabel('Mouse Clicks'), 16,0)
			self._create_event_doc('Right', 'Bring up the ion list', 17)


		self.setLayout(self.layout)

	def _create_event_doc(self, e_name='', e_doc='', row_num=0):
		name = QLabel('Key ' + e_name)
		doc = QLabel(e_doc)
		self.layout.addWidget(name, row_num, 0)
		self.layout.addWidget(doc, row_num, 1)