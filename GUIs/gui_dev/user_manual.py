from PyQt5.QtWidgets import QDialog, QLabel, QVBoxLayout
from PyQt5.QtCore import QSize
# This contains the complicated user manual
class UserManualDialog(QDialog):
	def __init__(self):
		super().__init__()

		self.setWindowTitle('GUI User Manual')
		content = QLabel(
		'''MAIN GUI HELP:
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
        'x':    Sets left x limit (xmin)
        'X':    Sets right x limit (xmax)
        ']':    Shifts canvas to the right
        '[':    Shifts canvas to the left

        'G':    Fit a Gaussian profile with 3 user-selected points
        'D':    Delete previous selected points
        ''')
		self.layout = QVBoxLayout()
		self.layout.addWidget(content)
		self.setLayout(self.layout)