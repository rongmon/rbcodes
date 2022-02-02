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
        'r':   resets/clears the axes and replots the spectra
        't':   Will restrict the ymax of the canvas to the users current mouse height
        'b':   Restricts ymin of the canvas to current mouse height
        's':   Smoothes the spectra
        'u':   Unsmooth spectra
        'x':   Sets left x limit (xmin)
        'X':   Sets right x limit (xmax)
        ']':   Shifts canvas to the right
        '[':   Shifts canvas to the left
        ''')
		self.layout = QVBoxLayout()
		self.layout.addWidget(content)
		self.setLayout(self.layout)