import sys
import pandas as pd

from PyQt5.QtWidgets import QWidget, QTextEdit, QVBoxLayout
from PyQt5.QtCore import pyqtSignal

from PyQt5 import QtCore
from PyQt5 import QtGui

class MessageBox(QWidget):

	def __init__(self):
		super().__init__()

		self.message = ''

		self.te = QTextEdit()
		self.te.setPlaceholderText('This meesage box will display important messages.')
		self.te.setReadOnly(True)

		self.setFixedHeight(200)
		self.setFixedWidth(250)

		layout = QVBoxLayout()
		layout.addWidget(self.te)
		layout.setAlignment(QtCore.Qt.AlignCenter)

		self.setLayout(layout)

	def on_sent_message(self, sent_message, hexColor='#000000'):
		#Receiving the signal
		#self.message += sent_message
		prefix = f'<span style=\"color:{hexColor};\" >'
		suffix = f'</span>'
		message = prefix + f'{sent_message}' + suffix
		self.te.setText(message)