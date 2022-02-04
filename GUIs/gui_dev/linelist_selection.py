import sys
import pandas as pd

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QLabel, QComboBox, QLineEdit

from PyQt5 import QtCore
from PyQt5 import QtGui

class lineListWidget(QWidget):
	def __init__(self, menubar):
		super().__init__()

		self.linelist_name = ''
		self.linelist = []

		# Widget column names
		layout = QGridLayout()
		layout.addWidget(QLabel('LineList Name'), 0, 0)
		layout.addWidget(QLabel('Ion Name'), 0, 1)
		layout.addWidget(QLabel('Estimated Z'), 0, 2)
		layout.addWidget(QLabel('Confidence'), 0, 3)
		layout.addWidget(QLabel('Flag'), 0, 4)

		# Selected line-list name
		self.l_lln = QLineEdit()
		self.l_lln.setReadOnly(True)
		layout.addWidget(self.l_lln, 1, 0)
		menubar.send_filename.connect(self.on_linelist_name_slot)

		# Ion Names in this selected line-list
		# Note: 'ALL' is in index 0
		self.l_combobox = QComboBox()
		layout.addWidget(self.l_combobox, 1, 1)
		menubar.send_linelist.connect(self.on_linelist_slot)
		self.l_combobox.currentIndexChanged.connect(self._index_changed)
		self.l_combobox.currentTextChanged.connect(self._text_changed)

		# 3 textedit box
		self.estZ = QLineEdit()
		self.estZ.setPlaceholderText('Guess redshift')
		self.conf = QLineEdit()
		self.conf.setPlaceholderText('% Confident?')
		self.flag = QLineEdit()
		self.flag.setPlaceholderText('Additional Info?')
		layout.addWidget(self.estZ)
		layout.addWidget(self.conf)
		layout.addWidget(self.flag)
		

		
		layout.setAlignment(QtCore.Qt.AlignLeft)
		self.setLayout(layout)

	# data receiption
	def on_linelist_name_slot(self, sent_linelist_name):
		self.l_lln.setText(sent_linelist_name)
	def on_linelist_slot(self, sent_linelist):
		self.linelist = sent_linelist
		self.l_combobox.addItems(['NONE', 'ALL'] + self.linelist['name'].tolist())
		#print(self.linelist)

	# combobox events
	def _index_changed(self, i): # i is an int
		print(type(i))
		print(i)
		if i < 2:
			print(self.linelist)
		else:
			print(self.linelist.loc[i-2])
	def _text_changed(self, s): # s is a str
		print(type(s))
		print(s)
		tmp_df = self.linelist.set_index('name')
		if s not in ['NONE', 'ALL'] :
			print(tmp_df.loc[s])
