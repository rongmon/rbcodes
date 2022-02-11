import sys
import pandas as pd

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QLabel, QComboBox, QLineEdit, QPushButton
from PyQt5.QtCore import pyqtSignal

from PyQt5 import QtCore
from PyQt5 import QtGui

from IGM.rb_setline import read_line_list

class LineListWidget(QWidget):
	send_lineindex = pyqtSignal(int)
	send_linelist = pyqtSignal(object)
	send_data = pyqtSignal(object)
	send_gauss_num = pyqtSignal(int)

	def __init__(self):
		super().__init__()

		self.linelist_name = ''
		self.linelist = []
		self.filename = ''

		# Widget column names
		layout = QGridLayout()
		layout.addWidget(QLabel('LineList Name'), 0, 0)
		layout.addWidget(QLabel('Ion Name'), 0, 1)
		layout.addWidget(QLabel('#Gauss'), 0, 2)
		layout.addWidget(QLabel('Estimated z'), 0, 3)
		layout.addWidget(QLabel('Confidence'), 0, 4)
		layout.addWidget(QLabel('Flag'), 0, 5)

		# Selected line-list name
		#self.l_lln = QLineEdit()
		#self.l_lln.setReadOnly(True)
		#self.l_lln.setMaximumWidth(120)

		self.l_lln = QComboBox()
		self.l_lln.setFixedWidth(120)
		self.l_lln.addItems(['NONE', 'LBG', 'Gal', 'LLS', 'LLS Small', 'DLA', 'atom'])
		layout.addWidget(self.l_lln, 1, 0)
		self.l_lln.currentTextChanged.connect(self._linelist_changed)
		#menubar.send_filename.connect(self.on_linelist_name_slot)

		# Ion Names in this selected line-list
		# Note: 'ALL' is in index 0
		self.l_combobox = QComboBox()
		self.l_combobox.setFixedWidth(100)
		layout.addWidget(self.l_combobox, 1, 1)
		#menubar.send_linelist.connect(self.on_linelist_slot)
		self.l_combobox.currentIndexChanged.connect(self._index_changed)
		self.l_combobox.currentTextChanged.connect(self._text_changed)

		self.gauss_num = QComboBox()
		self.gauss_num.setFixedWidth(50)
		self.gauss_num.addItems(['1', '2'])
		self.gauss_num.setCurrentIndex(0)
		self.gauss_num.activated.connect(self._on_gauss_num_acticated)
		layout.addWidget(self.gauss_num, 1,2)

		# 3 textedit box
		self.estZ = QLineEdit()
		self.estZ.setPlaceholderText('Guess redshift')
		self.estZ.setMaximumWidth(100)
		self.conf = QLineEdit()
		self.conf.setPlaceholderText('[0, 1.]')
		self.conf.setMaximumWidth(150)
		self.flag = QLineEdit()
		self.flag.setPlaceholderText('Additional Info?')
		button = QPushButton('Add to Table below')
		button.clicked.connect(self._on_button_clicked)
		layout.addWidget(self.estZ, 1,3)
		layout.addWidget(self.conf, 1,4)
		layout.addWidget(self.flag, 1,5)
		layout.addWidget(button, 1,6)

		

		
		layout.setAlignment(QtCore.Qt.AlignLeft)
		self.setLayout(layout)

	def zprint(self):
		print(self.estZ.text())

	# data receiption
	def on_linelist_name_slot(self, sent_linelist_name):
		self.l_lln.setText(sent_linelist_name)
	def on_linelist_slot(self, sent_linelist):
		self.linelist = sent_linelist
		#self.l_combobox.addItems(['NONE', 'ALL'] + self.linelist['name'].tolist())
		#print(self.linelist)

	# combobox events
	def _index_changed(self, i): # i is an int
		#print(type(i))
		#print(i)
		#if i < 2:
		#	print(self.linelist)
		#else:
		#	print(self.linelist.loc[i-2])
		self.send_lineindex.emit(i)

	def _text_changed(self, s): # s is a str
		#print(type(s))
		#print(s)
		tmp_df = self.linelist.set_index('name')
		#if s not in ['NONE', 'ALL'] :
			#print(tmp_df.loc[s])

	def _linelist_changed(self, s):
		llist = pd.DataFrame(columns=['wave', 'name'])
		if s in 'NONE':
			pass
		else:
			tmp = read_line_list(s)
		for li in tmp:
			newrow = {'wave': li['wrest'], 'name': li['ion']}
			llist = llist.append(newrow, ignore_index=True)
		self.linelist = llist
		self.l_combobox.clear()
		self.l_combobox.addItems(['NONE', 'ALL'] + self.linelist['name'].tolist())
		#self.l_combobox.setCurrentIndex(0)
		#self.l_combobox.setCurrentIndex(1)
		#print(self.linelist)
		self.send_linelist.emit(self.linelist)

	def _on_estZ_changed(self, newz):
		self.estZ.setText(str(newz))
		print(newz)

	def _on_sent_filename(self, sent_filename):
		self.filename = sent_filename

	def _on_button_clicked(self, sfilename):
		data = {'Name': self.filename,
				'z': float(self.estZ.text()),
				'Confidence': float(self.conf.text()),
				'Linelist': self.l_lln.currentText(),
				'Flag': self.flag.text()}
		#print(data)
		self.send_data.emit(data)

	def _on_gauss_num_acticated(self):
		self.send_gauss_num.emit(int(self.gauss_num.currentText()))
