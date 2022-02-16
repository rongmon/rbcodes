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
		layout.addWidget(QLabel('z error'), 0, 4)
		layout.addWidget(QLabel('Confidence'), 0, 5)
		layout.addWidget(QLabel('Flag'), 0, 6)

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
		self.l_combobox.setFixedWidth(150)
		layout.addWidget(self.l_combobox, 1, 1)
		self.l_combobox.addItem('NONE')
		self.l_combobox.setCurrentIndex(0)
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
		self.estZstd = QLineEdit()
		self.estZstd.setPlaceholderText('z Error')
		self.estZstd.setMaximumWidth(100)
		self.conf = QLineEdit()
		self.conf.setPlaceholderText('[0, 1.]')
		self.conf.setMaximumWidth(150)
		self.flag = QLineEdit()
		self.flag.setPlaceholderText('Additional Info?')
		button = QPushButton('Add to Table below')
		button.clicked.connect(self._on_button_clicked)
		layout.addWidget(self.estZ, 1,3)
		layout.addWidget(self.estZstd, 1, 4)
		layout.addWidget(self.conf, 1,5)
		layout.addWidget(self.flag, 1,6)
		layout.addWidget(button, 1,7)

		

		
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

	# line combobox events
	def _index_changed(self, i): # i is an int
		self.send_lineindex.emit(i)

	def _text_changed(self, s): # s is a str
		tmp_df = self.linelist.set_index('name')

	def _linelist_changed(self, s):
		llist = pd.DataFrame(columns=['wave', 'name'])
		if s in 'NONE':
			pass
		else:
			tmp = read_line_list(s)

		#need a line to append wrest to name if it doesn't have one
		if any(map(str.isdigit, tmp[1]['ion'])):
			# if name column has wrest
			for li in tmp:
				newrow = {'wave': li['wrest'], 'name': li['ion']}
				llist = llist.append(newrow, ignore_index=True)
		else:
			# if name column doesn't have wrest, need to append
			for li in tmp:
				newrow = {'wave': li['wrest'], 'name': li['ion']+' '+str(round(li['wrest']))}
				llist = llist.append(newrow, ignore_index=True)
		self.linelist = llist

		self.l_combobox.addItems(['ALL'] + self.linelist['name'].tolist())
		self.send_linelist.emit(self.linelist)

	def _on_estZ_changed(self, newz):
		show_prec = 5
		self.estZ.setText(str(round(newz[0], show_prec)))
		self.estZstd.setText(str(round(newz[1], show_prec)))

	def _on_sent_filename(self, sent_filename):
		self.filename = sent_filename

	def _on_button_clicked(self, sfilename):
		if len(self.estZ.text().strip()) < 1:
			self.estZ.setText('0')
		if len(self.estZstd.text().strip()) < 1:
			self.estZstd.setText('0')
		if len(self.conf.text().strip()) < 1:
			self.conf.setText('0')
		data = {'Name': self.filename,
				'z': float(self.estZ.text()),
				'z_err': float(self.estZstd.text()),
				'Confidence': float(self.conf.text()),
				'Linelist': self.l_lln.currentText(),
				'Flag': self.flag.text()}
		self.send_data.emit(data)

	def _on_gauss_num_acticated(self):
		self.send_gauss_num.emit(int(self.gauss_num.currentText()))
