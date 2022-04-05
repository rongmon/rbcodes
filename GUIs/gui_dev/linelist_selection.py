import sys
import pandas as pd
from numpy import floor, log10, isnan

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox, QLineEdit, QPushButton, QCheckBox
from PyQt5.QtCore import pyqtSignal

from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5.QtGui import QDoubleValidator

from IGM.rb_setline import read_line_list

class LineListWidget(QWidget):
	# Linelist constant
	# only need to update this one
	LINELISTS = ['NONE', 'Eiger_Strong', 'LBG', 'Gal', 'LLS', 'LLS Small', 'DLA', 'atom']
	# check function _get_linelist_df if any error showed up

	# sending out data
	send_lineindex = pyqtSignal(int)
	send_linelist = pyqtSignal(object)
	send_more_linelist = pyqtSignal(object)
	send_more_linelist_z = pyqtSignal(object)
	send_data = pyqtSignal(object)
	send_gauss_num = pyqtSignal(int)
	send_message = pyqtSignal(str)
	send_z_returnPressed = pyqtSignal(float)

	def __init__(self):
		super().__init__()

		#self.linelist_name = self.LINELISTS
		self.linelist = []
		self.filename = ''
		self.filenames = []
		self.newz = []

		glayout = QVBoxLayout()

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
		self.l_lln.addItems(self.LINELISTS)
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
		self.l_combobox.currentIndexChanged.connect(self._index_changed)
		#self.l_combobox.currentTextChanged.connect(self._text_changed)

		self.gauss_num = QComboBox()
		self.gauss_num.setFixedWidth(50)
		self.gauss_num.addItems(['1', '2', '3'])
		self.gauss_num.setCurrentIndex(0)
		self.gauss_num.activated.connect(self._on_gauss_num_activated)
		layout.addWidget(self.gauss_num, 1,2)

		# 3 textedit box
		self.onlyFloat = QDoubleValidator()
		self.estZ = QLineEdit()
		self.estZ.setPlaceholderText('Guess redshift')
		self.estZ.setMaximumWidth(100)
		self.estZ.setValidator(self.onlyFloat)
		self.estZ.returnPressed.connect(self._on_z_return_pressed)
		self.estZstd = QLineEdit()
		self.estZstd.setPlaceholderText('z Error')
		self.estZstd.setMaximumWidth(100)
		self.estZstd.setValidator(self.onlyFloat)
		self.conf = QLineEdit()
		self.conf.setPlaceholderText('[0, 1.]')
		self.conf.setMaximumWidth(150)
		self.conf_onlyFloat = QDoubleValidator(bottom=0., 
												top=1., 
												decimals=3,
												notation=QDoubleValidator.StandardNotation)
		self.conf.setValidator(self.conf_onlyFloat)

		self.flag = QLineEdit()
		self.flag.setPlaceholderText('Additional Info?')
		button = QPushButton('Add to Table below')
		button.clicked.connect(self._on_button_clicked)
		layout.addWidget(self.estZ, 1,3)
		layout.addWidget(self.estZstd, 1, 4)
		layout.addWidget(self.conf, 1,5)
		layout.addWidget(self.flag, 1,6)
		layout.addWidget(button, 1,7)

		# Secondary linelists
		l_checkbox = QCheckBox('Add More Linelists to Examine..')
		l_checkbox.stateChanged.connect(self._intialize_more_linelist)
		l_hlayout = QHBoxLayout()

		num_llists = 6
		self.llists_2 = []
		for i in range(num_llists):
			self.llists_2.append(self.add_linelists())
			#indices for llists_2:
			# 0-layout, 1-linelist combobox, 2-z lineedit
			l_hlayout.addLayout(self.llists_2[i][0])


		glayout.addLayout(layout)
		glayout.addWidget(l_checkbox)
		glayout.addLayout(l_hlayout)
		
		glayout.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
		self.setLayout(glayout)

	def add_linelists(self):
		ll_layout = QGridLayout()
		l_combox = QComboBox()

		l_combox.setFixedWidth(80)
		z_edit = QLineEdit()
		z_edit.setPlaceholderText('Guess z')
		z_edit.setReadOnly(True)
		z_edit.setMaximumWidth(60)
		#ll_layout.addWidget(QLabel('Linelist'), 0,0)
		#ll_layout.addWidget(QLabel('Guess z'), 0,1)
		ll_layout.addWidget(l_combox, 0,0)
		ll_layout.addWidget(z_edit, 0,1)
		ll_layout.setAlignment(QtCore.Qt.AlignLeft)

		return ll_layout, l_combox, z_edit

	def _intialize_more_linelist(self, s):
		if s == QtCore.Qt.Checked:
			# initialize more linelists for plotting
			colors = ['#A52A2A', '#FF7F50', '#40E0D0', '#DAA520', '#008000', '#4B0082']
			for i in range(len(self.llists_2)):
				self.llists_2[i][1].addItems(self.LINELISTS)
				t_color = 'QComboBox {color:' + colors[i] + '}'
				self.llists_2[i][1].setStyleSheet(t_color)
				# 2 parameters need to be passed: 1.selected linelist and 2.index of linelist widget changed
				self.llists_2[i][1].currentTextChanged.connect(lambda s, idx=i: self._addtional_linelist(idx, s))
				
				self.llists_2[i][2].setReadOnly(False)
				# only the index of linelist widget passed
				self.llists_2[i][2].returnPressed.connect(lambda idx=i: self._guess_z_return_pressed(idx))


		else:
			# grey out all things
			for i in range(len(self.llists_2)):
				self.llists_2[i][1].clear()
				self.llists_2[i][2].clear()
				self.llists_2[i][2].setReadOnly(True)
			

	def _on_z_return_pressed(self):
		self.send_z_returnPressed.emit(float(self.estZ.text()))

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
		if s in 'NONE':
			self.send_linelist.emit(s)
			self.l_combobox.clear()
			self.l_combobox.addItem('NONE')
			self.l_combobox.setCurrentIndex(0)
		else:
			llist = self._get_linelist_df(s)
			self.linelist = llist

			self.l_combobox.addItems(['ALL'] + self.linelist['name'].tolist())
			self.send_linelist.emit(self.linelist)
			self.l_combobox.setCurrentIndex(1)

	def _on_estZ_changed(self, newz):
		show_sigfig = 5
		self.newz = newz
		self.estZ.setText(str(self.round_to_sigfig(newz[0], show_sigfig)))
		self.estZstd.setText(str(self.round_to_sigfig(newz[1], show_sigfig)))

	def _on_sent_filename(self, sent_filename):
		self.filename = sent_filename

	def _on_sent_filenames(self, sent_filenames):
		self.filenames = sent_filenames

	def _on_sent_fitsobj(self, sent_fitsobj):
		self.fitsobj = sent_fitsobj

	def _on_button_clicked(self, sfilename):
		if len(self.estZ.text().strip()) < 1:
			self.estZ.setText('0')
		if len(self.estZstd.text().strip()) < 1:
			self.estZstd.setText('0')
		if len(self.conf.text().strip()) < 1:
			self.conf.setText('0')
		if len(self.flag.text().strip()) < 1:
			self.flag.setText('No comments')
		data = {'Name': self.filename,
				'z': self.newz[0], #float(self.estZ.text()),
				'z_err': self.newz[1], #float(self.estZstd.text()),
				'Confidence': float(self.conf.text()),
				'Linelist': self.l_lln.currentText(),
				'Flag': self.flag.text()}
		if self.fitsobj.ra is not None:
			data.update({'RA': self.fitsobj.ra,
							'DEC': self.fitsobj.dec})

		self.send_data.emit(data)

	def _on_sent_dictdata(self, sent_dict):
		#print(self.filename)
		#print(sent_dict)
		if len(sent_dict) > 0:
			if not isnan(sent_dict['z']):
				# extract z_estimated from z column
				if len(self.newz) < 2:
					self.newz.append(float(sent_dict['z']))
					self.newz.append(float(sent_dict['z_err']))
				else:
					self.newz[0] = float(sent_dict['z'])
					self.newz[1] = float(sent_dict['z_err'])

				show_sigfig = 5
				self.estZ.setText(str(self.round_to_sigfig(self.newz[0], show_sigfig)))
				self.estZstd.setText(str(self.round_to_sigfig(self.newz[1], show_sigfig)))
				self.send_message.emit('Found estimated z in table!')

			elif not isnan(sent_dict['z_guess']):
				# extract z_estimated from z_guess column
				if len(self.newz) < 2:
					self.newz.append(float(sent_dict['z_guess']))
					self.newz.append(0.)
					
				else:
					self.newz[0] = float(sent_dict['z_guess'])
					self.newz[1] = 0.				
				self.estZ.setText(str(self.newz[0]))
				self.estZstd.setText(str(self.newz[1]))
			
			self.conf.setText(str(sent_dict['Confidence']))
			self.flag.setText(sent_dict['Flag'])
			self.l_lln.setCurrentText(sent_dict['Linelist'])
			self.send_message.emit("Can't find z in table! Use z_guess now..")
		else:
			self.estZ.clear()
			self.estZstd.clear()
			self.newz = [0., 0.] # reset est_z and est_z_std back to zero

			self.conf.clear()
			self.flag.clear()
			self.l_lln.setCurrentIndex(0)

	def _on_gauss_num_activated(self):
		self.send_gauss_num.emit(int(self.gauss_num.currentText()))

	def round_to_sigfig(self, num=0., sigfig=1):
		return round(num, sigfig - int(floor(log10(abs(num)))) - 1)

	def _addtional_linelist(self, i, s):
		# i = index of the widget passing linelist
		# s = name of the linelist
		llist = pd.DataFrame(columns=['wave', 'name'])
		if s in 'NONE':
			self.send_more_linelist.emit({i:s})
		else:
			llist = self._get_linelist_df(s)

			self.send_more_linelist.emit({i:llist})

	def _guess_z_return_pressed(self, i):		
		llist = self._get_linelist_df(self.llists_2[i][1].currentText())
		z_guess = float(self.llists_2[i][2].text())
		self.send_more_linelist_z.emit([{i:llist}, z_guess])

		

	def _get_linelist_df(self, linelist_name):
		llist = pd.DataFrame(columns=['wave', 'name'])
		tmp = read_line_list(linelist_name)

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

		return llist