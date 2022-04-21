import sys
import os
import pandas as pd
from numpy import floor, log10, isnan, nan

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox, QLineEdit, QPushButton, QCheckBox
from PyQt5.QtCore import pyqtSignal

from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5.QtGui import QDoubleValidator

from IGM.rb_setline import read_line_list

LINELIST_DIR = os.path.dirname(os.path.abspath(__file__))

# This widget includes the primary linelist, redshift estimation
# and secondary linelists.
# Since this class has complicated structure in terms of widget composition,
# only import useful modules/functions to avoid slow-down running time
class LineListWidget(QWidget):
	# Linelist constant
	# only need to update this one
	LINELISTS = ['NONE']
	with open(LINELIST_DIR+'/gui_linelists.ascii') as f:
		next(f)
		for line in f:
			LINELISTS.append(line.strip())

	# check function _get_linelist_df if any error showed up

	# Exporting signals
	send_lineindex = pyqtSignal(int)
	send_linelist = pyqtSignal(object)
	send_more_linelist = pyqtSignal(object)
	send_more_linelist_z = pyqtSignal(object)
	send_data = pyqtSignal(object)
	send_gauss_num = pyqtSignal(int)
	send_message = pyqtSignal(str)
	send_z_returnPressed = pyqtSignal(float)
	send_linelists2multiG = pyqtSignal(list)

	# initialization - widget layout
	def __init__(self):
		super().__init__()

		#internal values
		self.linelist = []
		self.filename = ''
		self.filenames = []
		self.newz = []

		# Main(Grand) widget layout
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

		# linelist combobox
		self.l_lln = QComboBox()
		self.l_lln.setFixedWidth(120)
		self.l_lln.addItems(self.LINELISTS)
		layout.addWidget(self.l_lln, 1, 0)
		# selecting a linelist in linelists box triggers another action
		self.l_lln.currentTextChanged.connect(self._linelist_changed)

		# Ion Names in this selected line-list
		# Note: 'ALL' is in index 0
		self.l_combobox = QComboBox()
		self.l_combobox.setFixedWidth(150)
		layout.addWidget(self.l_combobox, 1, 1)
		self.l_combobox.addItem('NONE')
		self.l_combobox.setCurrentIndex(0)
		# selecting an ion in a linelist triggers another action
		self.l_combobox.currentIndexChanged.connect(self._index_changed)
		#self.l_combobox.currentTextChanged.connect(self._text_changed)

		# Number of Gaussian specified
		self.gauss_num = QComboBox()
		self.gauss_num.setFixedWidth(50)
		self.gauss_num.addItems(['0', '1', '2', '3'])
		self.gauss_num.setCurrentIndex(1)
		# selecting a number of the box triggers another action
		self.gauss_num.activated.connect(self._on_gauss_num_activated)
		layout.addWidget(self.gauss_num, 1,2)

		# User-input textboxes
		# this validator only allows user to type numbers
		self.onlyFloat = QDoubleValidator()
		# Estimated redshift
		self.estZ = QLineEdit()
		self.estZ.setPlaceholderText('Guess redshift')
		self.estZ.setMaximumWidth(100)
		self.estZ.setValidator(self.onlyFloat)
		# pressing return button triggers another action
		self.estZ.returnPressed.connect(self._on_z_return_pressed)
		# Errors in estimated redshift
		self.estZstd = QLineEdit()
		self.estZstd.setPlaceholderText('z Error')
		self.estZstd.setMaximumWidth(100)
		self.estZstd.setValidator(self.onlyFloat)
		self.estZstd.setReadOnly(True)
		# Confidence level
		self.conf = QLineEdit()
		self.conf.setPlaceholderText('[0, 1.]')
		self.conf.setMaximumWidth(150)
		self.conf_onlyFloat = QDoubleValidator(bottom=0., 
												top=1., 
												decimals=3,
												notation=QDoubleValidator.StandardNotation)
		self.conf.setValidator(self.conf_onlyFloat)
		# Flags/Comments on this
		self.flag = QLineEdit()
		self.flag.setPlaceholderText('Additional Info?')

		# Button to save current estimation to table/database
		button = QPushButton('Add to Table below')
		# clicking button triggers another action
		button.clicked.connect(self._on_button_clicked)


		layout.addWidget(self.estZ, 1,3)
		layout.addWidget(self.estZstd, 1, 4)
		layout.addWidget(self.conf, 1,5)
		layout.addWidget(self.flag, 1,6)
		layout.addWidget(button, 1,7)

		# Secondary linelists
		l_checkbox = QCheckBox('Add More Linelists to Examine..')
		# only triggers to set up more secondary linelist when checked
		l_checkbox.stateChanged.connect(self._intialize_more_linelist)
		l_hlayout = QHBoxLayout()

		# Number of secondary linelists needed
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
		self.setFixedHeight(150)

	# utility function to set up secondary linelist
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

	# intialize secondary linelist sytles
	def _intialize_more_linelist(self, s):
		# if the checkbox is checked, initialize secondary linelists
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
			# if checkbox is unchecked, grey out all things
			for i in range(len(self.llists_2)):
				self.llists_2[i][1].clear()
				self.llists_2[i][2].clear()
				self.llists_2[i][2].setReadOnly(True)
			
	# action to press return button on estZ LineEdit
	def _on_z_return_pressed(self):
		self.send_z_returnPressed.emit(float(self.estZ.text()))
		if self.gauss_num.currentIndex() < 1:
			self.estZstd.setText('nan')
			# this will not changed the default self.newz values

	# importing signal(linelist name) to slot
	def on_linelist_name_slot(self, sent_linelist_name):
		self.l_lln.setText(sent_linelist_name)
	def on_linelist_slot(self, sent_linelist):
		self.linelist = sent_linelist
		#self.l_combobox.addItems(['NONE', 'ALL'] + self.linelist['name'].tolist())
		#print(self.linelist)

	# ion line combobox events
	def _index_changed(self, i): # i is an int
		self.send_lineindex.emit(i)

	def _text_changed(self, s): # s is a str
		tmp_df = self.linelist.set_index('name')

	# Display all available ions in a selected linelist
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

	# display estimated redshifts and error
	# with specified significant figures
	def _on_estZ_changed(self, newz):
		show_sigfig = 5
		self.newz = newz
		if not isnan(float(self.newz[0])):
			self.estZ.setText(str(self.round_to_sigfig(newz[0], show_sigfig)))
			if self.gauss_num.currentIndex() > 0:
				# Except for manually guessing z (i.e., gauss_num=0),
				# sending out estZ to plot automatically without return pressed
				self.send_z_returnPressed.emit(self.newz[0])
		else:
			self.estZ.setText('nan')
		if not isnan(float(self.newz[1])):
			self.estZstd.setText(str(self.round_to_sigfig(newz[1], show_sigfig)))
		else:
			self.estZstd.setText('nan')

	# importing signal(filename) to slot
	def _on_sent_filename(self, sent_filename):
		self.filename = sent_filename

	# importing signal(filenames) to slot
	def _on_sent_filenames(self, sent_filenames):
		self.filenames = sent_filenames

	# importing signal(FitsObj) to slot
	def _on_sent_fitsobj(self, sent_fitsobj):
		self.fitsobj = sent_fitsobj
		if self.fitsobj.z_guess is not None:
			# replace estZ if z_guess is available
			self.estZ.setText(str(self.round_to_sigfig(self.fitsobj.z_guess, 3)))
			self.send_message.emit('Redshift posterior is found in the FITS file!')

	# action to "Add to Table below"
	# log available values to DataFrame
	def _on_button_clicked(self, sfilename):
		if len(self.estZ.text().strip()) < 1:
			self.estZ.setText('nan')
		if len(self.estZstd.text().strip()) < 1:
			self.estZstd.setText('nan')
		if len(self.conf.text().strip()) < 1:
			self.conf.setText('0')
		if len(self.flag.text().strip()) < 1:
			self.flag.setText('No comments')

		# prepare exporting data
		data = {'Name': self.filename,
				'z': self.newz[0], #float(self.estZ.text()),
				'z_err': self.newz[1], #float(self.estZstd.text()),
				'Confidence': float(self.conf.text()),
				'Linelist': self.l_lln.currentText(),
				'Flag': self.flag.text()}
		# add coordiantes if they are available
		if self.fitsobj.ra is not None:
			data.update({'RA': self.fitsobj.ra,
							'DEC': self.fitsobj.dec})
		# add z_guess if it is available
		if self.fitsobj.z_guess is not None:
			data.update({'z_guess': self.fitsobj.z_guess})

		# export data to DataFrame table
		self.send_data.emit(data)

	# importing data from table to slot
	def _on_sent_dictdata(self, sent_dict):
		#print(self.filename)
		#print(sent_dict)

		# if received data is non-empty
		# add data back to corresponding LineEdit widgets
		if len(sent_dict) > 0:
			#print(sent_dict['z'])
			if not isnan(float(sent_dict['z'])):
				# extract z_estimated from z column
				if len(self.newz) < 2:
					self.newz.append(float(sent_dict['z']))
					self.newz.append(float(sent_dict['z_err']))
				else:
					self.newz[0] = float(sent_dict['z'])
					self.newz[1] = float(sent_dict['z_err'])

				show_sigfig = 5
				self.estZ.setText(str(self.round_to_sigfig(self.newz[0], show_sigfig)))
				if not isnan(float(self.newz[1])):
					self.estZstd.setText(str(self.round_to_sigfig(self.newz[1], show_sigfig)))
				else:
					self.estZstd.setText('nan')
				self.send_message.emit('Found estimated z in table!')

			elif not isnan(float(sent_dict['z_guess'])):
				# extract z_estimated from z_guess column
				if len(self.newz) < 2:
					self.newz.append(float(sent_dict['z_guess']))
					self.newz.append(0.)
					
				else:
					self.newz[0] = float(sent_dict['z_guess'])
					self.newz[1] = nan				
				self.estZ.setText(str(self.newz[0]))
				self.estZstd.setText(str(self.newz[1]))
			
			self.conf.setText(str(sent_dict['Confidence']))
			self.flag.setText(sent_dict['Flag'])
			self.l_lln.setCurrentText(sent_dict['Linelist'])
			self.send_message.emit("Can't find z in table! Use z_guess now..")
		else:
			# sent_dict data is empy, reset everythin
			self.estZ.clear()
			self.estZstd.clear()
			self.newz = [nan, nan] # reset est_z and est_z_std back to nans

			self.conf.clear()
			self.flag.clear()
			self.l_lln.setCurrentIndex(0)

	# action to number of Gaussians selected
	def _on_gauss_num_activated(self):
		# exporting signals
		self.send_gauss_num.emit(int(self.gauss_num.currentText()))
		self.send_linelists2multiG.emit(self.LINELISTS)

	# utility function to round values to desired significant figures
	def round_to_sigfig(self, num=0., sigfig=1):
		if num is not None:
			return round(num, sigfig - int(floor(log10(abs(num)))) - 1)
		else:
			return None

	# load full content in all secondary linelist
	def _addtional_linelist(self, i, s):
		# i = index of the widget passing linelist
		# s = name of the linelist
		llist = pd.DataFrame(columns=['wave', 'name'])
		if s in 'NONE':
			self.send_more_linelist.emit({i:s})
		else:
			llist = self._get_linelist_df(s)

			self.send_more_linelist.emit({i:llist})

	# action to press return on secondary z_guess LineEdit widgets
	def _guess_z_return_pressed(self, i):		
		llist = self._get_linelist_df(self.llists_2[i][1].currentText())
		z_guess = float(self.llists_2[i][2].text())
		self.send_more_linelist_z.emit([{i:llist}, z_guess])

	# action to read linelist as dataframe
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