import sys
import os
import pandas as pd
import numpy as np

from astropy.io import fits, ascii
from astropy.table import Table

from PyQt5.QtCore import Qt, QSize, QUrl, pyqtSignal
from PyQt5.QtWidgets import QAction, QToolBar, QStatusBar, QMenuBar, QFileDialog, QComboBox, QLineEdit, QLabel
from PyQt5.QtGui import QKeySequence, QDesktopServices

from user_manual import UserManualDialog
from gui_io import LoadSpec
from utils import FitsObj

WORKING_DIR = os.path.abspath(os.getcwd()) + './example-data'
LINELIST_DIR = os.path.dirname(os.path.abspath(__file__)) + '/lines'

class Custom_ToolBar(QToolBar):
	send_fitsobj = pyqtSignal(object)
	send_filename = pyqtSignal(str)
	send_filenames = pyqtSignal(list)
	#Our personalized custom Toolbar
	def __init__(self, mainWindow):
		super().__init__()
		self.mW = mainWindow
		self.fitsobj = mainWindow.fitsobj
		self.filepaths = []
		self.filenames = []
		self.filename = ''
		self.scale2d = False
		self.scale = 0

		self.setWindowTitle('Customizable TooBar')
		#self.setFixedSize(QSize(200, 50))

		btn_loadtxt = self._create_button('Read All', 'Read all fits files from a TXT file. Make sure you run "ls *.fits > filenames.txt" before clicking this!')
		self.addAction(btn_loadtxt)
		btn_loadtxt.triggered.connect(self._load_from_txt)

		btn_loadf = self._create_button('Read FITS', 'Read a fits file and plot the spectrum')
		self.addAction(btn_loadf)
		#btn_loadf.triggered.connect(self._load_spec)
		btn_loadf.triggered.connect(self._load_more_specs)

		btn_savef = self._create_button('Save FITS', 'Save the fits file to unified format')
		self.addAction(btn_savef)
		btn_savef.triggered.connect(self._save_spec)
		self.addSeparator()

		btn_help = self._create_button('Help', 'Open the user manual for further help')
		self.addAction(btn_help)
		btn_help.triggered.connect(self._open_user_manual)	
		self.addSeparator()

		# file dropbox
		self.f_combobox = QComboBox()
		self.f_combobox.setFixedWidth(200)
		self.f_combobox.addItem('No FITS File')
		self.f_combobox.setCurrentIndex(0)
		self.f_combobox.currentIndexChanged.connect(self._read_selected_fits)
		self.addWidget(self.f_combobox)
		self.addSeparator()

		# Scaling label
		s_label = QLabel('Scale:')
		self.addWidget(s_label)
		# scaling dropbox
		self.s_combobox = QComboBox()
		self.s_combobox.setFixedWidth(100)
		self.addWidget(self.s_combobox)

		# normalization dropbox
		self.n_combobox = QComboBox()
		self.n_combobox.setFixedWidth(100)
		self.addWidget(self.n_combobox)

		# normalization range
		self.min_range = QLineEdit()
		self.min_range.setFixedWidth(100)
		self.addWidget(self.min_range)
		self.min_range.setPlaceholderText('Min')
		self.min_range.returnPressed.connect(self._return_pressed)
		self.min_range.setReadOnly(True)
		self.max_range = QLineEdit()
		self.max_range.setFixedWidth(100)
		self.max_range.setPlaceholderText('Max')
		self.max_range.returnPressed.connect(self._return_pressed)
		self.max_range.setReadOnly(True)
		self.addWidget(self.max_range)



	def _create_button(self, buttonName='', buttonTip=''):
		#Create buttons wrapper
		btn = QAction(buttonName, self.mW)
		btn.setStatusTip(buttonTip)

		return btn

	def _open_user_manual(self):
		#Open User manual to help
		manual = UserManualDialog(method=1)
		manual.exec_()

	def _load_spec(self):
		#Read a single fits file
		filepath, check = QFileDialog.getOpenFileName(None,
			'Load 1 spectrum FITS file',
			WORKING_DIR,
			'Fits Files (*.fits)')
		if check:
			#print(type(file), file)

			# read fits file
			fitsfile = fits.open(filepath)
			# find wavelength, flux, error
			if 'FLUX' in fitsfile:
				self.fitsobj.wave = fitsfile['WAVELENGTH'].data
				self.fitsobj.flux = fitsfile['FLUX'].data
				self.fitsobj.error = fitsfile['ERROR'].data 
			else:
				for i in range(len(fitsfile)):
					search_list = np.array(fitsfile[i].header.cards)
					if 'flux' in search_list:
						self.fitsobj.flux = fitsfile[i].data['flux']
					elif 'FLUX' in search_list:
						self.fitsobj.flux = fitsfile[i].data['FLUX']

					if 'loglam' in search_list:
						self.fitsobj.wave = 10**(fitsfile[i].data['loglam'])
					elif 'WAVELENGTH' in search_list:
						self.fitsobj.wave = fitsfile[i].data['WAVELENGTH']

					if 'ivar' in search_list:
						self.fitsobj.error = 1. / np.sqrt(fitsfile[i].data['ivar'])
					elif 'ERROR' in search_list:
						self.fitsobj.error = fitsfile[i].data['ERROR']

			

			filename = self._get_filename(filepath, extension=False)
			#print(filename)

			self.mW.sc.plot_spec(self.fitsobj.wave, 
							self.fitsobj.flux, 
							self.fitsobj.error,
							filename)
			self.send_fitsobj.emit(self.fitsobj)
			self.send_filename.emit(filename)

	def _load_more_specs(self):
		#Read multiple spec fits file
		#	only fetch filepath strings, do not open fits files here
		filepaths, check = QFileDialog.getOpenFileNames(None,
			'Load multiple FITS files',
			WORKING_DIR,
			'Fits Files (*.fits)')
		if check:
			filenames = [self._get_filename(fp, extension=False) for fp in filepaths]
			if len(self.filepaths)>0:
				newfiles = [fi for fi in filenames if fi not in self.filenames]
				self.filenames.extend(newfiles)
				newpaths = [fpi for fpi in filepaths if fpi not in self.filepaths]
				self.filepaths.extend(newpaths)
				self.f_combobox.addItems(newfiles)
				self.f_combobox.setCurrentIndex(len(self.filenames)-len(newfiles)+1)
			else:
				self.filenames = filenames
				self.filepaths = filepaths
				self.f_combobox.addItems(self.filenames)
				self.f_combobox.setCurrentIndex(1)

	def _load_from_txt(self):
		# read all fits files saved in a txt file
		# the txt file should be located within the same folder as fits
		txtpath, check = QFileDialog.getOpenFileName(None,
			'Read a TXT file containing all fits filenames in the current folder',
			WORKING_DIR,
			'Text Files (*.txt)')
		if check:
			with open(txtpath) as f:
				lines = f.readlines()
			txt_dirname = os.path.dirname(txtpath)
			filepaths = [txt_dirname + '/' + line.strip() for line in lines]
			filenames = [self._get_filename(fp, extension=False) for fp in filepaths]

			if len(self.filepaths)>0:
				newfiles = [fi for fi in filenames if fi not in self.filenames]
				self.filenames.extend(newfiles)
				newpaths = [fpi for fpi in filepaths if fpi not in self.filepaths]
				self.filepaths.extend(newpaths)
				self.f_combobox.addItems(newfiles)
				self.f_combobox.setCurrentIndex(len(self.filenames)-len(newfiles)+1)
			else:
				self.filenames = filenames
				self.filepaths = filepaths
				self.f_combobox.addItems(self.filenames)
				self.f_combobox.setCurrentIndex(1)


	def _save_spec(self):
		#Save spec fits file with our own fits format

		filepath, check = QFileDialog.getSaveFileName(None,
			'Save 1 spectrum FITS file',
			'',
			'Fits Files (*.fits')
		if check:
			table = Table()
			table['WAVELENGTH'] = self.fitsobj.wave
			table['FLUX'] = self.fitsobj.flux
			table['ERROR'] = self.fitsobj.error
			print('Saving a fits file to [{}]'.format(filepath))
			table.write(filepath, format='fits')

	def _get_filename(self, filepath, extension=False):
		# return the filename and ready to pass to other widgets
		base = os.path.basename(filepath)
		if extension:
			return base
		else:
			return os.path.splitext(base)[0]

	# combobox event
	def _read_selected_fits(self, i):
		# default is 0 ==> No fits file
		# filepaths index = i - 1
		if i < 1:
			pass
		else:
			loadspec = LoadSpec(self.filepaths[i-1])
			filename = self.filenames[i-1]
			self.filename = filename
			self.fitsobj = loadspec._load_spec()
			if len(self.fitsobj.flux.shape) > 1:
				# new plot 2d
				self.mW.sc.plot_spec2d(self.fitsobj.wave,
									self.fitsobj.flux,
									self.fitsobj.error,
									filename)
				self._add_scale2d()
			else:
				self.mW.sc.plot_spec(self.fitsobj.wave, 
								self.fitsobj.flux, 
								self.fitsobj.error,
								filename)
				self.s_combobox.clear()
			
			self.send_filename.emit(filename)
			self.send_fitsobj.emit(self.fitsobj)		

	def _add_scale2d(self):
		self.s_combobox.setMaxVisibleItems(4)
		self.s_combobox.addItems(['Linear', 'Log', 'Sqrt', 'Square'])
		self.s_combobox.setCurrentIndex(0)
		self.s_combobox.currentIndexChanged.connect(self._scaling_changed)

		self.n_combobox.addItems(['None','MinMax', '99.5%', '99%', '98%', '97%', '96%', '95%', '92.5%', '90%', 'Z-Score', 'Manual'])
		self.n_combobox.setCurrentIndex(0)
		self.n_combobox.currentIndexChanged.connect(self._normalization_changed)

		self.min_range.setReadOnly(False)
		self.max_range.setReadOnly(False)

	def _scaling_changed(self, i):
		if len(self.fitsobj.flux.shape) > 1:
			self.mW.sc.plot_spec2d(self.fitsobj.wave,
								self.fitsobj.flux,
								self.fitsobj.error,
								self.filename, scale=i)
			self.scale = i
		else:
			pass

	def _normalization_changed(self, i):
		if len(self.fitsobj.flux.shape) > 1:
			if i < 11:
				self.mW.sc.plot_spec2d(self.fitsobj.wave,
									self.fitsobj.flux,
									self.fitsobj.error,
									self.filename, 
									scale=self.scale,
									normalization=i)
		else:
			pass

	def _on_scale_limits_slot(self, sent_limits):
		self.min_range.setText(str(sent_limits[0]))
		self.max_range.setText(str(sent_limits[1]))

	def _return_pressed(self):
		manual_range = [float(self.min_range.text()), float(self.max_range.text())]
		self.n_combobox.setCurrentIndex(11)
		self.mW.sc.plot_spec2d(self.fitsobj.wave,
							self.fitsobj.flux,
							self.fitsobj.error,
							self.filename,
							scale=self.scale,
							normalization=manual_range)


#----------------------------- Menu bar ---------------------------
class Custom_MenuBar(QMenuBar):

	send_fitsobj = pyqtSignal(object)
	send_linelist = pyqtSignal(object)
	send_newlinelist = pyqtSignal(object)
	send_z_est = pyqtSignal(object)
	send_filename = pyqtSignal(str)

	#Create our Custom Menu bar
	def __init__(self, mainWindow):
		super().__init__()
		self.mW = mainWindow
		self.fitsobj = mainWindow.fitsobj
		self.linelist = mainWindow.linelist
		self.newlinelist = mainWindow.newlinelist
		self.z_est = mainWindow.z_est

		#------------------- Necessary Variable -------------------
		self.newlinelist = []


		#------------------- Menu Layout --------------------------

		# 1st level menu
		filemenu = self.addMenu('&File')
		editmenu = self.addMenu('&Edit')
		helpmenu = self.addMenu('&Help')

		# File - 2nd level menu
		loadspec_act = self._create_action_in_menu('Load Spectrum',
												   'Load a spectrum fits file from local folder')
		loadspec_act.triggered.connect(self._load_spec)

		savespec_act = self._create_action_in_menu('Save Spectrum',
												   'Save a spectrum fits file to local folder')
		savespec_act.triggered.connect(self._save_spec)

		loadz_act = self._create_action_in_menu('Load redshift',
												'Load redshift list file from local folder')
		loadz_act.triggered.connect(self._load_z)
		savez_act = self._create_action_in_menu('Save redshift',
												'Save redshift list file to local folder')


		loadlinelist_act = self._create_action_in_menu('Load LineList',
													   'Load LineList file to local folder')
		loadlinelist_act.triggered.connect(self._load_linelist)

		savelinelist_act = self._create_action_in_menu('Save LineList',
													   'Save LineList file to local folder')
		savelinelist_act.triggered.connect(self._save_linelist)

		self._add_actions_to_menu(filemenu, [loadspec_act, savespec_act, loadz_act, savez_act, loadlinelist_act, savelinelist_act])

		# Edit - 2nd level menu
		# do not know what goes here so far

		# Help - 2nd level menu
		about_act = self._create_action_in_menu('User Manual',
												'Introduction to Keystrokes and Mouse-clicks in the GUI')
		about_act.triggered.connect(self._open_user_manual)

		contact_act = self._create_action_in_menu('GitHub Update',
												  'Link to GitHub for Updates')
		contact_act.triggered.connect(lambda: QDesktopServices.openUrl(QUrl('https://github.com/rongmon/rbcodes')))
		self._add_actions_to_menu(helpmenu, [about_act, contact_act])
		
	def _create_action_in_menu(self, actionName='', actionTip=''):
		#Create buttons wrapper
		action = QAction(actionName, self.mW)
		action.setStatusTip(actionTip)
		return action

	def _add_actions_to_menu(self, parentMenu, actions):
		#Add multiple actions to specific menu
		for i in range(len(actions)):
			parentMenu.addAction(actions[i])	

	def _open_user_manual(self):
		#Open User manual to help
		manual = UserManualDialog()
		manual.exec_()

	def _load_spec(self):
		#Read spec fits file
		filepath, check = QFileDialog.getOpenFileName(None,
			'Load 1 spectrum FITS file',
			WORKING_DIR,
			'Fits Files (*.fits)')
		if check:
			#print(type(file), file)

			# read fits file
			fitsfile = fits.open(filepath)
			# find wavelength, flux, error
			self.fitsobj.wave = fitsfile['WAVELENGTH'].data
			self.fitsobj.flux = fitsfile['FLUX'].data
			self.fitsobj.error = fitsfile['ERROR'].data 

			self.send_fitsobj.emit(self.fitsobj)

			filename = self._get_filename(filepath, extension=False)
			#print(filename)

			self.mW.sc.plot_spec(self.fitsobj.wave, 
							self.fitsobj.flux, 
							self.fitsobj.error,
							filename)

	def _save_spec(self):
		#Save spec fits file with our own fits format

		filepath, check = QFileDialog.getSaveFileName(None,
			'Save 1 spectrum FITS file',
			'',
			'Fits Files (*.fits')
		if check:
			table = Table()
			table['WAVELENGTH'] = self.fitsobj.wave
			table['FLUX'] = self.fitsobj.flux
			table['ERROR'] = self.fitsobj.error
			print(filepath)
			table.write(filepath, format='fits')
			'''Output fits format
			table[1].data['WAVELENGTH']
			table[1].data['FLUX']
			table[1].data['ERROR']
			'''
	def _load_linelist(self):
		#Load linelist from lines folder
		filepath, check = QFileDialog.getOpenFileName(None,
			'Load 1 linelist',
			LINELIST_DIR,
			'ASCII Files (*.ascii);; LST Files (*.lst)')
		if check:
			if filepath.endswith('.ascii'):
				rawdata = ascii.read(filepath)
				self.linelist = rawdata.to_pandas()
			elif filepath.endswith('.lst'):
				if filepath.endswith('gal_vac.lst'):
					with open(filepath, 'r') as f:
						colnames = f.readline().replace(' ', '')[:-1].split(',')
						#self.linelist = pd.DataFrame(columns=colnames)
						lines = f.readlines()
					self.linelist = pd.DataFrame(columns=['wave', 'ID', 'name'])
					for line in lines:
						tmp = line.replace(' ', '')[:-1].split(',')
						row = [float(tmp[0]), int(tmp[1]), tmp[2]+'_'+tmp[3]]
						self.linelist = self.linelist.append({'wave': row[0],
															  'ID': row[1],
															  'name': row[2]},
															  ignore_index=True)
				elif filepath.endswith('lbg.lst'):
					with open(filepath, 'r') as f:
						colnames = f.readline().split()
						#self.linelist = pd.DataFrame(columns=colnames)
						lines = f.readlines()
					self.linelist = pd.DataFrame(columns=['wave', 'ID', 'name'])
					for line in lines:
						tmp = line.split()
						row = [float(tmp[0]), int(tmp[1]), tmp[2]+'_'+tmp[3]]
						self.linelist = self.linelist.append({'wave': row[0],
															  'ID': row[1],
															  'name': row[2]},
															  ignore_index=True)
					

			self.send_linelist.emit(self.linelist)

			filename = self._get_filename(filepath, extension=False)
			self.send_filename.emit(filename)

	def _save_linelist(self):
		#Save linelist to lines folder
		filepath, check = QFileDialog.getSaveFileName(None,
			'Save the current linelist',
			LINELIST_DIR,
			'ASCII FIles (*.ascii')
		if check:
			linetable = Table.from_pandas(self.newlinelist)
			print(filepath)
			linetable.write(filepath, format='ascii')

	def _load_z(self):
		#Load estimated redshift working file
		filepath, check = QFileDialog.getOpenFileName(None,
			'Load estimated redshifts',
			'',
			'TEXT Files (*.txt)')
		if check:
			self.z_est = pd.read_csv(filepath, sep=',')
			
			self.send_z_est.emit(self.z_est)

			filename = self._get_filename(filepath, extension=False)
			print(filename)

	def _save_z(self):
		#Save current estimated redshifts so far
		filepath, check = QFileDialog.getSaveFileName(None,
			'Save the current estimated redshifts',
			'',
			'TEXT Files (*.txt')
		if check:
			self.z_est.to_csv(filepath, index=False)

	def _get_filename(self, filepath, extension=False):
		# return the filename and ready to pass to other widgets
		base = os.path.basename(filepath)
		if extension:
			return base
		else:
			return os.path.splitext(base)[0]
