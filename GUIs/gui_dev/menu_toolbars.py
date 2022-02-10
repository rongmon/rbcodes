import sys
import os
import pandas as pd

from astropy.io import fits, ascii
from astropy.table import Table

from PyQt5.QtCore import Qt, QSize, QUrl, pyqtSignal
from PyQt5.QtWidgets import QAction, QToolBar, QStatusBar, QMenuBar, QFileDialog
from PyQt5.QtGui import QKeySequence, QDesktopServices

from user_manual import UserManualDialog
from utils import FitsObj

WORKING_DIR = os.path.abspath(os.getcwd()) + './example-data'
LINELIST_DIR = os.path.dirname(os.path.abspath(__file__)) + '/lines'

class Custom_ToolBar(QToolBar):
	send_fitsobj = pyqtSignal(object)
	send_filename = pyqtSignal(str)
	#Our personalized custom Toolbar
	def __init__(self, mainWindow):
		super().__init__()
		self.mW = mainWindow
		self.fitsobj = mainWindow.fitsobj

		self.setWindowTitle('Customizable TooBar')
		#self.setFixedSize(QSize(200, 50))

		btn_loadf = self._create_button('Read FITS', 'Read a fits file and plot the spectrum')
		self.addAction(btn_loadf)
		btn_loadf.triggered.connect(self._load_spec)

		btn_savef = self._create_button('Save FITS', 'Save the fits file to unified format')
		self.addAction(btn_savef)
		self.addSeparator()

		btn_help = self._create_button('Help', 'Open the user manual for further help')
		self.addAction(btn_help)
		btn_help.triggered.connect(self._open_user_manual)	

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
		self.send_filename.emit(filename)

	def _get_filename(self, filepath, extension=False):
		# return the filename and ready to pass to other widgets
		base = os.path.basename(filepath)
		if extension:
			return base
		else:
			return os.path.splitext(base)[0]


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
