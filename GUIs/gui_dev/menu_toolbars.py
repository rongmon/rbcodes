import sys

from astropy.io import fits
from astropy.table import Table

from PyQt5.QtCore import Qt, QSize, QUrl
from PyQt5.QtWidgets import QAction, QToolBar, QStatusBar, QMenuBar, QFileDialog
from PyQt5.QtGui import QKeySequence, QDesktopServices

from user_manual import UserManualDialog
from utils import FitsObj


class Custom_ToolBar(QToolBar):
	'''Our personalized custom Toolbar
	'''
	def __init__(self, mainWindow):
		super().__init__()
		self.mW = mainWindow

		self.setWindowTitle('Customizable TooBar')
		#self.setFixedSize(QSize(200, 50))

		btn1 = self._create_button('Btn1', 'First trial')
		btn2 = self._create_button('Btn2', 'Second trial')

		self.addAction(btn1)
		self.addSeparator()
		self.addAction(btn2)

	def _create_button(self, buttonName='', buttonTip=''):
		'''Create buttons wrapper
		Parameters:
		----------
		buttonName (str):		Name of the action
		buttonTip (str):		Brief notes for this button

		Return:
		------
		btn (QAction):			implemented button ready to be added
		'''
		btn = QAction(buttonName, self.mW)
		btn.setStatusTip(buttonTip)

		return btn

class Custom_MenuBar(QMenuBar):
	'''Create our Custom Menu bar
	'''
	def __init__(self, mainWindow):
		super().__init__()
		self.mW = mainWindow

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
		savez_act = self._create_action_in_menu('Save redshift',
												'Save redshift list file to local folder')
		loadlinelist_act = self._create_action_in_menu('Load LineList',
													   'Load LineList file to local folder')
		savelinelist_act = self._create_action_in_menu('Save LineList',
													   'Save LineList file to local folder')
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
		'''Create buttons wrapper
		Parameters:
		----------
		parentMenu (QMenu):		parent QMenu object which attaches the btn
		buttonName (str):		Name of the action
		buttonTip (str):		Brief notes for this button

		Return:
		------
		action (QAction):			implemented button ready to be added
		'''
		action = QAction(actionName, self.mW)
		action.setStatusTip(actionTip)
		return action

	def _add_actions_to_menu(self, parentMenu, actions):
		'''Add multiple actions to specific menu
		Parameters:
		----------
		parentMenu (QMenu):			parent QMenu consisting all corresponding actions
		action (list of QAction):	list of QAction ready to be added to parentMenu

		Return:
			None
		'''
		for i in range(len(actions)):
			parentMenu.addAction(actions[i])	

	def _open_user_manual(self):
		'''Open User manual to help
		'''
		manual = UserManualDialog()
		manual.exec_()

	def _load_spec(self):
		'''Read spec fits file
		'''
		file, check = QFileDialog.getOpenFileName(None,
			'Load 1 spectrum FITS file',
			'',
			'Fits Files (*.fits)')
		if check:
			#print(type(file), file)

			# read fits file
			fitsfile = fits.open(file)
			# find wavelength, flux, error
			self.mW.fitsobj.wave = fitsfile['WAVELENGTH'].data
			self.mW.fitsobj.flux = fitsfile['FLUX'].data
			self.mW.fitsobj.error = fitsfile['ERROR'].data 

			self.mW.sc.plot(self.mW.fitsobj.wave, 
							self.mW.fitsobj.flux, 
							self.mW.fitsobj.error)

	def _save_spec(self):
		'''Save spec fits file with our own fits format
		'''
		filename, check = QFileDialog.getSaveFileName(None,
			'Save 1 spectrum FITS file',
			'',
			'Fits Files (*.fits')
		if check:
			table = Table()
			table['WAVELENGTH'] = self.mW.fitsobj.wave
			table['FLUX'] = self.mW.fitsobj.flux
			table['ERROR'] = self.mW.fitsobj.error
			print(filename)
			table.write(filename, format='fits')
			'''Output fits format
			table[1].data['WAVELENGTH']
			table[1].data['FLUX']
			table[1].data['ERROR']
			'''
