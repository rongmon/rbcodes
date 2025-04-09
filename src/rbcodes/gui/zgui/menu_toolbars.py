import sys
import os
import pandas as pd
import numpy as np

from astropy.io import fits, ascii
from astropy.table import Table

from PyQt5.QtCore import Qt, QSize, QUrl, pyqtSignal
from PyQt5.QtWidgets import QAction, QToolBar, QStatusBar, QMenuBar, QFileDialog, QComboBox, QLineEdit, QLabel, QPushButton
from PyQt5.QtGui import QKeySequence, QDesktopServices

from user_manual import UserManualDialog

from utils import FitsObj, Fits_2dAux
from spec_advanced2d import ShowAdvanced, ZGuessPosterior

WORKING_DIR = os.path.abspath(os.getcwd()) + './example-data'
# This is not used anymore
#LINELIST_DIR = os.path.dirname(os.path.abspath(__file__)) + '/lines'

class Custom_ToolBar(QToolBar):
	send_fitsobj = pyqtSignal(object)
	send_filename = pyqtSignal(str)
	send_filenames = pyqtSignal(list)
	send_message = pyqtSignal(str)

	#Our personalized custom Toolbar
	def __init__(self, mainWindow):
		super().__init__()
		#------- internal variables -------
		self.mW = mainWindow
		self.loadspec = None
		self.fitsobj = mainWindow.fitsobj
		self.fits_2daux = Fits_2dAux()
		self.filepaths = []
		self.filenames = []
		self.filename = ''
		self.display2d = None # displayed 2d spectrum
		self.display2d_err = None # displayed 2d error spectrum
		self.scale = 0
		self.normalization = 0
		self.manual = None # No user manual yet
		self.extract1d = None # wait for sent_extract1d
		self.maxfile_lim = 50 # max displayed files in f_combobox at a time
		self.manual_range = [0,1] # manual scaling range preset
		self.num_length = 10 # number of total digits displayed in min/max boxes

		self.setWindowTitle('Customizable TooBar')
		#self.setFixedSize(QSize(200, 50))

		# "Read All" button
		btn_loadtxt = self._create_button('Read All', 'Read all fits files from a TXT file. Make sure you run "ls *.fits > filenames.txt" before clicking this!')
		self.addAction(btn_loadtxt)
		btn_loadtxt.triggered.connect(self._load_from_txt)

		# "Read FITS" button
		btn_loadf = self._create_button('Read FITS', 'Read a fits file and plot the spectrum')
		self.addAction(btn_loadf)
		btn_loadf.triggered.connect(self._load_specs)

		# "Save Extract1D" button
		btn_savef = self._create_button('Save Extract1D', 'Save the extracted 1D FITS file')
		self.addAction(btn_savef)
		btn_savef.triggered.connect(self._save_spec)
		self.addSeparator()

		# "PREV" buttons for file dropbox
		btn_prevf = self._create_button('PREV', 'Swtich to previous fits file loaded in backend.')
		self.addAction(btn_prevf)
		btn_prevf.triggered.connect(self._prev_fitsfile)
		# File dropbox
		self.f_combobox = QComboBox()
		self.f_combobox.setMinimumWidth(200)
		self.f_combobox.addItem('No FITS File')
		self.f_combobox.setCurrentIndex(0)
		# set the max visible items for one page
		self.f_combobox.setMaxVisibleItems(self.maxfile_lim)
		self.f_combobox.currentIndexChanged.connect(self._read_selected_fits)
		self.addWidget(self.f_combobox)
		# "NEXT" buttons for file dropbox
		btn_nextf = self._create_button('NEXT', 'Switch to next fits file loaded in backend.')
		self.addAction(btn_nextf)
		btn_nextf.triggered.connect(self._next_fitsfile)

		self.addSeparator()

		# Frame dropbox section
		if self.mW.toggle_frames:
			self.frames = None
			# '-' button for frame dropbox
			btn_prevframe = self._create_button('-', 'Swtich to previous available frame.')
			self.addAction(btn_prevframe)
			btn_prevframe.triggered.connect(self._prev_frame)
			# main frame dropbox
			self.frame_combobox = QComboBox()
			self.frame_combobox.setMinimumWidth(130)
			self.frame_combobox.addItem('NONE')
			self.frame_combobox.setCurrentIndex(0)
			self.frame_combobox.currentTextChanged.connect(self._select_frames)
			self.addWidget(self.frame_combobox)
			# '+' button for frame dropbox
			btn_nextframe = self._create_button('+', 'Switch to next available frame.')
			self.addAction(btn_nextframe)
			btn_nextframe.triggered.connect(self._next_frame)

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

		self.addSeparator()

		# "Advanced" Option
		btn_adv = self._create_button('Advanced', 'More 2D inspection')
		self.addAction(btn_adv)
		btn_adv.triggered.connect(self._on_advanced_option)

		self.addSeparator()
		# "Help" button
		btn_help = self._create_button('Help', 'Open the user manual for further help')
		self.addAction(btn_help)
		btn_help.triggered.connect(self._open_user_manual)	
		



	def _create_button(self, buttonName='', buttonTip=''):
		#Create buttons wrapper
		btn = QAction(buttonName, self.mW)
		btn.setStatusTip(buttonTip)

		return btn

	def _open_user_manual(self, checked):
		#Open User manual to help
		if self.manual is None:
			self.manual = UserManualDialog(method=1)
		self.manual.show()


	def _load_specs(self):
		# This function corresponds to the action of "Read FITS" button
		# Read multiple spec fits file
		#	only fetch filepath strings, do not open fits files here
		if self.mW.xspecio:
			filepaths, check = QFileDialog.getOpenFileNames(None,
				'Load multiple FITS/XSpec files',
				WORKING_DIR,
				'Fits Files (*.fits);; XSpec Files (*.xspec)')
		else:
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
				#if len(self.filenames) > self.maxfile_lim:
				#	self.send_message.emit('Max displayed files are reached! Some files are folded in the dropdown box.')
			else:
				self.filenames = filenames
				self.filepaths = filepaths
				self.f_combobox.addItems(self.filenames)
				self.f_combobox.setCurrentIndex(1)

	def _load_from_txt(self):
		# This function corresponds to the action of "Read All" button
		# Read all fits files saved in a txt file
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
				#if len(self.filenames) > self.maxfile_lim:
				#	self.send_message.emit('Max displayed files are reached! Some files are folded in the dropdown box.')
			else:
				self.filenames = filenames
				self.filepaths = filepaths
				self.f_combobox.addItems(self.filenames)
				self.f_combobox.setCurrentIndex(1)


	def _save_spec(self):
		# This function corresponds to the action of "Save Extract1D" button
		# Save extract1d (wave, flux1d, error1d) with a deepcopy of loaded 2D spec 
		if self.loadspec is not None:
			if self.extract1d is not None:
				newfilename = self.filename + '_ymin={}_ymax={}'.format(self.extract1d['YMIN'], self.extract1d['YMAX'])
				filepath, check = QFileDialog.getSaveFileName(None,
					'Save 1D Spectrum FITS file',
					newfilename,
					'Fits Files (*.fits)')
				print(filepath)
				if check:
					# create a deepcopy first
					hdul_copy = self.loadspec._save_copy()
					# and then append a BinTable with WAVELENGTH, FLUX1D, ERROR1D
					col1 = fits.Column(name='WAVELENGTH', format='D', array=self.extract1d['WAVELENGTH'])
					col2 = fits.Column(name='FLUX', format='D', array=self.extract1d['FLUX'])
					col3 = fits.Column(name='ERROR', format='D', array=self.extract1d['ERROR'])
					cols = fits.ColDefs([col1, col2, col3])
					tmp_hdu = fits.BinTableHDU.from_columns(cols)
					tmp_hdu.name = 'EXTRACT1D'
					hdul_copy.append(tmp_hdu)
					hdr = hdul_copy[-1].header
					# save ymin/ymax values in header as well
					hdr.set('YMIN', self.extract1d['YMIN'], 'ymin value of the extracted box')
					hdr.set('YMAX', self.extract1d['YMAX'], 'ymax value of the extracted box')

					# ymin/ymax are written in the final filename
					hdul_copy.writeto(filepath, overwrite=True)
					print('Saving a fits file to [{}]'.format(filepath))

			else:
				self.send_message.emit('Please select a new extraction box.')
				
				#self.extract1d = None # reset extraction object
		else:
			self.send_message.emit('Please read a FITS file first before saving!')

	def _get_filename(self, filepath, extension=False):
		# return the filename and ready to pass to other widgets
		base = os.path.basename(filepath)
		if extension:
			return base
		else:
			return os.path.splitext(base)[0]

	def _prev_fitsfile(self):
		# This function corresponds to the action of "Prev" button
		# Switch to the previous filename loaded in f_combobox
		# idx 0 is "No FITS File" placeholder
		idx_min = 1
		current = self.f_combobox.currentIndex()
		if current > idx_min:
			self.f_combobox.setCurrentIndex(current-1)
		else:
			pass

	def _next_fitsfile(self):
		# This function corresponds to the action of "Next" button
		# Switch to the next filename loaded in f_combobox
		current = self.f_combobox.currentIndex()
		idx_max = self.f_combobox.count() - 1
		if current < idx_max:
			self.f_combobox.setCurrentIndex(current+1)
		else:
			pass

	def _on_advanced_option(self):
		message = ''
		if self.fits_2daux.stamp is not None:
			img_sta = [self.fits_2daux.stamp]
			if self.fits_2daux.wcs is not None:
				img_sta.append(self.fits_2daux.wcs)
				# Scaling/Normalization only control the last img
				# make sure to swap positions for this one
				img_sta[0], img_sta[1] = img_sta[1], img_sta[0]
			self.img_sta = ShowAdvanced(img_sta, name='STAMP')
			self.img_sta.show()
			message += 'GUI found STAMP HDU.'
		else:
			message += 'GUI did not find STAMP HDU in the current fits file.'

		if self.fits_2daux.zpdf is not None:
			self.z_pdfplot = ZGuessPosterior(self.fits_2daux.zpdf)
			self.z_pdfplot.show()
			message += 'GUI found ZPosterior HDU.'
		else:
			message += 'GUI did not find ZPosterior HDU in the current fits file.'

		imgs, names = [], []
		if self.fitsobj.flux2d is not None:
			imgs.append(self.fitsobj.flux2d)
			names.append('FLUX')
		else:
			message += 'This fits file does not contain a 2D spectrum.'
		#if self.fitsobj.error2d is not None:
		#	imgs.append(self.fitsobj.error2d)
		#	names.append('ERROR')
		#else:
		#	message += 'This fits file does not contain a 2D error spectrum.'

		if self.fits_2daux.source is not None:
			imgs.append(self.fits_2daux.source)
			names.append('SOURCE')
			message += 'GUI found SRC HDU.'
		else:
			message += 'GUI did not find SRC HDU in the current fits file.'

		if self.fits_2daux.continuum is not None:
			imgs.append(self.fits_2daux.continuum)
			names.append('CONTINUUM')
			message += 'GUI found CONTINUUM HDU.'
		else:
			message += 'GUI did not find CONTINUUM HDU in the current fits file.'

		if self.fits_2daux.contamination is not None:
			imgs.append(self.fits_2daux.contamination)
			names.append('CONTAMINATION')
			message += 'GUI found CONTAMINATION HDU.'
		else:
			message += 'GUI did not find CONTAMINATION HDU in the current fits file.'
		if len(imgs) > 0:
			self.img_adv = ShowAdvanced(imgs, name=names)
			self.img_adv.show()
		else:
			message +='There is no additional images to inspect.'

		self.send_message.emit(message)

	# frame_combobox event
	def _select_frames(self, s):
		if self.fitsobj.flux2d is None:
			# only 1d spec exists
			pass
		else:
			# default is SCI
			if (s == 'SCI') | (s == ''):
				self.display2d = self.fitsobj.flux2d
				self.display2d_err = self.fitsobj.error2d
			else:
				self.display2d = self.frames[s]
				self.fitsobj.error2d = self.frames_err[s]
				# select individual frame



			# get current scaling and normalization variables
			s_cur = self.s_combobox.currentIndex()
			n_cur = self.n_combobox.currentIndex()
			self._scaling_changed(s_cur)
			self._normalization_changed(n_cur)

			'''
			# if we want to reset scaling and normalization
			# uncomment all of these
			self.mW.sc.plot_spec2d(self.fitsobj.wave,
										self.display2d,
										self.fitsobj.error2d,
										self.filename)
			# reset scale and normalization comboboxes
			#self.s_combobox.setCurrentIndex(0)
			#self.n_combobox.setCurrentIndex(0)
			'''

			if (s+'1d' in self.frames1d.keys()) and (self.frames1d[s+'1d'] is not None):
				#print('replotting...')
				self.mW.sc.replot(self.fitsobj.wave, 
								self.frames1d[s+'1d'], 
								self.frames1d['ERR1d'])
				self.send_message.emit('Optimal Extraction Found!')


	def _prev_frame(self):
		# This function corresponds to the action of "-" button
		# Switch to the previous frame loaded in frame_combobox
		idx_min = 0
		current = self.frame_combobox.currentIndex()
		if current > idx_min:
			self.frame_combobox.setCurrentIndex(current-1)
		else:
			self.frame_combobox.setCurrentIndex(self.frame_combobox.count() - 1)

	def _next_frame(self):
		# This function corresponds to the action of "-" button
		# Switch to the next frame loaded in frame_combobox
		current = self.frame_combobox.currentIndex()
		idx_max = self.frame_combobox.count() - 1
		if current < idx_max:
			self.frame_combobox.setCurrentIndex(current+1)
		else:
			# wrap around
			self.frame_combobox.setCurrentIndex(0)


	# combobox event
	def _read_selected_fits(self, i):
		# default is 0 ==> No fits file
		# filepaths index = i - 1
		if i < 1:
			pass
		else:
			self.filename = self.filenames[i-1]
			print(f'Currently, GUI reads FITS file from {self.filepaths[i-1]}')
			if self.mW.xspecio:
				# enable XSpectrum1D io
				#print('xspec io mode')
				from gui_io_xspec import LoadXSpec
				self.loadspec = LoadXSpec(self.filepaths[i-1])
				if len(self.loadspec.warning) > 0:
					# even XSpectrum1D cant read this file
					self.send_message.emit(self.loadspec.warning)
				else:
					self.fitsobj = self.loadspec._load_spec()
					# XSpectrum1D only reads 1D fits file
					self.mW.sc.plot_spec(self.fitsobj.wave, 
										self.fitsobj.flux, 
										self.fitsobj.error,
										self.filename)
					self.s_combobox.clear()
					self.send_filename.emit(self.filename)
					self.send_fitsobj.emit(self.fitsobj)

			else:
				# enable standalone IO class
				#print('default io mode')
				from gui_io import LoadSpec
				fnlist = self.filepaths[i-1].split('_')
				# searching for corresponding 2D for advanced display
				if '1D' in fnlist[-2]:
					fits2d_path = '_'.join(fnlist[:-2] + ['2D'] + [fnlist[-1]])
					if os.path.exists(fits2d_path):
						self.filepaths[i-1] = fits2d_path
						self.filename = fits2d_path.split('/')[-1].split('.')[0]

				self.loadspec = LoadSpec(self.filepaths[i-1])

				selfcheck = self.loadspec._load_spec()
				#print(selfcheck)
				if type(selfcheck) is str:
					print('bad fits format')
					self.send_message.emit(selfcheck)
				else:
					self.fitsobj = self.loadspec._load_spec()

					if (self.fitsobj.flux2d is not None) & (self.fitsobj.flux is None):
						# only 2d spec exists
						self.display2d = self.fitsobj.flux2d
						self.mW.sc.plot_spec2d(self.fitsobj.wave,
											self.display2d,
											self.fitsobj.error2d,
											self.filename)
						self._add_scale2d()

					elif (self.fitsobj.flux is not None) & (self.fitsobj.flux2d is None):
						# only 1d spec exists
						self.mW.sc.plot_spec(self.fitsobj.wave, 
										self.fitsobj.flux, 
										self.fitsobj.error,
										self.filename)
						self.s_combobox.clear()
						self.n_combobox.clear()
						self.min_range.clear()
						self.min_range.setReadOnly(True)
						self.max_range.clear()
						self.max_range.setReadOnly(True)

					elif (self.fitsobj.flux is not None) & (self.fitsobj.flux2d is not None):
						# both 1d and 2d specs exist

						# check if filename has keywords ymin/ymax for 2D fits
						if ('ymin' in self.filename) & ('ymax' in self.filename):
							flist = self.filename.split('.')[0].split('_')
							extraction_box = [int(flist[-2][5:]), int(flist[-1][5:])]
							self.display2d = self.fitsobj.flux2d
							self.mW.sc.plot_spec2d(self.fitsobj.wave,
												self.display2d,
												self.fitsobj.error2d,
												self.filename,
												prev_extraction=extraction_box)
							self.mW.sc.replot(self.fitsobj.wave, 
											self.fitsobj.flux, 
											self.fitsobj.error)
											
						else:
							# if not, set up extraction box as usual
							self.display2d = self.fitsobj.flux2d
							self.mW.sc.plot_spec2d(self.fitsobj.wave,
												self.display2d,
												self.fitsobj.error2d,
												self.filename,
												prev_extraction=None)
							self.mW.sc.replot(self.fitsobj.wave, 
											self.fitsobj.flux, 
											self.fitsobj.error)
						

						self._add_scale2d()

					self.send_filename.emit(self.filename)
					self.send_fitsobj.emit(self.fitsobj)

					self.fits_2daux = self.loadspec._check_advanced()

				# toggle_frames only available for default io
				if self.mW.toggle_frames:
					from gui_frame_io import ToggleFrames
					toggle_f = ToggleFrames(self.filepaths[i-1])
					self.frames, self.frames_err, self.frames1d = toggle_f._check_available_frames()

					self.frame_combobox.clear()
					for frame_name in self.frames.keys():
						if self.frames[frame_name] is not None:
							self.frame_combobox.addItem(str(frame_name))

			# reset scale and normalization comboboxes
			self.s_combobox.setCurrentIndex(0)
			self.n_combobox.setCurrentIndex(0)

					

	def _add_scale2d(self):
		self.s_combobox.setMaxCount(4)
		self.s_combobox.addItems(['Linear', 'Log', 'Sqrt', 'Square'])
		self.s_combobox.setCurrentIndex(0)
		#self.s_combobox.currentIndexChanged.connect(self._scaling_changed)
		self.s_combobox.activated.connect(self._scaling_changed)


		self.n_combobox.setMaxCount(12)
		self.n_combobox.addItems(['None','MinMax', '99.5%', '99%', '98%', '97%', '96%', '95%', '92.5%', '90%', 'Z-Score', 'Manual'])
		self.n_combobox.setCurrentIndex(0)
		#self.n_combobox.currentIndexChanged.connect(self._normalization_changed)
		self.n_combobox.activated.connect(self._normalization_changed)

		self.min_range.setReadOnly(False)
		self.max_range.setReadOnly(False)

	def _scaling_changed(self, i):
		if self.fitsobj.flux2d is not None:
			self.mW.sc.plot_spec2d(self.fitsobj.wave,
								self.display2d,
								self.fitsobj.error2d,
								self.filename, scale=i,
								normalization=self.normalization)
			
			self.scale = i
		else:
			pass

	def _normalization_changed(self, i):
		if self.fitsobj.flux2d is not None:
			if i < 11:
				self.mW.sc.plot_spec2d(self.fitsobj.wave,
									self.display2d,
									self.fitsobj.error2d,
									self.filename, 
									scale=self.scale,
									normalization=i)
				if i == 1:
					self.min_range.setText(str(0)[:self.num_length])
					self.max_range.setText(str(1)[:self.num_length])
				elif i == 10:
					self._on_scale_limits_slot([np.nan, np.nan])
				self.normalization = i
			else:
				self.n_combobox.setCurrentIndex(11)
				self.mW.sc.plot_spec2d(self.fitsobj.wave,
									self.display2d,
									self.fitsobj.error2d,
									self.filename,
									scale=self.scale,
									normalization=self.manual_range)
				self._on_scale_limits_slot(self.manual_range)


		else:
			pass

	def _on_scale_limits_slot(self, sent_limits):
		if np.isnan(sent_limits[0]) or np.isnan(sent_limits[-1]):
			self.min_range.setReadOnly(True)
			self.max_range.setReadOnly(True)
			self.min_range.setText('Mean')
			self.max_range.setText('STD')
		else:
			self.min_range.setReadOnly(False)
			self.max_range.setReadOnly(False)
			self.min_range.setText(str(sent_limits[0])[:self.num_length])
			self.max_range.setText(str(sent_limits[1])[:self.num_length])

	def _return_pressed(self):
		# min,max current values
		manual_range = [float(self.min_range.text()), float(self.max_range.text())]
		# sort and assign min,max values to avoid user errors
		manual_range.sort()
		self.min_range.setText(str(manual_range[0])[:self.num_length])
		self.max_range.setText(str(manual_range[-1])[:self.num_length])


		self.n_combobox.setCurrentIndex(11)
		self.mW.sc.plot_spec2d(self.fitsobj.wave,
							self.display2d,
							self.fitsobj.error2d,
							self.filename,
							scale=self.scale,
							normalization=manual_range)
		# save manual_range for other frames
		self.manual_range = manual_range

	def _on_sent_extract1d(self, sent_extract1d):
		self.extract1d = sent_extract1d


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
