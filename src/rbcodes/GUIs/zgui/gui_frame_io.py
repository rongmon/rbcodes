import sys
import os
import copy
from astropy.io import fits
import numpy as np

class ToggleFrames():
	'''Main class to search targeted frames'''
	#FRAMENAMES = ['SCI', 'MODEL', 'CNTM', 'EMLINE', 'CONT',
	#		'SCI_A', 'EMLINE_A', 'SCI_B', 'EMLINE_B']
	FRAMENAMES = ['SCI','SCI_A','SCI_B','SCIA','SCIB', 
				'EMLINE', 'EMLINE_A', 'EMLINEA', 'EMLINEB', 'EMLINE_B',
				'CONT','MODEL', 'CNTM']

	FRAMENAMES1D = ['SCI1d', 'ERR1d', 'MODEL1d', 'EMLINE1d', 'CONT1d', 'CNTM1d']
	OPT_EXT_NAMES = ['flux_opt_ext', 'flux_opt_ext_err', 'flux_opt_model',
					'flux_opt_emline', 'flux_opt_ext_cntm_removed', 'flux_opt_cntm']

	def __init__(self, filepath):
		# initialize internel variables
		self.filepath = filepath
		self.frames = {key:None for key in self.FRAMENAMES}
		self.frames_err = {key:None for key in self.FRAMENAMES}
		self.frames1d = {key:None for key in self.FRAMENAMES1D}

	def _check_available_frames(self):
		# check if the corresponding 2D fits file exists
		#fnlist = self.filepath.split('_')
		#fits2d = '_'.join(fnlist[:-2] + ['2D'] + [fnlist[-1]])
		#if not os.path.exists(fits2d):
		#	return self.frames, self.frames1d

		# read fits file from filepath
		hdul = fits.open(self.filepath)
		labels = [label.name for label in hdul]

		for label in labels:
			if label in self.FRAMENAMES:
				self.frames[label] = hdul[label].data
				# release internal memory
				#del hdul[label].data
			else:
				self.frames[label] = None

		# search for error frames
		for label in labels:
			if label in ['SCI', 'EMLINE','CONT','MODEL', 'CNTM']:
				self.frames_err[label] = hdul['ERR'].data
			elif label in ['SCI_B', 'SCIB', 'EMLINEB', 'EMLINE_B']:
				if 'ERRB' in labels:
					self.frames_err[label] = hdul['ERRB'].data
				elif 'ERR_B' in labels:
					self.frames_err[label] = hdul['ERR_B'].data
			elif label in ['SCI_A', 'SCIA', 'EMLINEA', 'EMLINE_A']:
				if 'ERRA' in labels:
					self.frames_err[label] = hdul['ERRA'].data
				elif 'ERR_A' in labels:
					self.frames_err[label] = hdul['ERR_A'].data



		# release internal memory
		hdul.close()
		del hdul
		return self.frames, self.frames_err, self.frames1d
