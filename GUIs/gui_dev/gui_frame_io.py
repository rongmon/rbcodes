import sys
import os
import copy
from astropy.io import fits
import numpy as np

class ToggleFrames():
	'''Main class to search targeted frames'''
	#FRAMENAMES = ['SCI', 'MODEL', 'CNTM', 'EMLINE', 'CONT',
	#		'SCI_A', 'EMLINE_A', 'SCI_B', 'EMLINE_B']
	FRAMENAMES = ['SCI','SCI_A','SCI_B','SCIA','SCIB', 'EMLINE',
			'EMLINE_A','EMLINEA','EMLINEB', 'EMLINE_B','CONT','MODEL', 'CNTM']
	FRAMENAMES1D = ['SCI1d', 'ERR1d', 'MODEL1d', 'EMLINE1d', 'CONT1d', 'CNTM1d']
	OPT_EXT_NAMES = ['flux_opt_ext', 'flux_opt_ext_err', 'flux_opt_model',
					'flux_opt_emline', 'flux_opt_ext_cntm_removed', 'flux_opt_cntm']

	def __init__(self, filepath):
		# initialize internel variables
		self.filepath = filepath
		self.frames = {key:None for key in self.FRAMENAMES}
		self.frames1d = {key:None for key in self.FRAMENAMES1D}

	def _check_available_frames(self):
		# check if the corresponding 2D fits file exists
		fnlist = self.filepath.split('_')
		fits2d = '_'.join(fnlist[:-2] + ['2D'] + [fnlist[-1]])
		if not os.path.exists(fits2d):
			return self.frames, self.frames1d

		# read fits file from filepath
		hdul = fits.open(fits2d)
		labels = [label.name for label in hdul]

		for label in labels:
			if label in self.FRAMENAMES:
				self.frames[label] = hdul[label].data
				# release internal memory
				del hdul[label].data
			else:
				self.frames[label] = None


		# search 1D file and open
		fnlist = hdul.filename().split('_')
		fits1d = '_'.join(fnlist[:-2] + ['1D'] + [fnlist[-1]])
		if os.path.exists(fits1d):
			hdul1d = fits.open(fits1d)
			labels1d = hdul1d[1].columns.names

			for i in range(len(self.OPT_EXT_NAMES)):
				if self.OPT_EXT_NAMES[i] in labels1d:
					self.frames1d[self.FRAMENAMES1D[i]] = hdul1d[1].data[self.OPT_EXT_NAMES[i]]
				else:
					self.frames1d[self.FRAMENAMES1D[i]] = None

			# release internal memory
			hdul1d.close()
			del hdul1d


		# release internal memory
		hdul.close()
		del hdul
		return self.frames, self.frames1d
