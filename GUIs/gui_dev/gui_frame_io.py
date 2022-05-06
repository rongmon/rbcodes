import sys
import os
import copy
from astropy.io import fits
import numpy as np

class ToggleFrames():
	'''Main class to search targeted frames'''
	FRAMENAMES = ['SCI','SCI_A','SCI_B', 'EMLINE',
			'EMLINE_A', 'EMLINE_B','CONT','MODEL', 'CNTM']

	def __init__(self, filepath):
		# initialize internel variables
		self.filepath = filepath
		self.frames = {key:None for key in self.FRAMENAMES}

	def _check_available_frames(self):
		# read fits file from filepath
		hdul = fits.open(self.filepath)
		labels = [label.name for label in hdul]

		for label in labels:
			if label in self.FRAMENAMES:
				self.frames[label] = hdul[label].data
				# release internal memory
				del hdul[label].data
			else:
				self.frames[label] = None
		# release internal memory
		hdul.close()
		del hdul

		return self.frames
