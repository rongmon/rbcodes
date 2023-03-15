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

		# pre-extract 1D from 2D, nanmean on y-axis
		for name2d,spec2d in self.frames.items():
			if spec2d is None:
				self.frames1d.update({name2d+'1d' : None})
			else:
				tmp_spec1d = np.nanmean(spec2d, axis=0)
				self.frames1d.update({name2d+'1d': tmp_spec1d})

			if name2d == 'ERR':
				self.frames1d.update({'ERR1d': np.nanmean(self.frames_err['SCI'], axis=0)})
			elif name2d in ['ERR_B', 'ERRB']:
				if 'ERRB' in labels:
					self.frames1d.update({'ERRB1d': np.nanmean(self.frames_err['SCIB'], axis=0)})
				elif 'ERR_B' in labels:
					self.frames1d.update({'ERR_B1d': np.nanmean(self.frames_err['SCI_B'], axis=0)})
			elif label in ['ERR_A', 'ERRA']:
				if 'ERRA' in labels:
					self.frames1d.update({'ERRA1d': np.nanmean(self.frames_err['SCIA'], axis=0)})
				elif 'ERR_A' in labels:
					self.frames_err['ERR_A1d'] = hdul['ERR_A'].data
					self.frames1d.update({'ERR_A1d': np.nanmean(self.frames_err['SCI_A'], axis=0)})

			#print(self.frames1d)

		# release internal memory
		hdul.close()
		del hdul
		return (self.frames, self.frames_err, self.frames1d)


'''PREVIOUS function used to find 2D/1D fits files and plot together
	def _check_available_frames(self):
		# check if the corresponding 2D fits file exists
		fnlist = self.filepath.split('_')
		fits2d = '_'.join(fnlist[:-2] + ['2D'] + [fnlist[-1]])
		print('TF 2D filename: ' + fits2d)

		# bug here...
		#if not os.path.exists(fits2d):
		#	return self.frames, self.frames1d

		# read fits file from filepath
		hdul = fits.open(fits2d)
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



		# search 1D file and open
		fnlist = hdul.filename().split('_')
		fits1d = '_'.join(fnlist[:-2] + ['1D'] + [fnlist[-1]])
		print('TF 1D filename: ' + fits1d)
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
		else:
			# if no 1D file available
			# pre-extract 1D from 2D, nanmean on y-axis
			for name2d,spec2d in self.frames.items():
				if spec2d is None:
					self.frames1d.update({name2d+'1d' : None})
				else:
					tmp_spec1d = np.nanmean(spec2d, axis=0)
					self.frames1d.update({name2d+'1d': tmp_spec1d})

				if name2d == 'ERR':
					self.frames1d.update({'ERR1d': np.nanmean(self.frames_err['SCI'], axis=0)})
				elif name2d in ['ERR_B', 'ERRB']:
					if 'ERRB' in labels:
						self.frames1d.update({'ERRB1d': np.nanmean(self.frames_err['SCIB'], axis=0)})
					elif 'ERR_B' in labels:
						self.frames1d.update({'ERR_B1d': np.nanmean(self.frames_err['SCI_B'], axis=0)})
				elif label in ['ERR_A', 'ERRA']:
					if 'ERRA' in labels:
						self.frames1d.update({'ERRA1d': np.nanmean(self.frames_err['SCIA'], axis=0)})
					elif 'ERR_A' in labels:
						self.frames_err['ERR_A1d'] = hdul['ERR_A'].data
						self.frames1d.update({'ERR_A1d': np.nanmean(self.frames_err['SCI_A'], axis=0)})

			#print(self.frames1d)

		# release internal memory
		hdul.close()
		del hdul
		return (self.frames, self.frames_err, self.frames1d)

'''