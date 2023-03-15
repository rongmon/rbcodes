import sys


# ------------- Utility classes ------------------------
# define as many classes we want that are not classified in other classes

class FitsObj():
	'''Main object to save and traverse data between widgets
	'''
	def __init__(self, wave=None, flux=None, error=None, 
				ra=None, dec=None, z_est=None, z_guess=None,
				flag=None, flux2d=None, error2d=None):
		# 1D wavelength array
		self.wave = wave
		# 1D spectrum
		self.flux = flux
		# 1D error spectrum
		self.error = error
		# RA
		self.ra = ra
		# DEC
		self.dec = dec
		# Estimated redshift
		self.z_est = z_est
		# Initially guessed redshift
		self.z_guess = z_guess
		# Any comments
		self.flag = flag

		# if fits file contains an image
		# 2D spectrum in SCI frame
		self.flux2d = flux2d
		# 2D error spectrum in SCI frame
		self.error2d = error2d

class Fits_2dAux():
	'''Auxiliary object to save and traverse non-important data
	'''
	def __init__(self, stamp=None, wcs=None, contamination=None, source=None, continuum=None, zpdf=None):
		# 2D stamp in STAMP frame
		self.stamp = stamp
		self.wcs = wcs
		# 2D contamination in CNTM frame
		self.contamination = contamination
		# 2D source in SRC_IMG or SRC frame
		self.source = source
		# 2D continuum in CONT frame
		self.continuum = continuum
		# Posterior redshift distribution for z_guess
		self.zpdf = zpdf