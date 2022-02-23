import sys
from astropy.io import fits
import numpy as np

from utils import FitsObj, FitsObj2d
# use test.fits from rbcodes/example-data as testing example
'''
file = fits.open('test.fits')
file.info()
Filename: test.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  FLUX          1 PrimaryHDU      10   (19663,)   float32
  1  ERROR         1 ImageHDU         7   (19663,)   float32
  2  WAVELENGTH    1 ImageHDU         7   (19663,)   float64
  3  CONTINUUM     1 ImageHDU         7   (19663,)   float32


NEED TO FIND A BETTER WAY TO TELL IF A FITS FILE HAS 1D(spec) OR 2D(image) DATA
'''

# This class is used to load spectrum from fits (unknown format) and save spectrum to fits file (our format)
class LoadSpec():
	def __init__(self, filepath):
		self.filepath = filepath
		self.fitsobj = FitsObj(wave=[], flux=[], error=[])
		#self.fitsobj2d = FitsObj2d(wave=[], flux2d=[])


	def _load_spec(self):
		# read fits file
		fitsfile = fits.open(self.filepath)

		# check if fits file contains an image
		if fitsfile[0].header['NAXIS'] > 1:
			self.fitsobj.flux = fitsfile[0].data
			wave0,wave1 = fitsfile[0].header['ADCWAVE0'], fitsfile[0].header['ADCWAVE1']
			self.fitsobj.wave = np.linspace(wave0, wave1, len(self.fitsobj.flux))
			# fake error 2d spectrum
			self.fitsobj.error = self.fitsobj.flux * 0.05

			return self.fitsobj

		# if not, find spectral info...
		# find wavelength, flux, error
		if 'FLUX' in fitsfile:
			'''FITS format 1:
				file['<variable>'].data
			'''
			self.fitsobj.wave = fitsfile['WAVELENGTH'].data
			self.fitsobj.flux = fitsfile['FLUX'].data
			self.fitsobj.error = fitsfile['ERROR'].data 
		else:
			'''FITS format 2:
				file[i].data['<variable>']
				for multiple i HDU table/image
			'''
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

		return self.fitsobj
