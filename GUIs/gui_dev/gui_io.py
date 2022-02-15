import sys
from astropy.io import fits
import numpy as np

from utils import FitsObj
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
'''

# This class is used to load spectrum from fits (unknown format) and save spectrum to fits file (our format)
class LoadSpec():
	def __init__(self, filepath):
		self.filepath = filepath
		self.fitsobj = FitsObj(wave=[], flux=[], error=[])


	def _load_spec(self):
		# read fits file
		fitsfile = fits.open(self.filepath)
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
