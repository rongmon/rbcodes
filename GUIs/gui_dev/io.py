import sys
from astropy.io import fits

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
	def __init__(self, filename):
		self.filename = filename
		file_ext = filename.split('.')[-1]

		if file_ext != 'fits':
			raise TypeError('Spectrum must be in fits format!')

		# read fits file
		file = fits.open(self.filename)
		# find wavelength, flux, error
		self.wave = file['WAVELENGTH'].data
		self.flux = file['FLUX'].data
		self.error = file['ERROR'].data 

