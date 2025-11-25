import sys
import os
import copy
from astropy.io import fits
import numpy as np

from rbcodes.GUIs.zgui.utils import FitsObj, Fits_2dAux
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
	# Initialization - requires filepath
	def __init__(self, filepath):
		# interval values
		self.filepath = filepath
		self.fitsobj = FitsObj(wave=[], flux=None, error=None)
		self.fits_2daux = Fits_2dAux()

	# Read FITS file with various format options
	def _load_spec(self):
		# read fits file

		fitsfile = fits.open(self.filepath)
		labels = [label.name for label in fitsfile]

		# Detect frame sources from FITS headers
		detected_frames = self._detect_frame_sources()
		self.fitsobj.frame_sources = detected_frames

		# ------------ FORMAT 1 --------------
		# check if filename has keywords ymin/ymax for 2D fits
		if 'ymin' in fitsfile.filename():
			fname = fitsfile.filename().split('.')[0]
			flist = fname.split('_')
			extraction_box = [int(flist[-2][5:]), int(flist[-1][5:])]

			self.fitsobj.flux2d = fitsfile['SCI'].data
			self.fitsobj.error2d = fitsfile['ERR'].data
			# Set RA DEC
			if 'RA' in fitsfile['SCI'].header:
				self.fitsobj.ra = fitsfile['SCI'].header['RA']
				self.fitsobj.dec = fitsfile['SCI'].header['DEC']

			# check if there is EXTRACT1D BinTableHDU
			if 'EXTRACT1D' in labels:
				# This is FITS format that GUI saved
				self.fitsobj.wave = fitsfile['EXTRACT1D'].data['WAVELENGTH']
				self.fitsobj.flux = fitsfile['EXTRACT1D'].data['FLUX']
				self.fitsobj.error = fitsfile['EXTRACT1D'].data['ERROR']

			fitsfile.close()
			return self.fitsobj



		# ------------ FORMAT 2 --------------
		# for a pair of fits files:
		# cal-2D, x1d-1D
		if fitsfile.filename().endswith('cal.fits'):
			self.fitsobj.flux2d = fitsfile['SCI'].data
			self.fitsobj.error2d = fitsfile['ERR'].data
			
			# search 1D file and open
			fits1d = fitsfile.filename()[:-8] + 'x1d.fits'
			if os.path.exists(fits1d):
				fitsfile1d = fits.open(fits1d)
				scale = self._scale_wave_unit(fitsfile1d['EXTRACT1D'].header)
				self.fitsobj.wave = fitsfile1d['EXTRACT1D'].data['WAVELENGTH'] * scale
				self.fitsobj.flux = fitsfile1d['EXTRACT1D'].data['FLUX']
				self.fitsobj.error = fitsfile1d['EXTRACT1D'].data['FLUX_ERROR']
				fitsfile1d.close()
			fitsfile.close()
			return self.fitsobj

		elif fitsfile.filename().endswith('x1d.fits'):
			scale = self._scale_wave_unit(fitsfile['EXTRACT1D'].header)
			self.fitsobj.wave = fitsfile['EXTRACT1D'].data['WAVELENGTH'] * scale
			self.fitsobj.flux = fitsfile['EXTRACT1D'].data['FLUX']
			self.fitsobj.error = fitsfile['EXTRACT1D'].data['FLUX_ERROR']

			# search 2D file and open
			fits2d = fitsfile.filename()[:-8] + 'cal.fits'
			#print('prepare 2D spec')
			if os.path.exists(fits2d):
				fitsfile2d = fits.open(fits2d)
				self.fitsobj.flux2d = fitsfile2d['SCI'].data
				self.fitsobj.error2d = fitsfile2d['ERR'].data
				fitsfile2d.close()
			fitsfile.close()
			return self.fitsobj



		# ------------ FORMAT 3 --------------
		# EIGER fits files 
		# for a pair of fits files:
		# XXX_2D_YYY.fits; XXX_1D_YYY.fits
		if '2D' in fitsfile.filename().split('_')[-2]:
			self.fitsobj.flux2d = fitsfile['SCI'].data
			self.fitsobj.error2d = fitsfile['ERR'].data
			self.fitsobj.wave = self._build_wave(fitsfile['SCI'].header)
			
			# Set RA DEC
			if 'RA' in fitsfile['SCI'].header:
				self.fitsobj.ra = fitsfile['SCI'].header['RA']
				self.fitsobj.dec = fitsfile['SCI'].header['DEC']

			# Check z_guess
			if 'EAZY_ZPDF' in labels:
				z = fitsfile['EAZY_ZPDF'].data['z']
				pdf = fitsfile['EAZY_ZPDF'].data['pdf']
				self.fitsobj.z_guess = z[np.argmax(pdf)]

			# search 1D file and open
			fnlist = fitsfile.filename().split('_')
			fits1d = '_'.join(fnlist[:-2] + ['1D'] + [fnlist[-1]])
			if os.path.exists(fits1d):
				fitsfile1d = fits.open(fits1d)
				self.fitsobj.flux = fitsfile1d[1].data['flux_opt_ext']
				self.fitsobj.error = fitsfile1d[1].data['flux_opt_ext_err']
				fitsfile1d.close()
			fitsfile.close()
			return self.fitsobj

		elif '1D' in fitsfile.filename().split('_')[-2]:
			self.fitsobj.flux = np.nan_to_num(fitsfile[1].data['flux_opt_ext'], nan=np.nan)
			self.fitsobj.error = np.nan_to_num(fitsfile[1].data['flux_opt_ext_err'], nan=np.nan)
			self.fitsobj.wave = fitsfile[1].data['wavelength']
			wscale = self._scale_wave_unit(fitsfile[1].header)
			self.fitsobj.wave *=wscale


			# search 2D file and open
			fnlist = fitsfile.filename().split('_')
			fits2d = '_'.join(fnlist[:-2] + ['2D'] + [fnlist[-1]])
			if os.path.exists(fits2d):
				fitsfile2d = fits.open(fits2d)
				self.fitsobj.flux2d = fitsfile2d['SCI'].data
				self.fitsobj.error2d = fitsfile2d['ERR'].data
				self.fitsobj.wave = self._build_wave(fitsfile2d['SCI'].header)
				
				# Set RA DEC
				if 'RA' in fitsfile2d['SCI'].header:
					self.fitsobj.ra = fitsfile2d['SCI'].header['RA']
					self.fitsobj.dec = fitsfile2d['SCI'].header['DEC']

				# Check z_guess
				if 'EAZY_ZPDF' in labels:
					z = fitsfile2d['EAZY_ZPDF'].data['z']
					pdf = fitsfile2d['EAZY_ZPDF'].data['pdf']
					self.fitsobj.z_guess = z[np.argmax(pdf)]
				fitsfile2d.close()
			fitsfile.close()
			return self.fitsobj




		# ------------ FORMAT 4 --------------
		# CHECK IF FITS FILE HAS AN IMAGE
		if (fitsfile[0].header['NAXIS'] > 1) & (len(fitsfile) < 2):
			# example file: long_radd.fits
			# delete this condition if long_radd.fits is no longer used
			if 'long_radd' in fitsfile.filename().split('.')[0]:
				self.fitsobj.flux2d = np.transpose(fitsfile[0].data)
			else:
				self.fitsobj.flux2d = fitsfile[0].data
			wave0,wave1 = fitsfile[0].header['ADCWAVE0'], fitsfile[0].header['ADCWAVE1']
			self.fitsobj.wave = np.linspace(wave0, wave1, len(self.fitsobj.flux2d[0]))
			# fake error 2d spectrum
			self.fitsobj.error2d = self.fitsobj.flux2d * 0.05

			fitsfile.close()
			return self.fitsobj

		# ------------ FORMAT 5 --------------
		elif 'SCI' in labels:
    	#Read in a specific format to account for EIGER emission line 2d spectrum
			# example file: spec2d_coadd_QSO_J0100_sID010242.fits
			self.fitsobj.flux2d = fitsfile['SCI'].data
			self.fitsobj.error2d = fitsfile['ERR'].data
			self.fitsobj.wave = self._build_wave(fitsfile['SCI'].header)
			
			# Set RA DEC
			if 'RA' in fitsfile['SCI'].header:
				self.fitsobj.ra = fitsfile['SCI'].header['RA']
				self.fitsobj.dec = fitsfile['SCI'].header['DEC']

			# Check z_guess
			if 'EAZY_ZPDF' in labels:
				z = fitsfile['EAZY_ZPDF'].data['z']
				pdf = fitsfile['EAZY_ZPDF'].data['pdf']
				self.fitsobj.z_guess = z[np.argmax(pdf)]

			fitsfile.close()
			return self.fitsobj


		# ------------ FORMAT 6 --------------
		elif 'COADD' in labels:
			#FITS format 2:
			#file[i].data['<variable>']
			#for multiple i HDU table/image
			#example file: SDSS_spec.fits
			#Note: sdss1.fits is okay; sdss2.fits has header PLUG_RA
			
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

				if ('RA' in search_list) & (self.fitsobj.ra is None):
					self.fitsobj.ra = fitsfile[i].header['RA']
				if ('DEC' in search_list) & (self.fitsobj.dec is None):
					self.fitsobj.dec = fitsfile[i].header['DEC']

			fitsfile.close()
			return self.fitsobj


		# ------------ FORMAT 7 --------------		
		# CHECK IF FITS FILE HAS 1D SPECTRAL INFO...
		# find wavelength, flux, error
		elif 'FLUX' in fitsfile:
			'''FITS format 1:
				file['<variable>'].data
			example file: test.fits
			'''
			self.fitsobj.wave = fitsfile['WAVELENGTH'].data
			self.fitsobj.flux = fitsfile['FLUX'].data
			self.fitsobj.error = fitsfile['ERROR'].data
	
			fitsfile.close()
			#print(type(self.fitsobj))
			return self.fitsobj
		
		# Soft warning.. No compatible fits format
		if self.fitsobj.flux is None:
			return f'This GUI currently is not compatible with {fitsfile.filename()}...'



	def _build_wave(self, header):
		'''Returns a NumPy array containing wavelength axis of a 2d specturm in Angstrom.
			Args:
				header (astropy.io.fits.Header): header that contains wavelength axis
				that is specified in 'CTYPE' keywords in NAXIS1 dimension.
			Returns:
				numpy.ndarray: Wavelength axis for this data.
		'''
		#Get keywords defining wavelength axis
		nwav = header['NAXIS1']
		wav0 = header['CRVAL1']
		dwav = header['CDELT1']
		pix0 = header['CRPIX1']
		wave=np.array([wav0 + (i - pix0) * dwav for i in range(nwav)])

		# Now check units to make sure everything is in angstrom
		card='CUNIT1'
		if not card in header:
			raise ValueError("Header must contain 'CUNIT1' keywords.")
		
		# check units
		wunit = header[card].rstrip().lstrip() # get rid of whitespace	
		#micrometer to Angstrom
		if wunit in ['um', 'micron', 'micrometer']:
			wave *=10000. 
		elif wunit in ['nm', 'nanometer']:
			wave *=10
		elif wunit in ['AA', 'Angstrom']:
		#elif header[card]=='AA':
			wave=wave
		else:
			raise ValueError("Predefined wavelength units are 'um','nm','AA'.")            
		return wave

	def _scale_wave_unit(self, header):
		# if wavelength array already existed,
		card = 'TUNIT1'
		if not card in header:
			raise ValueError("Header must contain 'TUNIT1' keywords.")
		
		# check units
		wunit = header[card].rstrip().lstrip() # get rid of whitespace
		#micrometer to Angstrom
		if wunit in ['um', 'micron', 'micrometer']:
			scale = 10000. 
		elif wunit in ['nm', 'nanometer']:
			scale = 10.
		elif wunit in ['AA', 'Angstrom']:
		#elif header[card]=='AA':
			scale = 1.
		else:
			raise ValueError("Predefined wavelength units are 'um','nm','AA'.")            
		return scale

	# Save the deep copy of current opened FITS file
	def _save_copy(self):
		fitsfile = fits.open(self.filepath)
		fitsfile_copy = copy.deepcopy(fitsfile)
		fitsfile.close()
		return fitsfile_copy

	# Extract frame sources (EMLINE* frames) from FITS headers
	def _detect_frame_sources(self):
		"""Detect available frame sources from FITS file headers.

		Returns:
			list: Uppercase frame names found in FITS file (e.g., ['EMLINEA', 'EMLINEB'])
		"""
		detected_frames = set()
		try:
			fitsfile = fits.open(self.filepath)

			# First priority: Check HDU names for EMLINE* frames (most reliable)
			# HDU names like EMLINEA, EMLINEB, EMLINEC are the frame sources
			labels = [label.name.upper() for label in fitsfile]
			for label in labels:
				# Match EMLINE followed by single letter (A, B, C, etc.)
				if label.startswith('EMLINE') and len(label) == 7:
					# This is EMLINEA, EMLINEB, etc.
					detected_frames.add(label)

			# Second priority: Check HDU headers for EMLINE* keywords
			# Only extract if we haven't found frames via HDU names
			if not detected_frames:
				for hdu in fitsfile:
					if hdu.header is not None:
						# Look for keywords that start with EMLINE followed by letter(s)
						for key in hdu.header:
							key_upper = key.upper()
							if key_upper.startswith('EMLINE') and len(key_upper) > 6:
								# Extract the frame designation (e.g., EMLINEA from EMLINEA_XYZ)
								# Get characters after EMLINE until underscore or end
								remainder = key_upper[6:]  # Skip 'EMLINE'
								if remainder[0].isalpha():  # Next char must be letter
									frame_base = 'EMLINE' + remainder[0]
									if frame_base not in ['EMLINE']:  # Skip invalid
										detected_frames.add(frame_base)

			fitsfile.close()
		except Exception as e:
			print(f"Warning: Could not detect frame sources from FITS header: {e}")

		# Return sorted list of uppercase frame names
		return sorted(list(detected_frames))

	# Extract addtional information from FITS file if available
	# Addtional infomation includes
	#													(stamp, contamination, posterior redshift estimation)
	def _check_advanced(self):
		# read fits file
		fitsfile = fits.open(self.filepath)
		labels = [label.name for label in fitsfile]

		# Check if STAMP exists
		if 'STAMP' in labels:
			self.fits_2daux.stamp = fitsfile['STAMP'].data
			from astropy.wcs import WCS
			try:
				self.fits_2daux.wcs = WCS(fitsfile['STAMP'].header)
			except AttributeError:
				self.fits_2daux.wcs = None
				print('Current FITS file does not have required WCS info.')
		elif 'SRC_IMG' in labels:
			self.fits_2daux.stamp = fitsfile['SRC_IMG'].data
			from astropy.wcs import WCS
			try:
				self.fits_2daux.wcs = WCS(fitsfile['SRC_IMG'].header)
			except AttributeError:
				self.fits_2daux.wcs = None
				print('Current FITS file does not have required WCS info.')
		else:
			self.fits_2daux.stamp = None

		# Check CONTAMINATION
		if 'CNTM' in labels:
			self.fits_2daux.contamination = fitsfile['CNTM'].data
		else:
			self.fits_2daux.contamination = None

		# Check CONTINUUM frame
		if 'CONT' in labels:
			self.fits_2daux.continuum = fitsfile['CONT'].data
		else:
			self.fits_2daux.continuum = None

		# Check SOURCE
		if 'SRC' in labels:
			self.fits_2daux.source = fitsfile['SRC'].data
		else:
			self.fits_2daux.source = None

		# check z_guess posterior
		if 'EAZY_ZPDF' in labels:
			self.fits_2daux.zpdf = fitsfile['EAZY_ZPDF'].data
		else:
			self.fits_2daux.zpdf = None

		fitsfile.close()
		return self.fits_2daux