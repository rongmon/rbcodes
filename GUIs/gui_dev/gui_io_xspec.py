import numpy as np
import astropy.units as u
from linetools.spectra.xspectrum1d import XSpectrum1D

from utils import FitsObj

class LoadXSpec():
	def __init__(self, filepath=''):
		self.fitsobj = FitsObj(wave=[])
		self.filepath = filepath
		self.warning = ''

		try:
			self.sp = XSpectrum1D.from_file(self.filepath)
		except (OSError, KeyError):
			self.warning += 'XSpectrum1D cannot read this file.'

	def _load_spec(self):
		if len(self.warning) < 1:
			# no error
			# need to convert astropy quantity class to numpy array
			self.fitsobj.flux = self.sp.flux.value
			self.fitsobj.error = self.sp.sig.value

			# need to convert wavelength unit to Angstrom
			self.fitsobj.wave = self.sp.wavelength.to(u.Angstrom).value
			

			return self.fitsobj



