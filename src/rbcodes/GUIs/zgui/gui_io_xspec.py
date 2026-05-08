import numpy as np
from astropy import units as u

from rbcodes.GUIs.zgui.utils import FitsObj
from rbcodes.utils.rb_spectrum import rb_spectrum


class LoadXSpec():
	def __init__(self, filepath=''):
		self.fitsobj = FitsObj(wave=[])
		self.filepath = filepath
		self.warning = ''

		try:
			self.sp = rb_spectrum.from_file(self.filepath)
			if self.sp._read_failed:
				self.warning = 'rb_spectrum could not read this file.'
		except Exception as e:
			self.warning = f'rb_spectrum could not read this file: {e}'

	def _load_spec(self):
		if len(self.warning) < 1:
			self.fitsobj.flux = self.sp.flux.value
			if self.sp.sig_is_set:
				self.fitsobj.error = self.sp.sig.value
			else:
				self.fitsobj.error = self.sp.flux.value * 0.001

			# wavelength already in Angstrom from rb_spectrum
			self.fitsobj.wave = self.sp.wavelength.to(u.Angstrom).value

			return self.fitsobj
