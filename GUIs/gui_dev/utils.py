import sys


# ------------- Utility classes ------------------------
# define as many classes we want that are not classified in other classes

class FitsObj():
	def __init__(self, wave, flux, error, 
				ra=None, dec=None, z_est=None, flag=None):
		self.wave = wave
		self.flux = flux
		self.error = error
		self.ra = ra
		self.dec = dec
		self.z_est = z_est
		self.flag = flag