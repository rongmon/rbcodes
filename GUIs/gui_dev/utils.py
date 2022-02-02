import sys


# ------------- Utility classes ------------------------
# define as many classes we want that are not classified in other classes

class FitsObj():
	def __init__(self, wave, flux, error):
		self.wave = wave
		self.flux = flux
		self.error = error