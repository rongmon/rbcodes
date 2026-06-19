import pandas as pd
from rbcodes.IGM.rb_setline import read_line_list


# ------------- Utility functions ----------------------

def clear_artists(ax, keep_lines=None):
	"""Clear texts and collections from a matplotlib axes.

	Parameters
	----------
	ax : matplotlib Axes
	keep_lines : int or None
	    If None (default), lines are not touched.
	    If int, lines beyond that index are removed (e.g. keep_lines=2
	    preserves the base flux+error lines; keep_lines=0 removes all).
	"""
	while ax.texts:
		ax.texts[-1].remove()
	while ax.collections:
		ax.collections[-1].remove()
	if keep_lines is not None:
		while len(ax.lines) > keep_lines:
			ax.lines[-1].remove()


def get_linelist_df(linelist_name):
	"""Load a linelist by name and return a DataFrame with columns ['wave', 'name'].

	Appends the rest wavelength to the ion name if it is not already present.
	"""
	llist = pd.DataFrame(columns=['wave', 'name'])
	tmp = read_line_list(linelist_name)

	if any(map(str.isdigit, tmp[1]['ion'])):
		# name column already contains the rest wavelength
		rows = [{'wave': li['wrest'], 'name': li['ion']} for li in tmp]
	else:
		# append rest wavelength to name
		rows = [{'wave': li['wrest'], 'name': li['ion'] + ' ' + str(round(li['wrest']))} for li in tmp]

	return pd.concat([llist, pd.DataFrame(rows)], ignore_index=True)


# ------------- Utility classes ------------------------
# define as many classes we want that are not classified in other classes

class FitsObj():
	'''Main object to save and traverse data between widgets
	'''
	def __init__(self, wave=None, flux=None, error=None,
				ra=None, dec=None, z_est=None, z_guess=None,
				flag=None, flux2d=None, error2d=None, frame_sources=None):
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

		# Detected frame sources from FITS headers (e.g., ['EMLINEA', 'EMLINEB'])
		self.frame_sources = frame_sources if frame_sources is not None else []

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