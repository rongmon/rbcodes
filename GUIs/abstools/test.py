
from linetools.spectra.xspectrum1d import XSpectrum1D
from rbcodes.GUIs.abstools import Absorber as A
from rbcodes.GUIs.abstools import Metal_Plot as M

# Read spectrum file
spectrum_file = '/Users/bordoloi/WORK/python/rbcodes/example-data/test.fits'
sp = XSpectrum1D.from_file(spectrum_file)
wave = sp.wavelength.value
flux = sp.flux.value
error = sp.sig.value

# Analysis parameters
z = 0.348
lines = [1031.93, 1037.62]
window_lim = [-2000.0, 2000.0]
mask_init = [-200.0, 200.0]

# Create absorber and prepare for analysis
absys = A.Absorber(z, wave, flux, error, lines=lines, window_lim=window_lim, mask_init=mask_init)
Abs = absys.ions
intervening_file = '/Users/bordoloi/WORK/python/rbcodes/GUIs/abstools/SpecPlot_Projects/Identified_Linelist.txt'
M.Transitions(Abs, intervening=intervening_file)
