'''
from astropy.io import fits
import matplotlib.pyplot as plt
file = fits.open(filename)
dat = file[1].data

plt.plot(dat['wave'], dat['flux'])
plt.plot(dat['wave'], dat['cont'])
plt.show()
'''


filename = input('filename??   ')
#linetools
import matplotlib
matplotlib.use('TkAgg') #python 3 only
from linetools.spectra.io import readspec

spec = readspec(filename+'.fits')
spec.fit_continuum()

# check
import matplotlib.pyplot as plt
plt.plot(spec.wavelength, spec.flux)
plt.plot(spec.wavelength, spec.co, 'o') # continuum
plt.show()

save_spec = input('Do you want to save the continuum?  ')
if save_spec == 'y':
    spec.write_to_fits(filename.split('_')[0]+'.fits')
'''
# once checked, save to fits file
# data saved as UPPERCASE letters...
# in order of: FLUX, ERROR, WAVELENGTH, CONTINUUM
spec.write_to_fits('output_filepath')

# read
from astropy.io import fits
import matplotlib.pyplot as plt
spec = fits.open('output_filepath')
flux = spec[0].data
error = spec[1].data
wave = spec[2].data
cont = spec[3].data
'''
