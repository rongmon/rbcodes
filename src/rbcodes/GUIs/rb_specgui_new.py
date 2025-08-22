import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
import pdb
import sys
import os
from astropy.table import Table

from linetools.spectra.xspectrum1d import XSpectrum1D  
# Deprecated: rb_plot_spec removed in v2.0.0
from GUIs import PlotSpec_Integrated as r

if __name__ == "__main__":
    # Get the filename of the spectrum from the command line, and plot it
    filename = sys.argv[1]
    wave = None
    flux = None
    error = None

    #Check if filetype is specified if not try to take what is given as extention to the file
    if (len(sys.argv) > 2):
        filetype = sys.argv[2]
    else:
        #Take File Extention and try to match it
        tt = os.path.splitext(filename)[1]
        if (tt == 'txt') or (tt == 'dat'):
            filetype = 'ascii'
        elif (tt == 'fits'):
            filetype = 'linetools'
        else:
            filetype = tt[1:len(tt)]
    
    cwd = os.getcwd()
    print(cwd + '/' + filename)

    efil = None
    if (len(sys.argv) > 3):
        efil = sys.argv[3]

    # Read in Files in differet formats
    if filetype == 'ascii':
        from astropy.io import ascii
        try:
            dat = ascii.read(cwd + '/' + filename)
            if isinstance(dat, Table):
                cols = dat.colnames
                if len(cols) >= 2:
                    wave = np.array(dat[cols[0]])
                    flux = np.array(dat[cols[1]])
                    if len(cols) >= 3:
                        error = np.array(dat[cols[2]])
        except Exception as e:
            print(f"Error reading ASCII file: {str(e)}")
            sys.exit(1)
            
    elif filetype == 'fits' or filetype == 'lt' or filetype == 'linetools':
        #Use linetools.io.readspec to read file
        try:
            if efil is not None:
                sp = XSpectrum1D.from_file(filename, efil=efil)
            else:
                sp = XSpectrum1D.from_file(filename)

            wave = sp.wavelength.value
            flux = sp.flux.value

            if not sp.sig_is_set:
                print('Assuming arbitrary 10% error on flux')
                error = 0.1 * flux
            else:
                error = sp.sig.value
        except Exception as e:
            print(f"Error reading FITS file: {str(e)}")
            sys.exit(1)

    elif filetype == 'lt_cont_norm':
        print('Using filetype:' + filetype)

        try:
            if efil is not None:
                sp = XSpectrum1D.from_file(filename, efil=efil)
            else:
                sp = XSpectrum1D.from_file(filename)

            wave = sp.wavelength.value

            if sp.co_is_set:
                print('Using provided continuum to normalize')
                flux = sp.flux.value / sp.co.value
            else:
                flux = sp.flux.value

            if not sp.sig_is_set:
                print('Assuming arbitrary 10% error on flux')
                error = 0.1 * flux
            else:
                if sp.co_is_set:
                    error = sp.sig.value / sp.co.value
                else:
                    error = sp.sig.value
        except Exception as e:
            print(f"Error reading continuum normalized file: {str(e)}")
            sys.exit(1)

    elif filetype == 'p':
        import pickle
        try:
            with open(cwd + '/' + filename, "rb") as f:
                dat = pickle.load(f)
            if isinstance(dat, dict):
                if 'wave' in dat and 'flux' in dat:
                    wave = np.array(dat['wave'])
                    flux = np.array(dat['flux'])
                    if 'error' in dat:
                        error = np.array(dat['error'])
        except Exception as e:
            print(f"Error reading pickle file: {str(e)}")
            sys.exit(1)

    if wave is None or flux is None:
        print("Error: Failed to load wave or flux data")
        sys.exit(1)

    # Handle flux normalization
    if np.all(np.isnan(flux)):
        print("Error: All flux values are NaN")
        sys.exit(1)

    cnt = np.nanmedian(flux)
    if cnt == 0:
        cnt = np.nanmean(flux)
        if cnt == 0:
            cnt = 1.0

    # Set default error if none provided
    if error is None:
        error = np.zeros_like(flux)

    r.rb_plotspec(wave, flux/cnt, error/cnt)
    plt.show() # show the window