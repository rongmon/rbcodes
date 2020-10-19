
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pdb
import sys
import os

from linetools.spectra.xspectrum1d import XSpectrum1D  
from GUIs import rb_plot_spec as r



if __name__ == "__main__":
    # Get the filename of the spectrum from the command line, and plot it
    filename = sys.argv[1]

    #Check if filetype is specified if not try to take what is given as extention to the file

    if (len(sys.argv) >2):
        filetype=sys.argv[2]
    else:
        #Take File Extention and try to match it
        tt=os.path.splitext(filename)[1]
        if (tt=='txt')| (tt=='dat'):
            filetype='ascii'
        elif (tt=='fits'):
        	filetype='linetools'
        else:
            filetype=tt[1:len(tt)]
    cwd=os.getcwd()
    print(cwd+'/'+filename)

    if (len(sys.argv) >3):
        efil=sys.argv[3]

    # Read in Files in differet formats
    if filetype=='ascii':
        from astropy.io import ascii
        dat=ascii.read(cwd+'/'+filename)
        tab=dat.keys()
        wave=np.array(dat[tab[0]])
        flux=np.array(dat[tab[1]])
        if (len(dat.keys())>=3):
            error=dat[tab[2]]
    elif filetype=='fits':
        from astropy.io import fits
        file=fits.open(cwd+'/'+filename)
        dat=file[1].data
        tab=dat.names
        wave=np.array(dat['wave'][0])
        flux=np.array(dat['flux'][0])
        if (len(tab)>=3):
            error=np.array(dat['error'][0])
    if filetype=='xfits':
        from astropy.io import fits
        hdu = fits.open(filename)
        wave = hdu['wavelength'].data
        flux = hdu['flux'].data
        error=hdu['error'].data
        hdu.close()
    elif filetype=='p':
        import pickle
        dat=pickle.load( open(cwd+'/'+filename, "rb" ))
        tab=dat.keys()
        wave=np.array(dat['wave'])
        flux=np.array(dat['flux'])
        if (len(tab)>=3):
            error=np.array(dat['error'])
    #Use linetools.io.readspec to read file
    elif filetype =='linetools':
        #from linetools.spectra import io as tio
        if (len(sys.argv) >3):
            sp=XSpectrum1D.from_file(filename,efil=efil)
        else:
            sp=XSpectrum1D.from_file(filename)

        #sp=tio.readspec(filename,inflg=None, efil=efil,**kwargs)
        wave=sp.wavelength.value
        flux=sp.flux.value

        if sp.sig_is_set == False:
            print('Assuiming arbiarbitrary 10% error on flux')
            error=0.1*flux
        else:
            error=sp.sig.value




    #sp=XSpectrum1D.from_file('PG0832+251_nbin3_coadd.fits') 
    r.rb_plot_spec(wave,flux/np.nanmedian(flux),error/np.nanmedian(flux))
    plt.show() # show the window
 