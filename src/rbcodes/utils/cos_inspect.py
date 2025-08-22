import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def cos_inspect(path='', filelist=None, xrange=None, yrange=None, pages=1):
    plt.rcParams['lines.linewidth'] = 0.5

    if not xrange:
        xrange = [1130, 1800]
    if not path:
        path = ''
    if not filelist:
        filelist = glob.glob(path + '*x1d*.fit*')
    else:
        filelist = [path + file for file in filelist]

    nfiles = len(filelist)
    if not pages:
        pages = 1

    plt.figure(figsize=(10, 6))

    for i, file in enumerate(filelist):
        # read header parameters
        print('Reading in: '+file)
        with fits.open(file, memmap=False) as hdulist:
            hdr0 = hdulist[0].header
            hdr1 = hdulist[1].header
            x1ddata=hdulist[1].data

            grating = hdr0['OPT_ELEM']
            cenwave = hdr0['CENWAVE']
            exptime = hdulist[1].data['EXPTIME']
            dateobs = hdr1['DATE-OBS']

            segments=hdulist[1].data['segment']

            try:
                x1d_fuva = x1ddata[x1ddata['segment'] == 'FUVA']
                wavelength = x1d_fuva['wavelength'].flatten()
                flux = x1d_fuva['flux'].flatten()
                err=x1d_fuva['error'].flatten()
                plt.plot(wavelength, flux,color='blue')
                plt.plot(wavelength, err,color='red')

            except:
                print("No FUVA segment found")

            try:
                x1d_fuvb = x1ddata[x1ddata['segment'] == 'FUVB']
                wavelength = x1d_fuvb['wavelength'].flatten()
                flux = x1d_fuvb['flux'].flatten()
                err=x1d_fuvb['error'].flatten()

                plt.plot(wavelength, flux,color='gray')
                plt.plot(wavelength, err,color='red')

            except:
                print("No FUVB segment found")







            if not yrange:
                yrange = [-0.01e-13, np.max(flux)]


            plt.xlim(xrange)
            plt.ylim(yrange)
            #plt.title(f"{hdr0['TARGNAME'].strip()[:2]} {hdr0['ROOTNAME'].strip()[:2]} {dateobs} {grating}{cenwave}{exptime:.2f} sec")

            # print summary information
            print('------------------------------------------------------------')
            print(f"{hdr0['TARGNAME']} Observation Summary")
        
            g130m = np.where(grating == 'G130M')[0]
            g160m = np.where(grating == 'G160M')[0]
            g140l = np.where(grating == 'G140L')[0]
        
            if len(g130m) > 0:
                print(f"  G130M: total exptime = {np.sum(exptime[g130m]):.2f}")
                print(f"         # exposures          {len(g130m)}")
                print("  grating positions "+grating)
            else:
                print("  G130M: NONE")
        
            if len(g160m) > 0:
                print(f"  G160M: total exptime = {np.sum(exptime[g160m]):.2f}")
                print(f"         # exposures          {len(g160m)}")
                print("  grating positions "+grating)
            else:
                print("  G160M: NONE")
        
            if len(g140l) > 0:
                print(f"  G140L: total exptime = {np.sum(exptime[g140l]):.2f}")
                print(f"         # exposures          {len(g140l)}")
                print("  grating positions "+grating)
            else:
                print("  G140L: NONE")
        
            print('------------------------------------------------------------')

    plt.tight_layout()
    plt.show()


