""" Read in HST/COS x1d header files and give some info"""
import glob
from astropy.io import fits
'''
This function will read all fits files in a folder and print out the header info for raw spectra.
It will print filename, exposure type, and object type entries from the header.

Written By: Rongmon Bordoloi April 2018
'''
def print_header(filenames):
    files=glob.glob(filenames)
    for i in range(0,len(files)):
        hdul = fits.open(files[i])
        objectID=hdul[0].header['TARGNAME']
        filename=hdul[0].header['FILENAME']
        filetype=hdul[0].header['FILETYPE']
        exptime= str(hdul[1].header['EXPTIME'])
        print(files[i]+' : '+filetype+' : '+objectID+ ' : '+exptime)

if __name__ == "__main__":
        files=glob.glob('*_x1d*.fits')

        for i in range(0,len(files)):
                hdul = fits.open(files[i])
                objectID=hdul[0].header['TARGNAME']
                filename=hdul[0].header['FILENAME']
                filetype=hdul[0].header['FILETYPE']
                exptime= str(hdul[1].header['EXPTIME'])
                print(files[i]+' : '+filetype+' : '+objectID+ ' : '+exptime)
