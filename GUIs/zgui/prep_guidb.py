import sys, os
#import astropy.units as u
#from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
import pandas as pd


# Prepare and initialize a database dataFrame for GUI
# Query all fits files in the curret working directory and loop through them to create the following db
# assuming the database is a 2D dataFrame:
# | filename | RA | DEC | z | z_err | Confidence | Flag | Linelist | z_guess |
# only fill out <filename>, <RA>, <DEC>, and <z_guess>
# assuming this is done in the current directory which contains all FITS files
class GUI_DataFrame():
    def __init__(self, dirpath):
        self.dirpath = os.path.abspath(dirpath)
        self.gui_df = pd.DataFrame(
            columns=['Name', 'RA', 'DEC', 'z', 'z_err', 'Confidence', 'Linelist', 'Flag', 'z_guess'])
        self.fitsfiles = []

    def prepare_database(self):
        if os.path.exists(self.dirpath):
            dirpaths, dirnames, filenames = [],[],[]
            for dp, dn, fn in os.walk(self.dirpath):
                dirpaths.append(dp)
                dirnames.append(dn)
                filenames.append(fn)

            # loop through all FITS files
            for i in range(len(dirpaths)):
                for file in filenames[i]:
                    if file.endswith('fits'):
                        # get rid of file extension
                        filename = file[:-5]
                        filepath = os.path.abspath(dirpaths[i]+'/'+file)
                        hdul = fits.open(filepath)
                        labels = [label.name for label in hdul]
                        if 'SCI' in labels:
                            # This is only compatible with EIGER emission line 2d spec FITS file
                            if 'RA' in hdul['SCI'].header:
                                # assume RA and DEC are always saved in pair
                                ra = hdul['SCI'].header['RA']
                                dec = hdul['SCI'].header['DEC']
                            else:
                                ra, dec = np.nan, np.nan
                            if 'EAZY_ZPDF' in labels:
                                z = hdul['EAZY_ZPDF'].data['z']
                                pdf = hdul['EAZY_ZPDF'].data['pdf']
                                z_guess = z[np.argmax(pdf)]
                            else:
                                z_guess = np.nan
                        hdul.close()
                        # Update GUI database
                        self.gui_df = self.gui_df.append(
                            {'Name': filename, 
                            'RA': ra, 
                            'DEC': dec, 
                            'z': np.nan, 
                            'z_err': np.nan, 
                            'Confidence': '', 
                            'Linelist': np.nan, 
                            'Flag': '', 
                            'z_guess': z_guess}, ignore_index=True)

                        # Update GUI FITS dataset
                        self.fitsfiles.append(file)

        else:
            print('Abort...')
            print('Given path does not exist.')

    def print_database(self):
        print(self.gui_df)

    def save_database(self, fname='z_guess_guidb'):
        save_path = self.dirpath+'/'+fname + '.csv'
        self.gui_df.to_csv(save_path,
                            sep=',',
                            header=True,
                            index=False)
        print('GUI database is saved at ' + save_path)
        
        # sort FITS filenames before saving as TXT file
        # before sorting
        #print(self.fitsfiles)
        # after sorting
        self.fitsfiles.sort()
        #print(self.fitsfiles)
        with open(self.dirpath+'/'+'FITS_files.txt', 'w') as f:
            f.write('\n'.join(self.fitsfiles))
        print('A TXT file containing all FITS files is also created within the same folder.')



#save the database

if __name__ == '__main__':
    current_path = os.getcwd()
    print('Prepare GUI database using path: ' + current_path)
    gui_database = GUI_DataFrame(current_path)
    gui_database.prepare_database()
    gui_database.print_database()
    confirm_db = input('Save the current database?(y/n)  ')
    if 'y' in confirm_db.lower():
        save_name = input('Name this database (NO file extension needed):\n')
        if len(save_name) > 0:
            gui_database.save_database(save_name)
        else:
            print('You did not name current database. GUI database is named with default name.')
            gui_database.save_database()
    else:
        print('Current GUI database is not saved.')