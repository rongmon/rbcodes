from astropy.coordinates import FK5, SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.io import fits

def convert_epoch(filename):
	''''--------------------------------------------
	Function to take Magellan FIRE spectrum and use the header information
	to transform the co-ordinates to J2000 epoch

	Input:  Magellan FIRE 1d spectrum filename in fits format

	Return:
	      astropy SkyCoord object in J2000 epoch



	 Calling Sequence:
	   from catalog.convert_FIRE_coordinates import convert_epoch
	   s=convert_epoch('J0100+28_F.fits')

	   #show results
	   s.to_string('hmsdms')


	Note: Only for specific use of older Magellan FIRE QSO spectra.


	Written By: Rongmon Bordoloi              May 2022
       ---------------------------------------------
	'''
	hdu=fits.open(filename)
	equinox     =  hdu[0].header['EQUINOX']
	ra_deg1      =  hdu[0].header['RA']
	dec_deg1     =  hdu[0].header['DEC']

	fk5_obs = FK5(equinox=Time(equinox, format='jyear'))
	fk5_2000 = FK5(equinox=Time(2000.0, format='jyear'))
	c = SkyCoord(ra_deg1,dec_deg1,frame=fk5_obs,unit=(u.deg,u.deg))
	c_j2000 = c.transform_to(fk5_2000)
	#ra_deg = c_j2000.ra.deg
	#dec_deg = c_j2000.dec.deg

	return c_j2000