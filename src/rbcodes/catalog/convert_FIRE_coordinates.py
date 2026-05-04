"""
Coordinate epoch conversion for Magellan/FIRE spectra.

Standalone module — not imported by other package modules; available for direct use.

Reads the EQUINOX keyword from a Magellan FIRE 1D spectrum FITS header and
transforms the stored (RA, Dec) to J2000 (FK5) coordinates using Astropy.

Note: Intended specifically for older Magellan FIRE QSO spectra that store
coordinates in a non-J2000 epoch.

Example
-------
    from rbcodes.catalog.convert_FIRE_coordinates import convert_epoch
    coord_j2000 = convert_epoch('J0100+28_F.fits')
    print(coord_j2000.to_string('hmsdms'))
"""
from astropy.coordinates import FK5, SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.io import fits

def convert_epoch(filename):
	"""Function to take Magellan FIRE spectrum and use the header information
	to transform the co-ordinates to J2000 epoch

	Parameters
	-----------

	 Magellan FIRE 1d spectrum filename in fits format

	Return
	------
	
	astropy SkyCoord object in J2000 epoch



	Example
	-------

	from catalog.convert_FIRE_coordinates import convert_epoch
	s=convert_epoch('J0100+28_F.fits')

	#show results
	s.to_string('hmsdms')


	Note: Only for specific use of older Magellan FIRE QSO spectra.


	Written By: Rongmon Bordoloi              May 2022
       ---------------------------------------------
	"""
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