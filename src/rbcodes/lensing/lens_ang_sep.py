"""
Gravitational lensing ray-tracing and angular separation utilities.

Standalone module — not imported by other package modules; available for direct use.

Provides tools to ray-trace image-plane coordinates to a source plane (or any
intermediate redshift plane) using deflection matrices from a lens model, along
with supporting cosmological distance calculations.

Key components:
- ``cosmic_D(w_m, w_l, z)``             — comoving, angular diameter, and luminosity distances
- ``ang_D12(w_m, w_l, z1, z2)``         — angular diameter distance between two redshifts
- ``ang_sep_D(xc1, yc1, xc2, yc2)``    — angular separation in radians
- ``model_delens``                      — class to load lens deflection matrices and ray-trace
  coordinates to a new source-plane redshift (nearest-pixel or bilinear lookup of the
  deflection; positions outside the maps give NaN)

Example
-------
    from rbcodes.lensing import lens_ang_sep as l
    s = l.model_delens(fits_x, fits_y, zl=0.392, zs_o=2.760)
    src_ra, src_dec = s.raytrace_new_z(ra, dec, zs_n=0.5)                    # nearest pixel
    src_ra, src_dec = s.raytrace_new_z(ra, dec, zs_n=0.5, interp='bilinear')

Originally written by Ahmed Shaban (2020); extended by Rongmon Bordoloi (2021).
"""
from __future__ import print_function, absolute_import, division, unicode_literals
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import map_coordinates
import matplotlib.pyplot as plt
import numpy as np

###    Compute the different cosmological distances based on Hogg 2000 paper
def cosmic_D(w_m, w_l, z):
    """
    The function takes the redshift and your deinfed cosmology and gives you the different cosmological distances in meters

    Input Parameters:
    w_m : Matter Density Parameter       (dimensionless)
    w_l : Dark Energy Density Parameter  (dimensionless)
    z   : is the redshift of the object

    Other Variables in the calculation:
    w_k : Curvature Parameter            (dimensionless)
    E_z : Defined function of the redshift
    Dc  : line of sight comoving distance
    Dm  : is the transverse comoving distance
    D_A : is the angular diameter distance of the object (will be given in meters)
    D_L : is the lumonsity distance

    Constants
    h  : Dimensionless Parameter
    H0 : Current Value of Hubble Parameter measured in km s^{−1} Mpc^{−1}
    DH : Hubble Distance in meters    # DH = c/H0 = 3000 * (1/h) # in Mpc
    """

    # Defining our constants
    h = 0.7                        # Dimensionless Parameter
    H0 = h * 100.                # Current Value of Hubble Parameter measured in km s^{−1} Mpc^{−1}
    DH = 9.26 * 10**(25) * (1/h)   # Hubble Distance in meters    # DH = c/H0 = 3000 * (1/h) # in Mpc


    w_k = 1. - w_m - w_l


    zd = np.linspace(0., z, 1000)
    E_z = np.sqrt( w_m * ((1. + zd)**3.) + w_k * ((1. + zd)**2)  + w_l)
    Dc = DH * np.trapezoid( (1. / E_z) , x=zd )     # line of sight comoving distance

    # Defining the trnasverse comoving distance
    #print(w_k)
    if w_k > 0.0:
        Dm = (DH / np.sqrt(w_k) * np.sinh( (np.sqrt(w_k) * Dc )/ DH ) )   # for w_k > 0.0
    elif w_k == 0.0:
        Dm = Dc     # for w_k = 0.0
    else:
        Dm = (DH / np.sqrt( np.abs(w_k)) * np.sinh( (np.sqrt(np.abs(w_k)) * Dc )/ DH ) )  # for w_k < 0.0


    D_A = Dm / (1. + z)
    
    D_L = (1. + z) * Dm
    
    output = {} 
    output['D_M'] = Dm 
    output['D_C'] = Dc 
    output['D_A'] = D_A
    output['D_L'] = D_L

    return output


def ang_D12(w_m, w_l, z1, z2):
    """
    w_m: Matter density parameter. (Baryonic+Dark Matters)
    w_l: Dark energy density parameter. 
    z_1: redshift of the near object. (e.g.: lens redshift)
    z_2: redshift of the far object. (e.g.: background source redshift)
    """
    # Defining our constants
    h = 0.7                        # Dimensionless Parameter
    H0 = h * 100.                # Current Value of Hubble Parameter measured in km s^{−1} Mpc^{−1}
    D_H = 9.26 * 10**(25) * (1/h)   # Hubble Distance in meters    # DH = c/H0 = 3000 * (1/h) # in Mpc

    
    w_k = 1. -  w_m - w_l
    out1 = cosmic_D(w_m, w_l, z1)
    out2 = cosmic_D(w_m, w_l, z2)
    
    D_M1 = out1['D_M']
    D_M2 = out2['D_M']
    
    D_A12 = (1. / (1. + z2)) * ( D_M2 * np.sqrt( 1. + w_k * (D_M1/D_H)**2. ) - D_M1 * np.sqrt( 1. + w_k * (D_M2/D_H)**2. ) )
    
    return D_A12

def ang_sep_D(xc1, yc1, xc2, yc2):
    """
    xc1: RA of the first point.
    yc1: Dec of the first point.
    
    xc2: RA of the second point.
    yc2: Dec of the second point.
    
    Ouput:
    The output is the angular separtion between the two points in radians
    """
    
    # Vincenty formula: accurate at all separations (the arccos form loses
    # precision below ~1 arcsec)
    d1, d2 = np.radians(yc1), np.radians(yc2)
    dra = np.radians(xc2 - xc1)
    num = np.hypot(np.cos(d2) * np.sin(dra),
                   np.cos(d1) * np.sin(d2) - np.sin(d1) * np.cos(d2) * np.cos(dra))
    den = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(dra)
    d = np.arctan2(num, den)
    # d is given in radians
    return d


class model_delens(object):

    def __init__(self,fits_x,fits_y, zl= 0.564,zs_o=1.701,**kwargs):
        r"""
        ---------------------------------------------------------------------
        Inititate the delensing object to load the required deflection matrices
        *****
        Note that this code uses a default 707 \Lambda CDM cosmology. 
        If your deflection matrices were created using a different comsology, 
        update using self.update_cosmology function. See example below.
        *****

        fits_x: Filepath+name of the x deflection matrix file
        fits_y: Filepath+name of the y deflection matrix file
        zs_o: Redshift of the lensed background source with which lens model was created.
        zl  : Lens redshift.
        ---------------------------------------------------------------------
        #EXAMPLE:

        #Set path to your deflection matrices
        path='/Users/bordoloi/Dropbox/Research/KCWI/SGAS1527/Lens_Model/'
        fits_x=path+'s1527_dplx_2.761.fits'
        fits_y=path+'s1527_dply_2.761.fits'


        # Array of coordinates from image plane that will be transported to source plane
        import numpy as np
        ra=np.array([231.937925,231.937803])
        dec=np.array([6.872142,6.872197])


        zBG_obj= 2.760  # The background lensed object redshift that was used to create the deflection matrices
        zl= 0.392 # Lens galaxy redshift

        zabs=0.5 Foreground redshift at which the ray tracing will be performed
        
        #Load the module
        from lensing import lens_ang_sep as l
        
        #Initialize the model
        s=l.model_delens(fits_x,fits_y,zl=zl,zs_o=zBG_obj)
        
        #Ray Trace to source plane at a given redshift
        s.raytrace_new_z(ra,dec,zs_n=zabs)

        #Print new coordinates
        print(s.src_ra, s.src_dec)

        # Print angular separation between source plane and image plane coordinates in arc sec

        print(l.ang_sep_D(ra,dec,s.src_ra,s.src_dec) * (180./np.pi) * 3600.)

        #------------- if you want to update cosmology 
        #Initialize the model
        s=l.model_delens(fits_x,fits_y,zl=zl,zs_o=zBG_obj)
        #Update cosmology 
        s.update_cosmology(Omega_m=0.25,Omega_l=0.75)

        #----------------------------------------------------

        Written By:
        
        Originally written By: Ahmed Shaban. 2020
        
        Edit: Rongmon Bordoloi.  March 2021 : Modified in a class format with flexible cosmology and file input formats








        """

        #For now assume a default 707 LambdaCDM cosmology
        self.w_m = 0.3; self.w_l = 0.7  # Lambda-CDM cosmology
        self.zl=zl
        self.zs_o=zs_o
    
        # read the data now so the files can be closed
        with fits.open(fits_x, memmap=False) as h:
            self.hdl_x = h[0]
            self.hdl_x.data
        with fits.open(fits_y, memmap=False) as h:
            self.hdl_y = h[0]
            self.hdl_y.data
    
        self.x_def = self.hdl_x.data    # x coordinate deflection matrix
        self.y_def = self.hdl_y.data    # y coordinate deflection matrix
    
        self.wcs_x = WCS(self.hdl_x.header).celestial
        self.wcs_y = WCS(self.hdl_y.header).celestial
    
        #Define pixel scales from the deflection matrix (degrees per pixel along
        #each axis; works for both CDELT and CD-matrix headers)
        self.pix_scale_x, self.pix_scale_y = proj_plane_pixel_scales(self.wcs_x)

    def update_cosmology(self,Omega_m=0.3,Omega_l=0.7):
        # Update the default cosmology to a different Omega_Lambda and Omega_matter value
        self.w_m=Omega_m
        self.w_l=Omega_l

    def deflection_scale(self, zs_n):
        """
        Factor by which the deflections (computed for zs_o) are rescaled for a plane
        at redshift zs_n: (D_ls/D_s)(zs_n) / (D_ls/D_s)(zs_o). It is 1 at zs_n = zs_o
        and 0 at the lens redshift.
        """
        D_so = cosmic_D(self.w_m, self.w_l, self.zs_o)['D_A']
        D_sn = cosmic_D(self.w_m, self.w_l, zs_n)['D_A']
        D_A_ls_o = ang_D12(self.w_m, self.w_l, self.zl, self.zs_o) # Angular diameter distance between the lens and the source
        D_A_ls_n = ang_D12(self.w_m, self.w_l, self.zl, zs_n) # Angular diameter distance between the lens and the new redshift plane
        return (D_A_ls_n/D_sn) / (D_A_ls_o/D_so)

    def deflection_at(self, theta_x, theta_y, interp='nearest'):
        """
        Deflection (in deflection-map pixels) at 0-based pixel positions of the
        deflection maps.

        interp: 'nearest' (value of the pixel containing the position) or
                'bilinear' (interpolated between the four nearest pixel centres).
        Positions outside the maps give NaN.
        """
        theta_x = np.asarray(theta_x, dtype=float)
        theta_y = np.asarray(theta_y, dtype=float)
        ny, nx = self.x_def.shape
        a_x = np.full(theta_x.shape, np.nan)
        a_y = np.full(theta_x.shape, np.nan)

        if interp == 'nearest':
            ix = np.rint(np.nan_to_num(theta_x, nan=-1.)).astype(int)
            iy = np.rint(np.nan_to_num(theta_y, nan=-1.)).astype(int)
            ok = np.isfinite(theta_x) & np.isfinite(theta_y) & \
                (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            a_x[ok] = self.x_def[iy[ok], ix[ok]]
            a_y[ok] = self.y_def[iy[ok], ix[ok]]
        elif interp == 'bilinear':
            ok = np.isfinite(theta_x) & np.isfinite(theta_y) & \
                (theta_x >= 0) & (theta_x <= nx - 1) & (theta_y >= 0) & (theta_y <= ny - 1)
            coords = np.vstack([theta_y[ok], theta_x[ok]])
            a_x[ok] = map_coordinates(self.x_def, coords, order=1, mode='nearest')
            a_y[ok] = map_coordinates(self.y_def, coords, order=1, mode='nearest')
        else:
            raise ValueError(f"interp must be 'nearest' or 'bilinear', got {interp!r}")

        # deflection maps are in arcsec; convert to pixels
        return a_x / 3600. / self.pix_scale_x, a_y / 3600. / self.pix_scale_y

    def raytrace_new_z(self,ralist,declist,zs_n=1.701,interp='nearest'):

        """
           Raytrace a set of input RA DECs to their corresponding source plane at the given redshift
           ralist: list of right ascension of the object at the image plane
           declist: list of declination of the object at the image plane.
           zs_n: redshift of the new source or plane at which you want to calculate the separation
           interp: 'nearest' (default) or 'bilinear' lookup of the deflection maps

           Output:
            returns the RA and Dec of the source plane location or the given redshift plane
            (also stored as self.src_ra, self.src_dec). Positions outside the
            deflection maps give NaN.
        """
        ra=np.asarray(ralist, dtype=float)
        dec=np.asarray(declist, dtype=float)

        # theta_x, theta_y: are the pixel positions in the image plane
        theta_x, theta_y = self.wcs_x.world_to_pixel_values(ra, dec)
        theta_x = np.asarray(theta_x); theta_y = np.asarray(theta_y)

        # alpha_x, alpha_y: are the deflection in the x or y direction (deflection in pixels)
        alpha_x, alpha_y = self.deflection_at(theta_x, theta_y, interp=interp)

        scale = self.deflection_scale(zs_n)
        src_px = theta_x - alpha_x * scale
        src_py = theta_y - alpha_y * scale

        self.src_ra, self.src_dec = self.wcs_x.pixel_to_world_values(src_px, src_py)
        return self.src_ra, self.src_dec
