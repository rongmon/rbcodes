""" Image generation and manipulation routines"""
from __future__ import print_function, absolute_import, division, unicode_literals
from astropy.io import fits
from astropy.wcs import WCS
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
    Dc = DH * np.trapz( (1. / E_z) , x=zd )     # line of sight comoving distance

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
    
    d  = np.arccos(np.sin(np.radians(yc1)) * np.sin(np.radians(yc2)) + np.cos(np.radians(yc1)) * np.cos(np.radians(yc2)) * np.cos(np.radians( xc1 - xc2 )))
    # d is given in radians
    return d


class model_delens(object):

    def __init__(self,fits_x,fits_y, zl= 0.564,zs_o=1.701,**kwargs):
        """
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
    
        self.hdl_x = fits.open(fits_x)[0]
        self.hdl_y = fits.open(fits_y)[0]
    
        self.x_def = self.hdl_x.data    # x coordinate deflection matrix
        self.y_def = self.hdl_y.data    # y coordinate deflection matrix
    
        self.wcs_x = WCS(self.hdl_x.header).celestial
        self.wcs_y = WCS(self.hdl_y.header).celestial
    
        #Define pixel scales from the deflection matrix
    
        pix_scale_arcsec_x=np.abs(self.hdl_x.header['CDELT1'])*3600.
        pix_scale_arcsec_y=np.abs(self.hdl_y.header['CDELT1'])*3600.
        
        self.pix_scale_x = pix_scale_arcsec_x * (1.0/3600.0)
        self.pix_scale_y = pix_scale_arcsec_y * (1.0/3600.0)

    def update_cosmology(self,Omega_m=0.3,Omega_l=0.7):
        # Update the default cosmology to a different Omega_Lambda and Omega_matter value
        self.w_m=Omega_m
        self.w_l=Omega_l

    def raytrace_new_z(self,ralist,declist,zs_n=1.701):

        """
           Raytrace a set of input RA DECs to their corresponding source plane at the given redshift
           ralist: list of right ascension of the object at the image plane
           declist: list of declination of the object at the image plane.
           zs_n: redshift of the new source or plane at which you want to calculate the separation
    
           Output:
            returns the RA and Dec of the source plane location or the given redshift plane
        """
        ra=np.array(ralist)
        dec=np.array(declist)

        out_so = cosmic_D(self.w_m, self.w_l, self.zs_o); D_so = out_so['D_A']
        out_sn = cosmic_D(self.w_m, self.w_l, zs_n); D_sn = out_sn['D_A']

        D_A_ls_o = ang_D12(self.w_m, self.w_l, self.zl, self.zs_o) # Angular diameter distance between the lens and the source
        D_A_ls_n = ang_D12(self.w_m, self.w_l, self.zl, zs_n) # Angular diameter distance between the lens and the new redshift plane

        theta_x, theta_y = self.wcs_x.world_to_pixel_values(ra, dec)

        # theta_x, theta_y: are the pixel poitions in the image plane
        #print(theta_x)
        #print(len(theta_x))
        theta_x = np.asarray(theta_x); theta_y = np.asarray(theta_y)

        # alpha_x, alpha_y: are the deflection in the x or y direction (deflection in pixels)
        alpha_x, alpha_y = [], []

        for k in range(len(theta_x)):
            a_x = self.x_def[int(theta_y[k]), int(theta_x[k])] * (1.0/3600.0)
            a_y = self.y_def[int(theta_y[k]), int(theta_x[k])] * (1.0/3600.0)

            apix_x = a_x / self.pix_scale_x;  alpha_x.append(apix_x)
            apix_y = a_y / self.pix_scale_y;  alpha_y.append(apix_y)

        #a_pix_x = np.asarray(a_pix_x); a_pix_y = np.asarray(a_pix_y)
        alpha_x = np.asarray(alpha_x); alpha_y = np.asarray(alpha_y)

        src_px = theta_x - alpha_x * ( (D_A_ls_n/D_sn) / (D_A_ls_o/D_so) )
        src_py = theta_y - alpha_y * ( (D_A_ls_n/D_sn) / (D_A_ls_o/D_so) )

        self.src_ra, self.src_dec = self.wcs_x.pixel_to_world_values(src_px, src_py)















    




