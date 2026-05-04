# Project Documentation
[Back to Main Page](../main_readme.md)

*Auto-generated documentation from docstrings*

## Modules

### lens_ang_sep

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
  coordinates to a new source-plane redshift

Example
-------
    from rbcodes.lensing import lens_ang_sep as l
    s = l.model_delens(fits_x, fits_y, zl=0.392, zs_o=2.760)
    s.raytrace_new_z(ra, dec, zs_n=0.5)
    print(s.src_ra, s.src_dec)

Originally written by Ahmed Shaban (2020); extended by Rongmon Bordoloi (2021).

## Functions

### cosmic_D() (`lens_ang_sep`)

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

### ang_D12() (`lens_ang_sep`)

w_m: Matter density parameter. (Baryonic+Dark Matters)
    w_l: Dark energy density parameter. 
    z_1: redshift of the near object. (e.g.: lens redshift)
    z_2: redshift of the far object. (e.g.: background source redshift)

### ang_sep_D() (`lens_ang_sep`)

xc1: RA of the first point.
    yc1: Dec of the first point.
    
    xc2: RA of the second point.
    yc2: Dec of the second point.
    
    Ouput:
    The output is the angular separtion between the two points in radians

### __init__() (`lens_ang_sep`)

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

### raytrace_new_z() (`lens_ang_sep`)

Raytrace a set of input RA DECs to their corresponding source plane at the given redshift
           ralist: list of right ascension of the object at the image plane
           declist: list of declination of the object at the image plane.
           zs_n: redshift of the new source or plane at which you want to calculate the separation
    
           Output:
            returns the RA and Dec of the source plane location or the given redshift plane

