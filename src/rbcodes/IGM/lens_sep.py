"""
Example code to plot to plot distances for differetn lens separations
"""
from astropy.cosmology import Planck18 as cosmo
import numpy as np
import matplotlib.pyplot as plt
import  astropy.units as u


def lens_sep(zlist=np.arange(0.8,3.,.05)):
    """
    Example code to plot to plot distances for differetn lens separations

    Parameters

        zlist= list of redshifts

    Returns

        Plot of physical separation vs redshift
    """
    


    arcsec2kpc= cosmo.arcsec_per_kpc_proper(zlist)
    
    
    z_L=0.493
    z_s = 2.7434
    z_c =1.62650
    
    
    theta_obs=1./3600.
    
    Ds= cosmo.comoving_distance(z_s)
    DL= cosmo.comoving_distance(z_L)
    Dc= cosmo.comoving_distance(z_c)
    
    dist= np.deg2rad(theta_obs)*DL *(Ds -Dc)/((1.+z_c)* (Ds-DL))
    
    print(dist.to(u.kpc))
    
    zabs=np.arange(.3,2.8,.01)
    plt.plot(zabs,1./cosmo.arcsec_per_kpc_proper(zabs),'k.')
    cosmo.arcsec_per_kpc_proper(zabs)
    for i in range(0,len(zlist)):
        z_L=zlist[i]
        z_s = max(zlist)
        q=np.where((zabs>z_L))
        z_c =zabs[q]
    
        if len(z_c)>0:
            Ds= cosmo.comoving_distance(z_s)
            DL= cosmo.comoving_distance(z_L)
            Dc= cosmo.comoving_distance(z_c)
    
            dist= np.deg2rad(theta_obs)*DL *(Ds -Dc)/((1.+z_c)* (Ds-DL))
    
    
            plt.plot(z_c,dist.to(u.kpc),'-')
    
    plt.show() 