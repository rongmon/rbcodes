from astropy.cosmology import Planck18_arXiv_v2 as cosmo
import numpy as np
import matplotlib.pyplot as plt
import  astropy.units as u


def lens_sep_to_kpc(delta_arcsec,zabs_list,z_lens,z_source):
        #------------------------------------------------------------------------------------------
    #   Function to compute physical separation between sightlines in a lensed quasar system.
    #   Input:- 
    #           delta_arcsec      :- Angular separation between two quasars in arcsecond
    #           zabs_list         :- List of absorber redshifts for which we want to compute physical separation
    #           z_lens            :- Lens galaxy redshift
    #           z_source          :- background Quasar redshift
    #
    #
    #   Output:- 
    #           distlist          :- numpy array with physical Separation for each absorber redshift in kpc
    #
    #   Example:- 
    #                >from IGM import lens_sep_to_kpc as l
    #                > delta_arcsec=1. # 1 arcsecond separation
    #                > zabs_list =[.2, .5,1.2,2.]
    #                >z_lens =0.55
    #                >z_source =3.5
    #                > out = l.lens_sep_to_kpc(delta_arcsec,zabs_list,z_lens,z_source)
    #
    #   Written :- Rongmon Bordoloi                           March 2 2021
    #
    #   Equation used is the equation (5) from Cooke et al 2010.
    #   [Cooke, R., Pettini, M., Steidel, C. C., et al. 2010, MNRAS, 409, 679]
    #   Uses Planck 2018 \LambdaCDM cosmology
    #
    #
    #------------------------------------------------------------------------------------------
   



    zabs_list=np.array(zabs_list)
    # Convert angular separation to radian
    theta_obs= np.deg2rad(delta_arcsec/3600.)
    arcsec2kpc= cosmo.arcsec_per_kpc_proper(zabs_list)


    Ds= cosmo.comoving_distance(z_source)
    DL= cosmo.comoving_distance(z_lens)
    Dc= cosmo.comoving_distance(zabs_list)

    dist= (theta_obs)*DL *(Ds -Dc)/((1.+z_source)* (Ds-DL))
    distlist=dist.to(u.kpc).value

    q=np.where((zabs_list <= z_lens))
    distlist[q]=delta_arcsec/arcsec2kpc.value[q]

    return distlist

