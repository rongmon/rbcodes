from astropy.cosmology import Planck18 as cosmo
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from typing import Union, List, Optional

def lens_sep_to_kpc(delta_arcsec, zabs_list, z_lens, z_source, 
                   custom_cosmo=None, return_with_units=False):
    """
    Compute physical separation between sightlines in a lensed quasar system.
    
    Parameters
    ----------
    delta_arcsec : float
        Angular separation between two quasars in arcsecond
    zabs_list : array-like
        List or array of absorber redshifts for which to compute physical separation
    z_lens : float
        Lens galaxy redshift
    z_source : float
        Background Quasar redshift
    custom_cosmo : astropy.cosmology object, optional
        Custom cosmology to use. Default is Planck18.
    return_with_units : bool, optional
        If True, returns result with astropy units attached. Default is False.
    
    Returns
    -------
    distlist : numpy.ndarray or astropy.units.Quantity
        Physical Separation for each absorber redshift in kpc
    
    Notes
    -----
    Equation used is equation (5) from Cooke et al. 2010.
    [Cooke, R., Pettini, M., Steidel, C. C., et al. 2010, MNRAS, 409, 679]
    Uses Planck 2018 LambdaCDM cosmology by default.
    
    Examples
    --------
    >>> from rbcodes.IGM import lens_sep_to_kpc as l
    >>> delta_arcsec = 1.  # 1 arcsecond separation
    >>> zabs_list = [.2, .5, 1.2, 2.]
    >>> z_lens = 0.55
    >>> z_source = 3.5
    >>> out = l.lens_sep_to_kpc(delta_arcsec, zabs_list, z_lens, z_source)
    ------------
    Written :- 
        Rongmon Bordoloi                               March 2 2021
        Updated for better error handling and units    April 2025 

    """
    # Input validation
    if not np.isscalar(delta_arcsec) or delta_arcsec <= 0:
        raise ValueError("delta_arcsec must be a positive number")
    
    if z_lens < 0 or z_source < 0:
        raise ValueError("Redshifts must be non-negative")
        
    if z_source <= z_lens:
        raise ValueError("Source redshift must be greater than lens redshift")
    
    # Use provided cosmology or default to Planck18
    cosmology = custom_cosmo if custom_cosmo is not None else cosmo
    
    # Ensure zabs_list is a numpy array for vectorized operations
    zabs_list = np.asarray(zabs_list)
    
    # Validate absorber redshifts
    if np.any(zabs_list < 0):
        raise ValueError("Absorber redshifts must be non-negative")
    
    # Convert angular separation to radian
    theta_obs = np.deg2rad(delta_arcsec/3600.)
    arcsec2kpc = cosmology.arcsec_per_kpc_proper(zabs_list)

    # Get comoving distances
    Ds = cosmology.comoving_distance(z_source)
    DL = cosmology.comoving_distance(z_lens)
    Dc = cosmology.comoving_distance(zabs_list)

    # Apply equation from Cooke et al. 2010
    dist = (theta_obs)*DL*(Ds-Dc)/((1.+zabs_list)*(Ds-DL))
    
    # Convert to kpc and remove units for backward compatibility
    distlist = dist.to(u.kpc)
    
    # For absorbers at redshifts <= lens redshift, use angular diameter distance
    # Vectorized operation instead of where loop
    mask = (zabs_list <= z_lens)
    if np.any(mask):
        dist_values = distlist.value  # Extract values for modification
        dist_values[mask] = delta_arcsec/arcsec2kpc.value[mask]
        distlist = dist_values * u.kpc  # Reattach units
    
    # Return with or without units based on parameter
    if return_with_units:
        return distlist
    else:
        return distlist.value