# Project Documentation
[Back to Main Page](../main_readme.md)

*Auto-generated documentation from docstrings*

## Modules

### mstar2mhalo

Stellar-mass to halo-mass conversion (Moster et al. 2010).

Standalone module — not imported by other package modules; available for direct use.

Converts a galaxy stellar mass to a halo mass and NFW virial radius using
the redshift-dependent SMHM relation of Moster et al. (2010), ApJ 710, 903.
Valid for redshifts 0 < z < 3.5.

Example
-------
    from rbcodes.halo.mstar2mhalo import mstar2mhalo
    import numpy as np
    mhalo, rvir = mstar2mhalo(mstar=10**10.4, z=0.5)

### halos

Stellar-mass to halo-mass conversions via the Behroozi UniverseMachine relation.

Standalone module — not imported by other package modules; available for direct use.

Provides:
- ``haloMassAtStellarMass(z, stellar_mass, population=None, centrals=True)`` —
  log halo mass at fixed stellar mass, read directly from the UniverseMachine DR1
  mock catalogues (``um_dr1_hm_averages.npz``).  ``population`` selects all
  (default), star-forming (``'sf'``) or quenched/passive (``'q'``) galaxies.
  Recommended for assigning halo masses to observed galaxies.
- ``stellarToHaloMass(z, stellar_mass, population=None, centrals=False)`` — the
  inverted relation: inverts the Behroozi et al. (2019) fit to median M* at fixed
  halo mass (``smhm_med*_params.txt``, copied unchanged from UniverseMachine DR1).
  Biased high at log M* >~ 10.8.
- ``R_200(Mh, z)``                       — compute the virial radius R_200 from a
  log halo mass and redshift (uses Planck18 cosmology).
- ``V_200(Mh, z)``, ``concentration_200c(Mh, z)``, ``nfw_vesc(Mh, z, R)`` — circular
  velocity, Dutton & Maccio (2014) concentration, and NFW escape velocity at R.
- ``halo_properties(z, stellar_mass, method='direct', population=None, centrals=True,
  R=None)`` — log_Mh, R200, V200, c200 and (given R) R/R200 and v_esc in one call.

Note: scatter in the SMHM relation is ~0.1 dex (range 0.04–0.13 dex).

Example
-------
    from rbcodes.halo import halos as h
    import numpy as np
    mh_sf = h.haloMassAtStellarMass(0.43, 11.0, population='sf')   # direct
    mh_inv = h.stellarToHaloMass(0.43, 11.0)                       # inverted
    p = h.halo_properties(0.43, 11.0, population='sf', R=100.)
    p['log_Mh'], p['R200'], p['V200'], p['R_over_R200'], p['v_esc']

### rb_nfw

NFW dark matter halo profile calculator.

Standalone module — not imported by other package modules; available for direct use.

Computes NFW halo profiles including density, circular velocity, escape velocity,
and velocity dispersion profiles given a halo mass, NFW concentration parameter,
and redshift.  Uses the Hoeft, Mucket & Gottlober (2004) ApJ 602 velocity
dispersion prescription.

Example
-------
    from rbcodes.halo.rb_nfw import rb_nfw
    nfw = rb_nfw(m200=1e12, c=10, z=0.5)
    # nfw keys: m200, c, r200, v200, maxvcirc, maxvcircr, r, rho, vcirc, M_r, sig_v, vesc

### halo_profile

NFW halo escape velocity and virial-radius/mass conversion utilities.

Standalone module — not imported by other package modules; available for direct use.

Provides:
- ``NFW_escape_vel``         — escape velocity at radius r for a (possibly truncated) NFW profile
- ``Deltavir``               — Bryan & Norman (1998) virial overdensity
- ``rvirmvir``               — convert between Rvir and Mvir
- ``mvir_to_cvir``           — concentration–mass relation (Dutton & Maccio 2014)
- ``NFW_property_from_Mvir`` — convenience wrapper returning Rvir, cvir, and vesc

Example
-------
    import astropy.units as u
    from rbcodes.halo.halo_profile import NFW_property_from_Mvir
    result = NFW_property_from_Mvir(r=10.*u.kpc, Mvir=3e11*u.solMass, z=0.5)
    print(result['vesc'], result['rvir'])

## Functions

### mstar2mhalo() (`mstar2mhalo`)

Function to compute the halo mass and virial radius of a galaxy.
     We use the Moster 2010 relation to compute these values.
     Moster et al 2010. The Astrophysical Journal, 710:903?923, 2010 February 20
    
    
    Input:-     mstar:- stellar mass of the galaxy under consideration
                                                   (e.g.  10^(10.4) )
                   z :-  Redshift of the galaxy
    
    Returns:- 
        m200:-   halo mass according to the paper
       r200:-   halo virial radius according to NFW profile with c=.2
       Written :- R.B., 2018
    ------------------------------------------------------------------------

### rb_nfw() (`rb_nfw`)

Function to compute a NFW profile.
    Velocity Dispersion equation taken from Hoeft M.; Mucket J. P. & Gottlober, S, 2004, ApJ 602,1
    http://adsabs.harvard.edu/cgi-bin/bib_query?2004ApJ...602..162H


    Input :-  
       m200 :-   Halo mass
       c    :-   NFW concentration paramter
       z    :- redshift


    Returns :-

       A bunch of stuff

### NFW_escape_vel() (`halo_profile`)

NFW profile escape velocity

    Parameters
    ----------
    r : Quantity w/ length units
        Radial distance at which to compute the escape velocity
    Mvir : Quantity w/ mass units
        Virial Mass
    CvirorRs : Quantity w/ dimensionless or distance units
        (Virial) Concentration parameter (if dimensionless), 
        or halo scale radius (if length units)
    Rvir : Quantity w/ length units
        Virial radius
    truncated : bool or float
        False for infinite-size NFW or a number to cut off the 
        halo  at this many times Rvir

### Deltavir() (`halo_profile`)

Standard Delta-vir from Bryan&Norman 98 (*not* Delta-c)

### rvirmvir() (`halo_profile`)

Convert between Rvir and Mvir
    
    Parameters
    ----------
    rvirormvir : Quantity w/ mass or length units
        Either Rvir or Mvir, depending on the input units
    cosmo : astropy cosmology
        The cosmology to assume
    z : float
        The redshift to assume for the conversion
        
    Returns
    -------
    mvirorrvir : Quantity w/ mass or length units
         Whichever ``rvirormvir`` is *not*

### mvir_to_cvir() (`halo_profile`)

Power-law fit to the c_vir-M_vir relation from
    Equations 12 & 13 of Dutton & Maccio 2014, arXiv:1402.7073.

### NFW_property_from_Mvir() (`halo_profile`)

Compute NFW halo properties 
    Rvir from Standard Delta-vir from Bryan&Norman 98
    cvir from Dutton & Maccio 2014
    
    Parameters
    ----------
    r : Quantity w/ length units
        Radial distance at which to compute the escape velocity
    Mvir : Quantity w/ mass units
        Virial Mass
    cosmo: astropy cosmology input
    truncated : bool or float
        False for infinite-size NFW or a number to cut off the 
        halo  at this many times Rvir
    z : float
        The redshift to assume for the conversion


   Returns 
   -------
   Returns a dictionary with following attributes
   
   Rvir : Quantity w/ length units
        Virial radius

   vesc: Quantity w km/s units
         Escape velocity

    CvirorRs : Quantity w/ dimensionless or distance units
        (Virial) Concentration parameter (if dimensionless), 
        or halo scale radius (if length units)
    -----------------------------------------------------------

    Example:
    
    import astropy.units as u
    from astropy import cosmology
    from halo import halo_profile as h
    r=10.*u.kpc
    mvir=3e11*u.solMass
    z=6.3
    o=h.NFW_property_from_Mvir(r,mvir,z=z)

    -----------------------------------------------------------

