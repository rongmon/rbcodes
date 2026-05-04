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
- ``stellarToHaloMass(z, stellar_mass)`` — convert log stellar mass to log halo mass
  using the Behroozi et al. UniverseMachine SMHM relation.  Requires the parameter
  file ``smhm_med_params.txt`` to be present in the working directory.
- ``R_200(Mh, z)``                       — compute the virial radius R_200 from a
  log halo mass and redshift (uses Planck18 cosmology).

Note: scatter in the SMHM relation is ~0.1 dex (range 0.04–0.13 dex).

Example
-------
    from rbcodes.halo import halos as h
    import numpy as np
    sm = np.arange(8, 10, 0.2)          # log stellar mass
    mh = [h.stellarToHaloMass(0.5, s) for s in sm]
    rv = [h.R_200(m, 0.5) for m in mh]

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

