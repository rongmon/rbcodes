#Load modules

import numpy as np
from importlib import reload
import astropy.units as u
from astropy import cosmology
from astropy import constants as cnst

from astropy.table import Table



def NFW_escape_vel(r, Mvir, Rvir, CvirorRs, truncated=False):
    """
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
    """
    CvirorRs = u.Quantity(CvirorRs)
    if CvirorRs.unit.is_equivalent(u.m):
        Cvir = Rvir/CvirorRs
    elif CvirorRs.unit.is_equivalent(u.one):
        Cvir = CvirorRs
    else:
        raise TypeError('CvirorRs must be length or dimensionless')
        
    a = Rvir / Cvir
    
    #"f-function" from the NFW literature (e.g. klypin 02) evaluated at Cvir
    fofC = np.log(1 + Cvir) - Cvir / (1 + Cvir)
    
    # value of the NFW potential at that point
    potential = (-cnst.G * Mvir / fofC) * np.log(1 + (r / a)) / r

    if truncated:
        rtrunc = Rvir * float(truncated)
        Ctrunc = rtrunc / a

        mtrunc = Mvir * (np.log(1 + Ctrunc) - Ctrunc / (1 + Ctrunc)) / fofC

        outer = r >= rtrunc
        potential[outer] = - Gkpc * mtrunc / r[outer]
        potential[~outer] = potential[~outer] + (Gkpc * Mvir / a) / (Ctrunc + 1) / fofC

    vesc = (2 * np.abs(potential)) ** 0.5
    return vesc.to(u.km/u.s)
def Deltavir(cosmo, z=0):
    """
    Standard Delta-vir from Bryan&Norman 98 (*not* Delta-c)
    """
    x = cosmo.Om(z) - 1
    return (18*np.pi**2 + 82*x - 39*x**2)/(x+1)
def rvirmvir(rvirormvir, cosmo, z=0):
    """
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
    """
    rhs = Deltavir(cosmo=cosmo, z=z) * cosmo.Om(z)*cosmo.H(z)**2 / (2*cnst.G)
    
    if rvirormvir.unit.is_equivalent(u.solMass):
        mvir = rvirormvir
        return ((mvir / rhs)**(1/3)).to(u.kpc)
    elif rvirormvir.unit.is_equivalent(u.kpc):
        rvir = rvirormvir
        return (rhs * rvir**3).to(u.solMass)
    else:
        raise ValueError('invalid input unit {}'.format(rvirormvir))
        
def mvir_to_cvir(mvir, z=0):
    """ Power-law fit to the c_vir-M_vir relation from
    Equations 12 & 13 of Dutton & Maccio 2014, arXiv:1402.7073.
    """
    a = 0.537 + (1.025 - 0.537) * np.exp(-0.718 * z**1.08)
    b = -0.097 + 0.024 * z
    m0 = 1e12 * u.solMass

    logc = a + b * np.log10(mvir / m0)
    return 10**logc

def NFW_property_from_Mvir(r, Mvir, z=0,
                             cosmo=cosmology.Planck18,
                             truncated=False):
    """
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
    """
    cvir = mvir_to_cvir(Mvir, z)
    rvir = rvirmvir(Mvir, cosmo, z)
    vesc=NFW_escape_vel(r, Mvir=Mvir, 
                          CvirorRs=cvir, 
                          Rvir=rvir, 
                          truncated=truncated)

    out={}
    out['rvir']=rvir
    out['cvir']=cvir
    out['vesc']=vesc
    out['Mvir']=Mvir 
    out['r']=r
    out['z']=z
    
    return out