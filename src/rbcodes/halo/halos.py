#!/usr/bin/python
"""
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
  log halo mass and redshift (Planck18 by default; any astropy cosmology via
  ``cosmo=``).
- ``V_200(Mh, z)``, ``concentration_200c(Mh, z)``, ``nfw_vesc(Mh, z, R)`` — circular
  velocity, Dutton & Maccio (2014) concentration, and NFW escape velocity at R.
- ``halo_properties(z, stellar_mass, method='direct', population=None, centrals=True,
  R=None)`` — log_Mh, R200, V200, c200 and (given R) R/R200 and v_esc in one call.

Note: scatter in the SMHM relation is ~0.1 dex (range 0.04–0.13 dex).

Which halo mass to use
----------------------
'direct'   (haloMassAtStellarMass) : halo mass AT FIXED stellar mass, from the
             UniverseMachine mocks.  Use this to assign halo masses to observed
             galaxies.
'inverted' (stellarToHaloMass)     : inverts the median M*(Mh) fit.  Matches older
             work, but is biased high where the relation is steep: at log M* = 11,
             z = 0.43 it gives log Mh = 13.17 vs 12.66 (all centrals) or 12.49
             (star-forming centrals) from the direct method.
Halo masses are treated as M200c when computing R200, V200, c200 and v_esc.

Examples
--------
    from rbcodes.halo import halos as h

    # 1. Halo mass of one galaxy (log M* = 11, z = 0.43), direct method
    h.haloMassAtStellarMass(0.43, 11.0)                        # all centrals -> 12.66
    h.haloMassAtStellarMass(0.43, 11.0, population='sf')       # star-forming -> 12.49
    h.haloMassAtStellarMass(0.43, 11.0, population='passive')  # quiescent    -> 12.77
    h.haloMassAtStellarMass(0.43, 11.0, centrals=False)        # incl. sats   -> 12.63

    # 2. Several stellar masses at the same redshift
    h.haloMassAtStellarMass(0.43, [10.5, 11.0, 11.5], population='sf')
    # -> array([12.06, 12.49, 13.18])

    # 3. The inverted relation (original method in this module)
    h.stellarToHaloMass(0.43, 11.0)                    # -> 13.17
    h.stellarToHaloMass(0.43, [10.5, 11.0])            # arrays work too

    # 4. Halo properties from a log halo mass
    h.R_200(12.5, 0.43)                  # -> 264.9 kpc
    h.V_200(12.5, 0.43)                  # -> 226.6 km/s
    h.concentration_200c(12.5, 0.43)     # -> 6.29
    h.nfw_vesc(12.5, 0.43, 100.)         # -> 542.6 km/s at R = 100 kpc
    h.nfw_vesc(12.5, 0.43, [30, 100, 200])   # -> [658.8, 542.6, 460.1] km/s

    # 5. Everything for one galaxy + sightline in one call
    p = h.halo_properties(0.43, 11.0, population='sf', R=100.)
    p['log_Mh'], p['R200'], p['V200'], p['c200'], p['R_over_R200'], p['v_esc']
    # -> 12.49, 262.5 kpc, 224.6 km/s, 6.30, 0.38, 536.8 km/s

    # 6. Same galaxy with the inverted relation, for comparison with older work
    p_old = h.halo_properties(0.43, 11.0, method='inverted', centrals=False, R=100.)
    # -> log_Mh 13.17, R200 443.4 kpc, R/R200 0.23, v_esc 1001.8 km/s

    # 7. A galaxy table (one redshift per galaxy -> loop), classified by sSFR
    import numpy as np
    z    = [0.4580, 0.3772, 0.2895]
    logM = [11.0,   11.1,   10.8]
    ssfr = [-8.8,  -12.2,  -10.3]          # log sSFR [1/yr]
    R    = [97.3,   68.9,  132.8]          # impact parameter [kpc]
    v_c  = [[-575.8, -205.6, 119.9],       # component velocities [km/s]
            [-156.1, -46.3],
            [-109.1, 46.4, 193.1, 840.6]]
    for zi, mi, si, Ri, vi in zip(z, logM, ssfr, R, v_c):
        pop = 'q' if si < -11 else 'sf'    # UniverseMachine quenched cut
        p = h.halo_properties(zi, mi, population=pop, R=Ri)
        print(pop, round(p['log_Mh'], 2), round(p['R_over_R200'], 2),
              np.round(np.array(vi) / p['v_esc'].value, 2))   # v / v_esc
    # sf 12.48 0.38 [-1.07 -0.38  0.22]
    # q  12.9  0.19 [-0.19 -0.06]
    # sf 12.32 0.54 [-0.26  0.11  0.46  2.  ]
"""
import sys
import math
import re
import warnings
from importlib.resources import files
import numpy as np
from scipy.interpolate import interp1d
from astropy.cosmology import Planck18
from astropy.constants import G
import astropy.units as u
from numpy import log10
_POPULATION_ALIASES = {
    None: '', 'all': '',
    'sf': '_sf', 'starforming': '_sf', 'star-forming': '_sf', 'active': '_sf',
    'q': '_q', 'quenched': '_q', 'quiescent': '_q', 'passive': '_q',
}

def _smhm_param_file(population=None, centrals=False):
    """Return the UniverseMachine DR1 parameter file name for a galaxy selection."""
    key = population.lower() if isinstance(population, str) else population
    if key not in _POPULATION_ALIASES:
        raise ValueError("population must be None/'all', 'sf' or 'q' "
                         "(got %r)" % (population,))
    return 'smhm_med%s%s_params.txt' % ('_cen' if centrals else '',
                                         _POPULATION_ALIASES[key])

def stellarToHaloMass(z,stellar_mass,population=None,centrals=False):
    """
    Convert log stellar mass to log peak halo mass (Behroozi et al. 2019,
    Appendix J; UniverseMachine DR1 observed-median fits).

    population : None/'all' (default, all galaxies), 'sf' (star-forming;
                 aliases 'starforming', 'active') or 'q' (quenched; aliases
                 'quiescent', 'passive').
    centrals   : if True, use the fit to central galaxies only.

    Note: the fits give the median M* at fixed Mpeak; this function inverts
    that relation, which is not the same as the median Mpeak at fixed M*.
    Where the relation is steep (log M* >~ 10.8) the inversion biases halo
    masses high (~0.5 dex at log M* = 11, z ~ 0.4).  Use
    ``haloMassAtStellarMass`` for halo mass at fixed stellar mass.

    Examples
    --------
    >>> stellarToHaloMass(0.43, 11.0)                    # all galaxies
    13.17
    >>> stellarToHaloMass(0.43, 11.0, population='sf')   # SF fit (same here)
    13.17
    >>> stellarToHaloMass(0.43, [10.5, 11.0, 11.5])      # array input
    array([12.1 , 13.17, 14.11])
    """

    #Load params
    param_file = files('rbcodes.halo').joinpath(
        _smhm_param_file(population, centrals)).open('r')
    param_list = []
    allparams = []
    for line in param_file:
        param_list.append(float((line.split(" "))[1]))
        allparams.append(line.split(" "))
        
    if (len(param_list) != 20):
        print("Parameter file not correct length.  (Expected 20 lines, got %d)." % len(param_list))
        quit()

    stellar_mass = np.asarray(stellar_mass, dtype=float)
    if np.any(stellar_mass > 20):
        stellar_mass = log10(stellar_mass)
        
    names = "EFF_0 EFF_0_A EFF_0_A2 EFF_0_Z M_1 M_1_A M_1_A2 M_1_Z ALPHA ALPHA_A ALPHA_A2 ALPHA_Z BETA BETA_A BETA_Z DELTA GAMMA GAMMA_A GAMMA_Z CHI2".split(" ");
    params = dict(zip(names, param_list))

    #Print SMHM relation
    a = 1.0/(1.0+z)
    a1 = a - 1.0
    lna = math.log(a)
    zparams = {}
    zparams['m_1'] = params['M_1'] + a1*params['M_1_A'] - lna*params['M_1_A2'] \
        + z*params['M_1_Z']
    zparams['sm_0'] = zparams['m_1'] + params['EFF_0'] + a1*params['EFF_0_A'] \
        - lna*params['EFF_0_A2'] + z*params['EFF_0_Z']
    zparams['alpha'] = params['ALPHA'] + a1*params['ALPHA_A'] \
        - lna*params['ALPHA_A2'] + z*params['ALPHA_Z']
    zparams['beta'] = params['BETA'] + a1*params['BETA_A'] + z*params['BETA_Z']
    zparams['delta'] = params['DELTA']
    zparams['gamma'] = 10**(params['GAMMA'] + a1*params['GAMMA_A'] + z*params['GAMMA_Z'])
    
    smhm_max = 14.5-0.35*z
    if (params['CHI2']>200):
        print('#Warning: chi^2 > 200 implies that not all features are well fit.')
        print('Comparison with the raw data (in data/smhm/median_raw/) is crucial.')

    Mhalo = [x*0.05 for x in range(int(10.5*20),int(smhm_max*20+1),1)]
    Mstar = []
        
    for m in [x*0.05 for x in range(int(10.5*20),int(smhm_max*20+1),1)]:
        dm = m-zparams['m_1'];
        dm2 = dm/zparams['delta'];
        sm = zparams['sm_0'] - math.log10(10**(-zparams['alpha']*dm) + \
                                          10**(-zparams['beta']*dm)) + \
                                          zparams['gamma']*math.exp(-0.5*(dm2*dm2));
        Mstar.append(sm)
        # print("%.2f %.6f %.6f" % (m,sm,sm-m))

    f = interp1d(Mstar, Mhalo, fill_value='extrapolate')
    Mh = f(stellar_mass).flatten()

    if (len(Mh) == 1):
        return(Mh[0])
    else:
        return(Mh)

##########################################
#
#  Takes a linear (not log) halo mass and returns Virial radius

def R_200(Mh, z, cosmo=Planck18):
    """
    R_200c (kpc) for a log10 halo mass Mh (Msun) at redshift z.
    cosmo is any astropy cosmology (default Planck18).

    Examples
    --------
    >>> R_200(12.5, 0.43)
    <Quantity 264.87 kpc>
    >>> R_200([12.0, 12.5, 13.0], 0.43)
    <Quantity [180.5, 264.9, 388.8] kpc>
    """

    #if (Mh > 20):
    #    Mh = log10(Mh)
    
    Mh = np.asarray(Mh, dtype=float)
    H = cosmo.H(z)
    r200 =(((G * (10**Mh * u.M_sun).to(u.kg) / (100 * H.cgs**2))**(1/3)).to(u.kpc)).flatten()

    if (len(r200) == 1):
        return(r200[0])
    else:
        return(r200)

##########################################
#
#  Direct halo mass at fixed stellar mass, from the UniverseMachine mocks

_DIRECT_GROUPS = ('all', 'cen', 'cen_sf', 'cen_q', 'sat', 'sf', 'q')
_hm_tables = None

def _load_hm_tables():
    global _hm_tables
    if _hm_tables is None:
        with files('rbcodes.halo').joinpath('um_dr1_hm_averages.npz').open('rb') as fh:
            d = np.load(fh)
            _hm_tables = {k: d[k] for k in d.files}
    return _hm_tables

def _direct_group(population, centrals):
    key = population.lower() if isinstance(population, str) else population
    if key not in _POPULATION_ALIASES:
        raise ValueError("population must be None/'all', 'sf' or 'q' "
                         "(got %r)" % (population,))
    pop = {'': 'all', '_sf': 'sf', '_q': 'q'}[_POPULATION_ALIASES[key]]
    if centrals:
        pop = 'cen' if pop == 'all' else 'cen_' + pop
    return _DIRECT_GROUPS.index(pop)

def _hm_at_snapshot(t, isnap, g, logsm, statistic, min_count):
    y = t['log_mh_wl' if statistic == 'wl' else 'log_mh'][isnap, :, g]
    ok = (t['counts'][isnap, :, g] >= min_count) & (y > 0)
    x, y = t['log_sm'][ok], y[ok]
    out = np.interp(logsm, x, y)
    lo, hi = logsm < x[0], logsm > x[-1]
    if lo.any() or hi.any():
        warnings.warn("log M* outside the well-sampled UniverseMachine range "
                      "[%.1f, %.1f] at a=%.3f; extrapolating linearly."
                      % (x[0], x[-1], t['scale_factor'][isnap]))
        out[lo] = y[0] + (logsm[lo] - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
        out[hi] = y[-1] + (logsm[hi] - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
    return out

def haloMassAtStellarMass(z, stellar_mass, population=None, centrals=True,
                          statistic='log', min_count=10):
    """
    Log peak halo mass at fixed observed stellar mass, read directly from the
    UniverseMachine DR1 mock catalogues (Behroozi et al. 2019,
    data/smhm/averages/hm_averages_*.dat).  Unlike ``stellarToHaloMass``, this
    is not an inversion of the median M*(Mh) relation, so it is not biased high
    where that relation is steep (log M* >~ 10.8).

    z            : redshift (scalar).  Linearly interpolated in scale factor
                   between the bracketing snapshots.
    stellar_mass : log10(M*/Msun) (scalar or array).  Interpolated between the
                   0.2 dex bin centres.
    population   : None/'all' (default), 'sf' (star-forming; aliases
                   'starforming', 'active') or 'q' (quenched; aliases
                   'quiescent', 'passive').  UniverseMachine classifies
                   galaxies as quenched if sSFR < 1e-11 /yr.
    centrals     : if True (default), use central galaxies only.
    statistic    : 'log' (default) -> <log10 Mpeak>; 'wl' -> weak-lensing
                   average log10(<Mpeak^(2/3)>^(3/2)).
    min_count    : ignore bins with fewer mock galaxies than this; outside
                   the remaining range the relation is extrapolated (with a
                   warning).

    Halo masses are peak virial masses (Bryan & Norman 1998), in Msun.

    Examples
    --------
    >>> haloMassAtStellarMass(0.43, 11.0)                        # all centrals
    12.66
    >>> haloMassAtStellarMass(0.43, 11.0, population='sf')       # SF centrals
    12.49
    >>> haloMassAtStellarMass(0.43, 11.0, population='q')        # quenched centrals
    12.77
    >>> haloMassAtStellarMass(0.43, [10.5, 11.0, 11.5], population='sf')
    array([12.06, 12.49, 13.18])
    >>> haloMassAtStellarMass(0.43, 11.0, population='sf', statistic='wl')
    12.58                                   # lensing-weighted average
    >>> haloMassAtStellarMass(0.6503, 11.6, population='sf', min_count=20)
    12.99                                   # stricter bin cut -> extrapolates
    """
    t = _load_hm_tables()
    g = _direct_group(population, centrals)
    logsm = np.atleast_1d(np.asarray(stellar_mass, dtype=float))
    a = 1.0 / (1.0 + z)
    A = t['scale_factor']
    i1 = int(np.clip(np.searchsorted(A, a), 1, len(A) - 1))
    i0 = i1 - 1
    w = np.clip((a - A[i0]) / (A[i1] - A[i0]), 0, 1)
    Mh = (1 - w) * _hm_at_snapshot(t, i0, g, logsm, statistic, min_count) \
        + w * _hm_at_snapshot(t, i1, g, logsm, statistic, min_count)

    if (len(Mh) == 1):
        return(Mh[0])
    else:
        return(Mh)

##########################################
#
#  Halo properties for a log M200c halo mass

def V_200(Mh, z, cosmo=Planck18):
    """
    Circular velocity at R_200c, sqrt(G M200c / R200c), in km/s.

    Examples
    --------
    >>> V_200(12.5, 0.43)
    <Quantity 226.6 km / s>
    """
    return np.sqrt(G * (10**np.asarray(Mh) * u.M_sun) / R_200(Mh, z, cosmo)).to(u.km / u.s)

def concentration_200c(Mh, z, cosmo=Planck18):
    """
    NFW concentration c200c(M, z) of Dutton & Maccio (2014), calibrated on a
    Planck cosmology; cosmo only sets h for the M200c [Msun/h] conversion.

    Examples
    --------
    >>> concentration_200c(12.5, 0.43)
    6.29
    """
    a = 0.520 + (0.905 - 0.520) * np.exp(-0.617 * z**1.21)
    b = -0.101 + 0.026 * z
    return 10**(a + b * (np.asarray(Mh) + log10(cosmo.h) - 12.0))

def nfw_vesc(Mh, z, R, c=None, cosmo=Planck18):
    """
    Escape velocity (to infinity) at radius R from an NFW halo of log mass
    M200c, in km/s.  R in kpc (or an astropy length).  c defaults to the
    Dutton & Maccio (2014) concentration.  cosmo: astropy cosmology (default
    Planck18).

    Evaluated at the projected impact parameter this is an upper limit on the
    local escape velocity of the absorbing gas, which lies at r >= R.  Compare
    line-of-sight velocities with v_esc/sqrt(3) to allow for projection of
    isotropic motions.

    Examples
    --------
    >>> nfw_vesc(12.5, 0.43, 100.)                  # R = 100 kpc
    <Quantity 542.6 km / s>
    >>> nfw_vesc(12.5, 0.43, [30, 100, 200])        # escape-velocity profile
    <Quantity [658.8, 542.6, 460.1] km / s>
    >>> import astropy.units as u
    >>> nfw_vesc(12.5, 0.43, 0.1 * u.Mpc)           # any length unit
    <Quantity 542.6 km / s>
    >>> nfw_vesc(12.5, 0.43, 100., c=10)            # fixed concentration
    <Quantity 534.5 km / s>
    >>> nfw_vesc(haloMassAtStellarMass(0.43, 11.0, population='sf'), 0.43, 100.)
    <Quantity 536.8 km / s>                       # from a stellar mass
    """
    if c is None:
        c = concentration_200c(Mh, z, cosmo)
    x = (u.Quantity(R, u.kpc) / R_200(Mh, z, cosmo)).decompose().value
    mc = np.log(1 + c) - c / (1 + c)
    return V_200(Mh, z, cosmo) * np.sqrt(2 * np.log(1 + c * x) / (x * mc))

def halo_properties(z, stellar_mass, method='direct', population=None,
                    centrals=True, R=None, cosmo=Planck18):
    """
    Halo mass and derived properties for a galaxy of log stellar mass at z.

    method     : 'direct' (default) -> haloMassAtStellarMass (halo mass at
                 fixed M*); 'inverted' -> stellarToHaloMass (inverted median
                 M*(Mh) fit, biased high at log M* >~ 10.8).
    population, centrals : as in haloMassAtStellarMass / stellarToHaloMass.
    R          : optional impact parameter (kpc) for R/R200 and v_esc.
    cosmo      : astropy cosmology for R200, V200, c200, v_esc (default Planck18).

    The halo mass is treated as M200c.  Returns a dict with log_Mh, R200
    (kpc), V200 (km/s), c200 and, if R is given, R_over_R200 and v_esc (km/s).
    z must be a scalar; stellar_mass and R may be arrays.

    Examples
    --------
    >>> p = halo_properties(0.43, 11.0, population='sf', R=100.)
    >>> p['log_Mh'], p['R200'], p['V200'], p['c200'], p['R_over_R200'], p['v_esc']
    (12.49, <Quantity 262.5 kpc>, <Quantity 224.6 km / s>, 6.30, 0.38,
     <Quantity 536.8 km / s>)

    Inverted relation, all galaxies (as in older work):
    >>> p = halo_properties(0.43, 11.0, method='inverted', centrals=False, R=100.)
    >>> p['log_Mh'], p['R200'], p['R_over_R200'], p['v_esc']
    (13.17, <Quantity 443.4 kpc>, 0.23, <Quantity 1001.8 km / s>)

    Several galaxies at the same redshift:
    >>> p = halo_properties(0.43, [10.5, 11.0, 11.5], population='sf', R=100.)
    >>> p['log_Mh'], p['R_over_R200']
    (array([12.06, 12.49, 13.18]), array([0.53, 0.38, 0.22]))

    v / v_esc for absorber components (km/s):
    >>> v_c = np.array([-575.8, -205.6, 119.9])
    >>> v_c / halo_properties(0.458, 11.0, population='sf', R=97.3)['v_esc'].value
    array([-1.07, -0.38,  0.22])
    """
    if method == 'direct':
        Mh = haloMassAtStellarMass(z, stellar_mass, population=population,
                                   centrals=centrals)
    elif method == 'inverted':
        Mh = stellarToHaloMass(z, stellar_mass, population=population,
                               centrals=centrals)
    else:
        raise ValueError("method must be 'direct' or 'inverted' (got %r)" % (method,))

    out = {'log_Mh': Mh, 'R200': R_200(Mh, z, cosmo), 'V200': V_200(Mh, z, cosmo),
           'c200': concentration_200c(Mh, z, cosmo)}
    if R is not None:
        R = u.Quantity(R, u.kpc)
        out['R_over_R200'] = (R / out['R200']).decompose().value
        out['v_esc'] = nfw_vesc(Mh, z, R, cosmo=cosmo)
    return out




    
# Osterbrock 2006
# L_Ha = 1.37e-12 * Qion
# SFR_Ha = 5.5e-42 * L_Ha (Calzetti et al 2012)
# L(Ha) = 2.86 * L(Hb) (Case B)

# Qion = 7.3e11 * L(Ha) = 2.04e12 * L(Hb)
# SFR(Hb) = 1.57e-41 * L(Hb)

