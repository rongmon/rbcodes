"""
Catalog cone-search utility.

Standalone module — not imported by other package modules; available for direct use.

Provides a simple angular cone search around a given (RA, Dec) pointing against
an input catalogue of coordinates, and an impact-parameter search that uses
catalogue redshifts to cut on physical separation, optionally keeping only
objects in front of / behind the QSO.

Example
-------
    from rbcodes.catalog.rb_search import cone_search, ip_search, ip_match
    result = cone_search(ra_center, dec_center, ra_list, dec_list, angular_scale_arcsec)
    out = ip_search(ra_qso, dec_qso, ra_gal, dec_gal, z_gal, max_ip=300.)
    pairs = ip_match(ra_qso, dec_qso, ra_gal, dec_gal, z_gal, max_ip=300.,
                     z_qso=z_qso, side='foreground', dv_buffer=3000.)
"""
import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.cosmology import Planck18
def cone_search(ra_center,dec_center,ra_list,dec_list,angular_scale):
    """
    Function to do a cone search around any (ra,dec) pointing with respect
    to a list of ra,dec entries.

    Parameters
    ----------
    
        ra_center     = RA of center (Degrees)
        dec_center    = DEC of center (Degrees) 
        ra_list       = RA list of input catalogue (Degrees)
        dec_list      = DEC list of input catalogue (Degrees) 
        angular_scale = angular search radius in arcsec
               
    
    Returns
    -------
        out:- struture containing RA,DEC and logical operator identifying the objects
    
    Written by R.B.  Mar 11 2013
    -------------------------------------------------------
    """
    deg2rad=np.pi/180
    rad2deg=180./np.pi

    ra_center=deg2rad*np.array(ra_center)
    dec_center=deg2rad*np.array(dec_center)
    ra_list=deg2rad*np.array(ra_list)
    dec_list=deg2rad*np.array(dec_list)

    theta =np.arccos(np.sin(dec_center)*np.sin(dec_list) + np.cos(dec_center)*np.cos(dec_list)*np.cos(ra_center-ra_list))
 
    sq= np.where(theta*rad2deg*3600<=angular_scale)
    
    if sq[0].size>0:
        outstr={'ra':ra_list[sq],'dec':dec_list[sq],'index':sq,'angle':theta[sq]*rad2deg*3600}
    else:
        outstr={'ra':-99,'dec':-99,'index':sq,'angle':-99}
    return outstr


def ip_search(ra_center, dec_center, ra_list, dec_list, zlist, max_ip,
              cosmo=None, comoving=False, z_qso=None, side='both', dv_buffer=0.):
    """
    Find catalogue objects within a physical impact parameter of a sightline.
    The impact parameter is always computed at each object's own redshift.

    Parameters
    ----------
        ra_center  = RA of center (Degrees)
        dec_center = DEC of center (Degrees)
        ra_list    = RA list of input catalogue (Degrees)
        dec_list   = DEC list of input catalogue (Degrees)
        zlist      = redshift list of input catalogue
        max_ip     = maximum impact parameter in kpc
        cosmo      = astropy cosmology [default: Planck18]
        comoving   = if True use comoving kpc, else proper kpc [default False]
        z_qso      = optional redshift of the sightline (QSO)
        side       = 'foreground', 'background' or 'both' w.r.t. z_qso
        dv_buffer  = exclude objects within this velocity (km/s) of z_qso
                     [default 0, i.e. strict z < z_qso / z > z_qso]

    Returns
    -------
        out:- dict with 'index' (into input lists), 'ra', 'dec', 'z',
              'angle' (arcsec), 'ip' (kpc) and 'vel' (km/s relative to
              z_qso, if z_qso given), sorted by impact parameter.
    """
    ra_list = np.asarray(ra_list, dtype=float)
    dec_list = np.asarray(dec_list, dtype=float)
    zlist = np.asarray(zlist, dtype=float)
    as_per_kpc = _arcsec_per_kpc(zlist, cosmo, comoving)
    return _ip_cut(ra_center, dec_center, ra_list, dec_list, zlist,
                   as_per_kpc, max_ip, z_qso, side, dv_buffer)


def ip_match(ra_qso, dec_qso, ra_gal, dec_gal, z_gal, max_ip,
             cosmo=None, comoving=False, z_qso=None, side='both', dv_buffer=0.,
             nbins=20, progress=True):
    """
    Catalogue-to-catalogue impact-parameter match: equivalent to running
    ip_search for every sightline in a QSO catalogue against a galaxy
    catalogue, but done with a KD-tree (astropy search_around_sky).

    Galaxies are split into bins of arcsec-per-kpc (i.e. redshift) so that
    each bin is searched with its own tight angular radius, instead of the
    huge radius set by the lowest-redshift galaxy.

    Parameters
    ----------
        ra_qso, dec_qso = RA, DEC lists of the QSO catalogue (Degrees)
        ra_gal, dec_gal = RA, DEC lists of the galaxy catalogue (Degrees)
        z_gal           = redshift list of the galaxy catalogue
        max_ip          = maximum impact parameter in kpc
        cosmo           = astropy cosmology [default: Planck18]
        comoving        = if True use comoving kpc, else proper kpc
        z_qso           = optional QSO redshift list (one per QSO)
        side            = 'foreground', 'background' or 'both' w.r.t. z_qso
        dv_buffer       = exclude galaxies within this velocity (km/s) of z_qso
        nbins           = number of redshift bins for the search [default 20]
        progress        = show a progress bar over the bins [default True]

    Returns
    -------
        out:- dict of flat arrays, one entry per QSO-galaxy pair:
              'qso_index', 'gal_index', 'ra', 'dec', 'z', 'angle' (arcsec),
              'ip' (kpc) and 'vel' (km/s relative to z_qso, if given),
              sorted by qso_index then impact parameter.
    """
    from astropy.coordinates import SkyCoord, search_around_sky

    if side not in ('foreground', 'background', 'both'):
        raise ValueError("side must be 'foreground', 'background' or 'both'")

    ra_qso = np.atleast_1d(np.asarray(ra_qso, dtype=float))
    dec_qso = np.atleast_1d(np.asarray(dec_qso, dtype=float))
    ra_gal = np.asarray(ra_gal, dtype=float)
    dec_gal = np.asarray(dec_gal, dtype=float)
    z_gal = np.asarray(z_gal, dtype=float)
    zq = None if z_qso is None else np.atleast_1d(np.asarray(z_qso, dtype=float))

    # cosmology evaluated once for the whole galaxy catalogue
    as_per_kpc = _arcsec_per_kpc(z_gal, cosmo, comoving)

    # Bin galaxies in arcsec/kpc (quantiles -> roughly equal numbers per bin)
    good = np.where(np.isfinite(as_per_kpc))[0]
    edges = np.unique(np.quantile(as_per_kpc[good], np.linspace(0, 1, nbins + 1)))
    gbin = np.clip(np.searchsorted(edges, as_per_kpc[good], side='right') - 1,
                   0, len(edges) - 2)

    qcoord = SkyCoord(ra_qso, dec_qso, unit='deg')
    bins = range(len(edges) - 1)
    if progress:
        try:
            from tqdm.auto import tqdm
            bins = tqdm(bins, desc='ip_match', unit='zbin')
        except ImportError:
            pass

    qi_all, gi_all, ang_all = [], [], []
    for b in bins:
        gsel = good[gbin == b]
        if gsel.size == 0:
            continue
        radius = max_ip * as_per_kpc[gsel].max()
        gcoord = SkyCoord(ra_gal[gsel], dec_gal[gsel], unit='deg')
        qi, gi, sep, _ = search_around_sky(qcoord, gcoord, radius * u.arcsec)
        qi_all.append(qi)
        gi_all.append(gsel[gi])
        ang_all.append(sep.to_value(u.arcsec))

    qi = np.concatenate(qi_all) if qi_all else np.array([], dtype=int)
    gi = np.concatenate(gi_all) if gi_all else np.array([], dtype=int)
    angle = np.concatenate(ang_all) if ang_all else np.array([])

    # Convert to impact parameter (at each galaxy's redshift) and cut
    ip = angle / as_per_kpc[gi]
    keep = ip <= max_ip

    # Foreground / background selection relative to each QSO
    vel = None
    if zq is not None:
        zqi = zq[qi]
        vel = c.to(u.km / u.s).value * (z_gal[gi] - zqi) / (1. + zqi)
        if side == 'foreground':
            keep &= vel < -dv_buffer
        elif side == 'background':
            keep &= vel > dv_buffer
        else:
            keep &= np.abs(vel) >= dv_buffer

    qi, gi, angle, ip = qi[keep], gi[keep], angle[keep], ip[keep]
    order = np.lexsort((ip, qi))
    qi, gi, angle, ip = qi[order], gi[order], angle[order], ip[order]

    out = {'gal_index': gi, 'ra': ra_gal[gi], 'dec': dec_gal[gi],
           'z': z_gal[gi], 'angle': angle, 'ip': ip}
    if vel is not None:
        out['vel'] = vel[keep][order]
    out['qso_index'] = qi
    return out


def _arcsec_per_kpc(zlist, cosmo, comoving):
    """arcsec per kpc at each redshift (NaN for z <= 0)."""
    if cosmo is None:
        cosmo = Planck18
    good = zlist > 0
    as_per_kpc = np.full(zlist.shape, np.nan)
    if comoving:
        as_per_kpc[good] = cosmo.arcsec_per_kpc_comoving(zlist[good]).value
    else:
        as_per_kpc[good] = cosmo.arcsec_per_kpc_proper(zlist[good]).value
    return as_per_kpc


def _ip_cut(ra_center, dec_center, ra_list, dec_list, zlist, as_per_kpc,
            max_ip, z_qso, side, dv_buffer):
    """Cone search one sightline, convert to impact parameter and cut."""
    if side not in ('foreground', 'background', 'both'):
        raise ValueError("side must be 'foreground', 'background' or 'both'")

    # Cone search with the largest angle any object could subtend at max_ip
    radius = max_ip * np.nanmax(as_per_kpc)
    cs = cone_search(ra_center, dec_center, ra_list, dec_list, radius)
    idx = cs['index'][0]
    angle = np.atleast_1d(cs['angle']) if idx.size > 0 else np.array([])

    # Convert to impact parameter (at each object's redshift) and cut
    ip = angle / as_per_kpc[idx]
    keep = ip <= max_ip

    # Foreground / background selection relative to the QSO
    vel = None
    if z_qso is not None:
        vel = c.to(u.km / u.s).value * (zlist[idx] - z_qso) / (1. + z_qso)
        if side == 'foreground':
            keep &= vel < -dv_buffer
        elif side == 'background':
            keep &= vel > dv_buffer
        else:
            keep &= np.abs(vel) >= dv_buffer

    idx, angle, ip = idx[keep], angle[keep], ip[keep]
    order = np.argsort(ip)
    idx, angle, ip = idx[order], angle[order], ip[order]

    out = {'index': idx, 'ra': ra_list[idx], 'dec': dec_list[idx],
           'z': zlist[idx], 'angle': angle, 'ip': ip}
    if vel is not None:
        out['vel'] = vel[keep][order]
    return out
