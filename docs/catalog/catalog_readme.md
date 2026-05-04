# Project Documentation
[Back to Main Page](../main_readme.md)

*Auto-generated documentation from docstrings*

## Modules

### rb_search

Catalog cone-search utility.

Standalone module — not imported by other package modules; available for direct use.

Provides a simple angular cone search around a given (RA, Dec) pointing against
an input catalogue of coordinates.

Example
-------
    from rbcodes.catalog.rb_search import cone_search
    result = cone_search(ra_center, dec_center, ra_list, dec_list, angular_scale_arcsec)

### GalaxyGroupFinder_slow

Friend-of-Friends galaxy group finder — reference (slow) implementation.

Standalone module — not imported by other package modules; available for direct use.

This is the original, explicit reference implementation of the FoF group finder.
It computes pairwise separations in a straightforward loop and is easier to
read and verify than the optimised version in ``galaxy_group_finder.py``, but
does not scale to large catalogues.

Prefer ``galaxy_group_finder.GalaxyGroupFinder`` for production use.

Example
-------
    from rbcodes.catalog.GalaxyGroupFinder_slow import GalaxyGroupFinder
    gf = GalaxyGroupFinder()
    gf.load_catalog(ra, dec, redshift)
    gf.find_groups(linking_length_kpc=200, velocity_gap_kms=500)

### convert_FIRE_coordinates

Coordinate epoch conversion for Magellan/FIRE spectra.

Standalone module — not imported by other package modules; available for direct use.

Reads the EQUINOX keyword from a Magellan FIRE 1D spectrum FITS header and
transforms the stored (RA, Dec) to J2000 (FK5) coordinates using Astropy.

Note: Intended specifically for older Magellan FIRE QSO spectra that store
coordinates in a non-J2000 epoch.

Example
-------
    from rbcodes.catalog.convert_FIRE_coordinates import convert_epoch
    coord_j2000 = convert_epoch('J0100+28_F.fits')
    print(coord_j2000.to_string('hmsdms'))

### galaxy_group_finder

Friend-of-Friends galaxy group finder.

Standalone module — not imported by other package modules; available for direct use.

Provides two group-finder classes:

- ``AngularGroupFinder``   — angular-only FoF clustering (no redshift); suitable for
  photometric catalogues or quick angular associations.
- ``GalaxyGroupFinder``    — full FoF finder using pairwise angular, comoving, and
  physical separations plus velocity linking (requires redshifts).

Both classes follow the same interface: load a catalogue, run the group finder,
and retrieve a group catalogue and member catalogue.

See also: ``GalaxyGroupFinder_slow.py`` — a reference (slower) implementation.

Example
-------
    from rbcodes.catalog.galaxy_group_finder import GalaxyGroupFinder
    gf = GalaxyGroupFinder()
    gf.load_catalog(ra, dec, redshift)
    gf.find_groups(linking_length_kpc=200, velocity_gap_kms=500)
    groups = gf.group_catalog

## Classes

### GalaxyGroupFinder (`GalaxyGroupFinder_slow`)

Friend-of-Friends galaxy group finder with flexible linking criteria.
    
    Computes pairwise angular, comoving, and physical distances between galaxies,
    then applies spatial and velocity linking to identify groups.

### AngularGroupFinder (`galaxy_group_finder`)

Angular-only friend-of-friends group finder.
    
    Simple clustering based purely on angular separation without redshift information.
    Useful for photometric catalogs, star clustering, or quick angular associations.

### GalaxyGroupFinder (`galaxy_group_finder`)

Friend-of-Friends galaxy group finder with flexible linking criteria.
    
    Computes pairwise angular, comoving, and physical distances between galaxies,
    then applies spatial and velocity linking to identify groups.
    Uses AngularGroupFinder internally for angular calculations.

## Functions

### cone_search() (`rb_search`)

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

### __init__() (`GalaxyGroupFinder_slow`)

Initialize the group finder.
        
        Parameters:
        -----------
        cosmology : astropy.cosmology object, optional
            Cosmological model. Default is Planck18.

### load_catalog() (`GalaxyGroupFinder_slow`)

Load galaxy catalog.
        
        Parameters:
        -----------
        ra : array-like
            Right ascension in degrees
        dec : array-like  
            Declination in degrees
        redshift : array-like
            Redshift values
        galaxy_ids : array-like, optional
            Galaxy IDs. If None, uses sequential integers.

### compute_distance_matrix() (`GalaxyGroupFinder_slow`)

Compute the requested distance matrix efficiently.
        Always computes angular separation first, then converts if needed.
        
        Parameters:
        -----------
        distance_type : str
            'angular_arcsec', 'comoving_mpc', or 'physical_mpc'

### compute_velocity_matrix() (`GalaxyGroupFinder_slow`)

Compute velocity separation matrix using δv = c × Δz/(1+z₁).

### find_groups() (`GalaxyGroupFinder_slow`)

Find galaxy groups using friend-of-friends algorithm.
        
        Parameters:
        -----------
        spatial_linking_length : float
            Spatial linking length threshold
        velocity_linking_length : float  
            Velocity linking threshold in km/s
        distance_type : str
            Type of distance to use: 'angular_arcsec', 'comoving_mpc', 'physical_mpc'
        linking_method : str
            'chain' for chain linking or 'full' for fully connected groups
        full_linking_strategy : str
            For full linking only: 'merge', 'closest', or 'largest'
            - 'merge': merge groups if all members satisfy criteria (original behavior)
            - 'closest': assign galaxy to closest group (by mean distance)
            - 'largest': assign galaxy to largest group

### _chain_linking() (`GalaxyGroupFinder_slow`)

Chain linking: galaxies linked if connected through any path.

### _full_linking() (`GalaxyGroupFinder_slow`)

Full linking: all group members must be connected to all others.
        
        Parameters:
        -----------
        adjacency : ndarray
            Boolean adjacency matrix
        distance_matrix : ndarray
            Distance matrix for closest group calculations
        strategy : str
            'merge', 'closest', or 'largest'

### _full_linking_merge() (`GalaxyGroupFinder_slow`)

Original full linking with group merging.

### _full_linking_assign() (`GalaxyGroupFinder_slow`)

Full linking with assignment to closest or largest group.
        
        Strategy:
        1. Find all maximal cliques (fully connected components)
        2. For galaxies that could join multiple groups, assign based on strategy

### _find_maximal_clique() (`GalaxyGroupFinder_slow`)

Find maximal clique starting from start_node using greedy approach.

### _create_catalogs() (`GalaxyGroupFinder_slow`)

Create group and member catalogs.

### save_catalogs() (`GalaxyGroupFinder_slow`)

Save group and member catalogs.
        
        Parameters:
        -----------
        output_prefix : str
            Prefix for output files
        format : str
            'csv' or 'fits'

### print_summary() (`GalaxyGroupFinder_slow`)

Print summary of group finding results.

### convert_epoch() (`convert_FIRE_coordinates`)

Function to take Magellan FIRE spectrum and use the header information
	to transform the co-ordinates to J2000 epoch

	Parameters
	-----------

	 Magellan FIRE 1d spectrum filename in fits format

	Return
	------
	
	astropy SkyCoord object in J2000 epoch



	Example
	-------

	from catalog.convert_FIRE_coordinates import convert_epoch
	s=convert_epoch('J0100+28_F.fits')

	#show results
	s.to_string('hmsdms')


	Note: Only for specific use of older Magellan FIRE QSO spectra.


	Written By: Rongmon Bordoloi              May 2022
       ---------------------------------------------

### __init__() (`galaxy_group_finder`)

Initialize the group finder.
        
        Parameters:
        -----------
        cosmology : astropy.cosmology object, optional
            Cosmological model. Default is Planck18.

### load_catalog() (`galaxy_group_finder`)

Load galaxy catalog.
        
        Parameters:
        -----------
        ra : array-like
            Right ascension in degrees
        dec : array-like  
            Declination in degrees
        redshift : array-like
            Redshift values
        galaxy_ids : array-like, optional
            Galaxy IDs. If None, uses sequential integers.

### _get_angular_candidate_pairs() (`galaxy_group_finder`)

Pre-filter object pairs based on angular separation.

### compute_angular_matrix() (`galaxy_group_finder`)

Compute angular separation matrix with optional pre-filtering.
        
        Parameters:
        -----------
        max_angular_separation : float, optional
            Maximum angular separation in arcsec for pre-filtering.
            If provided, only computes separations up to this limit.

### find_groups() (`galaxy_group_finder`)

Find galaxy groups using friend-of-friends algorithm.
        
        Parameters:
        -----------
        spatial_linking_length : float
            Spatial linking length threshold
        velocity_linking_length : float  
            Velocity linking threshold in km/s
        distance_type : str
            Type of distance to use: 'angular_arcsec', 'comoving_mpc', 'physical_mpc'
        linking_method : str
            'chain' for chain linking or 'full' for fully connected groups
        full_linking_strategy : str
            For full linking only: 'merge', 'closest', or 'largest'
            - 'merge': merge groups if all members satisfy criteria (original behavior)
            - 'closest': assign galaxy to closest group (by mean distance)
            - 'largest': assign galaxy to largest group

### _chain_linking() (`galaxy_group_finder`)

Chain linking: galaxies linked if connected through any path.

### _full_linking() (`galaxy_group_finder`)

Full linking: all group members must be connected to all others.
        
        Parameters:
        -----------
        adjacency : ndarray
            Boolean adjacency matrix
        distance_matrix : ndarray
            Distance matrix for closest group calculations
        strategy : str
            'merge', 'closest', or 'largest'

### _full_linking_merge() (`galaxy_group_finder`)

Original full linking with group merging.

### _full_linking_assign() (`galaxy_group_finder`)

Full linking with assignment to closest or largest group.
        
        Strategy:
        1. Find all maximal cliques (fully connected components)
        2. For galaxies that could join multiple groups, assign based on strategy

### _find_maximal_clique() (`galaxy_group_finder`)

Find maximal clique starting from start_node using greedy approach.

### _create_catalogs() (`galaxy_group_finder`)

Create group and member catalogs.

### save_catalogs() (`galaxy_group_finder`)

Save group and member catalogs.
        
        Parameters:
        -----------
        output_prefix : str
            Prefix for output files
        format : str
            'csv' or 'fits'

### print_summary() (`galaxy_group_finder`)

Print summary of group finding results.

### _estimate_max_angular_separation() (`galaxy_group_finder`)

Estimate maximum angular separation for pre-filtering based on linking length.
        Uses the closest galaxies to get conservative upper bound.

### _get_candidate_pairs() (`galaxy_group_finder`)

Pre-filter galaxy pairs based on angular separation and velocity criteria.
        Returns indices of pairs that could potentially be linked.

### compute_distance_matrix() (`galaxy_group_finder`)

Compute the requested distance matrix efficiently with pre-filtering.
        Always computes angular separation first, then converts if needed.
        
        Parameters:
        -----------
        distance_type : str
            'angular_arcsec', 'comoving_mpc', or 'physical_mpc'
        max_velocity_linking : float, optional
            Maximum velocity separation for pre-filtering (km/s).
            If provided, only computes distances for galaxy pairs that could be linked.

### compute_velocity_matrix() (`galaxy_group_finder`)

Compute velocity separation matrix using δv = c × Δz/(1+z₁).

