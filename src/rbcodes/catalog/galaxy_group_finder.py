import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import warnings
from typing import Optional, Tuple, Union
from scipy.spatial.distance import pdist, squareform


class AngularGroupFinder:
    """
    Angular-only friend-of-friends group finder.
    
    Simple clustering based purely on angular separation without redshift information.
    Useful for photometric catalogs, star clustering, or quick angular associations.
    """
    
    def __init__(self):
        """Initialize the angular group finder."""
        self.galaxies = None
        self.angular_matrix = None
        self.groups = None
        self.group_catalog = None
        self.member_catalog = None
    
    def load_catalog(self, ra, dec, object_ids=None):
        """
        Load object catalog (galaxies, stars, etc.).
        
        Parameters:
        -----------
        ra : array-like
            Right ascension in degrees
        dec : array-like  
            Declination in degrees
        object_ids : array-like, optional
            Object IDs. If None, uses sequential integers.
        """
        ra = np.asarray(ra)
        dec = np.asarray(dec)
        
        if object_ids is None:
            object_ids = np.arange(len(ra))
        else:
            object_ids = np.asarray(object_ids)
            
        # Validate inputs
        if not (len(ra) == len(dec) == len(object_ids)):
            raise ValueError("All input arrays must have the same length")
            
        self.galaxies = pd.DataFrame({
            'object_id': object_ids,
            'ra': ra,
            'dec': dec
        })
        
        print(f"Loaded {len(self.galaxies)} objects")
    
    def _get_angular_candidate_pairs(self, angular_matrix, max_angular_sep):
        """
        Pre-filter object pairs based on angular separation.
        """
        n_obj = len(self.galaxies)
        
        # Pre-filter by angular separation
        candidate_mask = angular_matrix < max_angular_sep
        np.fill_diagonal(candidate_mask, False)  # Remove self-pairs
        
        # Get indices of candidate pairs
        candidate_pairs = np.where(candidate_mask)
        n_candidates = len(candidate_pairs[0])
        n_total = n_obj * (n_obj - 1)
        
        print(f"Angular pre-filtering: {n_candidates:,} candidate pairs out of {n_total:,} total ({100*n_candidates/n_total:.1f}%)")
        
        return candidate_pairs, candidate_mask
    
    def compute_angular_matrix(self, max_angular_separation=None):
        """
        Compute angular separation matrix with optional pre-filtering.
        
        Parameters:
        -----------
        max_angular_separation : float, optional
            Maximum angular separation in arcsec for pre-filtering.
            If provided, only computes separations up to this limit.
        """
        if self.galaxies is None:
            raise ValueError("No catalog loaded. Call load_catalog() first.")
            
        # Return cached matrix if already computed
        if self.angular_matrix is not None:
            return self.angular_matrix
            
        n_obj = len(self.galaxies)
        print(f"Computing angular separation matrix for {n_obj} objects...")
        
        # Vectorized computation of all pairwise angular separations
        print("Computing angular separations (vectorized)...")
        coords = SkyCoord(
            ra=self.galaxies['ra'].values * u.deg,
            dec=self.galaxies['dec'].values * u.deg
        )
        
        # Vectorized computation
        separations = coords[:, None].separation(coords[None, :])
        angular_matrix = separations.arcsec
        
        # Optional: Zero out distances beyond threshold for memory/speed
        if max_angular_separation is not None:
            angular_matrix[angular_matrix > max_angular_separation] = np.inf
        
        self.angular_matrix = angular_matrix
        print("Angular separation matrix computed successfully!")
        return self.angular_matrix
    
    def find_groups(self, angular_linking_length, linking_method='chain', 
                   full_linking_strategy='merge'):
        """
        Find angular groups using friend-of-friends algorithm.
        
        Parameters:
        -----------
        angular_linking_length : float
            Angular linking length threshold in arcsec
        linking_method : str
            'chain' for chain linking or 'full' for fully connected groups
        full_linking_strategy : str
            For full linking only: 'merge', 'closest', or 'largest'
        """
        print(f"Finding groups with angular separation < {angular_linking_length} arcsec")
        print(f"Using {linking_method} linking method")
        
        n_obj = len(self.galaxies)
        
        # Compute angular separation matrix (with pre-filtering optimization)
        angular_matrix = self.compute_angular_matrix(angular_linking_length * 1.1)  # Small buffer
        
        # Create adjacency matrix based on linking criteria
        adjacency = angular_matrix < angular_linking_length
        np.fill_diagonal(adjacency, False)  # Remove self-links
        
        if linking_method == 'chain':
            groups = self._chain_linking(adjacency)
        elif linking_method == 'full':
            groups = self._full_linking(adjacency, angular_matrix, full_linking_strategy)
        else:
            raise ValueError("linking_method must be 'chain' or 'full'")
        
        self.groups = groups
        self._create_catalogs()
        
        n_groups = len(set(groups))
        n_singles = sum(1 for g in set(groups) if list(groups).count(g) == 1)
        print(f"Found {n_groups} groups ({n_groups - n_singles} multi-member, {n_singles} singles)")
    
    def _chain_linking(self, adjacency):
        """Chain linking: objects linked if connected through any path."""
        n_obj = len(adjacency)
        groups = [-1] * n_obj  # -1 means ungrouped
        current_group_id = 0
        
        for i in range(n_obj):
            if groups[i] == -1:  # Not yet assigned
                # Start new group with BFS
                group_members = []
                queue = [i]
                visited = set([i])
                
                while queue:
                    current = queue.pop(0)
                    group_members.append(current)
                    
                    # Find all neighbors
                    neighbors = np.where(adjacency[current])[0]
                    for neighbor in neighbors:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
                
                # Assign group ID to all members
                for member in group_members:
                    groups[member] = current_group_id
                
                current_group_id += 1
        
        return groups
    
    def _full_linking(self, adjacency, distance_matrix, strategy='merge'):
        """Full linking with different strategies."""
        n_obj = len(adjacency)
        
        if strategy == 'merge':
            return self._full_linking_merge(adjacency)
        elif strategy in ['closest', 'largest']:
            return self._full_linking_assign(adjacency, distance_matrix, strategy)
        else:
            raise ValueError("strategy must be 'merge', 'closest', or 'largest'")
    
    def _full_linking_merge(self, adjacency):
        """Original full linking with group merging."""
        n_obj = len(adjacency)
        groups = list(range(n_obj))  # Start with each object as its own group
        
        merged = True
        while merged:
            merged = False
            for i in range(n_obj):
                for j in range(i+1, n_obj):
                    if groups[i] != groups[j] and adjacency[i, j]:
                        # Check if merging these groups maintains full connectivity
                        group_i_members = [k for k in range(n_obj) if groups[k] == groups[i]]
                        group_j_members = [k for k in range(n_obj) if groups[k] == groups[j]]
                        
                        # Check if all members of both groups are connected to each other
                        can_merge = True
                        for gi in group_i_members:
                            for gj in group_j_members:
                                if not adjacency[gi, gj]:
                                    can_merge = False
                                    break
                            if not can_merge:
                                break
                        
                        if can_merge:
                            # Merge groups
                            old_group_j = groups[j]
                            for k in range(n_obj):
                                if groups[k] == old_group_j:
                                    groups[k] = groups[i]
                            merged = True
        
        # Renumber groups sequentially
        unique_groups = list(set(groups))
        group_mapping = {old: new for new, old in enumerate(unique_groups)}
        groups = [group_mapping[g] for g in groups]
        
        return groups
    
    def _full_linking_assign(self, adjacency, distance_matrix, strategy):
        """Full linking with assignment to closest or largest group."""
        n_obj = len(adjacency)
        groups = [-1] * n_obj  # -1 means unassigned
        current_group_id = 0
        
        # First pass: find initial fully connected groups (cliques)
        processed = [False] * n_obj
        
        for i in range(n_obj):
            if processed[i]:
                continue
                
            # Find maximal clique starting from object i
            clique = self._find_maximal_clique(adjacency, i, processed)
            
            if len(clique) > 1:
                # Assign group ID to clique members
                for member in clique:
                    groups[member] = current_group_id
                    processed[member] = True
                current_group_id += 1
        
        # Second pass: assign remaining objects to groups or as singles
        for i in range(n_obj):
            if groups[i] == -1:  # Unassigned
                # Find which existing groups this object can join
                candidate_groups = []
                
                for group_id in range(current_group_id):
                    group_members = [j for j in range(n_obj) if groups[j] == group_id]
                    
                    # Check if object i can join this group (connected to all members)
                    can_join = True
                    for member in group_members:
                        if not adjacency[i, member]:
                            can_join = False
                            break
                    
                    if can_join:
                        if strategy == 'closest':
                            # Calculate mean distance to group
                            mean_dist = np.mean([distance_matrix[i, member] for member in group_members])
                            candidate_groups.append((group_id, mean_dist, len(group_members)))
                        elif strategy == 'largest':
                            candidate_groups.append((group_id, 0, len(group_members)))
                
                if candidate_groups:
                    if strategy == 'closest':
                        # Sort by distance (ascending)
                        candidate_groups.sort(key=lambda x: x[1])
                        chosen_group = candidate_groups[0][0]
                    elif strategy == 'largest':
                        # Sort by size (descending)
                        candidate_groups.sort(key=lambda x: x[2], reverse=True)
                        chosen_group = candidate_groups[0][0]
                    
                    groups[i] = chosen_group
                else:
                    # No group to join, create new single-member group
                    groups[i] = current_group_id
                    current_group_id += 1
        
        return groups
    
    def _find_maximal_clique(self, adjacency, start_node, processed):
        """Find maximal clique starting from start_node using greedy approach."""
        n_obj = len(adjacency)
        clique = [start_node]
        
        # Find all neighbors of start_node
        candidates = [j for j in range(n_obj) 
                     if j != start_node and adjacency[start_node, j] and not processed[j]]
        
        # Greedily add nodes that are connected to all current clique members
        for candidate in candidates:
            # Check if candidate is connected to all clique members
            connected_to_all = True
            for clique_member in clique:
                if not adjacency[candidate, clique_member]:
                    connected_to_all = False
                    break
            
            if connected_to_all:
                clique.append(candidate)
        
        return clique
    
    def _create_catalogs(self):
        """Create group and member catalogs."""
        if self.groups is None:
            raise ValueError("No groups found. Run find_groups() first.")
        
        # Member catalog
        self.member_catalog = self.galaxies.copy()
        self.member_catalog['group_id'] = self.groups
        
        # Group catalog
        group_stats = []
        for group_id in set(self.groups):
            group_members = self.member_catalog[self.member_catalog['group_id'] == group_id]
            
            n_members = len(group_members)
            mean_ra = group_members['ra'].mean()
            mean_dec = group_members['dec'].mean()
            
            # Angular extent (maximum separation within group)
            if n_members > 1:
                group_coords = SkyCoord(
                    ra=group_members['ra'].values * u.deg,
                    dec=group_members['dec'].values * u.deg
                )
                # Calculate all pairwise separations within the group
                separations = []
                for i in range(len(group_coords)):
                    for j in range(i+1, len(group_coords)):
                        sep = group_coords[i].separation(group_coords[j]).arcsec
                        separations.append(sep)
                max_angular_extent = max(separations) if separations else 0.0
            else:
                max_angular_extent = 0.0
            
            group_stats.append({
                'group_id': group_id,
                'member_id': ','.join(map(str, group_members['object_id'].values)),
                'n_members': n_members,
                'mean_ra': mean_ra,
                'mean_dec': mean_dec,
                'max_angular_extent': max_angular_extent
            })
        
        self.group_catalog = pd.DataFrame(group_stats)
        self.group_catalog = self.group_catalog.sort_values('group_id').reset_index(drop=True)
    
    def save_catalogs(self, output_prefix, format='csv'):
        """Save group and member catalogs."""
        if self.group_catalog is None or self.member_catalog is None:
            raise ValueError("No catalogs to save. Run find_groups() first.")
        
        if format == 'csv':
            group_file = f"{output_prefix}_angular_groups.csv"
            member_file = f"{output_prefix}_angular_members.csv"
            
            self.group_catalog.to_csv(group_file, index=False)
            self.member_catalog.to_csv(member_file, index=False)
            
            print(f"Saved angular group catalog to {group_file}")
            print(f"Saved angular member catalog to {member_file}")
            
        elif format == 'fits':
            fits_file = f"{output_prefix}_angular_catalogs.fits"
            
            # Convert to astropy tables
            group_table = Table.from_pandas(self.group_catalog)
            member_table = Table.from_pandas(self.member_catalog)
            
            # Create HDU list
            primary_hdu = fits.PrimaryHDU()
            group_hdu = fits.BinTableHDU(group_table, name='ANGULAR_GROUPS')
            member_hdu = fits.BinTableHDU(member_table, name='ANGULAR_MEMBERS')
            
            hdul = fits.HDUList([primary_hdu, group_hdu, member_hdu])
            hdul.writeto(fits_file, overwrite=True)
            
            print(f"Saved angular catalogs to {fits_file}")
        
        else:
            raise ValueError("format must be 'csv' or 'fits'")
    
    def print_summary(self):
        """Print summary of group finding results."""
        if self.group_catalog is None:
            print("No groups found yet.")
            return
        
        total_groups = len(self.group_catalog)
        multi_member = len(self.group_catalog[self.group_catalog['n_members'] > 1])
        singles = total_groups - multi_member
        
        print(f"\n=== ANGULAR GROUP FINDING SUMMARY ===")
        print(f"Total objects: {len(self.galaxies)}")
        print(f"Total groups: {total_groups}")
        print(f"Multi-member groups: {multi_member}")
        print(f"Single objects: {singles}")
        
        if multi_member > 0:
            print(f"\nMulti-member group statistics:")
            multi_groups = self.group_catalog[self.group_catalog['n_members'] > 1]
            print(f"Largest group size: {multi_groups['n_members'].max()}")
            print(f"Mean group size: {multi_groups['n_members'].mean():.1f}")
            print(f"Mean angular extent: {multi_groups['max_angular_extent'].mean():.1f} arcsec")


class GalaxyGroupFinder:
    """
    Friend-of-Friends galaxy group finder with flexible linking criteria.
    
    Computes pairwise angular, comoving, and physical distances between galaxies,
    then applies spatial and velocity linking to identify groups.
    Uses AngularGroupFinder internally for angular calculations.
    """
    
    def __init__(self, cosmology=None):
        """
        Initialize the group finder.
        
        Parameters:
        -----------
        cosmology : astropy.cosmology object, optional
            Cosmological model. Default is Planck18.
        """
        self.cosmology = cosmology if cosmology is not None else Planck18
        self.galaxies = None
        self.distance_matrices = {}
        self.velocity_matrix = None
        self.groups = None
        self.group_catalog = None
        self.member_catalog = None
        
        # Internal angular group finder for angular calculations
        self._angular_finder = AngularGroupFinder()
    
    def load_catalog(self, ra, dec, redshift, galaxy_ids=None):
        """
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
        """
        ra = np.asarray(ra)
        dec = np.asarray(dec)
        redshift = np.asarray(redshift)
        
        if galaxy_ids is None:
            galaxy_ids = np.arange(len(ra))
        else:
            galaxy_ids = np.asarray(galaxy_ids)
            
        # Validate inputs
        if not (len(ra) == len(dec) == len(redshift) == len(galaxy_ids)):
            raise ValueError("All input arrays must have the same length")
            
        self.galaxies = pd.DataFrame({
            'galaxy_id': galaxy_ids,
            'ra': ra,
            'dec': dec,
            'redshift': redshift
        })
        
        # Load RA/Dec into internal angular finder
        self._angular_finder.load_catalog(ra, dec, galaxy_ids)
        
        print(f"Loaded {len(self.galaxies)} galaxies")
    
    def _estimate_max_angular_separation(self, spatial_linking_length, distance_type):
        """
        Estimate maximum angular separation for pre-filtering based on linking length.
        Uses the closest galaxies to get conservative upper bound.
        """
        if distance_type == 'angular_arcsec':
            return spatial_linking_length
        
        # For distance-based linking, convert to angular separation using closest galaxies
        min_redshift = self.galaxies['redshift'].min()
        
        if distance_type == 'comoving_mpc':
            # Convert comoving distance to angular separation at minimum redshift
            comoving_angular_diam_dist = self.cosmology.angular_diameter_distance(min_redshift).value * (1 + min_redshift)
            max_angular_rad = spatial_linking_length / comoving_angular_diam_dist
        elif distance_type == 'physical_mpc':
            # Convert physical distance to angular separation at minimum redshift  
            angular_diam_dist = self.cosmology.angular_diameter_distance(min_redshift).value
            max_angular_rad = spatial_linking_length / angular_diam_dist
        
        # Convert radians to arcsec and add safety factor
        max_angular_arcsec = max_angular_rad * (180/np.pi) * 3600 * 1.5  # 50% safety margin
        
        return max_angular_arcsec
    
    def _get_candidate_pairs(self, angular_matrix, max_angular_sep, max_velocity_sep):
        """
        Pre-filter galaxy pairs based on angular separation and velocity criteria.
        Returns indices of pairs that could potentially be linked.
        """
        n_gal = len(self.galaxies)
        redshifts = self.galaxies['redshift'].values
        c_kms = 299792.458
        
        # Pre-filter by angular separation
        angular_mask = angular_matrix < max_angular_sep
        
        # Pre-filter by velocity (symmetric check using minimum redshift)
        velocity_mask = np.zeros((n_gal, n_gal), dtype=bool)
        for i in range(n_gal):
            for j in range(n_gal):
                if angular_mask[i, j] and i != j:
                    # Use minimum redshift for conservative velocity estimate
                    z_min = min(redshifts[i], redshifts[j])
                    delta_z = abs(redshifts[j] - redshifts[i])
                    delta_v = c_kms * delta_z / (1 + z_min)
                    velocity_mask[i, j] = delta_v < max_velocity_sep
        
        # Combined mask
        candidate_mask = angular_mask & velocity_mask
        np.fill_diagonal(candidate_mask, False)  # Remove self-pairs
        
        # Get indices of candidate pairs
        candidate_pairs = np.where(candidate_mask)
        n_candidates = len(candidate_pairs[0])
        n_total = n_gal * (n_gal - 1)
        
        print(f"Pre-filtering: {n_candidates:,} candidate pairs out of {n_total:,} total ({100*n_candidates/n_total:.1f}%)")
        
        return candidate_pairs, candidate_mask
    
    def compute_distance_matrix(self, distance_type, max_velocity_linking=None):
        """
        Compute the requested distance matrix efficiently with pre-filtering.
        Always computes angular separation first, then converts if needed.
        
        Parameters:
        -----------
        distance_type : str
            'angular_arcsec', 'comoving_mpc', or 'physical_mpc'
        max_velocity_linking : float, optional
            Maximum velocity separation for pre-filtering (km/s).
            If provided, only computes distances for galaxy pairs that could be linked.
        """
        if self.galaxies is None:
            raise ValueError("No catalog loaded. Call load_catalog() first.")
            
        # Return cached matrix if already computed
        if distance_type in self.distance_matrices:
            return self.distance_matrices[distance_type]
            
        n_gal = len(self.galaxies)
        print(f"Computing {distance_type} distance matrix for {n_gal} galaxies...")
        
        # Always compute angular separation first (if not cached) - USE INTERNAL ANGULAR FINDER
        if 'angular_arcsec' not in self.distance_matrices:
            print("Computing angular separations (using AngularGroupFinder)...")
            # Estimate max angular separation for pre-filtering if velocity linking is provided
            max_angular_for_prefilter = None
            if max_velocity_linking is not None:
                max_angular_for_prefilter = self._estimate_max_angular_separation(1e6, distance_type)
            
            angular_matrix = self._angular_finder.compute_angular_matrix(max_angular_for_prefilter)
            self.distance_matrices['angular_arcsec'] = angular_matrix
        
        # If angular separation is what's needed, we're done
        if distance_type == 'angular_arcsec':
            print("Angular separation matrix computed successfully!")
            return self.distance_matrices['angular_arcsec']
        
        # For comoving/physical distances, convert from angular separation
        angular_matrix = self.distance_matrices['angular_arcsec']
        redshifts = self.galaxies['redshift'].values
        
        # Pre-filtering optimization
        if max_velocity_linking is not None:
            print("Applying pre-filtering optimization...")
            max_angular_sep = self._estimate_max_angular_separation(
                spatial_linking_length=1e6,  # Large value for pre-filtering
                distance_type=distance_type
            )
            candidate_pairs, candidate_mask = self._get_candidate_pairs(
                angular_matrix, max_angular_sep, max_velocity_linking
            )
        else:
            # No pre-filtering, compute all pairs
            candidate_pairs = None
            candidate_mask = None
        
        if distance_type == 'comoving_mpc':
            print("Converting to comoving distances (observer-centric)...")
            distance_matrix = np.zeros((n_gal, n_gal))
            
            if candidate_pairs is not None:
                # Only compute distances for candidate pairs
                for idx in range(len(candidate_pairs[0])):
                    i, j = candidate_pairs[0][idx], candidate_pairs[1][idx]
                    
                    # Observer-centric: use redshift of primary galaxy (i)
                    z_primary = redshifts[i]
                    
                    # Convert arcsec to comoving distance at primary galaxy's redshift
                    comoving_angular_diam_dist = self.cosmology.angular_diameter_distance(z_primary).value * (1 + z_primary)
                    theta_rad = angular_matrix[i, j] * (1/3600) * (np.pi/180)
                    comoving_sep = theta_rad * comoving_angular_diam_dist
                    
                    distance_matrix[i, j] = comoving_sep
            else:
                # Compute all pairs (original approach)
                for i in range(n_gal):
                    for j in range(n_gal):
                        if i != j:
                            z_primary = redshifts[i]
                            comoving_angular_diam_dist = self.cosmology.angular_diameter_distance(z_primary).value * (1 + z_primary)
                            theta_rad = angular_matrix[i, j] * (1/3600) * (np.pi/180)
                            comoving_sep = theta_rad * comoving_angular_diam_dist
                            distance_matrix[i, j] = comoving_sep
            
            self.distance_matrices['comoving_mpc'] = distance_matrix
            
        elif distance_type == 'physical_mpc':
            print("Converting to physical distances (observer-centric)...")
            distance_matrix = np.zeros((n_gal, n_gal))
            
            if candidate_pairs is not None:
                # Only compute distances for candidate pairs
                for idx in range(len(candidate_pairs[0])):
                    i, j = candidate_pairs[0][idx], candidate_pairs[1][idx]
                    
                    # Observer-centric: use redshift of primary galaxy (i)
                    z_primary = redshifts[i]
                    
                    # Convert arcsec to physical distance at primary galaxy's redshift
                    angular_diameter_dist = self.cosmology.angular_diameter_distance(z_primary).value
                    theta_rad = angular_matrix[i, j] * (1/3600) * (np.pi/180)
                    physical_sep = theta_rad * angular_diameter_dist
                    
                    distance_matrix[i, j] = physical_sep
            else:
                # Compute all pairs (original approach)
                for i in range(n_gal):
                    for j in range(n_gal):
                        if i != j:
                            z_primary = redshifts[i]
                            angular_diameter_dist = self.cosmology.angular_diameter_distance(z_primary).value
                            theta_rad = angular_matrix[i, j] * (1/3600) * (np.pi/180)
                            physical_sep = theta_rad * angular_diameter_dist
                            distance_matrix[i, j] = physical_sep
            
            self.distance_matrices['physical_mpc'] = distance_matrix
        
        else:
            raise ValueError("distance_type must be 'angular_arcsec', 'comoving_mpc', or 'physical_mpc'")
        
        print(f"{distance_type} matrix computed successfully!")
        return self.distance_matrices[distance_type]
    
    def compute_velocity_matrix(self):
        """
        Compute velocity separation matrix using δv = c × Δz/(1+z₁).
        """
        if self.galaxies is None:
            raise ValueError("No catalog loaded.")
            
        print("Computing velocity separations...")
        n_gal = len(self.galaxies)
        redshifts = self.galaxies['redshift'].values
        c_kms = 299792.458  # km/s
        
        velocity_matrix = np.zeros((n_gal, n_gal))
        for i in range(n_gal):
            for j in range(n_gal):
                if i != j:
                    delta_z = abs(redshifts[j] - redshifts[i])
                    z1 = redshifts[i]
                    delta_v = c_kms * delta_z / (1 + z1)
                    velocity_matrix[i, j] = delta_v
        
        self.velocity_matrix = velocity_matrix
        print("Velocity matrix computed successfully!")
    
    def find_groups(self, spatial_linking_length, velocity_linking_length, 
                   distance_type='comoving_mpc', linking_method='chain', 
                   full_linking_strategy='merge'):
        """
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
        """
        print(f"Finding groups with {distance_type} < {spatial_linking_length} and velocity < {velocity_linking_length} km/s")
        print(f"Using {linking_method} linking method")
        
        n_gal = len(self.galaxies)
        
        # Compute only the requested distance matrix (with optional pre-filtering)
        distance_matrix = self.compute_distance_matrix(distance_type, velocity_linking_length)
        
        # Compute velocity matrix if needed
        if self.velocity_matrix is None:
            self.compute_velocity_matrix()
        
        # Create adjacency matrix based on linking criteria
        spatial_links = distance_matrix < spatial_linking_length
        velocity_links = self.velocity_matrix < velocity_linking_length
        
        # Combined adjacency matrix (both spatial and velocity criteria must be satisfied)
        adjacency = spatial_links & velocity_links
        
        # Remove self-links
        np.fill_diagonal(adjacency, False)
        
        if linking_method == 'chain':
            groups = self._chain_linking(adjacency)
        elif linking_method == 'full':
            groups = self._full_linking(adjacency, distance_matrix, full_linking_strategy)
        else:
            raise ValueError("linking_method must be 'chain' or 'full'")
        
        self.groups = groups
        self._create_catalogs()
        
        n_groups = len(set(groups))
        n_singles = sum(1 for g in set(groups) if list(groups).count(g) == 1)
        print(f"Found {n_groups} groups ({n_groups - n_singles} multi-member, {n_singles} singles)")
    
    def _chain_linking(self, adjacency):
        """Chain linking: galaxies linked if connected through any path."""
        n_gal = len(adjacency)
        groups = [-1] * n_gal  # -1 means ungrouped
        current_group_id = 0
        
        for i in range(n_gal):
            if groups[i] == -1:  # Not yet assigned
                # Start new group with BFS
                group_members = []
                queue = [i]
                visited = set([i])
                
                while queue:
                    current = queue.pop(0)
                    group_members.append(current)
                    
                    # Find all neighbors
                    neighbors = np.where(adjacency[current])[0]
                    for neighbor in neighbors:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
                
                # Assign group ID to all members
                for member in group_members:
                    groups[member] = current_group_id
                
                current_group_id += 1
        
        return groups
    
    def _full_linking(self, adjacency, distance_matrix, strategy='merge'):
        """
        Full linking: all group members must be connected to all others.
        
        Parameters:
        -----------
        adjacency : ndarray
            Boolean adjacency matrix
        distance_matrix : ndarray
            Distance matrix for closest group calculations
        strategy : str
            'merge', 'closest', or 'largest'
        """
        n_gal = len(adjacency)
        
        if strategy == 'merge':
            return self._full_linking_merge(adjacency)
        elif strategy in ['closest', 'largest']:
            return self._full_linking_assign(adjacency, distance_matrix, strategy)
        else:
            raise ValueError("strategy must be 'merge', 'closest', or 'largest'")
    
    def _full_linking_merge(self, adjacency):
        """Original full linking with group merging."""
        n_gal = len(adjacency)
        groups = list(range(n_gal))  # Start with each galaxy as its own group
        
        # Find all cliques (fully connected subgraphs)
        merged = True
        while merged:
            merged = False
            for i in range(n_gal):
                for j in range(i+1, n_gal):
                    if groups[i] != groups[j] and adjacency[i, j]:
                        # Check if merging these groups maintains full connectivity
                        group_i_members = [k for k in range(n_gal) if groups[k] == groups[i]]
                        group_j_members = [k for k in range(n_gal) if groups[k] == groups[j]]
                        
                        # Check if all members of both groups are connected to each other
                        can_merge = True
                        for gi in group_i_members:
                            for gj in group_j_members:
                                if not adjacency[gi, gj]:
                                    can_merge = False
                                    break
                            if not can_merge:
                                break
                        
                        if can_merge:
                            # Merge groups
                            old_group_j = groups[j]
                            for k in range(n_gal):
                                if groups[k] == old_group_j:
                                    groups[k] = groups[i]
                            merged = True
        
        # Renumber groups sequentially
        unique_groups = list(set(groups))
        group_mapping = {old: new for new, old in enumerate(unique_groups)}
        groups = [group_mapping[g] for g in groups]
        
        return groups
    
    def _full_linking_assign(self, adjacency, distance_matrix, strategy):
        """
        Full linking with assignment to closest or largest group.
        
        Strategy:
        1. Find all maximal cliques (fully connected components)
        2. For galaxies that could join multiple groups, assign based on strategy
        """
        n_gal = len(adjacency)
        groups = [-1] * n_gal  # -1 means unassigned
        current_group_id = 0
        
        # First pass: find initial fully connected groups (cliques)
        processed = [False] * n_gal
        
        for i in range(n_gal):
            if processed[i]:
                continue
                
            # Find maximal clique starting from galaxy i
            clique = self._find_maximal_clique(adjacency, i, processed)
            
            if len(clique) > 1:
                # Assign group ID to clique members
                for member in clique:
                    groups[member] = current_group_id
                    processed[member] = True
                current_group_id += 1
        
        # Second pass: assign remaining galaxies to groups or as singles
        for i in range(n_gal):
            if groups[i] == -1:  # Unassigned
                # Find which existing groups this galaxy can join
                candidate_groups = []
                
                for group_id in range(current_group_id):
                    group_members = [j for j in range(n_gal) if groups[j] == group_id]
                    
                    # Check if galaxy i can join this group (connected to all members)
                    can_join = True
                    for member in group_members:
                        if not adjacency[i, member]:
                            can_join = False
                            break
                    
                    if can_join:
                        if strategy == 'closest':
                            # Calculate mean distance to group
                            mean_dist = np.mean([distance_matrix[i, member] for member in group_members])
                            candidate_groups.append((group_id, mean_dist, len(group_members)))
                        elif strategy == 'largest':
                            candidate_groups.append((group_id, 0, len(group_members)))
                
                if candidate_groups:
                    if strategy == 'closest':
                        # Sort by distance (ascending)
                        candidate_groups.sort(key=lambda x: x[1])
                        chosen_group = candidate_groups[0][0]
                    elif strategy == 'largest':
                        # Sort by size (descending)
                        candidate_groups.sort(key=lambda x: x[2], reverse=True)
                        chosen_group = candidate_groups[0][0]
                    
                    groups[i] = chosen_group
                else:
                    # No group to join, create new single-member group
                    groups[i] = current_group_id
                    current_group_id += 1
        
        return groups
    
    def _find_maximal_clique(self, adjacency, start_node, processed):
        """
        Find maximal clique starting from start_node using greedy approach.
        """
        n_gal = len(adjacency)
        clique = [start_node]
        
        # Find all neighbors of start_node
        candidates = [j for j in range(n_gal) 
                     if j != start_node and adjacency[start_node, j] and not processed[j]]
        
        # Greedily add nodes that are connected to all current clique members
        for candidate in candidates:
            # Check if candidate is connected to all clique members
            connected_to_all = True
            for clique_member in clique:
                if not adjacency[candidate, clique_member]:
                    connected_to_all = False
                    break
            
            if connected_to_all:
                clique.append(candidate)
        
        return clique
    
    def _create_catalogs(self):
        """Create group and member catalogs."""
        if self.groups is None:
            raise ValueError("No groups found. Run find_groups() first.")
        
        # Member catalog
        self.member_catalog = self.galaxies.copy()
        self.member_catalog['group_id'] = self.groups
        
        # Group catalog
        group_stats = []
        for group_id in set(self.groups):
            group_members = self.member_catalog[self.member_catalog['group_id'] == group_id]
            
            n_members = len(group_members)
            mean_z = group_members['redshift'].mean()
            mean_ra = group_members['ra'].mean()
            mean_dec = group_members['dec'].mean()
            
            # Velocity dispersion
            if n_members > 1:
                c_kms = 299792.458
                velocities = c_kms * (group_members['redshift'] - mean_z) / (1 + mean_z)
                vel_dispersion = velocities.std()
            else:
                vel_dispersion = 0.0
            
            group_stats.append({
                'group_id': group_id,
                'member_id': ','.join(map(str, group_members['galaxy_id'].values)),
                'n_members': n_members,
                'mean_z': mean_z,
                'mean_ra': mean_ra,
                'mean_dec': mean_dec,
                'velocity_dispersion': vel_dispersion
            })
        
        self.group_catalog = pd.DataFrame(group_stats)
        self.group_catalog = self.group_catalog.sort_values('group_id').reset_index(drop=True)
    
    def save_catalogs(self, output_prefix, format='csv'):
        """
        Save group and member catalogs.
        
        Parameters:
        -----------
        output_prefix : str
            Prefix for output files
        format : str
            'csv' or 'fits'
        """
        if self.group_catalog is None or self.member_catalog is None:
            raise ValueError("No catalogs to save. Run find_groups() first.")
        
        if format == 'csv':
            group_file = f"{output_prefix}_groups.csv"
            member_file = f"{output_prefix}_members.csv"
            
            self.group_catalog.to_csv(group_file, index=False)
            self.member_catalog.to_csv(member_file, index=False)
            
            print(f"Saved group catalog to {group_file}")
            print(f"Saved member catalog to {member_file}")
            
        elif format == 'fits':
            fits_file = f"{output_prefix}_catalogs.fits"
            
            # Convert to astropy tables
            group_table = Table.from_pandas(self.group_catalog)
            member_table = Table.from_pandas(self.member_catalog)
            
            # Create HDU list
            primary_hdu = fits.PrimaryHDU()
            group_hdu = fits.BinTableHDU(group_table, name='GROUPS')
            member_hdu = fits.BinTableHDU(member_table, name='MEMBERS')
            
            hdul = fits.HDUList([primary_hdu, group_hdu, member_hdu])
            hdul.writeto(fits_file, overwrite=True)
            
            print(f"Saved catalogs to {fits_file}")
        
        else:
            raise ValueError("format must be 'csv' or 'fits'")
    
    def print_summary(self):
        """Print summary of group finding results."""
        if self.group_catalog is None:
            print("No groups found yet.")
            return
        
        total_groups = len(self.group_catalog)
        multi_member = len(self.group_catalog[self.group_catalog['n_members'] > 1])
        singles = total_groups - multi_member
        
        print(f"\n=== GROUP FINDING SUMMARY ===")
        print(f"Total galaxies: {len(self.galaxies)}")
        print(f"Total groups: {total_groups}")
        print(f"Multi-member groups: {multi_member}")
        print(f"Single galaxies: {singles}")
        
        if multi_member > 0:
            print(f"\nMulti-member group statistics:")
            multi_groups = self.group_catalog[self.group_catalog['n_members'] > 1]
            print(f"Largest group size: {multi_groups['n_members'].max()}")
            print(f"Mean group size: {multi_groups['n_members'].mean():.1f}")
            print(f"Mean velocity dispersion: {multi_groups['velocity_dispersion'].mean():.1f} km/s")



# Example usage
if __name__ == "__main__":
    # Create sample data
    np.random.seed(42)
    n_galaxies = 100
    
    # Generate some clustered galaxies
    ra = np.random.uniform(0, 10, n_galaxies)  # degrees
    dec = np.random.uniform(-5, 5, n_galaxies)  # degrees
    redshift = np.random.uniform(0.01, 0.1, n_galaxies)
    
    # Add some close pairs
    for i in range(0, 10, 2):
        ra[i+1] = ra[i] + np.random.normal(0, 0.01)  # Close in RA
        dec[i+1] = dec[i] + np.random.normal(0, 0.01)  # Close in Dec
        redshift[i+1] = redshift[i] + np.random.normal(0, 0.001)  # Similar redshift
    
    # Example with angular-only group finding
    print("\n" + "="*50)
    print("ANGULAR-ONLY GROUP FINDING EXAMPLE")
    print("="*50)
    
    # Initialize angular group finder
    agf = AngularGroupFinder()
    
    # Load catalog (RA, Dec only)
    agf.load_catalog(ra, dec)
    
    # Find groups based on angular separation only
    agf.find_groups(
        angular_linking_length=600,  # 600 arcsec = 10 arcmin
        linking_method='chain'
    )
    
    # Print summary
    agf.print_summary()
    
    # Save catalogs
    #agf.save_catalogs("angular_groups", format='csv')
    
    print("\n" + "="*50)
    print("FULL COSMOLOGICAL GROUP FINDING EXAMPLE")
    print("="*50)
    
    # Initialize full group finder
    gf = GalaxyGroupFinder()
    
    # Load catalog
    gf.load_catalog(ra, dec, redshift)
    
    # Find groups with optimizations automatically enabled
    #gf.find_groups(
    #    spatial_linking_length=1.0,  # 1 Mpc comoving
    #    velocity_linking_length=500,  # 500 km/s
    #    distance_type='comoving_mpc',
    #    linking_method='chain'
    #)
    
    # Example with full linking and closest assignment
    gf.find_groups(
        spatial_linking_length=1.0,
        velocity_linking_length=500,
        distance_type='comoving_mpc',
        linking_method='full',
        full_linking_strategy='closest'  # or 'largest' or 'merge'
    )
    
    # Print summary
    gf.print_summary()
    
    # Save catalogs
    #gf.save_catalogs("galaxy_groups", format='csv')