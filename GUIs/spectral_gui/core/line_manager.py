#!/usr/bin/env python
"""
Line Manager - Component for handling spectral line catalogs and identifications.
"""
import os
import numpy as np
import pandas as pd
from pathlib import Path


class LineList:
    """
    A class to manage atomic line lists for different categories of absorbers.
    
    This provides a consistent interface for working with various line catalogs
    and accessing their properties for identification and analysis.
    """
    
    # Singleton pattern to ensure we don't load lists multiple times
    _instances = {}
    _line_lists = {}
    
    @classmethod
    def get_instance(cls, list_name):
        """
        Get an instance of the LineList for a specific list name.
        
        Parameters
        ----------
        list_name : str
            Name of the line list to retrieve
            
        Returns
        -------
        LineList
            Instance of the line list
        """
        if list_name not in cls._instances:
            cls._instances[list_name] = cls(list_name)
        return cls._instances[list_name]
    
    def __init__(self, list_name):
        """
        Initialize a LineList object.
        
        Parameters
        ----------
        list_name : str
            Name of the line list to load
        """
        self.name = list_name
        self._data = None
        self.load()
    
    def load(self):
        """
        Load line data from file or cache.
        
        Returns
        -------
        pandas.DataFrame
            DataFrame containing line data
        """
        if self.name in self._line_lists:
            self._data = self._line_lists[self.name]
            return self._data
        
        # Use the original line reading function from rb_setline module
        try:
            from IGM import rb_setline as line
            self._data = line.read_line_list(self.name)
            self._line_lists[self.name] = self._data
        except ImportError:
            # Fallback to local implementation
            self._data = self._read_line_list(self.name)
            self._line_lists[self.name] = self._data
        
        return self._data
    
    def _read_line_list(self, list_name):
        """
        Fallback method to read line lists if IGM module is not available.
        
        Parameters
        ----------
        list_name : str
            Name of the line list to read
            
        Returns
        -------
        pandas.DataFrame
            DataFrame containing line data
        """
        # Define file paths for common line lists
        line_lists = {
            'LLS': 'linelists/lls.lst',
            'LLS Small': 'linelists/lls_small.lst',
            'DLA': 'linelists/dla.lst',
            'LBG': 'linelists/lbg.lst',
            'Gal': 'linelists/gal.lst',
            'None': 'linelists/none.lst'
        }
        
        # Determine file path
        if list_name in line_lists:
            file_path = line_lists[list_name]
        else:
            raise ValueError(f"Unknown line list: {list_name}")
        
        # Check if we need to find the file in package directories
        if not os.path.exists(file_path):
            # Try to find the file in various locations
            package_dir = Path(__file__).parent.parent.parent
            file_path = package_dir / file_path
            
            if not file_path.exists():
                raise FileNotFoundError(f"Could not find line list file: {file_path}")
        
        # Read the line list file
        columns = ['ion', 'wrest', 'fval', 'gamma']
        try:
            data = pd.read_csv(file_path, sep='\s+', comment='#', names=columns)
            
            # Convert to format compatible with original code
            data_dict = []
            for _, row in data.iterrows():
                data_dict.append({
                    'ion': row['ion'],
                    'wrest': row['wrest'],
                    'f': row['fval'],
                    'gamma': row['gamma']
                })
            
            return data_dict
        except Exception as e:
            print(f"Error reading line list {list_name}: {e}")
            return []
    
    @property
    def data(self):
        """Get line list data."""
        return self._data
    
    def get_line_by_wavelength(self, wavelength, tolerance=0.1):
        """
        Find a line by wavelength.
        
        Parameters
        ----------
        wavelength : float
            Wavelength to search for
        tolerance : float, optional
            Matching tolerance, by default 0.1
            
        Returns
        -------
        dict or None
            Line data if found, None otherwise
        """
        for line in self._data:
            if abs(line['wrest'] - wavelength) <= tolerance:
                return line
        return None
    
    def get_line_by_name(self, ion_name):
        """
        Find a line by ion name.
        
        Parameters
        ----------
        ion_name : str
            Ion name to search for
            
        Returns
        -------
        dict or None
            Line data if found, None otherwise
        """
        for line in self._data:
            if line['ion'] == ion_name:
                return line
        return None
    
    def filter_by_wavelength_range(self, wmin, wmax):
        """
        Filter lines by wavelength range.
        
        Parameters
        ----------
        wmin : float
            Minimum wavelength
        wmax : float
            Maximum wavelength
            
        Returns
        -------
        list
            Filtered line data
        """
        return [line for line in self._data if wmin <= line['wrest'] <= wmax]
    
    def get_rest_wavelengths(self):
        """
        Get all rest wavelengths in the line list.
        
        Returns
        -------
        numpy.ndarray
            Array of rest wavelengths
        """
        return np.array([line['wrest'] for line in self._data])


class LineIdentifier:
    """
    A class to help identify spectral lines based on observed wavelengths.
    """
    
    # Common doublets and multiplets with rest wavelengths
    DOUBLETS = {
        'CIV': [1548.2049, 1550.77845],
        'MgII': [2796.3542699, 2803.5314853],
        'FeII': [2600.1724835, 2586.6495659, 2382.7641781],
        'OVI': [1031.9261, 1037.6167],
        'SiIV': [1393.76018, 1402.77291],
        'NeVIII': [770.409, 780.324],
        'Lyb': [1025.7223, 1215.6701],  # Lyman-beta and Lyman-alpha
        'Lya': [1215.6701, 1025.7223]   # Lyman-alpha and Lyman-beta
    }
    
    @staticmethod
    def identify_doublet(observed_wavelength, doublet_type):
        """
        Identify a potential absorber based on a doublet feature.
        
        Parameters
        ----------
        observed_wavelength : float
            Observed wavelength of the first line in the doublet
        doublet_type : str
            Type of doublet to identify
            
        Returns
        -------
        dict
            Information about the identified doublet
        """
        if doublet_type not in LineIdentifier.DOUBLETS:
            raise ValueError(f"Unknown doublet type: {doublet_type}")
        
        rest_wavelengths = LineIdentifier.DOUBLETS[doublet_type]
        primary_wavelength = rest_wavelengths[0]
        
        # Calculate redshift based on primary line
        redshift = observed_wavelength / primary_wavelength - 1.0
        
        # Calculate expected wavelengths for all lines in the doublet
        expected_wavelengths = [w * (1.0 + redshift) for w in rest_wavelengths]
        
        return {
            'type': doublet_type,
            'redshift': redshift,
            'rest_wavelengths': rest_wavelengths,
            'observed_wavelengths': expected_wavelengths
        }
    
    @staticmethod
    def identify_from_wavelength(observed_wavelength, rest_wavelength):
        """
        Calculate redshift from observed and rest wavelengths.
        
        Parameters
        ----------
        observed_wavelength : float
            Observed wavelength
        rest_wavelength : float
            Rest wavelength
            
        Returns
        -------
        float
            Calculated redshift
        """
        return observed_wavelength / rest_wavelength - 1.0
    
    @staticmethod
    def get_observed_wavelength(rest_wavelength, redshift):
        """
        Calculate observed wavelength from rest wavelength and redshift.
        
        Parameters
        ----------
        rest_wavelength : float
            Rest wavelength
        redshift : float
            Redshift
            
        Returns
        -------
        float
            Calculated observed wavelength
        """
        return rest_wavelength * (1.0 + redshift)
    
    @staticmethod
    def find_potential_matches(observed_wavelength, tolerance=2.0, line_lists=None):
        """
        Find potential line matches for an observed wavelength.
        
        Parameters
        ----------
        observed_wavelength : float
            Observed wavelength
        tolerance : float, optional
            Matching tolerance in Angstroms, by default 2.0
        line_lists : list, optional
            List of line list names to search, by default None (uses all available)
            
        Returns
        -------
        list
            List of potential matches with their redshifts
        """
        if line_lists is None:
            line_lists = ['LLS', 'DLA', 'LBG', 'Gal']
        
        matches = []
        
        for list_name in line_lists:
            try:
                line_list = LineList.get_instance(list_name)
                for line in line_list.data:
                    # Calculate what redshift would make this line match
                    redshift = observed_wavelength / line['wrest'] - 1.0
                    
                    # Skip negative redshifts
                    if redshift < 0:
                        continue
                    
                    # Check how well other lines in this list would match at this redshift
                    other_matches = 0
                    for other_line in line_list.data:
                        if other_line == line:
                            continue
                        
                        expected_wavelength = other_line['wrest'] * (1.0 + redshift)
                        # This is a simplified check - a real implementation would
                        # check against the actual spectrum
                        if expected_wavelength > 0:
                            other_matches += 1
                    
                    matches.append({
                        'ion': line['ion'],
                        'rest_wavelength': line['wrest'],
                        'redshift': redshift,
                        'line_list': list_name,
                        'other_potential_matches': other_matches
                    })
            except Exception as e:
                print(f"Error searching {list_name}: {e}")
        
        # Sort by number of potential matches (higher is better)
        matches.sort(key=lambda x: x['other_potential_matches'], reverse=True)
        
        return matches