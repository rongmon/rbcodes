#!/usr/bin/env python
"""
Absorber data model - Core component for representing spectral absorbers.
"""
import os
import numpy as np
import pandas as pd

class Absorber:
    """
    A class to represent an absorption system with multiple lines.
    
    Parameters
    ----------
    redshift : float
        Redshift of the absorber
    name : str, optional
        Name of the absorber, by default None
    color : str, optional
        Color for visualization, by default 'white'
    line_list : str, optional
        Line list identifier, by default 'None'
    """
    
    def __init__(self, redshift, name=None, color='white', line_list='None'):
        """Initialize an Absorber object."""
        self.redshift = float(redshift)
        self.name = name if name else f"z={redshift:.4f}"
        self.color = color
        self.line_list = line_list
        
        # Initialize empty list of lines
        self._lines = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs', 'Flag'])
    
    @property
    def lines(self):
        """Get DataFrame containing absorption lines."""
        return self._lines
    
    @property
    def line_count(self):
        """Get number of lines in this absorber."""
        return len(self._lines)
    
    def add_line(self, name, rest_wavelength, flag=1):
        """
        Add a line to this absorber.
        
        Parameters
        ----------
        name : str
            Line name or identifier
        rest_wavelength : float
            Rest wavelength of the line
        flag : int, optional
            Detection flag (0=non-detection, 1=detection, 2=blended, 3=low confidence),
            by default 1
            
        Returns
        -------
        int
            Index of the newly added line
        """
        # Calculate observed wavelength
        wave_obs = rest_wavelength * (1.0 + self.redshift)
        
        # Determine display name based on flag
        display_name = name
        if flag == 2:
            display_name = f"{name} [b]"
        elif flag == 3:
            display_name = f"{name} [p]"
        
        # Add to DataFrame
        new_row = pd.Series({
            'Name': display_name,
            'Wave_obs': wave_obs,
            'Zabs': self.redshift,
            'Flag': flag,
            'RestWave': rest_wavelength
        })
        
        self._lines = self._lines.append(new_row, ignore_index=True)
        return len(self._lines) - 1
    
    def remove_line(self, index):
        """
        Remove a line from this absorber.
        
        Parameters
        ----------
        index : int
            Index of the line to remove
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        if 0 <= index < len(self._lines):
            self._lines = self._lines.drop(index).reset_index(drop=True)
            return True
        return False
    
    def update_line(self, index, **kwargs):
        """
        Update properties of a line.
        
        Parameters
        ----------
        index : int
            Index of the line to update
        **kwargs
            Keyword arguments of properties to update
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        if 0 <= index < len(self._lines):
            # Update the row with new values
            for key, value in kwargs.items():
                if key in self._lines.columns:
                    self._lines.at[index, key] = value
            
            # Update display name if flag changed
            if 'Flag' in kwargs:
                name = self._lines.at[index, 'Name'].split(' [')[0]
                if kwargs['Flag'] == 2:
                    self._lines.at[index, 'Name'] = f"{name} [b]"
                elif kwargs['Flag'] == 3:
                    self._lines.at[index, 'Name'] = f"{name} [p]"
                elif kwargs['Flag'] == 1:
                    self._lines.at[index, 'Name'] = name
            
            return True
        return False
    
    def to_dict(self):
        """
        Convert absorber to dictionary for serialization.
        
        Returns
        -------
        dict
            Dictionary representation of absorber
        """
        return {
            'redshift': self.redshift,
            'name': self.name,
            'color': self.color,
            'line_list': self.line_list,
            'lines': self._lines.to_dict('records')
        }
    
    @classmethod
    def from_dict(cls, data):
        """
        Create absorber from dictionary.
        
        Parameters
        ----------
        data : dict
            Dictionary representation of absorber
            
        Returns
        -------
        Absorber
            New absorber object
        """
        absorber = cls(
            redshift=data['redshift'],
            name=data.get('name'),
            color=data.get('color', 'white'),
            line_list=data.get('line_list', 'None')
        )
        
        # Add lines
        if 'lines' in data:
            absorber._lines = pd.DataFrame(data['lines'])
        
        return absorber


class AbsorberManager:
    """
    A class to manage multiple absorbers.
    
    This class provides methods for adding, removing, and querying absorbers,
    as well as for importing/exporting them to various formats.
    """
    
    def __init__(self):
        """Initialize an empty absorber manager."""
        self._absorbers = []
        
        # Keep track of visualization objects
        self._line_plots = []
        self._text_objects = []
    
    @property
    def absorbers(self):
        """Get list of absorbers."""
        return self._absorbers
    
    @property
    def redshifts(self):
        """Get list of absorber redshifts."""
        return [abs.redshift for abs in self._absorbers]
    
    @property
    def count(self):
        """Get number of absorbers."""
        return len(self._absorbers)
    
    def add_absorber(self, absorber):
        """
        Add an absorber to the manager.
        
        Parameters
        ----------
        absorber : Absorber
            Absorber object to add
            
        Returns
        -------
        int
            Index of the newly added absorber
        """
        self._absorbers.append(absorber)
        # Initialize visualization storage for this absorber
        self._line_plots.append([])
        self._text_objects.append([])
        return len(self._absorbers) - 1
    
    def remove_absorber(self, index):
        """
        Remove an absorber from the manager.
        
        Parameters
        ----------
        index : int
            Index of the absorber to remove
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        if 0 <= index < len(self._absorbers):
            self._absorbers.pop(index)
            self._line_plots.pop(index)
            self._text_objects.pop(index)
            return True
        return False
    
    def get_absorber(self, index):
        """
        Get an absorber by index.
        
        Parameters
        ----------
        index : int
            Index of the absorber
            
        Returns
        -------
        Absorber or None
            Absorber object if found, None otherwise
        """
        if 0 <= index < len(self._absorbers):
            return self._absorbers[index]
        return None
    
    def get_absorber_by_redshift(self, redshift, tolerance=1e-4):
        """
        Find an absorber by redshift.
        
        Parameters
        ----------
        redshift : float
            Redshift to search for
        tolerance : float, optional
            Matching tolerance, by default 1e-4
            
        Returns
        -------
        tuple or None
            (index, absorber) if found, None otherwise
        """
        for i, abs in enumerate(self._absorbers):
            if abs.redshift >= redshift - tolerance and abs.redshift <= redshift + tolerance:
                return i, abs
        return None
    
    def export_to_dataframe(self):
        """
        Export absorbers to pandas DataFrame.
        
        Returns
        -------
        pandas.DataFrame
            DataFrame containing absorber data
        """
        data = []
        for abs in self._absorbers:
            data.append({
                'Zabs': abs.redshift,
                'list': abs.line_list,
                'color': abs.color
            })
        return pd.DataFrame(data)
    
    def export_lines_to_dataframe(self):
        """
        Export all lines from all absorbers to pandas DataFrame.
        
        Returns
        -------
        pandas.DataFrame
            DataFrame containing line data
        """
        all_lines = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs', 'Flag'])
        
        for abs in self._absorbers:
            if abs.line_count > 0:
                all_lines = all_lines.append(abs.lines, ignore_index=True)
        
        return all_lines
    
    def import_from_dataframe(self, df):
        """
        Import absorbers from pandas DataFrame.
        
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing absorber data
            
        Returns
        -------
        int
            Number of absorbers imported
        """
        # Clear existing absorbers
        self._absorbers = []
        self._line_plots = []
        self._text_objects = []
        
        # Import new absorbers
        for _, row in df.iterrows():
            abs = Absorber(
                redshift=row['Zabs'],
                color=row['color'],
                line_list=row['list']
            )
            self.add_absorber(abs)
        
        return len(self._absorbers)
    
    def import_lines(self, lines_df):
        """
        Import lines for absorbers from pandas DataFrame.
        
        Parameters
        ----------
        lines_df : pandas.DataFrame
            DataFrame containing line data
            
        Returns
        -------
        int
            Number of lines imported
        """
        count = 0
        
        # Group by redshift
        grouped = lines_df.groupby('Zabs')
        
        for zabs, group in grouped:
            # Find matching absorber
            match = self.get_absorber_by_redshift(zabs)
            
            if match:
                abs_idx, abs = match
                
                # Clear existing lines for this absorber
                abs._lines = pd.DataFrame(columns=['Name', 'Wave_obs', 'Zabs', 'Flag'])
                
                # Add lines
                for _, line in group.iterrows():
                    name = line['Name']
                    wave_obs = line['Wave_obs']
                    
                    # Determine flag based on name
                    flag = 1
                    if '[b]' in name:
                        flag = 2
                        name = name.replace(' [b]', '')
                    elif '[p]' in name:
                        flag = 3
                        name = name.replace(' [p]', '')
                    
                    # Calculate rest wavelength
                    rest_wave = wave_obs / (1.0 + zabs)
                    
                    # Add line
                    abs.add_line(name, rest_wave, flag)
                    count += 1
        
        return count
    
    def save(self, file_path):
        """
        Save absorbers to file.
        
        Parameters
        ----------
        file_path : str
            Path to save file
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        try:
            # Export absorbers
            abs_df = self.export_to_dataframe()
            abs_df.to_csv(file_path, index=False)
            
            # Export lines
            lines_df = self.export_lines_to_dataframe()
            lines_path = file_path.replace('.csv', '_lines.txt')
            if lines_df.shape[0] > 0:
                lines_df.to_csv(lines_path, sep=' ', index=False)
                
            return True
        except Exception as e:
            print(f"Error saving absorbers: {e}")
            return False
    
    def load(self, file_path):
        """
        Load absorbers from file.
        
        Parameters
        ----------
        file_path : str
            Path to load file
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        try:
            # Import absorbers
            abs_df = pd.read_csv(file_path)
            self.import_from_dataframe(abs_df)
            
            # Import lines if available
            lines_path = file_path.replace('.csv', '_lines.txt')
            if os.path.exists(lines_path):
                lines_df = pd.read_csv(lines_path, sep=' ')
                self.import_lines(lines_df)
                
            return True
        except Exception as e:
            print(f"Error loading absorbers: {e}")
            return False
    
    def store_visualization(self, absorber_index, line_plots, text_objects):
        """
        Store visualization objects for an absorber.
        
        Parameters
        ----------
        absorber_index : int
            Index of the absorber
        line_plots : list
            List of line plot objects
        text_objects : list
            List of text objects
        """
        if 0 <= absorber_index < len(self._absorbers):
            self._line_plots[absorber_index] = line_plots
            self._text_objects[absorber_index] = text_objects
    
    def get_visualization(self, absorber_index):
        """
        Get visualization objects for an absorber.
        
        Parameters
        ----------
        absorber_index : int
            Index of the absorber
            
        Returns
        -------
        tuple or None
            (line_plots, text_objects) if found, None otherwise
        """
        if 0 <= absorber_index < len(self._absorbers):
            return self._line_plots[absorber_index], self._text_objects[absorber_index]
        return None, None