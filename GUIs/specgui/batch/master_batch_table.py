# rbcodes/GUIs/specgui/batch/master_batch_table.py
import os
import pandas as pd
import numpy as np
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from PyQt5.QtCore import QObject, pyqtSignal


class MasterBatchTable(QObject):
    """Master table using pandas DataFrame as single source of truth."""
    
    # Signals
    item_added = pyqtSignal(int)  # row_index
    item_removed = pyqtSignal(int)  # row_index
    item_updated = pyqtSignal(int, str)  # row_index, what_changed
    table_cleared = pyqtSignal()
    validation_changed = pyqtSignal()
    
    def __init__(self):
        super().__init__()
        # Create empty DataFrame with all required columns
        self.df = self._create_empty_dataframe()
        self.metadata = {
            'version': '1.0',
            'created': datetime.now().isoformat(),
            'description': '',
            'software_version': '1.0'
        }
        self.rb_spec_objects = {}  # Dict[int, rb_spec] - keyed by row index
    
    def _create_empty_dataframe(self):
        """Create empty DataFrame with all required columns."""
        columns = [
            # Template parameters
            'filename', 'redshift', 'transition', 'transition_name',
            'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax', 'linelist', 'method',
            
            # Analysis state
            'continuum_method', 'continuum_order', 'continuum_masks', 'continuum_fit_params',
            'use_weights', 'optimize_cont', 'calculate_snr', 'binsize', 'processing_status', 'last_modified', 'error_message',
            
            # Results
            'W', 'W_e', 'N', 'N_e', 'logN', 'logN_e', 'vel_centroid', 'vel_disp', 'SNR', 'calculation_timestamp'
        ]
        
        df = pd.DataFrame(columns=columns)
        
        # Set default values for new rows
        df = df.astype({
            'redshift': 'float64',
            'transition': 'float64',
            'slice_vmin': 'int64',
            'slice_vmax': 'int64',
            'ew_vmin': 'int64',
            'ew_vmax': 'int64',
            'continuum_order': 'int64',
            'binsize': 'int64',
            'W': 'float64',
            'W_e': 'float64',
            'N': 'float64',
            'N_e': 'float64',
            'logN': 'float64',
            'logN_e': 'float64',
            'vel_centroid': 'float64',
            'vel_disp': 'float64',
            'SNR': 'float64'
        }, errors='ignore')
        
        return df
    
    def add_item(self, template_params) -> int:
        """Add a new batch item. Returns row index."""
        # Convert template_params to dict if it's an object
        if hasattr(template_params, '__dict__'):
            params = template_params.__dict__
        else:
            params = template_params
        
        # Create new row with defaults
        new_row = {
            # Template parameters
            'filename': params.get('filename', ''),
            'redshift': params.get('redshift', 0.0),
            'transition': params.get('transition', 0.0),
            'transition_name': params.get('transition_name', ''),
            'slice_vmin': params.get('slice_vmin', -1500),
            'slice_vmax': params.get('slice_vmax', 1500),
            'ew_vmin': params.get('ew_vmin', -200),
            'ew_vmax': params.get('ew_vmax', 200),
            'linelist': params.get('linelist', 'atom'),
            'method': params.get('method', 'closest'),
            
            # Analysis defaults
            'continuum_method': 'polynomial',
            'continuum_order': 3,
            'continuum_masks': '[]',  # Store as JSON string
            'continuum_fit_params': '{}',
            'use_weights': False,
            'optimize_cont': True,  # NEW: Default to True for optimized continuum fit
            'calculate_snr': True,
            'binsize': 3,
            'processing_status': 'ready',
            'last_modified': datetime.now().isoformat(),
            'error_message': '',
            
            # Results defaults
            'W': 0.0, 'W_e': 0.0, 'N': 0.0, 'N_e': 0.0,
            'logN': 0.0, 'logN_e': 0.0, 'vel_centroid': 0.0,
            'vel_disp': 0.0, 'SNR': 0.0, 'calculation_timestamp': ''
        }
        
        # Add row to DataFrame
        row_index = len(self.df)
        self.df.loc[row_index] = new_row
        
        self.item_added.emit(row_index)
        return row_index
    
    def remove_item(self, row_index: int) -> bool:
        """Remove a batch item by row index."""
        if 0 <= row_index < len(self.df):
            # Remove rb_spec object if it exists
            if row_index in self.rb_spec_objects:
                del self.rb_spec_objects[row_index]
            
            # Remove row from DataFrame
            self.df.drop(index=row_index, inplace=True)
            self.df.reset_index(drop=True, inplace=True)
            
            # Update rb_spec_objects keys (shift indices down)
            new_rb_spec_objects = {}
            for old_idx, spec_obj in self.rb_spec_objects.items():
                new_idx = old_idx if old_idx < row_index else old_idx - 1
                if new_idx >= 0:
                    new_rb_spec_objects[new_idx] = spec_obj
            self.rb_spec_objects = new_rb_spec_objects
            
            self.item_removed.emit(row_index)
            return True
        return False
    
    def get_item_count(self) -> int:
        """Get the total number of items."""
        return len(self.df)
    
    def get_all_items(self):
        """Get all items as a list of row dictionaries."""
        items = []
        for idx, row in self.df.iterrows():
            # Create a simple object with the row data
            item = type('BatchItem', (), {
                'id': str(idx),  # Use row index as ID
                'template': type('Template', (), row.to_dict())(),
                'analysis': type('Analysis', (), {
                    'processing_status': row['processing_status'],
                    'continuum_method': row['continuum_method'],
                    'continuum_order': row['continuum_order'],
                    'continuum_masks': eval(row['continuum_masks']) if row['continuum_masks'] else [],
                    'use_weights': row['use_weights'],
                    'optimize_cont': row['optimize_cont'],  # NEW field
                    'calculate_snr': row['calculate_snr'],
                    'binsize': row['binsize'],
                    'error_message': row['error_message'],
                    'last_modified': row['last_modified']
                })(),
                'results': type('Results', (), {
                    'W': row['W'], 'W_e': row['W_e'], 'N': row['N'], 'N_e': row['N_e'],
                    'logN': row['logN'], 'logN_e': row['logN_e'],
                    'vel_centroid': row['vel_centroid'], 'vel_disp': row['vel_disp'],
                    'SNR': row['SNR'], 'calculation_timestamp': row['calculation_timestamp']
                })()
            })()
            items.append(item)
        return items
    
    def get_item(self, row_index: int):
        """Get a single item by row index."""
        if 0 <= row_index < len(self.df):
            items = self.get_all_items()
            return items[row_index]
        return None
    
    def update_template(self, row_index: int, **kwargs) -> bool:
        """Update template parameters for an item."""
        if 0 <= row_index < len(self.df):
            for field, value in kwargs.items():
                if field in self.df.columns:
                    self.df.loc[row_index, field] = value
            
            self.df.loc[row_index, 'last_modified'] = datetime.now().isoformat()
            self.item_updated.emit(row_index, 'template_update')
            return True
        return False
    
    def update_analysis(self, row_index: int, **kwargs) -> bool:
        """Update analysis parameters for an item."""
        if 0 <= row_index < len(self.df):
            for field, value in kwargs.items():
                if field == 'continuum_masks':
                    # Store as JSON string
                    import json
                    self.df.loc[row_index, field] = json.dumps(value)
                elif field in self.df.columns:
                    self.df.loc[row_index, field] = value
            
            self.df.loc[row_index, 'last_modified'] = datetime.now().isoformat()
            self.item_updated.emit(row_index, 'analysis_update')
            return True
        return False
    
    def update_results(self, row_index: int, **kwargs) -> bool:
        """Update results for an item."""
        if 0 <= row_index < len(self.df):
            for field, value in kwargs.items():
                if field in self.df.columns:
                    self.df.loc[row_index, field] = value
            
            self.df.loc[row_index, 'calculation_timestamp'] = datetime.now().isoformat()
            self.item_updated.emit(row_index, 'results_update')
            return True
        return False
    
    def set_rb_spec_object(self, row_index: int, spec_object) -> bool:
        """Store the rb_spec object for an item."""
        if 0 <= row_index < len(self.df):
            self.rb_spec_objects[row_index] = spec_object
            return True
        return False
    
    def get_rb_spec_object(self, row_index: int):
        """Get the rb_spec object for an item."""
        return self.rb_spec_objects.get(row_index, None)
    
    def clear_all(self):
        """Clear all items."""
        self.df = self._create_empty_dataframe()
        self.rb_spec_objects = {}
        self.table_cleared.emit()
    
    def validate_all(self) -> Tuple[List[int], List[Tuple[int, str]]]:
        """Validate all items. Returns (valid_indices, invalid_items_with_errors)."""
        valid_indices = []
        invalid_items = []
        
        for idx, row in self.df.iterrows():
            try:
                # Check file exists
                if not os.path.exists(row['filename']):
                    invalid_items.append((idx, f"File does not exist: {row['filename']}"))
                    continue
                
                # Check redshift range
                if not (0 <= row['redshift'] <= 10):
                    invalid_items.append((idx, f"Redshift {row['redshift']} is outside reasonable range (0-10)"))
                    continue
                
                # Check transition wavelength
                if not (1 <= row['transition'] <= 10000):
                    invalid_items.append((idx, f"Transition wavelength {row['transition']} Å is outside reasonable range"))
                    continue
                
                # Check velocity ranges
                if row['slice_vmin'] >= row['slice_vmax']:
                    invalid_items.append((idx, f"Slice velocity range invalid: vmin >= vmax"))
                    continue
                
                if row['ew_vmin'] >= row['ew_vmax']:
                    invalid_items.append((idx, f"EW velocity range invalid: vmin >= vmax"))
                    continue
                
                valid_indices.append(idx)
                
            except Exception as e:
                invalid_items.append((idx, f"Validation error: {str(e)}"))
        
        return valid_indices, invalid_items
    
    def get_items_by_status(self, status: str) -> List[int]:
        """Get row indices by processing status."""
        return self.df[self.df['processing_status'] == status].index.tolist()
    
    def to_dict(self) -> Dict:
        """Convert entire table to dictionary for serialization."""
        return {
            'metadata': self.metadata,
            'dataframe': self.df.to_dict('records')  # Convert to list of row dicts
        }
    
    def from_dict(self, data: Dict):
        """Load table from dictionary."""
        self.clear_all()
        
        # Load metadata
        self.metadata = data.get('metadata', {})
        if 'version' not in self.metadata:
            self.metadata['version'] = '1.0'
        
        # Load DataFrame
        if 'dataframe' in data and data['dataframe']:
            self.df = pd.DataFrame(data['dataframe'])
            
            # Add missing columns with defaults for backward compatibility
            if 'optimize_cont' not in self.df.columns:
                self.df['optimize_cont'] = True  # Default to True for existing configurations
            
            # Fix data types after loading from JSON
            if not self.df.empty:
                # Convert numeric columns that might be loaded as strings/objects
                numeric_int_columns = ['slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax', 'continuum_order', 'binsize']
                numeric_float_columns = ['redshift', 'transition', 'W', 'W_e', 'N', 'N_e', 'logN', 'logN_e', 
                                       'vel_centroid', 'vel_disp', 'SNR']
                
                # Convert integer columns
                for col in numeric_int_columns:
                    if col in self.df.columns:
                        self.df[col] = pd.to_numeric(self.df[col], errors='coerce').fillna(0).astype('int64')
                
                # Convert float columns  
                for col in numeric_float_columns:
                    if col in self.df.columns:
                        self.df[col] = pd.to_numeric(self.df[col], errors='coerce').fillna(0.0).astype('float64')
                
                # Convert boolean columns
                bool_columns = ['use_weights', 'calculate_snr', 'optimize_cont']
                for col in bool_columns:
                    if col in self.df.columns:
                        self.df[col] = self.df[col].astype('bool')


    def import_template_csv(self, filepath: str) -> Tuple[bool, str]:
        """Import template from CSV file."""
        try:
            import pandas as pd
            
            # Read CSV file
            df = pd.read_csv(filepath)
            
            # Validate required columns
            required_columns = ['filename', 'redshift', 'transition', 'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax']
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                return False, f"Missing required columns: {missing_columns}"
            
            # Clear existing items
            self.clear_all()
            
            imported_count = 0
            for _, row in df.iterrows():
                try:
                    # Create template parameters from CSV row
                    template_params = {
                        'filename': str(row['filename']),
                        'redshift': float(row['redshift']),
                        'transition': float(row['transition']),
                        'transition_name': str(row.get('transition_name', f"λ {row['transition']:.2f}")),
                        'slice_vmin': int(row['slice_vmin']),
                        'slice_vmax': int(row['slice_vmax']),
                        'ew_vmin': int(row['ew_vmin']),
                        'ew_vmax': int(row['ew_vmax']),
                        'linelist': str(row.get('linelist', 'atom')),
                        'method': str(row.get('method', 'closest'))
                    }
                    
                    # Validate file exists
                    if not os.path.exists(template_params['filename']):
                        print(f"Warning: File does not exist: {template_params['filename']}")
                    
                    # Add item to master table
                    self.add_item(template_params)
                    imported_count += 1
                    
                except Exception as e:
                    print(f"Error importing row {imported_count + 1}: {str(e)}")
                    continue
            
            return True, f"Successfully imported {imported_count} items from CSV"
            
        except Exception as e:
            return False, f"Error reading CSV file: {str(e)}"
    
    def export_template_csv(self, filepath: str) -> bool:
        """Export just the template part as CSV."""
        try:
            import csv
            
            items = self.get_all_items()
            if not items:
                return False
            
            with open(filepath, 'w', newline='') as csvfile:
                fieldnames = ['filename', 'redshift', 'transition', 'transition_name', 
                             'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax', 'linelist', 'method']
                
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for item in items:
                    row = {
                        'filename': item.template.filename,
                        'redshift': item.template.redshift,
                        'transition': item.template.transition,
                        'transition_name': item.template.transition_name,
                        'slice_vmin': item.template.slice_vmin,
                        'slice_vmax': item.template.slice_vmax,
                        'ew_vmin': item.template.ew_vmin,
                        'ew_vmax': item.template.ew_vmax,
                        'linelist': item.template.linelist,
                        'method': item.template.method
                    }
                    writer.writerow(row)
            
            return True
        except Exception as e:
            print(f"Error exporting template CSV: {str(e)}")
            return False    