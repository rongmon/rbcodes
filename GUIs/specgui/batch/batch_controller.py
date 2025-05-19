# rbcodes/GUIs/specgui/batch/batch_controller.py
import os
import csv
import json
import numpy as np
from PyQt5.QtCore import QObject, pyqtSignal
import traceback  # Add this for better error reporting

from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object


class BatchController(QObject):
    """Controller for managing batch processing of absorption line spectra."""
    
    # Signals
    status_updated = pyqtSignal(str)  # For status updates
    batch_progress = pyqtSignal(int, int)  # current, total
    batch_item_completed = pyqtSignal(int, dict)  # index, results
    batch_item_failed = pyqtSignal(int, str)  # index, error message
    
    def __init__(self):
        super().__init__()
        self.batch_items = []  # List of batch items to process
        self.results = []  # List of processed results
        self.error_log = []  # List of errors
        
        # Default batch settings
        self.batch_settings = {
            'continuum_method': 'polynomial',  # 'polynomial', 'flat', 'interactive'
            'polynomial_order': 3,
            'use_weights': False,
            'calculate_snr': True,
            'binsize': 3,
            'save_individual_json': True,
            'output_directory': '',
        }
    
    def update_batch_items(self, items):
        """Update the list of batch items."""
        self.batch_items = items
        self.status_updated.emit(f"Batch configuration updated with {len(items)} items.")
    
    def update_batch_settings(self, settings):
        """Update the batch processing settings."""
        self.batch_settings.update(settings)
        self.status_updated.emit("Batch settings updated.")
    

    def process_batch(self, selected_indices=None):
        """
        Process all batch items or selected indices.
        
        Parameters:
        -----------
        selected_indices : list or None
            If provided, only process the items at these indices. 
            If None, process all items.
        
        Returns:
        --------
        success : bool
            True if batch processing completed, False if it failed or was cancelled.
        results : list
            List of processing results.
        """
        if not self.batch_items:
            self.status_updated.emit("No batch items to process.")
            return False, []
        
        # Clear previous results
        self.results = []
        self.error_log = []
        
        # Determine which items to process
        if selected_indices is not None:
            items_to_process = [self.batch_items[i] for i in selected_indices if i <     len(self.batch_items)]
            indices_to_process = selected_indices
        else:
            items_to_process = self.batch_items
            indices_to_process = range(len(self.batch_items))
        
        total_items = len(items_to_process)
        if total_items == 0:
            self.status_updated.emit("No valid items to process.")
            return False, []
        
        self.status_updated.emit(f"Starting batch processing of {total_items}     items...")
        
        # Create output directory if needed
        if self.batch_settings['save_individual_json'] and     self.batch_settings['output_directory']:
            output_dir = self.batch_settings['output_directory']
            if not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir)
                    self.status_updated.emit(f"Created output directory:     {output_dir}")
                except Exception as e:
                    self.status_updated.emit(f"Error creating output directory: {str(    e)}")
                    self.error_log.append({
                        'item_index': -1,  # Special value for general errors
                        'item': None,
                        'error': f"Failed to create output directory: {str(e)}"
                    })
                    return False, []
        
        # Process each item
        for i, (item, original_index) in enumerate(zip(items_to_process,     indices_to_process)):
            try:
                item_desc = self._get_item_description(item)
                self.status_updated.emit(f"Processing item {i+1}/{total_items}:     {item_desc}")
                self.batch_progress.emit(i+1, total_items)
                
                # Process the item based on its type
                if item['type'] == 'multiple_systems':
                    result = self._process_system(item)
                elif item['type'] == 'multiple_transitions':
                    result = self._process_transition(item)
                elif item['type'] == 'multiple_files':
                    result = self._process_file(item)
                else:
                    raise ValueError(f"Unknown batch item type: {item['type']}")
                
                # Add to results
                self.results.append(result)
                
                # Emit signal for item completion
                self.batch_item_completed.emit(original_index, result)
                
            except Exception as e:
                error_msg = f"Error processing item {i+1}/{total_items}: {str(e)}"
                self.status_updated.emit(error_msg)
                self.error_log.append({
                    'item_index': original_index,
                    'item': item,
                    'error': str(e)
                })
                self.batch_item_failed.emit(original_index, str(e))
                
                # Add a placeholder result with error status
                self.results.append({
                    'type': item['type'],
                    'filename': item.get('filename', 'Unknown'),
                    'redshift': item.get('redshift', 0),
                    'transition': item.get('transition', 0),
                    'transition_name': item.get('transition_name', 'Unknown'),
                    'vmin': item.get('vmin', 0),
                    'vmax': item.get('vmax', 0),
                    'W': 0.0,
                    'W_e': 0.0,
                    'logN': 0.0,
                    'logN_e': 0.0,
                    'status': 'error',
                    'error': str(e),
                    'spec_object': None
                })
        
        # Batch processing complete
        if self.error_log:
            self.status_updated.emit(f"Batch processing completed with {len(    self.error_log)} errors.")
        else:
            self.status_updated.emit("Batch processing completed successfully.")
        
        return len(self.error_log) < len(items_to_process), self.results    

    
    def _get_item_description(self, item):
        """Get a description string for a batch item."""
        if item['type'] == 'multiple_systems':
            return f"z={item['redshift']}, {item['transition']} Å"
        elif item['type'] == 'multiple_transitions':
            return f"z={item['redshift']}, {item['transition']} Å"
        elif item['type'] == 'multiple_files':
            return f"{os.path.basename(item['filename'])}, z={item['redshift']}, {item['transition']} Å"
        return "Unknown item"
    
# batch_controller.py - update the processing methods

    def _process_system(self, item):
        """
        Process a batch item with multiple absorption systems in one file.
        
        This processes one file at a specific redshift and transition.
        """
        try:
            # Extract parameters
            filename = item['filename']
            redshift = item['redshift']
            transition = item['transition']
            transition_name = item.get('transition_name', 'Unknown')
            vmin = item['vmin']
            vmax = item['vmax']
            linelist = item.get('linelist', 'atom')
            
            # Status update
            self.status_updated.emit(f"Processing {os.path.basename(filename)},     z={redshift:.6f}, {transition_name}")
            
            # Load the spectrum file
            if filename.lower().endswith('.json'):
                # Load saved rb_spec object
                from rbcodes.GUIs.rb_spec import load_rb_spec_object
                spec = load_rb_spec_object(filename)
            else:
                # Load spectrum from file
                from rbcodes.GUIs.rb_spec import rb_spec
                # Determine file type based on extension
                ext = os.path.splitext(filename)[1].lower()
                if ext == '.fits':
                    filetype ='linetools'# 'fits'
                elif ext in ['.txt', '.dat']:
                    filetype = 'ascii'
                else:
                    filetype = None  # Auto-detect
                    
                spec = rb_spec.from_file(filename, filetype=filetype)
            
            # Apply redshift
            spec.shift_spec(redshift)
            
            # Slice spectrum around transition
            spec.slice_spec(
                transition,
                vmin, vmax,
                use_vel=True,
                linelist=linelist
            )
            
            # Fit continuum based on settings
            method = self.batch_settings['continuum_method']
            if method == 'polynomial':
                # Use polynomial fitting
                spec.fit_continuum(
                    Legendre=self.batch_settings['polynomial_order'],
                    use_weights=self.batch_settings['use_weights'],
                    sigma_clip=True  # Use sigma clipping for robustness
                )
            elif method == 'flat':
                # Use flat continuum
                spec.cont = np.ones_like(spec.flux_slice)
                spec.fnorm = spec.flux_slice.copy()
                spec.enorm = spec.error_slice.copy()
            else:  # 'interactive'
                # Can't use interactive in batch mode, use polynomial instead
                self.status_updated.emit("Interactive mode not available in batch     processing. Using polynomial fit.")
                spec.fit_continuum(
                    Legendre=self.batch_settings['polynomial_order'],
                    use_weights=self.batch_settings['use_weights'],
                    sigma_clip=True
                )
            
            # Compute EW
            spec.compute_EW(
                transition,
                vmin=vmin,
                vmax=vmax,
                SNR=self.batch_settings['calculate_snr'],
                _binsize=self.batch_settings['binsize']
            )
            
            # Save results as JSON if enabled
            if self.batch_settings['save_individual_json'] and     self.batch_settings['output_directory']:
                # Create output filename
                output_dir = self.batch_settings['output_directory']
                basename = os.path.splitext(os.path.basename(filename))[0]
                output_file = f"{basename}_{transition_name}_z{redshift:.3f}.json"
                output_path = os.path.join(output_dir, output_file)
                
                # Save the rb_spec object
                spec.save_slice(output_path, file_format='json')
            
            # Create result dictionary
            result = {
                'type': 'multiple_systems',
                'filename': filename,
                'redshift': redshift,
                'transition': transition,
                'transition_name': transition_name,
                'vmin': vmin,
                'vmax': vmax,
                'W': spec.W,
                'W_e': spec.W_e,
                'logN': spec.logN,
                'logN_e': spec.logN_e,
                'vel_centroid': getattr(spec, 'vel_centroid', 0),
                'SNR': getattr(spec, 'SNR', 0),
                'status': 'success',
                'spec_object': spec  # Store the rb_spec object for later use
            }
            
            return result
        
        except Exception as e:
            # Log the error
            error_msg = f"Error processing {os.path.basename(item.get('filename',     'Unknown'))}: {str(e)}"
            self.status_updated.emit(error_msg)
            print(f"Detailed error: {e}")
            import traceback
            traceback.print_exc()
            
            # Return result with error status
            return {
                'type': item['type'],
                'filename': item.get('filename', 'Unknown'),
                'redshift': item.get('redshift', 0),
                'transition': item.get('transition', 0),
                'transition_name': item.get('transition_name', 'Unknown'),
                'vmin': item.get('vmin', 0),
                'vmax': item.get('vmax', 0),
                'W': 0.0,
                'W_e': 0.0,
                'logN': 0.0,
                'logN_e': 0.0,
                'status': 'error',
                'error': str(e),
                'spec_object': None
            }
    
    def _process_transition(self, item):
        """
        Process a batch item with multiple transitions at one redshift.
        
        This is very similar to _process_system since we're still processing
        one file, one redshift, and one transition at a time.
        """
        # This can reuse the same implementation as _process_system
        return self._process_system(item)
    
    def _process_file(self, item):
        """
        Process a batch item with a single file, redshift, and transition.
        
        This is also similar to _process_system.
        """
        # This can reuse the same implementation as _process_system
        return self._process_system(item)    
    


    def _save_item_json(self, result):
        """Save an individual item's rb_spec object to JSON."""
        if not result.get('spec_object') or not self.batch_settings['output_directory']:
            return False
        
        try:
            # Create output filename
            output_dir = self.batch_settings['output_directory']
            filename = result.get('filename', 'Unknown')
            basename = os.path.splitext(os.path.basename(filename))[0]
            transition_name = result.get('transition_name', 'unknown').replace(' ', '_')
            redshift = result.get('redshift', 0)
            
            output_file = f"{basename}_{transition_name}_z{redshift:.3f}.json"
            output_path = os.path.join(output_dir, output_file)
            
            # Save the rb_spec object
            spec_object = result['spec_object']
            spec_object.save_slice(output_path, file_format='json')
            
            return True
        except Exception as e:
            print(f"Error saving JSON: {str(e)}")
            return False
    
    def save_batch_configuration(self, filepath):
        """Save the current batch configuration to a file."""
        try:
            config = {
                'batch_items': self.batch_items,
                'batch_settings': self.batch_settings
            }
            
            with open(filepath, 'w') as f:
                json.dump(config, f, indent=2)
            
            self.status_updated.emit(f"Batch configuration saved to {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error saving batch configuration: {str(e)}")
            return False
    
    def load_batch_configuration(self, filepath):
        """Load a batch configuration from a file."""
        try:
            with open(filepath, 'r') as f:
                config = json.load(f)
            
            if 'batch_items' in config:
                self.batch_items = config['batch_items']
            
            if 'batch_settings' in config:
                self.batch_settings.update(config['batch_settings'])
            
            self.status_updated.emit(f"Batch configuration loaded from {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error loading batch configuration: {str(e)}")
            return False
    
    def export_results_csv(self, filepath):
        """Export batch results to a CSV file."""
        if not self.results:
            self.status_updated.emit("No results to export.")
            return False
        
        try:
            with open(filepath, 'w', newline='') as csvfile:
                fieldnames = ['filename', 'redshift', 'transition', 'transition_name', 
                             'vmin', 'vmax', 'W', 'W_e', 'logN', 'logN_e', 'status']
                
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for result in self.results:
                    # Create a dict with only the fields in fieldnames
                    row = {field: result.get(field, '') for field in fieldnames}
                    writer.writerow(row)
            
            self.status_updated.emit(f"Results exported to {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error exporting results: {str(e)}")
            return False
    
    def export_error_log(self, filepath):
        """Export the error log to a file."""
        if not self.error_log:
            self.status_updated.emit("No errors to export.")
            return False
        
        try:
            with open(filepath, 'w') as f:
                for i, error in enumerate(self.error_log):
                    f.write(f"Error #{i+1}:\n")
                    f.write(f"Item index: {error['item_index']}\n")
                    f.write(f"Item: {self._get_item_description(error['item'])}\n")
                    f.write(f"Error: {error['error']}\n")
                    f.write("\n")
            
            self.status_updated.emit(f"Error log exported to {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error exporting error log: {str(e)}")
            return False