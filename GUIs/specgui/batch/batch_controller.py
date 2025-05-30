# rbcodes/GUIs/specgui/batch/batch_controller.py
import os
import json
import numpy as np
from typing import Tuple
from PyQt5.QtCore import QObject, pyqtSignal

from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object
from rbcodes.GUIs.specgui.batch.master_batch_table import MasterBatchTable

class BatchController(QObject):
    """Controller for managing batch processing using master batch table."""
    
    # Signals
    status_updated = pyqtSignal(str)  # For status updates
    batch_progress = pyqtSignal(int, int)  # current, total
    batch_item_completed = pyqtSignal(int, dict)  # row_index, results
    batch_item_failed = pyqtSignal(int, str)  # row_index, error message
    table_changed = pyqtSignal()  # When table structure changes
    
    def __init__(self):
        super().__init__()
        self.master_table = MasterBatchTable()
        self.error_log = []  # List of errors during processing
        self.config_load_directory = None  # Track where config was loaded from
        
        # Connect master table signals
        self.master_table.item_added.connect(lambda: self.table_changed.emit())
        self.master_table.item_removed.connect(lambda: self.table_changed.emit())
        self.master_table.item_updated.connect(self._on_item_updated)
        self.master_table.table_cleared.connect(lambda: self.table_changed.emit())

    def load_batch_configuration(self, filepath: str) -> bool:
        """Load a batch configuration using hybrid format."""
        try:        
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            # Load master table
            self.master_table.from_dict(data)
            
            # Store the directory where config was loaded from
            self.config_load_directory = os.path.dirname(filepath)
            
            # RECREATE rb_spec objects for completed items
            self._recreate_rb_spec_objects()
            
            self.status_updated.emit(f"Batch configuration loaded from {filepath}")
            return True
            
        except Exception as e:
            error_msg = f"Error loading batch configuration: {str(e)}"
            self.status_updated.emit(error_msg)
            return False

    def _on_item_updated(self, row_index: int, what_changed: str):
        """Handle item updates from master table."""
        self.table_changed.emit()
        
        # Update status based on what changed
        item = self.master_table.get_item(row_index)
        if item:
            if what_changed in ['slice_range', 'major_change']:
                self.status_updated.emit(f"Item updated - will need reprocessing")
            elif what_changed in ['continuum_method', 'continuum_masks']:
                self.status_updated.emit(f"Continuum settings updated - will need refitting")
            elif what_changed == 'ew_range':
                self.status_updated.emit(f"EW range updated - will need recalculation")

    # Properties for backward compatibility
    @property
    def batch_items(self):
        """Get batch items in old format for backward compatibility."""
        items = []
        for item in self.master_table.get_all_items():
            # Convert to old format
            old_format = {
                'type': 'multiple_transitions',  # Default type
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
            items.append(old_format)
        return items
    
    def add_batch_item(self, template_params: dict) -> int:
        """Add a batch item from dictionary (for backward compatibility)."""
        return self.master_table.add_item(template_params)
    
    def remove_batch_item(self, row_index: int) -> bool:
        """Remove a batch item by row index."""
        return self.master_table.remove_item(row_index)
    
    def clear_all_items(self):
        """Clear all batch items."""
        self.master_table.clear_all()
        self.error_log.clear()
    
    def validate_all_items(self):
        """Validate all items in the master table."""
        valid_ids, invalid_items = self.master_table.validate_all()
        
        if invalid_items:
            error_messages = [f"Item {i+1}: {error}" for i, (item_id, error) in enumerate(invalid_items)]
            error_summary = f"Found {len(invalid_items)} validation errors:\n" + "\n".join(error_messages[:5])
            if len(invalid_items) > 5:
                error_summary += f"\n... and {len(invalid_items) - 5} more errors"
            self.status_updated.emit(error_summary)
            return False, error_summary
        
        self.status_updated.emit(f"All {len(valid_ids)} batch items are valid.")
        return True, f"All {len(valid_ids)} items valid"

    def process_batch(self, selected_row_indices=None):
        """
        Process batch items using the master table.
        
        Parameters:
        -----------
        selected_row_indices : list or None
            If provided, only process these row indices. If None, process all items.
        
        Returns:
        --------
        success : bool
            True if batch processing completed successfully.
        results : list
            List of processing results in old format for backward compatibility.
        """
        total_items = self.master_table.get_item_count()
        
        if total_items == 0:
            self.status_updated.emit("No batch items to process.")
            return False, []
        
        # Determine which items to process
        if selected_row_indices is not None:
            items_to_process = selected_row_indices
        else:
            items_to_process = list(range(total_items))
        
        if not items_to_process:
            self.status_updated.emit("No valid items to process.")
            return False, []
        
        self.status_updated.emit(f"Starting batch processing of {len(items_to_process)} items...")
        self.error_log.clear()
        
        # Process each item
        results_old_format = []  # For backward compatibility
        
        for i, row_index in enumerate(items_to_process):
            try:
                item = self.master_table.get_item(row_index)
                if not item:
                    continue
                    
                self.status_updated.emit(f"Processing item {i+1}/{len(items_to_process)}: {item.template.transition_name}")
                self.batch_progress.emit(i+1, len(items_to_process))
                
                # Process the item
                success = self._process_single_item(row_index)
                
                if success:
                    # Get updated item data
                    updated_item = self.master_table.get_item(row_index)
                    
                    # Create result in old format for backward compatibility
                    result_old_format = {
                        'type': 'multiple_transitions',
                        'filename': updated_item.template.filename,
                        'redshift': updated_item.template.redshift,
                        'transition': updated_item.template.transition,
                        'transition_name': updated_item.template.transition_name,
                        'slice_vmin': updated_item.template.slice_vmin,
                        'slice_vmax': updated_item.template.slice_vmax,
                        'ew_vmin': updated_item.template.ew_vmin,
                        'ew_vmax': updated_item.template.ew_vmax,
                        'W': updated_item.results.W,
                        'W_e': updated_item.results.W_e,
                        'logN': updated_item.results.logN,
                        'logN_e': updated_item.results.logN_e,
                        'vel_centroid': updated_item.results.vel_centroid,
                        'SNR': updated_item.results.SNR,
                        'status': 'success',
                        'spec_object': self.master_table.get_rb_spec_object(row_index)
                    }
                    results_old_format.append(result_old_format)
                    
                    # Emit signal for item completion
                    self.batch_item_completed.emit(row_index, result_old_format)
                    
                else:
                    # Get error info
                    item = self.master_table.get_item(row_index)
                    error_message = item.analysis.error_message if item else "Unknown error"
                    
                    # Create error result
                    error_result = {
                        'type': 'multiple_transitions',
                        'filename': item.template.filename if item else '',
                        'transition_name': item.template.transition_name if item else '',
                        'status': 'error',
                        'error': error_message,
                        'spec_object': None
                    }
                    results_old_format.append(error_result)
                    
                    # Emit signal for item failure
                    self.batch_item_failed.emit(row_index, error_message)
                
            except Exception as e:
                error_msg = f"Error processing item {i+1}/{len(items_to_process)}: {str(e)}"
                self.status_updated.emit(error_msg)
                self.error_log.append(error_msg)
                
                # Update item status
                self.master_table.update_analysis(row_index, 
                                                processing_status="error", 
                                                error_message=str(e))
                
                # Emit failure signal
                self.batch_item_failed.emit(row_index, str(e))
        
        # Batch processing complete
        if self.error_log:
            self.status_updated.emit(f"Batch processing completed with {len(self.error_log)} errors.")
        else:
            self.status_updated.emit("Batch processing completed successfully.")
        
        return len(self.error_log) == 0, results_old_format

    def _process_single_item(self, row_index: int) -> bool:
        """Process a single batch item using the master table."""
        try:
            # Get item data from DataFrame
            item = self.master_table.get_item(row_index)
            if not item:
                return False
            
            # Update status
            self.master_table.update_analysis(row_index, processing_status="processing")
            
            # Load the spectrum file
            if item.template.filename.lower().endswith('.json'):
                # Load saved rb_spec object
                spec = load_rb_spec_object(item.template.filename)
            else:
                # Load spectrum from file
                # Determine file type based on extension
                ext = os.path.splitext(item.template.filename)[1].lower()
                if ext == '.fits':
                    filetype = 'linetools'
                elif ext in ['.txt', '.dat']:
                    filetype = 'ascii'
                else:
                    filetype = None  # Auto-detect
                    
                spec = rb_spec.from_file(item.template.filename, filetype=filetype)
            
            # Apply redshift
            spec.shift_spec(item.template.redshift)
            
            # Slice spectrum around transition using SLICE velocity range
            spec.slice_spec(
                item.template.transition,
                item.template.slice_vmin, item.template.slice_vmax,
                use_vel=True,
                linelist=item.template.linelist
            )
            
            # Fit continuum based on settings from master table
            method = item.analysis.continuum_method
            if method == 'polynomial':
                # Build mask from continuum_masks
                mask = []
                for vmin, vmax in item.analysis.continuum_masks:
                    mask.extend([vmin, vmax])
                
                # Use polynomial fitting
                spec.fit_continuum(
                    mask=mask if mask else False,
                    Legendre=item.analysis.continuum_order,
                    use_weights=item.analysis.use_weights,
                    sigma_clip=True  # Use sigma clipping for robustness
                )
            elif method == 'flat':
                # Use flat continuum
                spec.cont = np.ones_like(spec.flux_slice)
                spec.fnorm = spec.flux_slice.copy()
                spec.enorm = spec.error_slice.copy()
            
            # Store the rb_spec object
            self.master_table.set_rb_spec_object(row_index, spec)
            
            # Compute EW using EW velocity range (different from slice range!)
            spec.compute_EW(
                item.template.transition,
                vmin=item.template.ew_vmin,
                vmax=item.template.ew_vmax,
                SNR=item.analysis.calculate_snr,
                _binsize=item.analysis.binsize
            )
            
            # Update results in master table
            self.master_table.update_results(row_index,
                W=spec.W,
                W_e=spec.W_e,
                N=spec.N,
                N_e=spec.N_e,
                logN=spec.logN,
                logN_e=spec.logN_e,
                vel_centroid=getattr(spec, 'vel_centroid', 0),
                vel_disp=getattr(spec, 'vel_disp', 0),
                SNR=getattr(spec, 'SNR', 0)
            )
            
            # Update analysis status
            self.master_table.update_analysis(row_index, processing_status="complete")
              
            return True
            
        except Exception as e:
            error_msg = str(e) if str(e) else f"Unknown error: {type(e).__name__}"
            
            # Update item with error status
            self.master_table.update_analysis(row_index, 
                                            processing_status="error", 
                                            error_message=error_msg)
            return False

    def save_batch_configuration(self, filepath: str) -> bool:
        """Save the current batch configuration using hybrid format."""
        try:
            data = self.master_table.to_dict()
            
            with open(filepath, 'w') as f:
                json.dump(data, f, indent=2)
            
            self.status_updated.emit(f"Batch configuration saved to {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error saving batch configuration: {str(e)}")
            return False

    def load_batch_configuration(self, filepath: str) -> bool:
        """Load a batch configuration using hybrid format."""
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            # Load master table
            self.master_table.from_dict(data)
            
            # Store the directory where config was loaded from
            self.config_load_directory = os.path.dirname(filepath)
            
            # RECREATE rb_spec objects for completed items
            self._recreate_rb_spec_objects()
            
            self.status_updated.emit(f"Batch configuration loaded from {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error loading batch configuration: {str(e)}")
            return False
    
    def _recreate_rb_spec_objects(self):
        """Recreate rb_spec objects for items that have been processed."""
        completed_indices = self.master_table.get_items_by_status("complete")
        
        recreated_count = 0
        for row_index in completed_indices:
            item = self.master_table.get_item(row_index)
            if item:
                # Check if rb_spec object already exists
                existing_spec = self.master_table.get_rb_spec_object(row_index)
                if not existing_spec:
                    try:
                        success = self._recreate_single_rb_spec(item, row_index)
                        if success:
                            recreated_count += 1
                            print(f"Recreated rb_spec object for {item.template.transition_name}")
                        else:
                            print(f"Failed to recreate rb_spec object for {item.template.transition_name}")
                    except Exception as e:
                        print(f"Error recreating rb_spec for {item.template.transition_name}: {e}")
        
        # EMIT SIGNAL TO REFRESH UI
        if recreated_count > 0:
            self.table_changed.emit()
    
    def _recreate_single_rb_spec(self, item, row_index) -> bool:
        """Recreate a single rb_spec object from saved parameters."""
        try:
            from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object
            
            # Load the spectrum file (same as in _process_single_item)
            if item.template.filename.lower().endswith('.json'):
                spec = load_rb_spec_object(item.template.filename)
            else:
                ext = os.path.splitext(item.template.filename)[1].lower()
                if ext == '.fits':
                    filetype = 'linetools'
                elif ext in ['.txt', '.dat']:
                    filetype = 'ascii'
                else:
                    filetype = None
                    
                spec = rb_spec.from_file(item.template.filename, filetype=filetype)
            
            # Apply redshift
            spec.shift_spec(item.template.redshift)
            
            # Slice spectrum
            spec.slice_spec(
                item.template.transition,
                item.template.slice_vmin, item.template.slice_vmax,
                use_vel=True,
                linelist=item.template.linelist
            )
            
            # Recreate continuum fit
            method = item.analysis.continuum_method
            if method == 'polynomial':
                # Parse continuum_masks from JSON string
                import json
                try:
                    continuum_masks = json.loads(item.analysis.continuum_masks) if item.analysis.continuum_masks else []
                except:
                    continuum_masks = []
                
                mask = []
                for vmin, vmax in continuum_masks:
                    mask.extend([vmin, vmax])
                
                spec.fit_continuum(
                    mask=mask if mask else False,
                    Legendre=item.analysis.continuum_order,
                    use_weights=item.analysis.use_weights,
                    sigma_clip=True
                )
            elif method == 'flat':
                spec.cont = np.ones_like(spec.flux_slice)
                spec.fnorm = spec.flux_slice.copy()
                spec.enorm = spec.error_slice.copy()
            
            # Restore the EW calculation results to the spec object
            spec.W = item.results.W
            spec.W_e = item.results.W_e
            spec.N = item.results.N
            spec.N_e = item.results.N_e
            spec.logN = item.results.logN
            spec.logN_e = item.results.logN_e
            spec.vel_centroid = item.results.vel_centroid
            spec.vel_disp = item.results.vel_disp
            spec.SNR = item.results.SNR
            spec.vmin = item.template.ew_vmin
            spec.vmax = item.template.ew_vmax
            spec.trans = item.template.transition_name
            spec.trans_wave = item.template.transition
            
            # Store the recreated object using row_index
            self.master_table.set_rb_spec_object(row_index, spec)
            
            return True
            
        except Exception as e:
            print(f"Error recreating rb_spec object: {e}")
            return False
    
    def export_template_csv(self, filepath: str) -> bool:
        """Export just the template part as CSV."""
        return self.master_table.export_template_csv(filepath)
    
    def import_template_csv(self, filepath: str) -> Tuple[bool, str]:
        """Import template from CSV."""
        return self.master_table.import_template_csv(filepath)

    def export_results_csv(self, filepath: str) -> bool:
        """Export batch results to a CSV file with enhanced fields including velocity centroid and dispersion."""
        try:
            import csv
            
            items = self.master_table.get_all_items()
            if not items:
                self.status_updated.emit("No results to export.")
                return False
            
            with open(filepath, 'w', newline='') as csvfile:
                fieldnames = ['filename', 'redshift', 'transition', 'transition_name', 
                             'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax',
                             'W', 'W_e', 'logN', 'logN_e', 'SNR', 
                             'vel_centroid', 'vel_disp', 'status']
                
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for item in items:
                    row = {
                        'filename': os.path.basename(item.template.filename),
                        'redshift': item.template.redshift,
                        'transition': item.template.transition,
                        'transition_name': item.template.transition_name,
                        'slice_vmin': item.template.slice_vmin,
                        'slice_vmax': item.template.slice_vmax,
                        'ew_vmin': item.template.ew_vmin,
                        'ew_vmax': item.template.ew_vmax,
                        'W': item.results.W,
                        'W_e': item.results.W_e,
                        'logN': item.results.logN,
                        'logN_e': item.results.logN_e,
                        'SNR': item.results.SNR,
                        'vel_centroid': item.results.vel_centroid,
                        'vel_disp': item.results.vel_disp,
                        'status': item.analysis.processing_status
                    }
                    writer.writerow(row)
            
            self.status_updated.emit(f"Enhanced CSV results exported to {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error exporting enhanced CSV results: {str(e)}")
            return False
    
    def export_error_log(self, filepath: str) -> bool:
        """Export the error log to a file."""
        if not self.error_log:
            self.status_updated.emit("No errors to export.")
            return False
        
        try:
            with open(filepath, 'w') as f:
                for i, error in enumerate(self.error_log):
                    f.write(f"Error #{i+1}: {error}\n")
            
            self.status_updated.emit(f"Error log exported to {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error exporting error log: {str(e)}")
            return False
    
    # Convenience methods for UI
    def get_item_count(self) -> int:
        """Get total number of items."""
        return self.master_table.get_item_count()
    
    def get_valid_item_count(self) -> int:
        """Get number of valid items."""
        valid_indices, _ = self.master_table.validate_all()
        return len(valid_indices)
    
    def get_completed_item_count(self) -> int:
        """Get number of completed items."""
        return len(self.master_table.get_items_by_status("complete"))
    
    def get_error_item_count(self) -> int:
        """Get number of items with errors."""
        return len(self.master_table.get_items_by_status("error"))