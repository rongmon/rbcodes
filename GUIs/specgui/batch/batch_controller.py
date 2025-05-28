# rbcodes/GUIs/specgui/batch/batch_controller.py
import os
import json
import numpy as np
from typing import Tuple
from PyQt5.QtCore import QObject, pyqtSignal

from rbcodes.GUIs.rb_spec import rb_spec, load_rb_spec_object
from .master_batch_table import MasterBatchTable, TemplateParams


class BatchController(QObject):
    """Controller for managing batch processing using master batch table."""
    
    # Signals
    status_updated = pyqtSignal(str)  # For status updates
    batch_progress = pyqtSignal(int, int)  # current, total
    batch_item_completed = pyqtSignal(str, dict)  # item_id, results
    batch_item_failed = pyqtSignal(str, str)  # item_id, error message
    table_changed = pyqtSignal()  # When table structure changes
    
    def __init__(self):
        super().__init__()
        self.master_table = MasterBatchTable()
        self.error_log = []  # List of errors during processing
        
        # Default batch settings
        self.batch_settings = {
            'continuum_method': 'polynomial',  # 'polynomial', 'flat'
            'polynomial_order': 3,
            'use_weights': False,
            'calculate_snr': True,
            'binsize': 3,
            'save_individual_json': True,
            'output_directory': '',
        }
        
        # Connect master table signals
        self.master_table.item_added.connect(lambda: self.table_changed.emit())
        self.master_table.item_removed.connect(lambda: self.table_changed.emit())
        self.master_table.item_updated.connect(self._on_item_updated)
        self.master_table.table_cleared.connect(lambda: self.table_changed.emit())
    
    def _on_item_updated(self, item_id: str, what_changed: str):
        """Handle item updates from master table."""
        self.table_changed.emit()
        
        # Update status based on what changed
        item = self.master_table.get_item(item_id)
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
    
    def add_batch_item(self, template_params: dict) -> str:
        """Add a batch item from dictionary (for backward compatibility)."""
        params = TemplateParams(
            filename=template_params['filename'],
            redshift=template_params['redshift'],
            transition=template_params['transition'],
            transition_name=template_params.get('transition_name', 'Unknown'),
            slice_vmin=template_params['slice_vmin'],
            slice_vmax=template_params['slice_vmax'],
            ew_vmin=template_params['ew_vmin'],
            ew_vmax=template_params['ew_vmax'],
            linelist=template_params.get('linelist', 'atom'),
            method=template_params.get('method', 'closest')
        )
        return self.master_table.add_item(params)
    
    def remove_batch_item(self, item_id: str) -> bool:
        """Remove a batch item."""
        return self.master_table.remove_item(item_id)
    
    def clear_all_items(self):
        """Clear all batch items."""
        self.master_table.clear_all()
        self.error_log.clear()
    
    def update_batch_settings(self, settings):
        """Update the batch processing settings."""
        self.batch_settings.update(settings)
        self.status_updated.emit("Batch settings updated.")
    
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
    
    def process_batch(self, selected_item_ids=None):
        """
        Process batch items using the master table.
        
        Parameters:
        -----------
        selected_item_ids : list or None
            If provided, only process these item IDs. If None, process all items.
        
        Returns:
        --------
        success : bool
            True if batch processing completed successfully.
        results : list
            List of processing results in old format for backward compatibility.
        """
        all_items = self.master_table.get_all_items()
        
        if not all_items:
            self.status_updated.emit("No batch items to process.")
            return False, []
        
        # Determine which items to process
        if selected_item_ids is not None:
            items_to_process = [item for item in all_items if item.id in selected_item_ids]
        else:
            items_to_process = all_items
        
        total_items = len(items_to_process)
        if total_items == 0:
            self.status_updated.emit("No valid items to process.")
            return False, []
        
        self.status_updated.emit(f"Starting batch processing of {total_items} items...")
        self.error_log.clear()
        
        # Create output directory if needed
        if self.batch_settings['save_individual_json'] and self.batch_settings['output_directory']:
            output_dir = self.batch_settings['output_directory']
            if not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir)
                    self.status_updated.emit(f"Created output directory: {output_dir}")
                except Exception as e:
                    error_msg = f"Failed to create output directory: {str(e)}"
                    self.status_updated.emit(error_msg)
                    self.error_log.append(error_msg)
                    return False, []
        
        # Process each item
        results_old_format = []  # For backward compatibility
        
        for i, item in enumerate(items_to_process):
            try:
                self.status_updated.emit(f"Processing item {i+1}/{total_items}: {item.template.transition_name}")
                self.batch_progress.emit(i+1, total_items)
                
                # Process the item
                success = self._process_single_item(item)
                
                if success:
                    # Create result in old format for backward compatibility
                    result_old_format = {
                        'type': 'multiple_transitions',
                        'filename': item.template.filename,
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
                        'vel_centroid': item.results.vel_centroid,
                        'SNR': item.results.SNR,
                        'status': 'success',
                        'spec_object': item.rb_spec_object
                    }
                    results_old_format.append(result_old_format)
                    
                    # Emit signal for item completion
                    self.batch_item_completed.emit(item.id, result_old_format)
                    
                else:
                    # Create error result
                    error_result = {
                        'type': 'multiple_transitions',
                        'filename': item.template.filename,
                        'redshift': item.template.redshift,
                        'transition': item.template.transition,
                        'transition_name': item.template.transition_name,
                        'slice_vmin': item.template.slice_vmin,
                        'slice_vmax': item.template.slice_vmax,
                        'ew_vmin': item.template.ew_vmin,
                        'ew_vmax': item.template.ew_vmax,
                        'W': 0.0,
                        'W_e': 0.0,
                        'logN': 0.0,
                        'logN_e': 0.0,
                        'status': 'error',
                        'error': item.analysis.error_message,
                        'spec_object': None
                    }
                    results_old_format.append(error_result)
                    
                    # Emit signal for item failure
                    self.batch_item_failed.emit(item.id, item.analysis.error_message)
                
            except Exception as e:
                error_msg = f"Error processing item {i+1}/{total_items}: {str(e)}"
                self.status_updated.emit(error_msg)
                self.error_log.append(error_msg)
                
                # Update item status
                self.master_table.update_analysis(item.id, 
                                                processing_status="error", 
                                                error_message=str(e))
                
                # Create error result for backward compatibility
                error_result = {
                    'type': 'multiple_transitions',
                    'filename': item.template.filename,
                    'redshift': item.template.redshift,
                    'transition': item.template.transition,
                    'transition_name': item.template.transition_name,
                    'slice_vmin': item.template.slice_vmin,
                    'slice_vmax': item.template.slice_vmax,
                    'ew_vmin': item.template.ew_vmin,
                    'ew_vmax': item.template.ew_vmax,
                    'W': 0.0,
                    'W_e': 0.0,
                    'logN': 0.0,
                    'logN_e': 0.0,
                    'status': 'error',
                    'error': str(e),
                    'spec_object': None
                }
                results_old_format.append(error_result)
                
                # Emit failure signal
                self.batch_item_failed.emit(item.id, str(e))
        
        # Batch processing complete
        if self.error_log:
            self.status_updated.emit(f"Batch processing completed with {len(self.error_log)} errors.")
        else:
            self.status_updated.emit("Batch processing completed successfully.")
        
        return len(self.error_log) == 0, results_old_format
    
    def _process_single_item(self, item) -> bool:
        """Process a single batch item using the master table."""
        try:
            # Update status
            self.master_table.update_analysis(item.id, processing_status="processing")
            
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
            
            # Fit continuum based on settings
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
            self.master_table.set_rb_spec_object(item.id, spec)
            
            # Compute EW using EW velocity range (different from slice range!)
            spec.compute_EW(
                item.template.transition,
                vmin=item.template.ew_vmin,
                vmax=item.template.ew_vmax,
                SNR=item.analysis.calculate_snr,
                _binsize=item.analysis.binsize
            )
            
            # Update results in master table
            self.master_table.update_results(item.id,
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
            self.master_table.update_analysis(item.id, processing_status="complete")
            
            # Save results as JSON if enabled
            if self.batch_settings['save_individual_json'] and self.batch_settings['output_directory']:
                output_dir = self.batch_settings['output_directory']
                basename = os.path.splitext(os.path.basename(item.template.filename))[0]
                
                # Create a more descriptive transition identifier
                transition_id = f"{item.template.transition_name}_{item.template.transition:.0f}"
                output_file = f"{basename}_{transition_id}_z{item.template.redshift:.3f}.json"
                output_path = os.path.join(output_dir, output_file)
                
                # Save the rb_spec object
                spec.save_slice(output_path, file_format='json')            
            return True
            
        except Exception as e:
            # Update item with error status
            self.master_table.update_analysis(item.id, 
                                            processing_status="error", 
                                            error_message=str(e))
            return False
    
    def recalculate_item(self, item_id: str, what_changed: str) -> bool:
        """Recalculate a specific item based on what changed."""


        print(f"DEBUG: recalculate_item called for {item_id}, changed: {what_changed}")
        
        item = self.master_table.get_item(item_id)
        if not item:
            print(f"DEBUG: Item not found!")
            return False
        
        if not item.rb_spec_object:
            print(f"DEBUG: rb_spec_object is None!")
            return False
        
        print(f"DEBUG: rb_spec_object exists, proceeding with recalculation")
            
        item = self.master_table.get_item(item_id)
        if not item or not item.rb_spec_object:
            return False
        
        try:
            spec = item.rb_spec_object
            
            recalc_needed = item.needs_recalculation(what_changed)
            
            if 'slice' in recalc_needed:
                # Need to re-slice spectrum
                spec.shift_spec(item.template.redshift)
                spec.slice_spec(
                    item.template.transition,
                    item.template.slice_vmin, item.template.slice_vmax,
                    use_vel=True,
                    linelist=item.template.linelist
                )
            
            if 'continuum' in recalc_needed:
                # Need to refit continuum
                method = item.analysis.continuum_method
                if method == 'polynomial':
                    # Build mask from continuum_masks
                    mask = []
                    for vmin, vmax in item.analysis.continuum_masks:
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
            
            if 'ew' in recalc_needed:
                # Need to recalculate EW
                spec.compute_EW(
                    item.template.transition,
                    vmin=item.template.ew_vmin,
                    vmax=item.template.ew_vmax,
                    SNR=item.analysis.calculate_snr,
                    _binsize=item.analysis.binsize
                )
                
                # Update results in master table
                self.master_table.update_results(item.id,
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
            
            # Update processing status
            self.master_table.update_analysis(item.id, processing_status="complete")
            
            return True
            
        except Exception as e:
            self.master_table.update_analysis(item.id, 
                                            processing_status="error", 
                                            error_message=str(e))
            return False
    
    def update_continuum_masks(self, item_id: str, masks: list) -> bool:
        """Update continuum masks for an item and trigger recalculation."""
        # Convert masks to list of tuples
        mask_tuples = []
        for i in range(0, len(masks), 2):
            if i + 1 < len(masks):
                mask_tuples.append((masks[i], masks[i+1]))
        
        # Update in master table
        success = self.master_table.update_analysis(item_id, continuum_masks=mask_tuples)
        
        if success:
            # Trigger recalculation
            return self.recalculate_item(item_id, 'continuum_masks')
        
        return False
    

    def update_ew_range(self, item_id: str, ew_vmin: int, ew_vmax: int) -> bool:
        print(f"DEBUG: update_ew_range called for {item_id}")
        
        success = self.master_table.update_template(item_id, ew_vmin=ew_vmin, ew_vmax=ew_vmax)
        print(f"DEBUG: master_table update_template returned: {success}")
        
        if success:
            result = self.recalculate_item(item_id, 'ew_range')
            print(f"DEBUG: recalculate_item returned: {result}")
            return result
        
        return False

    
    def save_batch_configuration(self, filepath: str) -> bool:
        """Save the current batch configuration using hybrid format."""
        try:
            data = self.master_table.to_dict()
            data['batch_settings'] = self.batch_settings
            
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
            
            # Load batch settings if available
            if 'batch_settings' in data:
                self.batch_settings.update(data['batch_settings'])
            
            # RECREATE rb_spec objects for completed items
            self._recreate_rb_spec_objects()
            
            self.status_updated.emit(f"Batch configuration loaded from {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error loading batch configuration: {str(e)}")
            return False
    
    def _recreate_rb_spec_objects(self):
        """Recreate rb_spec objects for items that have been processed."""
        completed_items = self.master_table.get_items_by_status("complete")
        
        recreated_count = 0
        for item_id in completed_items:
            item = self.master_table.get_item(item_id)
            if item and not item.rb_spec_object:
                try:
                    success = self._recreate_single_rb_spec(item)
                    if success:
                        recreated_count += 1
                        print(f"Recreated rb_spec object for {item.template.transition_name}")
                    else:
                        print(f"Failed to recreate rb_spec object for {item.template.transition_name}")
                except Exception as e:
                    print(f"Error recreating rb_spec for {item.template.transition_name}: {e}")
        
        # EMIT SIGNAL TO REFRESH UI - ADD THIS
        if recreated_count > 0:
            self.table_changed.emit()
            print(f"DEBUG: Recreated {recreated_count} rb_spec objects, emitted table_changed signal")
    
    def _recreate_single_rb_spec(self, item) -> bool:
        """Recreate a single rb_spec object from saved parameters."""
        try:
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
                mask = []
                for vmin, vmax in item.analysis.continuum_masks:
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
            
            # Store the recreated object
            self.master_table.set_rb_spec_object(item.id, spec)
            
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
        """Export batch results to a CSV file."""
        try:
            import csv
            
            items = self.master_table.get_all_items()
            if not items:
                self.status_updated.emit("No results to export.")
                return False
            
            with open(filepath, 'w', newline='') as csvfile:
                fieldnames = ['filename', 'redshift', 'transition', 'transition_name', 
                             'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax',
                             'W', 'W_e', 'logN', 'logN_e', 'SNR', 'status']
                
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
                        'W': item.results.W,
                        'W_e': item.results.W_e,
                        'logN': item.results.logN,
                        'logN_e': item.results.logN_e,
                        'SNR': item.results.SNR,
                        'status': item.analysis.processing_status
                    }
                    writer.writerow(row)
            
            self.status_updated.emit(f"Results exported to {filepath}")
            return True
        except Exception as e:
            self.status_updated.emit(f"Error exporting results: {str(e)}")
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
        return len(self.master_table.items)
    
    def get_valid_item_count(self) -> int:
        """Get number of valid items."""
        valid_ids, _ = self.master_table.validate_all()
        return len(valid_ids)
    
    def get_completed_item_count(self) -> int:
        """Get number of completed items."""
        return len(self.master_table.get_items_by_status("complete"))
    
    def get_error_item_count(self) -> int:
        """Get number of items with errors."""
        return len(self.master_table.get_items_by_status("error"))