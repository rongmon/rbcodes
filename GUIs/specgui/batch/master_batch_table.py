# rbcodes/GUIs/specgui/batch/master_batch_table.py
import os
import uuid
import json
import csv
import numpy as np
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
from PyQt5.QtCore import QObject, pyqtSignal


@dataclass
class TemplateParams:
    """Template parameters - the human-readable/editable part."""
    filename: str
    redshift: float
    transition: float
    transition_name: str
    slice_vmin: int
    slice_vmax: int
    ew_vmin: int
    ew_vmax: int
    linelist: str = "atom"
    method: str = "closest"


@dataclass
class AnalysisState:
    """Analysis state - created during processing."""
    continuum_method: str = "polynomial"
    continuum_order: int = 3
    continuum_masks: List[Tuple[int, int]] = None
    continuum_fit_params: Dict = None
    use_weights: bool = False
    calculate_snr: bool = True
    binsize: int = 3
    processing_status: str = "ready"  # ready, processing, complete, error
    last_modified: str = ""
    error_message: str = ""
    
    def __post_init__(self):
        if self.continuum_masks is None:
            self.continuum_masks = []
        if self.continuum_fit_params is None:
            self.continuum_fit_params = {}
        if not self.last_modified:
            self.last_modified = datetime.now().isoformat()


@dataclass
class Results:
    """Analysis results."""
    W: float = 0.0
    W_e: float = 0.0
    N: float = 0.0
    N_e: float = 0.0
    logN: float = 0.0
    logN_e: float = 0.0
    vel_centroid: float = 0.0
    vel_disp: float = 0.0
    SNR: float = 0.0
    calculation_timestamp: str = ""
    
    def __post_init__(self):
        if not self.calculation_timestamp:
            self.calculation_timestamp = datetime.now().isoformat()


class BatchItem:
    """Single batch item containing template, analysis state, and results."""
    
    def __init__(self, template_params: TemplateParams, item_id: str = None):
        self.id = item_id or str(uuid.uuid4())
        self.template = template_params
        self.analysis = AnalysisState()
        self.results = Results()
        self.rb_spec_object = None  # Will hold the actual rb_spec instance
        
    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        return {
            'id': self.id,
            'template': asdict(self.template),
            'analysis': asdict(self.analysis),
            'results': asdict(self.results),
            # Note: rb_spec_object is not serialized - it's recreated as needed
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'BatchItem':
        """Create from dictionary."""
        template = TemplateParams(**data['template'])
        item = cls(template, data['id'])
        
        # Restore analysis state
        analysis_data = data.get('analysis', {})
        item.analysis = AnalysisState(**analysis_data)
        
        # Restore results
        results_data = data.get('results', {})
        item.results = Results(**results_data)
        
        return item
    
    def validate(self) -> Tuple[bool, str]:
        """Validate this batch item."""
        # Check file exists
        if not os.path.exists(self.template.filename):
            return False, f"File does not exist: {self.template.filename}"
        
        # Check redshift range
        if not (0 <= self.template.redshift <= 10):
            return False, f"Redshift {self.template.redshift} is outside reasonable range (0-10)"
        
        # Check transition wavelength
        if not (1 <= self.template.transition <= 10000):
            return False, f"Transition wavelength {self.template.transition} Å is outside reasonable range (1-10000 Å)"
        
        # Check velocity ranges
        if self.template.slice_vmin >= self.template.slice_vmax:
            return False, f"Slice velocity range invalid: vmin ({self.template.slice_vmin}) >= vmax ({self.template.slice_vmax})"
        
        if self.template.ew_vmin >= self.template.ew_vmax:
            return False, f"EW velocity range invalid: vmin ({self.template.ew_vmin}) >= vmax ({self.template.ew_vmax})"
        
        # EW range should be within slice range (with tolerance)
        tolerance = 50  # km/s
        if (self.template.ew_vmin < (self.template.slice_vmin - tolerance) or 
            self.template.ew_vmax > (self.template.slice_vmax + tolerance)):
            return False, (f"EW range ({self.template.ew_vmin}, {self.template.ew_vmax}) should be within "
                          f"slice range ({self.template.slice_vmin}, {self.template.slice_vmax}) ± {tolerance} km/s")
        
        return True, ""
    
    def needs_recalculation(self, what_changed: str) -> List[str]:
        """Determine what needs to be recalculated based on what changed."""
        recalc_needed = []
        
        if what_changed in ['slice_range', 'redshift', 'transition', 'filename']:
            # Major changes - need to redo everything
            recalc_needed = ['slice', 'continuum', 'ew']
        elif what_changed in ['continuum_method', 'continuum_order', 'continuum_masks']:
            # Continuum changes - need to refit continuum and recalculate EW
            recalc_needed = ['continuum', 'ew']
        elif what_changed in ['ew_range', 'snr_settings']:
            # Just need to recalculate EW
            recalc_needed = ['ew']
        
        return recalc_needed
    
    def mark_for_recalculation(self, what_changed: str):
        """Mark item as needing recalculation."""
        if what_changed in ['slice_range', 'redshift', 'transition', 'filename', 'continuum_method', 'continuum_order', 'continuum_masks']:
            self.analysis.processing_status = "needs_processing"
        elif what_changed in ['ew_range', 'snr_settings']:
            self.analysis.processing_status = "needs_ew_recalc"
        
        self.analysis.last_modified = datetime.now().isoformat()


class MasterBatchTable(QObject):
    """Master table that serves as single source of truth for all batch processing."""
    
    # Signals
    item_added = pyqtSignal(str)  # item_id
    item_removed = pyqtSignal(str)  # item_id
    item_updated = pyqtSignal(str, str)  # item_id, what_changed
    table_cleared = pyqtSignal()
    validation_changed = pyqtSignal()
    
    def __init__(self):
        super().__init__()
        self.items: Dict[str, BatchItem] = {}
        self._item_order: List[str] = []  # Maintain order for display
        self.metadata = {
            'version': '1.0',
            'created': datetime.now().isoformat(),
            'description': '',
            'software_version': '1.0'
        }
    
    def add_item(self, template_params: TemplateParams) -> str:
        """Add a new batch item."""
        item = BatchItem(template_params)
        self.items[item.id] = item
        self._item_order.append(item.id)
        
        self.item_added.emit(item.id)
        return item.id
    
    def remove_item(self, item_id: str) -> bool:
        """Remove a batch item."""
        if item_id in self.items:
            del self.items[item_id]
            if item_id in self._item_order:
                self._item_order.remove(item_id)
            self.item_removed.emit(item_id)
            return True
        return False
    
    def get_item(self, item_id: str) -> Optional[BatchItem]:
        """Get a batch item by ID."""
        return self.items.get(item_id)

    def get_item_count(self) -> int:
        """Get the total number of items in the table."""
        return len(self.items)
    
    def get_all_items(self) -> List[BatchItem]:
        """Get all items in order."""
        return [self.items[item_id] for item_id in self._item_order if item_id in self.items]
    

    def update_template(self, item_id: str, **kwargs) -> bool:
        """Update template parameters for an item."""
        print(f"DEBUG: update_template called for item_id: {item_id}")
        print(f"DEBUG: kwargs: {kwargs}")
        
        if item_id not in self.items:
            print(f"DEBUG: item_id {item_id} not found in self.items")
            print(f"DEBUG: Available item_ids: {list(self.items.keys())}")
            return False
        
        item = self.items[item_id]
        print(f"DEBUG: Found item: {item.template.transition_name}")
        
        changed_fields = []
        
        # Update template fields
        for field, value in kwargs.items():
            print(f"DEBUG: Trying to update field '{field}' to value '{value}'")
            
            if hasattr(item.template, field):
                old_value = getattr(item.template, field)
                print(f"DEBUG: Field '{field}' exists, old value: {old_value}, new value: {value}")
                
                if old_value != value:
                    setattr(item.template, field, value)
                    changed_fields.append(field)
                    print(f"DEBUG: Updated field '{field}' from {old_value} to {value}")
                else:
                    print(f"DEBUG: Field '{field}' unchanged (same value)")
            else:
                print(f"DEBUG: Field '{field}' does not exist on template")
                print(f"DEBUG: Available template fields: {dir(item.template)}")
        
        print(f"DEBUG: Changed fields: {changed_fields}")
        
        # ... rest of the method (determine what needs recalculation)
        
        print(f"DEBUG: update_template returning True")
        return True

    
    def update_analysis(self, item_id: str, **kwargs) -> bool:
        """Update analysis parameters for an item."""
        if item_id not in self.items:
            return False
        
        item = self.items[item_id]
        changed_fields = []
        
        # Update analysis fields
        for field, value in kwargs.items():
            if hasattr(item.analysis, field):
                old_value = getattr(item.analysis, field)
                if old_value != value:
                    setattr(item.analysis, field, value)
                    changed_fields.append(field)
        
        # Mark last modified
        item.analysis.last_modified = datetime.now().isoformat()
        
        if changed_fields:
            if any(field in ['continuum_method', 'continuum_order'] for field in changed_fields):
                item.mark_for_recalculation('continuum_method')
                self.item_updated.emit(item_id, 'continuum_method')
            elif 'continuum_masks' in changed_fields:
                item.mark_for_recalculation('continuum_masks')
                self.item_updated.emit(item_id, 'continuum_masks')
            else:
                self.item_updated.emit(item_id, 'analysis_update')
        
        return True
    
    def update_results(self, item_id: str, **kwargs) -> bool:
        """Update results for an item."""
        if item_id not in self.items:
            return False
        
        item = self.items[item_id]
        
        # Update results fields
        for field, value in kwargs.items():
            if hasattr(item.results, field):
                setattr(item.results, field, value)
        
        # Mark calculation timestamp
        item.results.calculation_timestamp = datetime.now().isoformat()
        
        self.item_updated.emit(item_id, 'results_update')
        return True
    
    def set_rb_spec_object(self, item_id: str, spec_object) -> bool:
        """Store the rb_spec object for an item."""
        if item_id not in self.items:
            return False
        
        self.items[item_id].rb_spec_object = spec_object
        return True
    
    def get_rb_spec_object(self, item_id: str):
        """Get the rb_spec object for an item."""
        if item_id in self.items:
            return self.items[item_id].rb_spec_object
        return None
    
    def clear_all(self):
        """Clear all items."""
        self.items.clear()
        self._item_order.clear()
        self.table_cleared.emit()
    
    def validate_all(self) -> Tuple[List[str], List[Tuple[str, str]]]:
        """Validate all items. Returns (valid_ids, invalid_items_with_errors)."""
        valid_ids = []
        invalid_items = []
        
        for item_id in self._item_order:
            if item_id in self.items:
                is_valid, error_msg = self.items[item_id].validate()
                if is_valid:
                    valid_ids.append(item_id)
                else:
                    invalid_items.append((item_id, error_msg))
        
        return valid_ids, invalid_items
    
    def get_items_by_status(self, status: str) -> List[str]:
        """Get item IDs by processing status."""
        return [item_id for item_id, item in self.items.items() 
                if item.analysis.processing_status == status]
    
    def to_dict(self) -> Dict:
        """Convert entire table to dictionary for serialization."""
        return {
            'metadata': self.metadata,
            'template': {
                'version': self.metadata['version'],
                'description': self.metadata.get('description', ''),
                'items': [item.template.__dict__ for item in self.get_all_items()]
            },
            'analysis': {
                'created': self.metadata.get('created', ''),
                'software_version': self.metadata.get('software_version', '1.0'),
                'items': [
                    {
                        'id': item.id,
                        'analysis': asdict(item.analysis),
                        'results': asdict(item.results)
                    }
                    for item in self.get_all_items()
                ]
            }
        }
    
    def from_dict(self, data: Dict):
        """Load table from dictionary."""
        self.clear_all()
        
        # Load metadata
        self.metadata = data.get('metadata', {})
        if 'version' not in self.metadata:
            self.metadata['version'] = '1.0'
        
        # Load template items
        template_data = data.get('template', {})
        template_items = template_data.get('items', [])
        
        # Load analysis data if available
        analysis_data = data.get('analysis', {})
        analysis_items = analysis_data.get('items', [])
        analysis_by_index = {i: item for i, item in enumerate(analysis_items)}
        
        # Create items
        for i, template_item in enumerate(template_items):
            # Create template params
            template_params = TemplateParams(**template_item)
            
            # Add item
            item_id = self.add_item(template_params)
            item = self.items[item_id]
            
            # Restore analysis state if available
            if i in analysis_by_index:
                analysis_item = analysis_by_index[i]
                if 'analysis' in analysis_item:
                    item.analysis = AnalysisState(**analysis_item['analysis'])
                if 'results' in analysis_item:
                    item.results = Results(**analysis_item['results'])
                # Use the ID from the file if available
                if 'id' in analysis_item:
                    old_id = item.id
                    new_id = analysis_item['id']
                    self.items[new_id] = self.items.pop(old_id)
                    self._item_order[self._item_order.index(old_id)] = new_id
    
    def export_template_csv(self, filepath: str) -> bool:
        """Export just the template part as CSV."""
        try:
            with open(filepath, 'w', newline='') as csvfile:
                if not self.items:
                    # Create empty template with headers
                    fieldnames = ['filename', 'redshift', 'transition', 'transition_name', 
                                'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax', 'linelist']
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    return True
                
                # Get fieldnames from first item
                first_item = next(iter(self.items.values()))
                fieldnames = list(first_item.template.__dict__.keys())
                
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for item in self.get_all_items():
                    writer.writerow(asdict(item.template))
            
            return True
        except Exception as e:
            print(f"Error exporting CSV: {e}")
            return False
    
    def import_template_csv(self, filepath: str) -> Tuple[bool, str]:
        """Import template from CSV. Returns (success, error_message)."""
        try:
            items_added = 0
            with open(filepath, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                
                if not reader.fieldnames:
                    return False, "CSV file has no headers"
                
                # Check required fields
                required_fields = ['filename', 'redshift', 'transition', 'slice_vmin', 'slice_vmax', 'ew_vmin', 'ew_vmax']
                missing_fields = [field for field in required_fields if field not in reader.fieldnames]
                if missing_fields:
                    return False, f"Missing required fields: {', '.join(missing_fields)}"
                
                for row_num, row in enumerate(reader, 1):
                    try:
                        # Convert types
                        template_params = TemplateParams(
                            filename=row['filename'],
                            redshift=float(row['redshift']),
                            transition=float(row['transition']),
                            transition_name=row.get('transition_name', f"λ {float(row['transition']):.2f}"),
                            slice_vmin=int(row['slice_vmin']),
                            slice_vmax=int(row['slice_vmax']),
                            ew_vmin=int(row['ew_vmin']),
                            ew_vmax=int(row['ew_vmax']),
                            linelist=row.get('linelist', 'atom'),
                            method=row.get('method', 'closest')
                        )
                        
                        self.add_item(template_params)
                        items_added += 1
                        
                    except (ValueError, KeyError) as e:
                        print(f"Error parsing row {row_num}: {e}")
            
            return True, f"Imported {items_added} items successfully"
            
        except Exception as e:
            return False, f"Error reading CSV file: {str(e)}"