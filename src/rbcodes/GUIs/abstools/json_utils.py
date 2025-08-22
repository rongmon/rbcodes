"""
JSON saving and loading functionality for the absorption line analysis toolbox.
This module provides functions to save and load analysis data in JSON format.
"""

import json
import numpy as np
import os
from datetime import datetime

class NumPyJSONEncoder(json.JSONEncoder):
    """
    Custom JSON encoder that handles NumPy arrays and other special types.
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return {
                "_type": "ndarray",
                "data": obj.tolist(),
                "dtype": str(obj.dtype)
            }
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif hasattr(obj, '__dict__'):
            # For objects with a __dict__ attribute
            result = obj.__dict__.copy()
            result["_type"] = obj.__class__.__name__
            return result
        else:
            # Let the base class default method raise the TypeError
            return super().default(obj)

def json_object_hook(obj):
    """
    Custom object hook for JSON deserialization that handles special types.
    """
    if "_type" in obj:
        if obj["_type"] == "ndarray":
            return np.array(obj["data"], dtype=obj.get("dtype", None))
        # Additional types can be handled here if needed
    return obj

def save_to_json(data, filename):
    """
    Save the analysis data to a JSON file.
    
    Parameters:
    -----------
    data : dict
        The ions dictionary containing all analysis data
    filename : str
        Path to the output JSON file
        
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Remove matplotlib text objects that can't be serialized
        clean_data = {}
        for key, value in data.items():
            if key == 'Target':
                clean_data[key] = value
                continue
                
            clean_value = {}
            for k, v in value.items():
                if k == 'EW_text':
                    # Skip matplotlib text objects
                    continue
                elif k == 'pco':
                    # Convert polynomial object to coefficients
                    if v is not None:
                        clean_value[k] = {
                            "_type": "legendre",
                            "coef": v.coef.tolist(),
                            "domain": v.domain.tolist(),
                            "window": v.window.tolist()
                        }
                    else:
                        clean_value[k] = None
                else:
                    clean_value[k] = v
            clean_data[key] = clean_value
            
        # Add metadata
        metadata = {
            "version": "1.0",
            "saved_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "format": "AbsTools JSON"
        }
        
        output_data = {
            "metadata": metadata,
            "ions": clean_data
        }
        
        # Write to file
        with open(filename, 'w') as f:
            json.dump(output_data, f, cls=NumPyJSONEncoder, indent=2)
            
        return True
        
    except Exception as e:
        print(f"Error saving to JSON: {e}")
        return False

def load_from_json(filename):
    """
    Load analysis data from a JSON file.
    
    Parameters:
    -----------
    filename : str
        Path to the input JSON file
        
    Returns:
    --------
    dict or None
        The loaded ions dictionary, or None if loading failed
    """
    try:
        # Read from file
        with open(filename, 'r') as f:
            data = json.load(f, object_hook=json_object_hook)
            
        # Extract ions dictionary
        ions = data.get("ions", {})
        
        # Reconstruct polynomial objects
        import numpy.polynomial.legendre as L
        
        for key, value in ions.items():
            if key == 'Target':
                continue
                
            if 'pco' in value and value['pco'] is not None and isinstance(value['pco'], dict):
                pco_data = value['pco']
                if pco_data.get("_type") == "legendre":
                    # Create a new Legendre polynomial object
                    pco = L.Legendre(pco_data["coef"])
                    if "domain" in pco_data:
                        pco.domain = np.array(pco_data["domain"])
                    if "window" in pco_data:
                        pco.window = np.array(pco_data["window"])
                    value['pco'] = pco
            
            # Initialize EW_text to None
            value['EW_text'] = None
                
        return ions
        
    except Exception as e:
        print(f"Error loading from JSON: {e}")
        return None

# Example implementation of save and load functionality for Metal_Plot.py
def add_json_support_to_save_page():
    """
    Code to add to SavePage class in ui_components.py to support JSON saving
    """
    def onjson(self, parentvals):
        """Save analysis data as a JSON file"""
        try:
            from json_utils import save_to_json
            
            json_file = self.jsonline.text()
            success = save_to_json(parentvals.ions, json_file)
            
            if success:
                self.jsonsave.setStyleSheet('background-color : green')
                parentvals.save = True
            else:
                self.jsonsave.setStyleSheet('background-color : red')
        except Exception as e:
            print(f"Error saving JSON file: {e}")
            self.jsonsave.setStyleSheet('background-color : red')
            
    # Add the function to the SavePage class
    # SavePage.onjson = onjson
    
    # Also need to add UI elements to the setup_ui method of SavePage:
    """
    # Create JSON save UI elements
    jsonlabel = QLabel(r"Enter path and filename: (e.g. pathname\analysis.json)", self)
    jsonlabel.setGeometry(100, 325, 400, 30)
    
    self.jsonline = QLineEdit(self)
    self.jsonline.setText(f"Spectrum_Analysis_z_{str(parentvals.z)}.json")
    self.jsonline.setGeometry(100, 350, 300, 30)
    
    self.jsonsave = QPushButton(r"Save as JSON", self)
    self.jsonsave.setGeometry(410, 350, 200, 30)
    self.jsonsave.clicked.connect(lambda: self.onjson(parentvals))
    """

# Example implementation of loading from JSON in Metal_Plot.py
def load_from_json_example():
    """
    Example code to load analysis data from a JSON file in Metal_Plot.py
    """
    from json_utils import load_from_json
    
    # Assume we have a JSON file to load
    json_file = "Spectrum_Analysis_z_0.348.json"
    
    if os.path.exists(json_file):
        # Load the ions dictionary
        ions = load_from_json(json_file)
        
        if ions is not None:
            # Initialize the GUI with the loaded data
            from GUIs.abstools import Metal_Plot as M
            M.Transitions(ions)
        else:
            print(f"Failed to load {json_file}")
    else:
        print(f"File not found: {json_file}")