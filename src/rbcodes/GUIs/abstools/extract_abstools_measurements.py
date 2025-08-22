import json
import numpy as np
import pandas as pd

def extract_measurements_table(json_file):
    """
    Extract measurements from an AbsTools JSON file into a pandas DataFrame
    
    Parameters:
    -----------
    json_file : str
        Path to the JSON file saved from AbsTools
        
    Returns:
    --------
    df : pandas.DataFrame
        DataFrame containing all measurements
    """
    # Custom JSON decoder to handle special types
    def json_object_hook(obj):
        if "_type" in obj:
            if obj["_type"] == "ndarray":
                return np.array(obj["data"], dtype=obj.get("dtype", None))
        return obj
    
    # Read the JSON file
    with open(json_file, 'r') as f:
        data = json.load(f, object_hook=json_object_hook)
    
    # Extract ions dictionary
    if "ions" in data:
        ions = data["ions"]
    else:
        ions = data
    
    # Extract redshift
    redshift = ions['Target']['z']
    
    # Create a list to hold the measurements
    measurements = []
    
    # Loop through each ion (skip Target)
    for ion_name in ions:
        if ion_name == 'Target':
            continue
            
        ion = ions[ion_name]
        
        # Get the rest wavelength and oscillator strength
        rest_wavelength = ion['lam_0']
        f_value = ion['f']
        
        # Get EW measurement
        ew = ion['EW']
        ew_error = ion['EWsig']
        
        # Get velocity limits for the measurement
        vel_min = ion['EWlims'][0]
        vel_max = ion['EWlims'][1]
        
        # Get column density (if available)
        if ion['N'] is not None and ion['N'] > 0:
            log_N = np.log10(ion['N'])
            log_N_err = np.log10(ion['N'] + ion['Nsig']) - log_N
        else:
            log_N = None
            log_N_err = None
        
        # Get the flag (detection type)
        flag = ion['flag']  # 0=detection, 1=upper limit, 2=lower limit
        
        # Convert flag to string for clarity
        flag_str = {0: "Detection", 1: "Upper Limit", 2: "Lower Limit"}.get(flag, "Unknown")
        
        # Add to measurements list
        measurements.append({
            'Ion': ion_name,
            'Rest Wavelength (Å)': rest_wavelength,
            'f-value': f_value,
            'EW (mÅ)': ew,
            'EW Error (mÅ)': ew_error,
            'Velocity Min (km/s)': vel_min,
            'Velocity Max (km/s)': vel_max,
            'log N (cm^-2)': log_N,
            'log N Error': log_N_err,
            'Flag': flag,
            'Detection Type': flag_str
        })
    
    # Create DataFrame
    df = pd.DataFrame(measurements)
    
    # Add metadata as attributes
    df.attrs['redshift'] = redshift
    df.attrs['source_file'] = json_file
    
    return df

# Example usage
if __name__ == "__main__":
    # Replace with your JSON file path
    file_path = "Spectrum_Analysis_z_0.348.json"
    
    # Extract measurements to a DataFrame
    df = extract_measurements_table(file_path)
    
    # Print the table
    print(f"Measurements for spectrum at z = {df.attrs['redshift']}")
    print(df)
    
    # Export to CSV
    csv_file = file_path.replace('.json', '_measurements.csv')
    df.to_csv(csv_file, index=False)
    print(f"Measurements saved to {csv_file}")
    
    # You can also export to other formats like Excel
    excel_file = file_path.replace('.json', '_measurements.xlsx')
    df.to_excel(excel_file, index=False)
    print(f"Measurements saved to {excel_file}")