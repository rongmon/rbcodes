import os
import json
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def load_abstools_json(json_file):
    """Load an AbsTools JSON file and return the ions dictionary"""
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
        return data["ions"]
    else:
        return data

def batch_process_abstools_files(json_pattern):
    """
    Process multiple AbsTools JSON files matching the given pattern
    and compile their measurements
    
    Parameters:
    -----------
    json_pattern : str
        Glob pattern for JSON files, e.g., "data/*.json"
        
    Returns:
    --------
    all_measurements : pandas.DataFrame
        DataFrame containing measurements from all files
    """
    # Find all files matching the pattern
    file_list = glob.glob(json_pattern)
    
    if not file_list:
        print(f"No files found matching pattern: {json_pattern}")
        return None
    
    print(f"Found {len(file_list)} files to process")
    
    # List to hold all measurements
    all_measurements = []
    
    # Process each file
    for file_path in file_list:
        try:
            # Get the spectrum ID from the filename
            spectrum_id = os.path.basename(file_path).replace('.json', '')
            
            # Load the data
            ions = load_abstools_json(file_path)
            
            # Extract redshift
            redshift = ions['Target']['z']
            
            # Process each ion (skip Target)
            for ion_name in ions:
                if ion_name == 'Target':
                    continue
                    
                ion = ions[ion_name]
                
                # Get measurements
                ew = ion['EW']
                ew_error = ion['EWsig']
                
                # Check if column density is available
                if ion['N'] is not None and ion['N'] > 0:
                    log_N = np.log10(ion['N'])
                    log_N_err = np.log10(ion['N'] + ion['Nsig']) - log_N
                else:
                    log_N = None
                    log_N_err = None
                
                # Get the flag (detection type)
                flag = ion['flag']
                flag_str = {0: "Detection", 1: "Upper Limit", 2: "Lower Limit"}.get(flag, "Unknown")
                
                # Add row to measurements
                all_measurements.append({
                    'Spectrum ID': spectrum_id,
                    'Redshift': redshift,
                    'Ion': ion_name,
                    'Rest Wavelength (Å)': ion['lam_0'],
                    'f-value': ion['f'],
                    'EW (mÅ)': ew,
                    'EW Error (mÅ)': ew_error,
                    'log N (cm^-2)': log_N,
                    'log N Error': log_N_err,
                    'Flag': flag,
                    'Detection Type': flag_str
                })
            
            print(f"Processed {ion_name} from {spectrum_id}")
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    # Create DataFrame
    df = pd.DataFrame(all_measurements)
    return df

def plot_ew_comparison(df, ion_names=None):
    """
    Create a comparison plot of EW measurements for selected ions
    across different spectra
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with measurements from batch_process_abstools_files
    ion_names : list, optional
        List of ion names to include in the plot, if None use all
    """
    if ion_names is None:
        ion_names = df['Ion'].unique()
    
    # Filter the data to only include the specified ions
    plot_data = df[df['Ion'].isin(ion_names)]
    
    # Get unique spectrum IDs
    spectra = plot_data['Spectrum ID'].unique()
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Set up colors and markers
    colors = plt.cm.tab10(np.linspace(0, 1, len(ion_names)))
    markers = ['o', 's', '^', 'D', '*', 'p', 'h', 'v', '>', '<']
    
    # Plot each ion
    for i, ion in enumerate(ion_names):
        ion_data = plot_data[plot_data['Ion'] == ion]
        
        # Get x positions with a small offset for each ion
        x_positions = np.arange(len(spectra)) + (i - len(ion_names)/2) * 0.15
        
        # For each spectrum, find the EW for this ion
        y_values = []
        y_errors = []
        x_labels = []
        
        for j, spectrum in enumerate(spectra):
            spectrum_data = ion_data[ion_data['Spectrum ID'] == spectrum]
            
            if not spectrum_data.empty:
                ew = spectrum_data['EW (mÅ)'].values[0]
                ew_err = spectrum_data['EW Error (mÅ)'].values[0]
                flag = spectrum_data['Flag'].values[0]
                
                y_values.append(ew)
                y_errors.append(ew_err)
                x_labels.append(spectrum)
                
                # Different marker for upper/lower limits
                marker = markers[0]  # Default
                if flag == 1:  # Upper limit
                    marker = markers[1]
                elif flag == 2:  # Lower limit
                    marker = markers[2]
                
                # Plot the point
                ax.errorbar(x_positions[j], ew, yerr=ew_err, 
                            fmt=marker, color=colors[i], markersize=8,
                            label=ion if j == 0 else "")
    
    # Customize plot
    ax.set_xticks(np.arange(len(spectra)))
    ax.set_xticklabels(spectra, rotation=45, ha='right')
    ax.set_ylabel('Equivalent Width (mÅ)')
    ax.set_title('Equivalent Width Comparison Across Spectra')
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(title='Ion')
    
    plt.tight_layout()
    return fig, ax

# Example usage
if __name__ == "__main__":
    # Replace with your directory containing JSON files
    json_files = "data/Spectrum_Analysis_*.json"
    
    # Process all files
    measurements_df = batch_process_abstools_files(json_files)
    
    if measurements_df is not None and not measurements_df.empty:
        # Print summary
        print("\nSummary of all measurements:")
        print(f"Total measurements: {len(measurements_df)}")
        print(f"Unique spectra: {measurements_df['Spectrum ID'].nunique()}")
        print(f"Unique ions: {measurements_df['Ion'].nunique()}")
        
        # Save to a master CSV file
        measurements_df.to_csv("all_abstools_measurements.csv", index=False)
        print("All measurements saved to all_abstools_measurements.csv")
        
        # Create a comparison plot for some common ions
        common_ions = ['OVI 1031.93', 'OVI 1037.62', 'CIV 1548.20', 'CIV 1550.78']
        available_ions = [ion for ion in common_ions if ion in measurements_df['Ion'].values]
        
        if available_ions:
            fig, ax = plot_ew_comparison(measurements_df, available_ions)
            plt.savefig("ew_comparison.png", dpi=300)
            plt.show()
            print("Comparison plot saved to ew_comparison.png")