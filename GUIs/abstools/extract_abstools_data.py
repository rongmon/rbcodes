import json
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.legendre as L

def load_abstools_analysis(json_file):
    """
    Load an AbsTools JSON analysis file and return the data as a dictionary
    
    Parameters:
    -----------
    json_file : str
        Path to the JSON file saved from AbsTools
        
    Returns:
    --------
    ions : dict
        Dictionary containing all the analysis data
    """
    # Custom JSON decoder to handle special types
    def json_object_hook(obj):
        if "_type" in obj:
            if obj["_type"] == "ndarray":
                return np.array(obj["data"], dtype=obj.get("dtype", None))
            elif obj["_type"] == "legendre":
                # Reconstruct Legendre polynomial
                pco = L.Legendre(obj["coef"])
                if "domain" in obj:
                    pco.domain = np.array(obj["domain"])
                if "window" in obj:
                    pco.window = np.array(obj["window"])
                return pco
        return obj
    
    # Read the JSON file
    with open(json_file, 'r') as f:
        data = json.load(f, object_hook=json_object_hook)
    
    # Extract ions dictionary from the data
    if "ions" in data:
        return data["ions"]
    else:
        return data  # For older format files that might not have the "ions" wrapper

# Example usage
if __name__ == "__main__":
    # Replace with your JSON file path
    file_path = "Spectrum_Analysis_z_0.348.json"
    
    # Load the analysis data
    ions = load_abstools_analysis(file_path)
    
    # Print some metadata information
    print(f"Target redshift: {ions['Target']['z']}")
    print(f"Number of analyzed transitions: {len(ions) - 1}")  # -1 for the Target entry
    
    # List all transitions (ions) in the file
    print("\nTransitions in the analysis:")
    for key in ions.keys():
        if key != 'Target':
            print(f"- {key}")
    
    # Let's extract data for a specific transition (e.g., first ion in the file)
    first_ion = [key for key in ions.keys() if key != 'Target'][0]
    print(f"\nExtracting data for: {first_ion}")
    
    # Get the continuum-normalized spectrum data
    ion_data = ions[first_ion]
    
    # Wavelength array
    wavelength = ion_data['wave']
    print(f"Wavelength range: {wavelength.min():.2f} - {wavelength.max():.2f} Å")
    
    # Flux and error arrays
    flux = ion_data['flux']
    error = ion_data['error']
    
    # Continuum fit
    continuum = ion_data['cont']
    
    # Velocity array 
    velocity = ion_data['vel']
    print(f"Velocity range: {velocity.min():.2f} - {velocity.max():.2f} km/s")
    
    # Polynomial fit coefficients (can be used to regenerate continuum)
    poly_order = ion_data['order']
    poly_coeffs = ion_data['pco']
    print(f"Polynomial order used: {poly_order}")
    
    # Measured values (if available)
    if ion_data['EW'] is not None:
        print(f"Equivalent Width: {ion_data['EW']:.2f} ± {ion_data['EWsig']:.2f} mÅ")
    
    if ion_data['N'] is not None and ion_data['N'] > 0:
        log_N = np.log10(ion_data['N'])
        log_N_err = np.log10(ion_data['N'] + ion_data['Nsig']) - log_N
        print(f"Column Density: log N = {log_N:.2f} ± {log_N_err:.2f}")
    
    # Integration limits for the EW measurement
    ew_limits = ion_data['EWlims']
    print(f"EW integration limits: {ew_limits[0]:.2f} to {ew_limits[1]:.2f} km/s")
    
    # Plot the data
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Plot the raw spectrum and continuum fit
    ax1.step(velocity, flux, 'k', where='mid', label='Flux')
    ax1.step(velocity, error, 'r', where='mid', label='Error')
    ax1.plot(velocity, continuum, 'b-', label='Continuum Fit')
    ax1.set_ylabel('Flux')
    ax1.set_title(f'{first_ion} - Raw Spectrum and Continuum Fit')
    ax1.legend()
    
    # Plot the normalized spectrum
    ax2.step(velocity, flux/continuum, 'k', where='mid', label='Normalized Flux')
    ax2.step(velocity, error/continuum, 'r', where='mid', label='Normalized Error')
    ax2.axhline(y=1, ls='--', c='b', label='Continuum')
    
    # Add vertical lines for the EW integration limits if they exist
    if ew_limits[0] is not None:
        ax2.axvline(x=ew_limits[0], ls='--', c='g', label='EW Lower Limit')
    if ew_limits[1] is not None:
        ax2.axvline(x=ew_limits[1], ls='--', c='m', label='EW Upper Limit')
    
    ax2.set_xlabel('Velocity (km/s)')
    ax2.set_ylabel('Normalized Flux')
    ax2.set_title(f'{first_ion} - Normalized Spectrum')
    ax2.legend()
    
    plt.tight_layout()
    plt.show()

    # Extract data for all transitions
    print("\nExtract data to a dict of numpy arrays for all transitions:")
    result_data = {}
    
    for key in ions.keys():
        if key != 'Target':
            result_data[key] = {
                'wavelength': ions[key]['wave'],
                'flux': ions[key]['flux'],
                'error': ions[key]['error'],
                'continuum': ions[key]['cont'],
                'velocity': ions[key]['vel'],
                'normalized_flux': ions[key]['flux'] / ions[key]['cont'],
                'normalized_error': ions[key]['error'] / ions[key]['cont'],
                'rest_wavelength': ions[key]['lam_0'],
                'oscillator_strength': ions[key]['f'],
                'ew_value': ions[key]['EW'],
                'ew_error': ions[key]['EWsig'],
                'column_density': ions[key]['N'],
                'column_density_error': ions[key]['Nsig'],
                'ew_limits': ions[key]['EWlims'],
                'flag': ions[key]['flag']  # 0=detection, 1=upper limit, 2=lower limit
            }
    
    print("Data extracted for the following transitions:")
    for key in result_data:
        print(f"- {key}: {len(result_data[key]['wavelength'])} data points")
    
    # Now you can access any transition's data using result_data[ion_name]
    # For example: result_data['OVI 1031.93']['normalized_flux']