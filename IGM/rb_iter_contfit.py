import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
import warnings
import os
from matplotlib.gridspec import GridSpec

# Check if linetools is available, if not provide installation instructions
try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
    LINETOOLS_AVAILABLE = True
except ImportError:
    LINETOOLS_AVAILABLE = False
    print("linetools package not found. To install, run:")
    print("pip install linetools")
    print("For full installation instructions, visit: https://github.com/linetools/linetools")

# The cleaned version of the function
def rb_iter_contfit(wave, flux, error=None, **kwargs):
    """Iterative continuum fitter using Legendre polynomials
    
    Parameters
    ----------
    wave : ndarray
        Wavelength array
    flux : ndarray
        Flux array 
    error : ndarray, optional
        Error array. If None, constant error is estimated from flux.
    
    Optional Parameters
    ------------------
    maxiter : int
        Maximum iteration (default: 25)
    order : int
        Polynomial order of fit (default: 4)
    sigma : float
        Sigma clipping threshold (default: 3.0)
    use_weights : bool
        Whether to use error spectrum as weights (default: True)
    return_model : bool
        Whether to return the astropy model (default: False)
    
    Returns
    ---------
    fit_final : ndarray
        Final fitted continuum array
    resid_final : ndarray
        Residual error array (flux/continuum)
    fit_error : float
        Standard deviation of the residual
    fit_model : astropy.modeling.Model (optional)
        The fitted model object (returned if return_model=True)        
    fitter : astropy.modeling.fitting.LevMarLSQFitter (optional)
        The fitter object with fit information including covariance matrix
        (returned if return_model=True)
    ------------------------------
    Written by:  Rongmon Bordoloi
    Tested on Python 3.7  Sep 4 2019
    Highly modified for optimization April 16, 2025. 
    --------------------------
   
    """
    # Get optional parameters
    maxiter = kwargs.get('maxiter', 25)
    order = kwargs.get('order', 4)
    sigma_level = kwargs.get('sigma', 3.0)
    use_weights = kwargs.get('use_weights', True)
    return_model = kwargs.get('return_model', False)
    
    # Check if error array is provided, if not, estimate it
    if error is None:
        # Create a constant error array based on flux standard deviation
        # Use a robust estimator - median absolute deviation
        from astropy.stats import mad_std
        error_value = mad_std(flux)
        error = np.ones_like(flux) * error_value
        warnings.warn(
            f"No error spectrum provided. Using constant error of {error_value:.3e} "
            f"estimated from flux standard deviation.",
            UserWarning
        )
    
    # Initialize a mask
    mask = np.ones(np.size(wave))

    # Looking for chip gaps 
    chip_gap = np.where(((error == 0) & (flux == 0)) | (error - flux == 0))
    if np.size(chip_gap) > 0:
        mask[chip_gap] = 0

    # Now get rid of negative error values
    qq = np.where(error <= 0)
    if np.size(qq) > 0:
        median_error = np.median(error)
        error[qq] = median_error
        warnings.warn(
            f"Found {len(qq[0])} negative or zero error values. "
            f"Replaced with median error ({median_error:.3e}).",
            UserWarning
        )

    # Do a sanity check to avoid bad flux values
    q = np.where(flux <= 0)
    if np.size(q) > 0:
        flux[q] = error[q]
        warnings.warn(
            f"Found {len(q[0])} negative or zero flux values. "
            f"Replaced with corresponding error values.",
            UserWarning
        )

    # Calculate weights from error spectrum (1/σ²)
    w = 1/error**2 if use_weights else None

    # Find points outside chip gaps for statistics
    outside_chip_gap = np.where(((error != 0) & (flux != 0)) | (error - flux != 0))
    med_err = np.median(error[outside_chip_gap])
    med_flux = np.median(flux[outside_chip_gap])

    # Mask potential emission or absorption features (flux < median error)
    bd = np.where(flux < med_err)
    if len(bd[0]) > 0:
        mask[bd] = 0
    
    # Apply the mask
    qq = np.where(mask == 1)
    
    # Extract the data for fitting (excluding masked points)
    flux_new = flux[qq]
    wave_new = wave[qq]
    weights = w[qq] if use_weights else None
    
    # Fit the data using Legendre Polynomials
    g_init = models.Legendre1D(order)
    
    # Initialize fitters
    fit = fitting.LevMarLSQFitter()
    
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore') 
        
        # Use FittingWithOutlierRemoval for sigma clipping
        new_fit = fitting.FittingWithOutlierRemoval(
            fit, 
            sigma_clip, 
            niter=maxiter, 
            sigma=sigma_level
        )
        
        # Perform the fitting with or without weights
        if use_weights:
            filtered_fit, filtered_data = new_fit(g_init, wave_new, flux_new, weights=weights)
        else:
            filtered_fit, filtered_data = new_fit(g_init, wave_new, flux_new)

    # Apply the fit to the full wavelength array
    fit_final = filtered_fit(wave)
    
    # Calculate residuals
    resid_final = flux / fit_final
    
    # Calculate fit error from residual standard deviation
    # Only use non-masked points for calculating the error
    fit_error = np.std(resid_final[qq])
    
    # Get parameter errors from covariance matrix if using weights
    if use_weights and hasattr(fit, 'fit_info') and 'param_cov' in fit.fit_info:
        param_cov = fit.fit_info['param_cov']
        if param_cov is not None:
            # Extract the parameter errors from the diagonal of the covariance matrix
            param_errors = np.sqrt(np.diag(param_cov))
            # Could store or return these if needed
    
    # Return results
    if return_model:
        return fit_final, resid_final, fit_error, filtered_fit, fit
    else:
        return fit_final, resid_final, fit_error


# Function to extract spectral chunks
def extract_spectral_chunks(wave, flux, error, window_size=50, n_chunks=4, method='uniform'):
    """
    Extract small chunks from a spectrum for testing continuum fitting
    
    Parameters
    ----------
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    window_size : float
        Size of each chunk in Angstroms
    n_chunks : int
        Number of chunks to extract
    method : str
        Method for selecting chunks
        'uniform' - Uniformly spaced chunks across the spectrum
        'random' - Randomly selected chunks
        'features' - Chunks with significant features
        
    Returns
    -------
    chunks : list of tuples
        List of (wave_chunk, flux_chunk, error_chunk, start_idx, end_idx) for each chunk
    """
    # Get wavelength range
    wave_min = np.min(wave)
    wave_max = np.max(wave)
    full_range = wave_max - wave_min
    
    # Check if we have enough data for the requested window size and number of chunks
    if full_range < window_size * n_chunks:
        print(f"Warning: Spectrum range ({full_range:.1f} Å) is smaller than requested chunks ({window_size*n_chunks:.1f} Å)")
        print(f"Reducing window size to {full_range/n_chunks:.1f} Å")
        window_size = full_range / n_chunks
    
    chunks = []
    
    if method == 'uniform':
        # Select uniformly spaced chunks
        chunk_centers = np.linspace(wave_min + window_size/2, wave_max - window_size/2, n_chunks)
        
        for center in chunk_centers:
            start_wave = center - window_size/2
            end_wave = center + window_size/2
            
            # Find indices corresponding to this wavelength range
            start_idx = np.argmin(np.abs(wave - start_wave))
            end_idx = np.argmin(np.abs(wave - end_wave))
            
            # Ensure we have enough points
            if end_idx - start_idx < 10:
                print(f"Warning: Chunk at {center:.1f} has too few points, expanding window")
                # Expand window until we have at least 10 points
                while end_idx - start_idx < 10 and (start_idx > 0 or end_idx < len(wave)-1):
                    if start_idx > 0:
                        start_idx -= 1
                    if end_idx < len(wave)-1:
                        end_idx += 1
            
            # Extract the chunk
            wave_chunk = wave[start_idx:end_idx+1]
            flux_chunk = flux[start_idx:end_idx+1]
            error_chunk = error[start_idx:end_idx+1]
            
            chunks.append((wave_chunk, flux_chunk, error_chunk, start_idx, end_idx))
    
    elif method == 'random':
        # Select random chunks
        np.random.seed(42)  # For reproducibility
        
        for _ in range(n_chunks):
            # Select a random center point
            center_idx = np.random.randint(0, len(wave))
            center = wave[center_idx]
            
            # Define chunk boundaries
            start_wave = center - window_size/2
            end_wave = center + window_size/2
            
            # Find indices
            start_idx = np.argmin(np.abs(wave - start_wave))
            end_idx = np.argmin(np.abs(wave - end_wave))
            
            # Ensure boundaries are within range
            start_idx = max(0, start_idx)
            end_idx = min(len(wave)-1, end_idx)
            
            # Extract the chunk
            wave_chunk = wave[start_idx:end_idx+1]
            flux_chunk = flux[start_idx:end_idx+1]
            error_chunk = error[start_idx:end_idx+1]
            
            chunks.append((wave_chunk, flux_chunk, error_chunk, start_idx, end_idx))
    
    elif method == 'features':
        # Find regions with significant features
        
        # Calculate running standard deviation to find variable regions
        window = max(21, int(len(wave)/50))  # Use an odd-sized window
        window = window + 1 if window % 2 == 0 else window  # Ensure it's odd
        
        # Pad flux array for edge handling
        padded_flux = np.pad(flux, (window//2, window//2), mode='edge')
        
        # Calculate rolling standard deviation
        std_array = np.zeros_like(flux)
        for i in range(len(flux)):
            std_array[i] = np.std(padded_flux[i:i+window])
        
        # Find the regions with highest variance
        sorted_indices = np.argsort(std_array)[::-1]  # Sort in descending order
        
        # Select top regions, but ensure they don't overlap
        selected_centers = []
        min_separation = int(len(wave) / (n_chunks * 2))  # Minimum separation between chunks
        
        for idx in sorted_indices:
            # Check if this index is far enough from already selected centers
            if all(abs(idx - center) > min_separation for center in selected_centers):
                selected_centers.append(idx)
                
                # Break if we've found enough chunks
                if len(selected_centers) >= n_chunks:
                    break
        
        # Extract chunks around these centers
        for center_idx in selected_centers:
            center = wave[center_idx]
            
            # Define chunk boundaries
            start_wave = center - window_size/2
            end_wave = center + window_size/2
            
            # Find indices
            start_idx = np.argmin(np.abs(wave - start_wave))
            end_idx = np.argmin(np.abs(wave - end_wave))
            
            # Ensure boundaries are within range
            start_idx = max(0, start_idx)
            end_idx = min(len(wave)-1, end_idx)
            
            # Extract the chunk
            wave_chunk = wave[start_idx:end_idx+1]
            flux_chunk = flux[start_idx:end_idx+1]
            error_chunk = error[start_idx:end_idx+1]
            
            chunks.append((wave_chunk, flux_chunk, error_chunk, start_idx, end_idx))
    
    else:
        raise ValueError(f"Unknown chunk selection method: {method}")
    
    return chunks


# Test chunks with different polynomial orders
def test_chunks_with_orders(chunks, orders=[3, 4, 5]):
    """Test continuum fitting on each chunk with different polynomial orders"""
    
    all_results = []
    
    for i, (wave_chunk, flux_chunk, error_chunk, start_idx, end_idx) in enumerate(chunks):
        chunk_results = []
        
        # Try different polynomial orders
        for order in orders:
            # Run the continuum fitting
            fit_result, resid_final, fit_error = rb_iter_contfit(
                wave_chunk, 
                flux_chunk, 
                error=error_chunk, 
                order=order,
                use_weights=False
            )
            
            chunk_results.append((order, fit_result, resid_final, fit_error))
        
        all_results.append(chunk_results)
    
    return all_results


# Function to plot chunk results
def plot_chunk_results(chunks, all_results, chunk_info=None):
    """Plot the results of chunk testing with different polynomial orders"""
    
    n_chunks = len(chunks)
    n_orders = len(all_results[0])
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 4*n_chunks))
    gs = GridSpec(n_chunks, 2, figure=fig, width_ratios=[3, 1])
    
    for i, (chunk, chunk_results) in enumerate(zip(chunks, all_results)):
        wave_chunk, flux_chunk, error_chunk, start_idx, end_idx = chunk
        
        # Plot data and fits
        ax1 = fig.add_subplot(gs[i, 0])
        ax1.plot(wave_chunk, flux_chunk, 'k-', alpha=0.7, label='Observed flux')
        
        best_error = float('inf')
        best_order = None
        
        for j, (order, fit_result, resid_final, fit_error) in enumerate(chunk_results):
            ax1.plot(wave_chunk, fit_result, '-', label=f'Order {order} (std={fit_error:.4f})')
            
            if fit_error < best_error:
                best_error = fit_error
                best_order = order
        
        # Add title with chunk info if provided
        if chunk_info is not None and i < len(chunk_info):
            ax1.set_title(f"Chunk {i+1}: {chunk_info[i]} - Best Order: {best_order}")
        else:
            ax1.set_title(f"Chunk {i+1}: {wave_chunk[0]:.1f} - {wave_chunk[-1]:.1f} Å - Best Order: {best_order}")
        
        ax1.set_xlabel('Wavelength (Å)')
        ax1.set_ylabel('Flux')
        ax1.legend()
        
        # Plot residuals
        ax2 = fig.add_subplot(gs[i, 1])
        
        for j, (order, fit_result, resid_final, fit_error) in enumerate(chunk_results):
            # Offset residuals for clarity
            offset = j * 0.5
            ax2.plot(wave_chunk, resid_final + offset, '-', label=f'Order {order}')
            ax2.axhline(y=1.0 + offset, color='r', linestyle='--', alpha=0.5)
        
        ax2.set_title('Residuals')
        ax2.set_xlabel('Wavelength (Å)')
        ax2.set_ylabel('Flux / Continuum + offset')
        ax2.legend()
    
    plt.tight_layout()
    return fig


# Main function to test with spectral chunks
def test_spectral_chunks(filename, window_size=50, n_chunks=4, method='uniform', orders=[3, 4, 5, 6]):
    """
    Test continuum fitting on chunks of a spectrum with different polynomial orders
    
    Parameters
    ----------
    filename : str
        Path to the FITS file
    window_size : float
        Size of each chunk in Angstroms
    n_chunks : int
        Number of chunks to extract
    method : str
        Method for selecting chunks: 'uniform', 'random', or 'features'
    orders : list
        List of polynomial orders to test
    """
    if not LINETOOLS_AVAILABLE:
        print("Cannot run test: linetools package not available")
        return
    
    # Check if the file exists
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        print("Please ensure the file exists or provide the correct path.")
        return
    
    try:
        # Load the spectrum
        print(f"Loading spectrum from {filename}...")
        sp = XSpectrum1D.from_file(filename)
        
        # Extract wavelength, flux, and error arrays
        wave = sp.wavelength.value  # Convert from Quantity to numpy array
        flux = sp.flux.value
        error = sp.sig.value if sp.sig is not None else np.ones_like(flux) * 0.1 * np.median(flux)
        
        print(f"Loaded spectrum with {len(wave)} data points")
        print(f"Wavelength range: {wave.min():.2f} - {wave.max():.2f} Å")
        
        # Extract chunks
        print(f"Extracting {n_chunks} chunks of {window_size} Å using method '{method}'...")
        chunks = extract_spectral_chunks(wave, flux, error, window_size, n_chunks, method)
        
        # Test chunks with different polynomial orders
        print(f"Testing polynomial orders {orders} on each chunk...")
        all_results = test_chunks_with_orders(chunks, orders)
        
        # Generate chunk info for plots
        chunk_info = []
        for i, (wave_chunk, _, _, _, _) in enumerate(chunks):
            chunk_info.append(f"{wave_chunk[0]:.1f} - {wave_chunk[-1]:.1f} Å")
        
        # Plot results
        fig = plot_chunk_results(chunks, all_results, chunk_info)
        
        # Save and show the figure
        #plt.savefig(f'chunk_test_results_{method}_w{window_size}.png')
        plt.show()
        
        # Print summary of best orders for each chunk
        print("\nSummary of best polynomial orders for each chunk:")
        for i, (chunk_results, info) in enumerate(zip(all_results, chunk_info)):
            best_result = min(chunk_results, key=lambda x: x[3])  # Find lowest fit_error
            best_order, _, _, best_error = best_result
            print(f"Chunk {i+1} ({info}): Best order = {best_order} (std = {best_error:.4f})")
        
        return chunks, all_results
    
    except Exception as e:
        print(f"Error processing the spectrum: {str(e)}")
        import traceback
        traceback.print_exc()
        return None, None


# Plot full spectrum with chunk locations
def plot_full_spectrum_with_chunks(wave, flux, chunks):
    """Plot the full spectrum and highlight the extracted chunks"""
    
    plt.figure(figsize=(12, 6))
    plt.plot(wave,  flux, 'k-', alpha=0.5)
    
    colors = ['r', 'g', 'b', 'c', 'm', 'y']
    
    for i, (wave_chunk, flux_chunk, _, start_idx, end_idx) in enumerate(chunks):
        color = colors[i % len(colors)]
        plt.axvspan(wave_chunk[0], wave_chunk[-1], alpha=0.2, color=color)
        plt.plot(wave_chunk, flux_chunk, color=color, linewidth=2, label=f'Chunk {i+1}')
    
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux')
    plt.title('Full spectrum with extracted chunks')
    plt.legend()
    plt.tight_layout()
    #plt.savefig('full_spectrum_with_chunks.png')
    plt.show()


if __name__ == "__main__":
    # Default parameters
    default_filename = '/Users/bordoloi/WORK/python/rbcodes/example-data/test.fits'
    default_window_size = 50  # Angstroms
    default_n_chunks = 4
    default_method = 'uniform'  # 'uniform', 'random', or 'features'
    default_orders = [3, 4, 5, 6]
    
    # Check if the default file exists
    if not os.path.exists(default_filename):
        filename = input(f"Default file '{default_filename}' not found. Enter path to a FITS spectrum file: ")
    else:
        filename = default_filename
    
    # Get user input for parameters
    print("\nParameters for spectral chunk testing:")
    print(f"1. Window size (default: {default_window_size} Å)")
    print(f"2. Number of chunks (default: {default_n_chunks})")
    print(f"3. Selection method (default: '{default_method}')")
    print("   Options: 'uniform', 'random', 'features'")
    print(f"4. Polynomial orders to test (default: {default_orders})")
    print("\nPress Enter to use defaults or enter new values.")
    
    # Get window size
    window_size_input = input("Window size in Angstroms: ")
    window_size = float(window_size_input) if window_size_input.strip() else default_window_size
    
    # Get number of chunks
    n_chunks_input = input("Number of chunks: ")
    n_chunks = int(n_chunks_input) if n_chunks_input.strip() else default_n_chunks
    
    # Get selection method
    method_input = input("Selection method ('uniform', 'random', 'features'): ")
    method = method_input if method_input.strip() in ['uniform', 'random', 'features'] else default_method
    
    # Get polynomial orders
    orders_input = input("Polynomial orders (comma-separated integers, e.g., '3,4,5,6'): ")
    if orders_input.strip():
        try:
            orders = [int(x.strip()) for x in orders_input.split(',')]
        except ValueError:
            print("Invalid input for polynomial orders. Using default.")
            orders = default_orders
    else:
        orders = default_orders
    
    print(f"\nRunning test with: window_size={window_size}, n_chunks={n_chunks}, method='{method}', orders={orders}")
    
    # Load the spectrum to get wave and flux arrays for the full spectrum plot
    try:
        sp = XSpectrum1D.from_file(filename)
        wave = sp.wavelength.value
        flux = sp.flux.value
        
        # Run the test
        chunks, all_results = test_spectral_chunks(filename, window_size, n_chunks, method, orders)
        
        # Plot full spectrum with chunk locations
        if chunks is not None:
            plot_full_spectrum_with_chunks(wave, flux, chunks)
            
    except Exception as e:
        print(f"Error: {str(e)}")