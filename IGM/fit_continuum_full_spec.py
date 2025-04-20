"""
Full spectrum continuum fitting module that divides the spectrum into chunks, 
applies optimal polynomial fitting to each chunk, and then blends the results.

This module implements an advanced continuum fitting approach for quasar or galaxy spectra
by breaking the spectrum into overlapping chunks, finding the optimal polynomial order
for each chunk using BIC, and then blending the chunks together to create a smooth continuum.

Author: Rongmon Bordoloi
Contributors: Various
Last updated: April 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
from scipy.interpolate import interp1d
import warnings

# Try to import required modules from different possible locations
try:
    from rbcodes.IGM.rb_iter_contfit import rb_iter_contfit, fit_optimal_polynomial, calculate_bic
except ImportError:
    try:
        from IGM.rb_iter_contfit import rb_iter_contfit, fit_optimal_polynomial, calculate_bic
    except ImportError:
        print("Error: rb_iter_contfit module not found. Please ensure it's in your Python path.")
        print("You can install it with: pip install rbcodes")
        raise

# Check if linetools is available
try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
    LINETOOLS_AVAILABLE = True
except ImportError:
    LINETOOLS_AVAILABLE = False
    print("linetools package not found. To install, run:")
    print("pip install linetools")
    print("For full installation instructions, visit: https://github.com/linetools/linetools")


def generate_full_coverage_chunks(wave, flux, error, window_size=100, overlap=30, 
                                 method='uniform', wmin=None, wmax=None):
    """
    Generate chunks that ensure full coverage of the spectrum within the specified wavelength range.
    
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
    overlap : float
        Overlap between chunks in Angstroms
    method : str
        Method for selecting chunks: 'uniform' or 'features'
    wmin : float, optional
        Minimum wavelength to include (default: min of wave array)
    wmax : float, optional
        Maximum wavelength to include (default: max of wave array)
    
    Returns
    -------
    chunks : list of tuples
        List of (wave_chunk, flux_chunk, error_chunk, start_idx, end_idx) for each chunk
    """
    # Set default wavelength range if not specified
    if wmin is None:
        wmin = np.min(wave)
    if wmax is None:
        wmax = np.max(wave)
    
    # Ensure wmin and wmax are within the range of the wavelength array
    wmin = max(wmin, np.min(wave))
    wmax = min(wmax, np.max(wave))
    
    # Find indices corresponding to the wavelength range
    range_indices = np.where((wave >= wmin) & (wave <= wmax))[0]
    if len(range_indices) == 0:
        raise ValueError(f"No data points found in wavelength range {wmin}-{wmax}")
    
    # Extract the selected wavelength range
    range_wave = wave[range_indices]
    range_flux = flux[range_indices]
    range_error = error[range_indices]
    
    # Get wavelength range
    wave_min = np.min(range_wave)
    wave_max = np.max(range_wave)
    full_range = wave_max - wave_min
    
    # Effective window size (accounting for overlap)
    effective_window = window_size - overlap
    
    # Calculate number of chunks needed to cover the spectrum
    n_chunks = int(np.ceil(full_range / effective_window))
    
    # Ensure at least 3 chunks for better blending (unless the range is too small)
    min_chunks = min(3, max(1, int(full_range / 20)))  # At least 1 chunk, at most 3
    n_chunks = max(min_chunks, n_chunks)
    
    print(f"Creating {n_chunks} chunks with {window_size}Å size and {overlap}Å overlap")
    print(f"Wavelength range: {wave_min:.2f} - {wave_max:.2f} Å")
    
    chunks = []
    
    if method == 'uniform':
        # Calculate starting points for each chunk
        if n_chunks > 1:
            # Distribute chunks evenly across spectrum
            effective_range = full_range - window_size
            step = effective_range / (n_chunks - 1) if effective_range > 0 else 0
            chunk_starts = [wave_min + i * step for i in range(n_chunks)]
        else:
            # Just one chunk covering the whole spectrum
            chunk_starts = [wave_min]
        
        # Process each chunk
        for i, start_wave in enumerate(chunk_starts):
            # For the last chunk, ensure it reaches the end
            if i == n_chunks - 1:
                end_wave = wave_max
            else:
                end_wave = start_wave + window_size
            
            # Find indices corresponding to this wavelength range in the full array
            start_idx = np.max([range_indices[0], np.argmin(np.abs(wave - start_wave))])
            end_idx = np.min([range_indices[-1], np.argmin(np.abs(wave - end_wave))])
            
            # Ensure we have enough points for fitting
            if end_idx - start_idx < 10:
                print(f"Warning: Chunk at {start_wave:.1f} has too few points, expanding window")
                # Expand window until we have at least 10 points
                while end_idx - start_idx < 10 and (start_idx > range_indices[0] or end_idx < range_indices[-1]):
                    if start_idx > range_indices[0]:
                        start_idx -= 1
                    if end_idx < range_indices[-1]:
                        end_idx += 1
            
            # Extract the chunk
            wave_chunk = wave[start_idx:end_idx+1]
            flux_chunk = flux[start_idx:end_idx+1]
            error_chunk = error[start_idx:end_idx+1]
            
            chunks.append((wave_chunk, flux_chunk, error_chunk, start_idx, end_idx))
    
    elif method == 'features':
        # Identify features to guide chunk placement
        # Use a rolling standard deviation to identify variable regions
        window = max(21, int(len(range_wave)/50))  # Use an odd-sized window
        window = window + 1 if window % 2 == 0 else window  # Ensure it's odd
        
        # Pad flux array for edge handling
        padded_flux = np.pad(range_flux, (window//2, window//2), mode='edge')
        
        # Calculate rolling standard deviation
        std_array = np.zeros_like(range_flux)
        for i in range(len(range_flux)):
            std_array[i] = np.std(padded_flux[i:i+window])
        
        # Smooth the standard deviation array
        smooth_window = int(window/2)
        if smooth_window > 1:
            kernel = np.ones(smooth_window) / smooth_window
            std_array = np.convolve(std_array, kernel, mode='same')
        
        # Find significant features
        median_std = np.median(std_array)
        feature_threshold = median_std * 1.5
        feature_indices = np.where(std_array > feature_threshold)[0]
        
        # If no significant features found, fall back to uniform distribution
        if len(feature_indices) < n_chunks:
            print("No significant features found. Using uniform distribution.")
            return generate_full_coverage_chunks(wave, flux, error, window_size, overlap, 
                                              method='uniform', wmin=wmin, wmax=wmax)
        
        # Choose feature centers to guide chunk placement
        step = len(feature_indices) / n_chunks
        feature_centers = [feature_indices[int(i * step)] for i in range(n_chunks)]
        
        # Add first and last points to ensure coverage of edges
        if 0 not in feature_centers:
            feature_centers = [0] + feature_centers
        if len(range_wave)-1 not in feature_centers:
            feature_centers.append(len(range_wave)-1)
        
        # Sort indices
        feature_centers.sort()
        
        # Create chunks centered around these feature points
        for i, center_idx in enumerate(feature_centers):
            # Map the center index from range_wave to the full wave array
            center_wave = range_wave[center_idx]
            
            # Define chunk boundaries
            start_wave = center_wave - window_size/2
            end_wave = center_wave + window_size/2
            
            # Adjust boundaries to stay within spectrum
            start_wave = max(wave_min, start_wave)
            end_wave = min(wave_max, end_wave)
            
            # Find corresponding indices in the full array
            start_idx = np.argmin(np.abs(wave - start_wave))
            end_idx = np.argmin(np.abs(wave - end_wave))
            
            # Ensure indices are within the selected range
            start_idx = max(range_indices[0], start_idx)
            end_idx = min(range_indices[-1], end_idx)
            
            # Ensure we have enough points
            if end_idx - start_idx < 10:
                # Expand window until we have at least 10 points
                while end_idx - start_idx < 10 and (start_idx > range_indices[0] or end_idx < range_indices[-1]):
                    if start_idx > range_indices[0]:
                        start_idx -= 1
                    if end_idx < range_indices[-1]:
                        end_idx += 1
            
            # Extract the chunk
            wave_chunk = wave[start_idx:end_idx+1]
            flux_chunk = flux[start_idx:end_idx+1]
            error_chunk = error[start_idx:end_idx+1]
            
            chunks.append((wave_chunk, flux_chunk, error_chunk, start_idx, end_idx))
    
    else:
        raise ValueError(f"Unknown chunk selection method: {method}")
    
    # Verify full coverage of the selected wavelength range
    coverage = np.zeros_like(range_wave, dtype=bool)
    for _, _, _, start_idx, end_idx in chunks:
        # Map the chunk indices back to the range_wave indices
        chunk_start_in_range = max(0, start_idx - range_indices[0])
        chunk_end_in_range = min(len(range_wave)-1, end_idx - range_indices[0])
        if chunk_start_in_range <= chunk_end_in_range:
            coverage[chunk_start_in_range:chunk_end_in_range+1] = True
    
    # Check if entire selected range is covered
    if not np.all(coverage):
        uncovered = np.where(~coverage)[0]
        print(f"Warning: {len(uncovered)} points in selected range not covered by any chunk.")
        print(f"Wavelengths: {range_wave[uncovered[0]]:.2f} - {range_wave[uncovered[-1]]:.2f}")
        
        # Add extra chunks to cover any gaps
        while not np.all(coverage):
            # Find the largest gap
            gap_starts = []
            gap_ends = []
            
            # Identify gaps
            in_gap = False
            for i in range(len(coverage)):
                if not coverage[i] and not in_gap:
                    gap_starts.append(i)
                    in_gap = True
                elif coverage[i] and in_gap:
                    gap_ends.append(i-1)
                    in_gap = False
            
            # Handle the case where the last point is uncovered
            if in_gap:
                gap_ends.append(len(coverage)-1)
            
            # Find the largest gap
            gap_lengths = [gap_ends[i] - gap_starts[i] + 1 for i in range(len(gap_starts))]
            largest_gap_idx = np.argmax(gap_lengths)
            
            # Add a chunk to cover this gap
            gap_start_idx = gap_starts[largest_gap_idx]
            gap_end_idx = gap_ends[largest_gap_idx]
            
            # Map back to the full wave array indices
            start_idx = range_indices[0] + gap_start_idx
            end_idx = range_indices[0] + gap_end_idx
            
            # Expand chunk boundaries to include some overlap with existing chunks
            overlap_points = int(min(20, gap_lengths[largest_gap_idx] * 0.2))
            start_idx = max(range_indices[0], start_idx - overlap_points)
            end_idx = min(range_indices[-1], end_idx + overlap_points)
            
            # Extract the chunk
            wave_chunk = wave[start_idx:end_idx+1]
            flux_chunk = flux[start_idx:end_idx+1]
            error_chunk = error[start_idx:end_idx+1]
            
            chunks.append((wave_chunk, flux_chunk, error_chunk, start_idx, end_idx))
            
            # Update coverage
            chunk_start_in_range = max(0, start_idx - range_indices[0])
            chunk_end_in_range = min(len(range_wave)-1, end_idx - range_indices[0])
            coverage[chunk_start_in_range:chunk_end_in_range+1] = True
        
        print(f"Added extra chunks to ensure full coverage. Total chunks: {len(chunks)}")
    
    return chunks


def create_cosine_taper_weights(length):
    """
    Create a weight array with a cosine taper at both ends.
    
    Parameters
    ----------
    length : int
        Length of the array
    
    Returns
    -------
    array
        Weight array with values between 0 and 1
    """
    # Use a cosine taper (Tukey window)
    # Taper 30% at each end
    taper_fraction = 0.3
    taper_length = int(length * taper_fraction / 2)
    
    weights = np.ones(length)
    
    # Only apply taper if array is long enough
    if length > 4 and taper_length > 1:
        # Apply left taper
        for i in range(taper_length):
            weights[i] = 0.5 * (1 - np.cos(np.pi * i / taper_length))
        
        # Apply right taper
        for i in range(taper_length):
            weights[-(i+1)] = 0.5 * (1 - np.cos(np.pi * i / taper_length))
    
    return weights


def plot_continuum_results(wave, flux, error, continuum, chunks, chunk_results, wmin=None, wmax=None):
    """
    Plot the results of continuum fitting, showing how chunks contribute to the full continuum.
    
    Parameters
    ----------
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    continuum : array
        Fitted continuum array
    chunks : list of tuples
        List of chunks used for fitting
    chunk_results : list of dict
        Results for each chunk
    wmin : float, optional
        Minimum wavelength to display
    wmax : float, optional
        Maximum wavelength to display
    """
    # Set default wavelength range if not specified
    if wmin is None:
        wmin = np.min(wave)
    if wmax is None:
        wmax = np.max(wave)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(4, 1, figure=fig, height_ratios=[3, 1, 1, 1])
    
    # Filter indices within the wavelength range
    range_mask = (wave >= wmin) & (wave <= wmax)
    plot_wave = wave[range_mask]
    plot_flux = flux[range_mask]
    plot_continuum = continuum[range_mask]
    
    # Plot full spectrum with continuum
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(plot_wave, plot_flux, 'k-', alpha=0.7, label='Observed flux')
    ax1.plot(plot_wave, plot_continuum, 'r-', linewidth=2, label='Fitted continuum')
    
    # Highlight the chunks
    colors = plt.cm.tab10(np.linspace(0, 1, len(chunks)))
    
    # Plot chunk boundaries within the visible wavelength range
    for i, ((wave_chunk, _, _, _, _), color) in enumerate(zip(chunks, colors)):
        # Only show chunks that overlap with the displayed range
        if max(wave_chunk) >= wmin and min(wave_chunk) <= wmax:
            # Draw vertical lines at chunk boundaries
            chunk_min = max(min(wave_chunk), wmin)
            chunk_max = min(max(wave_chunk), wmax)
            
            ax1.axvline(x=chunk_min, color=color, linestyle='--', alpha=0.5)
            ax1.axvline(x=chunk_max, color=color, linestyle='--', alpha=0.5)
            
            # Add text label with polynomial order
            result = chunk_results[i]
            ypos = np.max(plot_flux) * (0.9 - 0.05 * (i % 5))  # Stagger labels
            ax1.text((chunk_min + chunk_max)/2, 
                     ypos,
                     f"Order {result['best_order']}", 
                     color=color, ha='center', fontsize=9)
    
    ax1.set_xlabel('Wavelength (Å)')
    ax1.set_ylabel('Flux')
    ax1.set_xlim(wmin, wmax)
    ax1.set_title('Spectrum with Fitted Continuum')
    ax1.legend()
    
    # Plot individual chunk fits
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    
    # Plot each chunk's fit with its own color
    for i, ((wave_chunk, flux_chunk, error_chunk, _, _), color) in enumerate(zip(chunks, colors)):
        # Skip chunks that don't overlap with the displayed range
        if max(wave_chunk) < wmin or min(wave_chunk) > wmax:
            continue
            
        # Get the continuum model from chunk results
        if 'continuum' in chunk_results[i]:
            chunk_continuum = chunk_results[i]['continuum']
            ax2.plot(wave_chunk, chunk_continuum, '-', color=color, 
                     label=f"Chunk {i+1} (Order {chunk_results[i]['best_order']})")
    
    ax2.set_ylabel('Flux')
    ax2.set_title('Individual Chunk Fits')
    ax2.legend(loc='upper right', fontsize=8, ncol=3)
    
    # Plot normalized flux
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax3.plot(plot_wave, plot_flux/plot_continuum, 'k-')
    ax3.axhline(y=1.0, color='r', linestyle='--')
    ax3.set_ylim([0, 2])
    ax3.set_ylabel('Normalized Flux')
    ax3.set_title('Continuum-Normalized Spectrum')
    
    # Plot BIC values for each chunk
    ax4 = fig.add_subplot(gs[3])
    
    # Only plot for chunks that overlap with the display range
    for i, result in enumerate(chunk_results):
        # Check if the chunk overlaps with the displayed range
        wave_chunk = chunks[i][0]
        if max(wave_chunk) < wmin or min(wave_chunk) > wmax:
            continue
            
        if 'bic_results' in result:
            bic_data = result['bic_results']
            orders = [order for order, _ in bic_data]
            bic_values = [bic for _, bic in bic_data]
            
            # Normalize BIC values for better comparison
            bic_min = min(bic_values)
            bic_values = [b - bic_min for b in bic_values]
            
            ax4.plot(orders, bic_values, 'o-', color=colors[i], 
                     label=f"Chunk {i+1} ({min(wave_chunk):.0f}-{max(wave_chunk):.0f}Å)")
    
    ax4.set_xlabel('Polynomial Order')
    ax4.set_ylabel('ΔBIC')
    ax4.set_title('BIC Model Selection (Lower is Better)')
    ax4.legend(loc='upper right', fontsize=8, ncol=3)
    
    plt.tight_layout()
    plt.show()


def fit_quasar_continuum(data, chunk_params=None, fitting_params=None, 
                        save_output=True, output_filename=None, plot=True):
    """
    Fit the continuum of a quasar or galaxy spectrum by dividing it into chunks,
    fitting each with an optimal polynomial order, and blending the results.
    
    This function implements a powerful approach to continuum fitting that works well
    for complex spectra with a mix of emission and absorption features. It breaks the
    spectrum into overlapping segments, finds the optimal polynomial order for each
    using the Bayesian Information Criterion (BIC), and then blends them together
    to create a smooth, continuous fit.
    
    Parameters
    ----------
    data : str or tuple
        Either a path to a FITS file or a tuple of (wave, flux, error) arrays.
        If tuple, error can be None and will be estimated from flux.
    
    chunk_params : dict, optional
        Parameters for dividing the spectrum:
        - window_size : float, size of each chunk in Angstroms (default: 100)
        - overlap_fraction : float, fraction of chunk size to overlap (default: 0.3)
        - method : str, method for selecting chunks: 'uniform' or 'features' (default: 'uniform')
        - wmin : float, minimum wavelength to include (default: min of wave array)
        - wmax : float, maximum wavelength to include (default: max of wave array)
    
    fitting_params : dict, optional
        Parameters for continuum fitting:
        - min_order : int, minimum polynomial order to test (default: 2)
        - max_order : int, maximum polynomial order to test (default: 6)
        - maxiter : int, maximum iterations for fitting (default: 25)
        - sigma : float, sigma clipping threshold (default: 3.0)
        - use_weights : bool, whether to use error spectrum as weights (default: True)
    
    save_output : bool, optional
        Whether to save the continuum-normalized spectrum (default: True)
    
    output_filename : str, optional
        Filename for the output normalized spectrum (default: input filename with _normalized suffix)
    
    plot : bool, optional
        Whether to plot the results (default: True)
    
    Returns
    -------
    dict
        A dictionary containing:
        - 'wave': wavelength array
        - 'flux': original flux array
        - 'error': error array
        - 'continuum': fitted continuum array
        - 'normalized_flux': flux/continuum
        - 'chunk_results': list of results for each chunk
        - 'output_file': path to saved output file (if save_output=True)
    
    Examples
    --------
    Example 1: Basic usage with a FITS file
    
    >>> from rbcodes.IGM import fit_continuum_full_spec as fc
    >>> results = fc.fit_quasar_continuum('quasar.fits')
    
    Example 2: Specify wavelength range and customized fitting parameters
    
    >>> results = fc.fit_quasar_continuum(
    ...     'quasar.fits',
    ...     chunk_params={
    ...         'window_size': 150,         # 150 Angstroms per chunk
    ...         'overlap_fraction': 0.4,    # 40% overlap between chunks
    ...         'method': 'uniform',        # Uniform chunk spacing
    ...         'wmin': 1200,               # Start at 1200 Angstroms
    ...         'wmax': 1800                # End at 1800 Angstroms
    ...     },
    ...     fitting_params={
    ...         'min_order': 2,             # Test polynomials starting at order 2
    ...         'max_order': 8,             # Maximum polynomial order to test
    ...         'maxiter': 20,              # Maximum fitting iterations
    ...         'sigma': 2.5,               # Sigma-clipping threshold
    ...         'use_weights': True         # Use error spectrum for weighting
    ...     },
    ...     save_output=True,
    ...     output_filename='normalized_quasar.fits'
    ... )
    
    Example 3: Provide wavelength, flux, and error arrays directly
    
    >>> import numpy as np
    >>> wave = np.linspace(1000, 2000, 1000)
    >>> flux = np.random.normal(1.0, 0.1, 1000)
    >>> error = np.ones_like(flux) * 0.1
    >>> results = fc.fit_quasar_continuum(
    ...     (wave, flux, error),
    ...     chunk_params={'wmin': 1200, 'wmax': 1800}
    ... )
    """
    # Set default parameters for chunks
    if chunk_params is None:
        chunk_params = {}
    
    # Set default parameters for fitting
    if fitting_params is None:
        fitting_params = {}
    
    # Check if data is a filename or arrays
    if isinstance(data, str):
        # It's a filename
        filename = data
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            return None
            
        if not LINETOOLS_AVAILABLE:
            print("Cannot run function: linetools package not available")
            return None
        
        # Load the spectrum
        print(f"Loading spectrum from {filename}...")
        sp = XSpectrum1D.from_file(filename)
        
        # Extract wavelength, flux, and error arrays
        wave = sp.wavelength.value
        flux = sp.flux.value
        error = sp.sig.value if sp.sig_is_set else np.ones_like(flux) * 0.1 * np.median(flux)
        
        # Set default output filename based on input file
        if output_filename is None and save_output:
            base, ext = os.path.splitext(filename)
            output_filename = f"{base}_normalized{ext}"
    
    elif isinstance(data, tuple) and len(data) >= 2:
        # It's a tuple of (wave, flux, error) arrays
        if len(data) >= 3:
            wave, flux, error = data
        else:
            wave, flux = data
            error = None
        
        # Validate arrays
        if wave is None or flux is None:
            print("Wavelength and flux arrays must be provided")
            return None
        
        if len(wave) != len(flux):
            print("Wavelength and flux arrays must have the same length")
            return None
        
        if error is not None and len(error) != len(flux):
            print("Error array must have the same length as flux array")
            return None
        
        # Ensure we have an error array
        if error is None:
            # Create a constant error array (10% of median flux)
            error = np.ones_like(flux) * 0.1 * np.median(flux)
            print("No error array provided. Using constant error (10% of median flux).")
        
        # Set a default output filename if none provided
        if output_filename is None and save_output:
            output_filename = "normalized_spectrum.fits"
    
    else:
        print("Data must be either a filename or a tuple of (wave, flux, error) arrays")
        return None
    
    # Extract wavelength range parameters
    wmin = chunk_params.get('wmin', None)
    wmax = chunk_params.get('wmax', None)
    
    # Filter data to selected wavelength range if specified
    if wmin is not None or wmax is not None:
        # Use actual array bounds if not specified
        if wmin is None:
            wmin = np.min(wave)
        if wmax is None:
            wmax = np.max(wave)
            
        # Ensure wmin and wmax are within the range of the wavelength array
        wmin = max(wmin, np.min(wave))
        wmax = min(wmax, np.max(wave))
        
        print(f"Selected wavelength range: {wmin:.2f} - {wmax:.2f} Å")
        
        # Get indices within the range
        range_mask = (wave >= wmin) & (wave <= wmax)
        range_indices = np.where(range_mask)[0]
        
        if len(range_indices) == 0:
            print(f"No data points found in wavelength range {wmin}-{wmax}")
            return None
    else:
        range_mask = np.ones_like(wave, dtype=bool)
        range_indices = np.arange(len(wave))
        wmin = np.min(wave)
        wmax = np.max(wave)
    
    print(f"Spectrum has {len(range_indices)} data points in selected range")
    
    try:
        # Extract parameters for chunking
        window_size = chunk_params.get('window_size', 100)
        overlap_fraction = chunk_params.get('overlap_fraction', 0.3)
        method = chunk_params.get('method', 'uniform')
        
        # Extract parameters for fitting
        min_order = fitting_params.get('min_order', 2)
        max_order = fitting_params.get('max_order', 6)
        maxiter = fitting_params.get('maxiter', 25)
        sigma = fitting_params.get('sigma', 3.0)
        use_weights = fitting_params.get('use_weights', True)
        
        # Calculate the overlap in Angstroms
        overlap = window_size * overlap_fraction
        
        # Generate chunks to cover the selected spectrum range
        chunks = generate_full_coverage_chunks(
            wave, flux, error, 
            window_size=window_size,
            overlap=overlap,
            method=method,
            wmin=wmin,
            wmax=wmax
        )
        
        print(f"Generated {len(chunks)} chunks to cover the spectrum")
        
        # Arrays to store full continuum and weights
        full_continuum = np.zeros_like(flux)
        weights = np.zeros_like(flux)
        chunk_results = []
        
        # Process each chunk using fit_optimal_polynomial
        for i, (wave_chunk, flux_chunk, error_chunk, start_idx, end_idx) in enumerate(chunks):
            print(f"Processing chunk {i+1}/{len(chunks)}: {wave_chunk[0]:.1f} - {wave_chunk[-1]:.1f} Å")
            
            # Call fit_optimal_polynomial for this chunk
            chunk_result = fit_optimal_polynomial(
                wave_chunk, 
                flux_chunk, 
                error=error_chunk, 
                min_order=min_order,
                max_order=max_order,
                maxiter=maxiter,
                sigma=sigma,
                use_weights=use_weights,
                plot=False  # Don't plot individual chunks
            )
            
            # Extract the best fit continuum and other data from the result
            best_fit = chunk_result['continuum']
            best_order = chunk_result['best_order']
            bic_results = chunk_result['bic_results']
            
            # Store results for this chunk
            chunk_results.append({
                'chunk_id': i+1,
                'wavelength_range': (wave_chunk[0], wave_chunk[-1]),
                'best_order': best_order,
                'bic_results': bic_results,
                'std_error': chunk_result['fit_error'],
                'indices': (start_idx, end_idx),
                'continuum': best_fit
            })
            
            # Create a smooth weight function for this chunk (cosine taper)
            chunk_weights = create_cosine_taper_weights(len(wave_chunk))
            
            # Add this chunk's continuum to the full continuum with weights
            full_continuum[start_idx:end_idx+1] += best_fit * chunk_weights
            weights[start_idx:end_idx+1] += chunk_weights
        
        # Verify complete coverage
        if np.any(weights == 0):
            zero_weight_indices = np.where(weights == 0)[0]
            print(f"Warning: {len(zero_weight_indices)} points have no continuum coverage!")
            print(f"Points at indices: {zero_weight_indices}")
            
            # Handle any gaps by interpolation
            if len(zero_weight_indices) > 0:
                print("Fixing gaps with interpolation...")
                # Find valid points
                valid_indices = np.where(weights > 0)[0]
                
                # Create interpolation function
                interp_func = interp1d(
                    wave[valid_indices], 
                    full_continuum[valid_indices], 
                    bounds_error=False, 
                    fill_value="extrapolate"
                )
                
                # Fill in missing values
                full_continuum[zero_weight_indices] = interp_func(wave[zero_weight_indices])
                weights[zero_weight_indices] = 1.0
        
        # Normalize by weights to get final continuum
        full_continuum = full_continuum / weights
        
        # Calculate normalized flux
        normalized_flux = flux / full_continuum
        
        # Save the normalized spectrum if requested
        output_file = None
        if save_output and LINETOOLS_AVAILABLE:
            if output_filename is None:
                # Create a default output filename
                if isinstance(data, str):
                    base, ext = os.path.splitext(data)
                    output_file = f"{base}_normalized{ext}"
                else:
                    output_file = "normalized_spectrum.fits"
            else:
                output_file = output_filename
            
            # Create the XSpectrum1D object with continuum
            normalized_spec = XSpectrum1D.from_tuple((wave, flux, error, full_continuum))
            normalized_spec.write_to_fits(output_file)
            print(f"Saved normalized spectrum to {output_file}")
        
        # Plot the results if requested
        if plot:
            plot_continuum_results(wave, flux, error, full_continuum, chunks, chunk_results, wmin=wmin, wmax=wmax)
        
        # Return the results
        return {
            'wave': wave,
            'flux': flux,
            'error': error,
            'continuum': full_continuum,
            'normalized_flux': normalized_flux,
            'chunk_results': chunk_results,
            'output_file': output_file
        }
        
    except Exception as e:
        print(f"Error processing the spectrum: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


# Example usage
if __name__ == "__main__":
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Fit a continuum across a full spectrum")
    parser.add_argument("filename", help="Path to spectrum FITS file")
    parser.add_argument("--window", type=float, default=100, help="Window size in Angstroms (default: 100)")
    parser.add_argument("--overlap", type=float, default=0.3, help="Overlap fraction (default: 0.3)")
    parser.add_argument("--method", choices=["uniform", "features"], default="uniform", help="Method for selecting chunks (default: uniform)")
    parser.add_argument("--min-order", type=int, default=2, help="Minimum polynomial order (default: 2)")
    parser.add_argument("--max-order", type=int, default=6, help="Maximum polynomial order (default: 6)")
    parser.add_argument("--maxiter", type=int, default=20, help="Maximum iterations (default: 20)")
    parser.add_argument("--sigma", type=float, default=2.5, help="Sigma clipping threshold (default: 2.5)")
    parser.add_argument("--no-weights", action="store_true", help="Don't use error spectrum as weights")
    parser.add_argument("--no-plot", action="store_true", help="Don't show plot")
    parser.add_argument("--no-save", action="store_true", help="Don't save output file")
    parser.add_argument("--output", help="Output filename")
    parser.add_argument("--wmin", type=float, help="Minimum wavelength to include")
    parser.add_argument("--wmax", type=float, help="Maximum wavelength to include")
    
    args = parser.parse_args()
    
    # Set parameters for the continuum fitting
    chunk_params = {
        'window_size': args.window,          # Angstroms
        'overlap_fraction': args.overlap,    # Fraction of window to overlap
        'method': args.method,               # Method for selecting chunks
        'wmin': args.wmin,                   # Minimum wavelength to include
        'wmax': args.wmax                    # Maximum wavelength to include
    }
    
    fitting_params = {
        'min_order': args.min_order,         # Minimum polynomial order to test
        'max_order': args.max_order,         # Maximum polynomial order to test
        'maxiter': args.maxiter,             # Maximum iterations
        'sigma': args.sigma,                 # Sigma clipping threshold
        'use_weights': not args.no_weights   # Use error as weights
    }
    
    # Run the function
    results = fit_quasar_continuum(
        args.filename, 
        chunk_params=chunk_params,
        fitting_params=fitting_params,
        save_output=not args.no_save,
        output_filename=args.output,
        plot=not args.no_plot
    )
    
    if results is not None:
        print("\nSummary of Chunk Results (Ordered by Wavelength):")
        for chunk in sorted(results['chunk_results'], key=lambda x: x['wavelength_range'][0]):
            print(f"Chunk {chunk['chunk_id']}: "
                  f"{chunk['wavelength_range'][0]:.1f} - {chunk['wavelength_range'][1]:.1f} Å, "
                  f"Best Order: {chunk['best_order']}, "
                  f"Std Error: {chunk['std_error']:.4f}")
        
        if not args.no_save:
            print(f"\nFull continuum fitting complete. "
                  f"Normalized spectrum saved to {results['output_file']}")
        else:
            print("\nFull continuum fitting complete. Output not saved (--no-save option used).")