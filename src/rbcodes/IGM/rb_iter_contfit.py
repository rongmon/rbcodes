import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
import warnings
import os
from matplotlib.gridspec import GridSpec
# Read in the 1D spectrum to be analyzed
from pkg_resources import resource_filename
# Check if linetools is available, if not provide installation instructions
try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
    LINETOOLS_AVAILABLE = True
except ImportError:
    LINETOOLS_AVAILABLE = False
    print("linetools package not found. To install, run:")
    print("pip install linetools")
    print("We only use linetools to read a fits file in the example.")

def rb_iter_contfit(wave, flux, error=None, **kwargs):
    """Iterative continuum fitter using Legendre polynomials with sigma clipping
    
    This function fits a continuum to spectral data using an iterative process with 
    Legendre polynomials and sigma clipping to exclude outliers (such as absorption 
    or emission features). It uses the astropy FittingWithOutlierRemoval and 
    sigma_clip tools for robust fitting.
    
    Parameters
    ----------
    wave : ndarray
        Wavelength array
    flux : ndarray
        Flux array 
    error : ndarray, optional
        Error array. If None, constant error is estimated from flux using robust
        statistical methods (median absolute deviation).
    
    Optional Parameters
    ------------------
    maxiter : int
        Maximum number of fitting iterations (default: 25)
    order : int
        Polynomial order of Legendre fit (default: 4)
    sigma : float
        Sigma clipping threshold for outlier rejection (default: 3.0)
    use_weights : bool
        Whether to use error spectrum as weights (1/σ²) in the fitting process (default: True)
    include_model : bool
        Whether to include the astropy model and fitter objects in the return dictionary (default: False)
        This replaces the previous return_model parameter for consistency.
    silent : bool
        If True, suppress all print statements and plots (default: False)
    
    Returns
    ---------
    dict
        A dictionary containing the following keys:
        - 'continuum': Final fitted continuum array evaluated at input wavelength points
        - 'residuals': Residual array (flux/continuum) - useful for normalization
        - 'fit_error': Standard deviation of the residuals (measure of fit quality)
        - 'model': The fitted Legendre polynomial model object (included if include_model=True)
        - 'fitter': The fitter object with fit information including covariance matrix
                   (included if include_model=True)
        - 'wave': Input wavelength array
        - 'flux': Input flux array
        - 'error': Error array (input or estimated)
    
    Notes
    -----
    The function automatically handles:
    - Chip gaps (where flux and error are both zero)
    - Negative or zero error values (replaced with median error)
    - Negative or zero flux values (replaced with corresponding error values)
    - Potential emission or absorption features during initial fitting
    
    If `use_weights=True`, the fitting uses inverse variance weighting (1/σ²)
    to give more weight to points with smaller errors.
    
    Examples
    --------
    Basic usage with error array:
    
    >>> result = rb_iter_contfit(wavelength, flux, error)
    >>> cont = result['continuum']
    >>> normalized_flux = flux / cont
    
    With model return and custom parameters:
    
    >>> result = rb_iter_contfit(
    ...     wavelength, flux, error, 
    ...     order=5, 
    ...     sigma=2.5, 
    ...     include_model=True
    ... )
    >>> model = result['model']
    >>> new_wave = np.linspace(min(wavelength), max(wavelength), 1000)
    >>> new_cont = model(new_wave)
    
    Without error array (will be estimated):
    
    >>> result = rb_iter_contfit(wavelength, flux)
    >>> cont = result['continuum']
    
    Notes
    -----
    This function uses Legendre polynomials which are more numerically stable
    than simple polynomials for fitting across wide wavelength ranges.
    
    The sigma clipping approach progressively removes outliers that don't fit
    the polynomial model, which is particularly useful for spectra with
    absorption or emission features.
   
    ------------------------------
    Written by:  Rongmon Bordoloi
    Tested on Python 3.7  Sep 4 2019
    Modified April 2025 - Standardized return format to use dictionary
    Modified May 2025 - Added silent parameter to suppress output
    --------------------------
    """
    # Get optional parameters
    maxiter = kwargs.get('maxiter', 25)
    order = kwargs.get('order', 4)
    sigma_level = kwargs.get('sigma', 3.0)
    use_weights = kwargs.get('use_weights', False)
    silent = kwargs.get('silent', False)
    
    # For backward compatibility, check for both include_model and return_model
    include_model = kwargs.get('include_model', False) or kwargs.get('return_model', False)
    
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
    param_errors = None
    if use_weights and hasattr(fit, 'fit_info') and 'param_cov' in fit.fit_info:
        param_cov = fit.fit_info['param_cov']
        if param_cov is not None:
            # Extract the parameter errors from the diagonal of the covariance matrix
            param_errors = np.sqrt(np.diag(param_cov))
    
    # Prepare the result dictionary
    result = {
        'continuum': fit_final,
        'residuals': resid_final,
        'fit_error': fit_error,
        'wave': wave,
        'flux': flux,
        'error': error
    }
    
    # Include model objects if requested
    if include_model:
        result.update({
            'model': filtered_fit,
            'fitter': fit
        })
        
        # Include parameter errors if available
        if param_errors is not None:
            result['param_errors'] = param_errors
    
    return result

def calculate_bic(residuals, n_params, n_points, error=None):
    """
    Calculate the Bayesian Information Criterion (BIC) for a model fit.
    
    The BIC is defined as:
    
    BIC = k * ln(n) + n * ln(RSS/n)
    
    where:
    - k is the number of parameters in the model
    - n is the number of data points
    - RSS is the residual sum of squares (Σ(y_i - f_i)²)
    
    If error weights are provided, the weighted BIC is calculated as:
    
    BIC = k * ln(n) + n * ln(WRSS/n)
    
    where:
    - WRSS is the weighted residual sum of squares (Σ((y_i - f_i)/σ_i)²)
    - σ_i are the error values
    
    BIC penalizes models with more parameters to prevent overfitting.
    A lower BIC indicates a better fit, balancing goodness-of-fit with model complexity.
    
    Parameters
    ----------
    residuals : array
        Residuals from the model fit (observed - predicted)
    n_params : int
        Number of parameters in the model (polynomial order + 1)
    n_points : int
        Number of data points
    error : array, optional
        Error array for weighted BIC calculation. If provided, the residuals
        are weighted by 1/error² in the calculation.
    
    Returns
    -------
    float
        BIC value (lower is better)
    """
    if error is None:
        # Calculate the residual sum of squares
        rss = np.sum(residuals**2)
    else:
        rss = np.sum((residuals/error)**2)

    
    # Calculate BIC
    # BIC = k * ln(n) + n * ln(RSS/n)
    bic = n_params * np.log(n_points) + n_points * np.log(rss/n_points)
    
    return bic

def fit_optimal_polynomial(wave, flux, error=None, min_order=1, max_order=6, 
                           maxiter=20, sigma=3.0, use_weights=False, include_model=True, 
                           plot=True, save_plot=False, plot_filename=None, silent=False, **kwargs):
    """
    Fit a spectral region with the optimal polynomial order determined by Bayesian Information Criterion.
    
    This function tests multiple polynomial orders and selects the optimal one based on the
    Bayesian Information Criterion (BIC), which balances goodness-of-fit with model complexity.
    For each polynomial order, it uses the rb_iter_contfit function to perform the actual fitting
    with sigma clipping.
    
    Parameters
    ----------
    wave : ndarray
        Wavelength array
    flux : ndarray
        Flux array
    error : ndarray, optional
        Error array. If None, a constant error is estimated (10% of median flux).
    min_order : int, optional
        Minimum polynomial order to test (default: 1)
    max_order : int, optional
        Maximum polynomial order to test (default: 6)
    maxiter : int, optional
        Maximum iterations for sigma-clipping in rb_iter_contfit (default: 20)
    sigma : float, optional
        Sigma clipping threshold for outlier rejection (default: 3.0)
    use_weights : bool, optional
        Whether to use error spectrum as weights (default: False)
    include_model : bool, optional
        Whether to include the model objects in the result dictionary (default: True)
        Note: During testing of polynomial orders, include_model is always set to True internally.
    plot : bool, optional
        Whether to plot the results, including all tested polynomial fits and BIC values (default: True)
    save_plot : bool, optional
        Whether to save the plot to a file (default: False)
    plot_filename : str, optional
        Filename for the saved plot (default: "polynomial_fit.png")
    silent : bool, optional
        If True, suppress all print statements and plots (default: False)
    **kwargs : 
        Additional parameters passed directly to rb_iter_contfit
    
    Returns
    -------
    dict
        A dictionary containing the following keys:
        - 'wave': Input wavelength array
        - 'flux': Input flux array
        - 'error': Error array (input or estimated)
        - 'continuum': Best-fit continuum based on optimal polynomial order
        - 'normalized_flux': Flux divided by continuum
        - 'residuals': Residual array (flux/continuum)
        - 'best_order': Optimal polynomial order determined by BIC
        - 'fit_error': Standard deviation of residuals for best fit
        - 'bic_results': List of (order, BIC) tuples for all tested orders
        - 'model': The astropy model for the best fit (if include_model=True)
        - 'fitter': The astropy fitter object (if include_model=True)
        - 'param_errors': Parameter errors if available (if include_model=True)
        - 'plot_file': path to saved plot file (if save_plot=True)

    
    Notes
    -----
    The Bayesian Information Criterion (BIC) is calculated as:
        BIC = k * ln(n) + n * ln(RSS/n)
    where k is the number of parameters, n is the number of data points,
    and RSS is the residual sum of squares.
    
    If weights are used, the weighted residuals are used in the BIC calculation.
    
    The function automatically adjusts max_order if there are too few data points
    for the requested polynomial order.
    
    Examples
    --------
    Basic usage:
    
    >>> result = fit_optimal_polynomial(wavelength, flux, error)
    >>> normalized_flux = result['normalized_flux']
    >>> best_order = result['best_order']
    
    Fitting without plotting:
    
    >>> result = fit_optimal_polynomial(wavelength, flux, error, plot=False)
    
    Testing higher-order polynomials:
    
    >>> result = fit_optimal_polynomial(wavelength, flux, error, 
    ...                                min_order=3, max_order=10)
    
    Accessing the model (for predicting at new wavelengths):
    
    >>> result = fit_optimal_polynomial(wavelength, flux, error, include_model=True)
    >>> new_wave = np.linspace(wavelength.min(), wavelength.max(), 1000)
    >>> new_cont = result['model'](new_wave)
    
    Running in silent mode (no output or plots):
    
    >>> result = fit_optimal_polynomial(wavelength, flux, error, silent=True)
    """
    # Check input arrays
    if wave is None or flux is None:
        raise ValueError("Wavelength and flux arrays must be provided")
    
    if len(wave) != len(flux):
        raise ValueError("Wavelength and flux arrays must have the same length")
    
    # Ensure we have an error array
    if error is None:
        # Create a constant error array (10% of median flux)
        error = np.ones_like(flux) * 0.1 * np.median(flux)
        if not silent:
            print("No error array provided. Using constant error (10% of median flux).")
    
    # Ensure the number of points is sufficient for fitting
    n_points = len(wave)
    if n_points < max_order + 2:
        if not silent:
            print(f"Too few points ({n_points}) for requested max_order ({max_order})")
        max_order = n_points - 2
        if not silent:
            print(f"Reduced max_order to {max_order}")
    
    # Define the range of polynomial orders to test
    test_orders = range(min_order, max_order + 1)
    
    # Track results for each order
    all_fits = []
    bic_values = []
    
    if not silent:
        print(f"Testing polynomial orders from {min_order} to {max_order}...")
    
    # Test different polynomial orders
    for order in test_orders:
        # Call rb_iter_contfit for each order with dictionary return format
        result = rb_iter_contfit(
            wave, 
            flux, 
            error=error, 
            order=order,
            maxiter=maxiter,
            sigma=sigma,
            use_weights=use_weights,
            include_model=True,  # Always include model for BIC calculations
            silent=silent,  # Pass silent parameter to rb_iter_contfit
            **kwargs 
        )
        
        # Extract results from dictionary
        fit_result = result['continuum']
        residuals = result['residuals']
        std_error = result['fit_error']
        fit_model = result['model']
        fitter = result['fitter']
        
        # Calculate raw residuals (observed - predicted)
        raw_residuals = flux - fit_result
        
        # Calculate BIC (n_params = order + 1)
        if use_weights:
            bic_value = calculate_bic(raw_residuals, order + 1, n_points, error=error)
        else:
            bic_value = calculate_bic(raw_residuals, order + 1, n_points)
        
        # Store results
        all_fits.append((order, fit_result, residuals, std_error, bic_value, fit_model, fitter))
        bic_values.append((order, bic_value))
        
        if not silent:
            print(f"Order {order}: BIC = {bic_value:.2f}, StdDev = {std_error:.4f}")
    
    # Find the best order (minimum BIC)
    best_idx = np.argmin([bic for _, bic in bic_values])
    best_order, best_fit, best_residuals, best_std, best_bic, best_fit_model, best_fitter = all_fits[best_idx]
    
    if not silent:
        print(f"\nOptimal polynomial order: {best_order}")
        print(f"BIC value: {best_bic:.2f}")
        print(f"Standard deviation of residuals: {best_std:.4f}")
    
    # Calculate normalized flux
    normalized_flux = flux / best_fit
    
    # Prepare the result dictionary
    result_dict = {
        'wave': wave,
        'flux': flux,
        'error': error,
        'continuum': best_fit,
        'normalized_flux': normalized_flux,
        'residuals': best_residuals,
        'best_order': best_order,
        'fit_error': best_std,
        'bic_results': bic_values,
        'all_fits': all_fits  # Include all fits for plotting
    }
    
    # Only include model objects if requested
    if include_model:
        model_info = {
            'model': best_fit_model,
            'fitter': best_fitter
        }
        
        # Get parameter errors if available and weights were used
        if use_weights:
            param_errors = None
            if hasattr(best_fitter, 'fit_info') and 'param_cov' in best_fitter.fit_info:
                param_cov = best_fitter.fit_info['param_cov']
                if param_cov is not None:
                    param_errors = np.sqrt(np.diag(param_cov))
            
            if param_errors is not None:
                model_info['param_errors'] = param_errors
        
        # Update the result dictionary with model information
        result_dict.update(model_info)
 
    # Plot the results if requested and not in silent mode
    if plot and not silent:
        fig = plot_polynomial_fit_results(result_dict)
        
        # Save the plot if requested
        if save_plot:
            if plot_filename is None:
                # Create default plot filename
                plot_filename = "polynomial_fit.png"
            
            fig.savefig(plot_filename, dpi=300, bbox_inches='tight')
            if not silent:
                print(f"Saved plot to {plot_filename}")
            result_dict['plot_file'] = plot_filename
        
        plt.show()
    else:
        # If not showing plot but saving, and not in silent mode
        if save_plot and not silent:
            fig = plot_polynomial_fit_results(result_dict)
            
            if plot_filename is None:
                plot_filename = "polynomial_fit.png"
            
            fig.savefig(plot_filename, dpi=300, bbox_inches='tight')
            if not silent:
                print(f"Saved plot to {plot_filename}")
            result_dict['plot_file'] = plot_filename
            plt.close(fig)    
    # Remove 'all_fits' from the result dictionary if not needed for return
    # This avoids cluttering the return value with large data
    if 'all_fits' in result_dict and not plot:
        del result_dict['all_fits']

    return result_dict

def plot_polynomial_fit_results(result_dict, figure_size=(10, 8)):
    """
    Plot the results of polynomial fitting with BIC comparison.
    
    Creates a three-panel figure showing:
    1. Original spectrum with all polynomial fits
    2. Normalized flux (flux/continuum)
    3. BIC values for different polynomial orders
    
    Parameters
    ----------
    result_dict : dict
        Dictionary containing fitting results with the following keys:
        - 'wave': Wavelength array
        - 'flux': Flux array
        - 'continuum': Best-fit continuum
        - 'normalized_flux': Flux divided by continuum
        - 'best_order': Optimal polynomial order
        - 'bic_results': List of (order, BIC) tuples for all tested orders
        - 'all_fits': List of (order, fit, residuals, std_error, bic) for all orders
    figure_size : tuple, optional
        Size of the figure in inches (width, height) (default: (10, 8))
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure object, which can be further customized or saved
    
    Notes
    -----
    This function requires matplotlib to be installed.
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    
    # Extract data from result dictionary
    wave = result_dict['wave']
    flux = result_dict['flux']
    best_fit = result_dict['continuum']
    normalized_flux = result_dict['normalized_flux']
    best_order = result_dict['best_order']
    all_fits = result_dict['all_fits']
    bic_values = result_dict['bic_results']
    
    # Create figure with GridSpec for multiple panels
    fig = plt.figure(figsize=figure_size)
    gs = fig.add_gridspec(3, 1, height_ratios=[2, 1, 1], hspace=0.3)
    
    # Main plot - spectrum and fits
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(wave, flux, 'k-', label='Observed flux')
    
    # Get best fit BIC value
    best_bic = next(bic for order, bic in bic_values if order == best_order)
    
    # Plot best fit with emphasis
    ax1.plot(wave, best_fit, 'r-', linewidth=2, 
             label=f'Order {best_order} fit (BIC: {best_bic:.2f})')
    
    # Plot other fits for comparison
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_fits)))
    for i, (order, fit, _, _, bic, _, _) in enumerate(all_fits):
        if order != best_order:  # Skip the best fit (already plotted)
            ax1.plot(wave, fit, '-', color=colors[i], alpha=0.5,
                    label=f'Order {order} (BIC: {bic:.2f})')
    
    ax1.set_ylabel('Flux')
    
    # Adjust the figure to make room for the legend
    plt.subplots_adjust(right=0.75)
    
    # Add the legend outside the plot
    legend = ax1.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), 
                        title="Polynomial Fits", frameon=True, fancybox=True, shadow=True)
    
    # Second panel - normalized flux
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax2.plot(wave, normalized_flux, 'k-')
    ax2.axhline(y=1.0, color='r', linestyle='--')
    ax2.set_ylabel('Normalized Flux')
    
    # Third panel - BIC values
    ax3 = fig.add_subplot(gs[2])
    bic_values_array = np.array(bic_values)
    orders = bic_values_array[:, 0]
    bics = bic_values_array[:, 1]
    
    # Normalize BIC values for clearer visualization
    min_bic = np.min(bics)
    bics_normalized = bics - min_bic
    
    ax3.plot(orders, bics_normalized, 'o-', color='blue')
    ax3.scatter([best_order], [0], color='red', s=80, zorder=5)  # Best order is at 0 after normalization
    
    ax3.set_xlabel('Polynomial Order')
    ax3.set_ylabel('ΔBIC')
    ax3.set_title('BIC Values (lower is better)')
    ax3.grid(True, linestyle='--', alpha=0.5)
    
    # Only mark best order
    for i, (order, bic_norm) in enumerate(zip(orders, bics_normalized)):
        if order == best_order:
            ax3.annotate(f'Best Order: {int(order)}', 
                         xy=(order, bic_norm),
                         xytext=(10, +20),
                         textcoords='offset points',
                         arrowprops=dict(arrowstyle='->'))
    
    # Adjust x limits for first two subplots
    for ax in [ax1, ax2]:
        ax.set_xlim(wave.min(), wave.max())
    
    # Only show x-axis labels on the bottom plots
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    # Add a common title at the top of all subplots
    fig.suptitle('Polynomial Continuum Fitting Results', fontsize=14)
    
    # Ensure the layout adjusts properly
    plt.tight_layout(rect=[0, 0, 0.75, 0.95])  # Adjust for the legend
    
    return fig

def fit_spectral_region(filename, lam_min, lam_max, min_order=2, max_order=6, 
                        maxiter=20, sigma=3.0, use_weights=True, plot=True,
                        save_output=False, output_filename=None,
                        save_plot=False, plot_filename=None, silent=False):
    """
    Load a spectrum and fit a specific wavelength region with optimal polynomial order.
    
    Parameters
    ----------
    filename : str
        Path to the FITS file containing the spectrum
    lam_min : float
        Minimum wavelength of the region to fit
    lam_max : float
        Maximum wavelength of the region to fit
    min_order, max_order, maxiter, sigma, use_weights, plot :
        Parameters passed to fit_optimal_polynomial
    save_output : bool, optional
        Whether to save the normalized spectrum region (default: False)
    output_filename : str, optional
        Filename for the output normalized spectrum region
    silent : bool, optional
        If True, suppress all print statements and plots (default: False)
    
    Returns
    -------
    dict
        Results from fit_optimal_polynomial
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    if not silent:
        print(f"Loading spectrum from {filename}...")
    sp = XSpectrum1D.from_file(filename)
    
    # Extract wavelength, flux, and error arrays
    wave = sp.wavelength.value
    flux = sp.flux.value
    error = sp.sig.value if sp.sig_is_set else None
    
    # Select the wavelength region
    region_mask = (wave >= lam_min) & (wave <= lam_max)
    
    if not np.any(region_mask):
        raise ValueError(f"No data points found in wavelength range {lam_min}-{lam_max}")
    
    region_wave = wave[region_mask]
    region_flux = flux[region_mask]
    region_error = error[region_mask] if error is not None else None
    
    if not silent:
        print(f"Selected region: {region_wave[0]:.2f} - {region_wave[-1]:.2f} Å")
        print(f"Number of data points: {len(region_wave)}")
    
    # Fit the region with optimal polynomial order
    result = fit_optimal_polynomial(
        region_wave, 
        region_flux, 
        error=region_error,
        min_order=min_order,
        max_order=max_order,
        maxiter=maxiter,
        sigma=sigma,
        use_weights=use_weights,
        plot=plot,
        save_plot=save_plot,
        plot_filename=plot_filename,
        silent=silent  # Pass silent parameter to fit_optimal_polynomial
        )
    
    # Save the normalized region if requested
    if save_output and output_filename is not None:
        # Create a new XSpectrum1D object with the continuum
        normalized_spec = XSpectrum1D.from_tuple(
            (region_wave, region_flux, region_error, result['continuum'])
        )
        normalized_spec.write_to_fits(output_filename)
        if not silent:
            print(f"Saved normalized spectrum region to {output_filename}")
        result['output_file'] = output_filename
    
    return result

# Example usage
if __name__ == "__main__":
    # Path to the example FITS file
    filename = resource_filename('rbcodes', 'example-data/test.fits')

    
    # Define wavelength region to fit
    lam_min = 1390  # Ångstroms
    lam_max = 1530  # Ångstroms
    
    # Parameters for fitting
    min_order = 0
    max_order = 8
    
    # Fit the spectral region
    result = fit_spectral_region(
        filename,
        lam_min, 
        lam_max,
        min_order=min_order,
        max_order=max_order,
        use_weights=False,
        plot=True,
        save_output=False,
        output_filename='normalized_region.fits'
    )
    
    print(f"\nContinuum fitting complete for region {lam_min}-{lam_max} Å")
    print(f"Best polynomial order: {result['best_order']}")