# Project Documentation

*Auto-generated documentation from docstrings*

## Modules

### rb_setline

Read in atomic line information for a given or approximate rest frame  wavelength.

### ransac_contfit

A continuum fitter class. This reads in a 1D spectrum and allows continuum fitting.

### rb_specbin

Rebin 1D spectrum to new pixel scale.

### lens_sep

Example code to plot to plot distances for differetn lens separations

## Classes

### cont_fitter (`ransac_contfit`)

A continuum fitter class. This reads in a 1D spectrum and allows continuum fitting.
   The initial continuum is fitted using a random sample consensus model.
   Then user has the option to tweak the fitted continuum by hand usling interactive 
   linetools API and/or save the spectrum object.


    Parameters
    -----------
    
    filename=input spectrum filename
    efil = input error spectrum name [If it exists, otherwise None]
    window=    default [149], smoothing window
    mednorm = False [default], if set normalizes the spectrum with median value


    Returns
    --------

    This gives a cont_fitter object with following attributes:

    self.wave= wavelength

    self.flux= Flux
    self.error= error
    self.cont= Fitted continuum 


    Written : Rongmon Bordoloi      August 2020
    Edit    : Rongmon Bordoloi      March 2022: Added more input options

    Based on RANSAC continuum fitter written by Bin Liu Summer 2020.
    --------------------------------------------------------------------------------------------
    
    Example
    -------

               from IGM import ransac_contfit as c 

            Two ways to read in spectrum, from file: 
                 #efil = optional if error spectrum is defined in another file
               sp=c.cont_fitter.from_file(fluxfilename,efil=errorfilename)

            or from input wave,flux,error array. 
               sp=c.cont_fitter.from_data(wave,flux,error=error)
            

            Now fit continuum
               sp.fit_continuum(window=149)

               #AND YOU ARE DONE.

               #OPTIONAL:

               #Tweak the fitted continuum 
               sp.tweak_continuum()

               #Show new continuum
               sp.plot_spectrum()

               #Save it as a fits file
               sp.save_spectrum(outputfilename)

            # If the user wants to tweak an old fitted continuum, do the following sequence with linetools...
            ---->
            from linetools.spectra.xspectrum1d import XSpectrum1D
            sp=XSpectrum1D.from_file(filename)
            sp.fit_continuum()



    --------------------------------------------------------------------------------------------

### LLSFitter (`LLSFitter`)

Class for measuring the column density of Lyman Limit Systems (LLS)
    using both curve_fit and MCMC methods.

## Functions

### compute_EW() (`compute_EW`)

Function to compute the equivalent width (EW) within a given velocity window.

    Enhanced version with improved error handling, robustness, and dynamic saturation detection.

    Parameters:
        lam (array): Observed wavelength vector (in Angstroms).
        flx (array): Flux vector (same length as lam, preferably continuum normalized).
        wrest (float): Rest frame wavelength of the absorption line (in Angstroms).
        lmts (list): Velocity limits [vmin, vmax] in km/s.
        flx_err (array): Error spectrum (same length as flx).
        plot (bool): If True, will plot the spectrum and equivalent width. Default is False.

    Optional Parameters:
        zabs (float): Absorber redshift. Default is 0.
        f0 (float): Oscillator strength of the transition.
        sat_limit (float or str): Limit for saturation. Default is 'auto' which uses median error.
                                  Can also be a float value or None to disable saturation handling.
        verbose (bool): If True, will print detailed output. Default is False.
        SNR (bool): If True, computes the Signal-to-Noise Ratio. Default is False.
        normalization (str): Method of flux normalization. Options: 'median', 'mean', 'none'. Default is 'none'.
        _binsize (int): Binning size for SNR calculation. Default is 1.

    Returns:
        dict: A dictionary containing equivalent width measurements and related information.
            - 'ew_tot': Total rest frame equivalent width (in Angstroms).
            - 'err_ew_tot': Error on the total equivalent width.
            - 'vel_disp': 1-sigma velocity dispersion.
            - 'vel50_err': Error on the velocity centroid.
            - 'line_saturation': Boolean flag indicating if line is saturated.
            - 'saturation_fraction': Fraction of integration window that is saturated.
            - 'col': AOD column density.
            - 'colerr': Error on the AOD column density.
            - 'Tau_a': Apparent optical depth.
            - 'med_vel': Velocity centroid (EW-weighted velocity within velocity limits).
    
    Written:
        - Rongmon Bordoloi, 2nd November 2016
        - Translated from Matlab code `compute_EW.m`, which in turn is based on Chris Thom's `eqwrange.pro`.
        - Tested with COS-Halos/Dwarfs data.
        
    Edits:
        - RB, July 5, 2017: Output is a dictionary; edited minor dictionary arrangement.
        - RB, July 25, 2019: Added `med_vel` to the output.
        - RB, April 28, 2021: Modified `med_vel` to be weighted by `EW` and `vel_disp`.
        - RB, February 21, 2025: rewritten for clarity, added SNR keyword to compute signal-to-noise of the spectrum.
        - RB, April 8, 2025: Major improvements.
        - RB, April 9, 2025: Added dynamic saturation detection based on median error.

    Improvements:
    - Dynamic saturation detection based on median error in the integration window
    - Saturation flags and metrics in the output dictionary
    - More robust input validation
    - Better handling of NaN and inf values
    - Flexible normalization options
    - Improved error handling
    - Consistent error propagation

### lens_sep_to_kpc() (`lens_sep_to_kpc`)

------------------------------------------------------------------------------------------
       Function to compute physical separation between sightlines in a lensed quasar system.
       Input:- 
               delta_arcsec      :- Angular separation between two quasars in arcsecond
               zabs_list         :- List of absorber redshifts for which we want to compute physical separation
               z_lens            :- Lens galaxy redshift
               z_source          :- background Quasar redshift
    
    
       Output:- 
               distlist          :- numpy array with physical Separation for each absorber redshift in kpc
    
       Example:- 
                    >from IGM import lens_sep_to_kpc as l
                    > delta_arcsec=1. # 1 arcsecond separation
                    > zabs_list =[.2, .5,1.2,2.]
                    >z_lens =0.55
                    >z_source =3.5
                    > out = l.lens_sep_to_kpc(delta_arcsec,zabs_list,z_lens,z_source)
    
       Written :- Rongmon Bordoloi                           March 2 2021
    
       Equation used is the equation (5) from Cooke et al 2010.
       [Cooke, R., Pettini, M., Steidel, C. C., et al. 2010, MNRAS, 409, 679]
       Uses Planck 2018 \LambdaCDM cosmology
    
    
    ------------------------------------------------------------------------------------------

### rb_setline() (`rb_setline`)

Function to read in atomic line information for a given rest frame  wavelength.
                           Or 
    For the line matching the closest wavelength. 

    Parameters
    ----------
    lambda_rest :-  Rest Frame wavelength (in \AA) of the line to match
    method     :-   'closest' ->  If set will match the closest line.
                    'Exact'  ->  If set will match the exact wavelength.
                    'Name'   -> Match by name, USE WITH CARE. MUST INPUT OPTIONAL NAMELIST
 
    Returns
    ----------
    
    dic :- Dictionary with fval,lambda and species name.

    Example
    -------

       str=rb_setline(2796.3,'closest')


    Written By: Rongmon Bordoloi                Jan 2018, Python 2.7
    Edit:       Rongmon Bordoloi                            Sep 2018, Depreciated kwargs to be compatible with python 3

### read_line_list() (`rb_setline`)

Module to read a linelist defined by the label

    Parameters
    ----------
    lable : Label string [e.g. atom, LLS, LLS Small, LBG, Gal, Eiger_Strong]
      Must include redshift

    Returns
    ----------
    a dictionary with wrest, ion name and fvalues

### from_file() (`ransac_contfit`)

Read spectrum from filename given.

### from_data() (`ransac_contfit`)

read spectrum from input wave,flux,error array.

### rb_specbin() (`rb_specbin`)

This function bins up 1D spectra in integer pixels. The routine returns a
       structure of flux and wavelength and variance that has been rebinned.
    
    Parameters
    -----------
     
        fx       - Flux
        nbin     - Number of pixels to bin on
        VAR=  -- Input variance array [Optional]
        WAV=  -- Input wavelength array [Optional]
    
    Returns
    --------
        bin       - Structure of data
      
    Example
    --------
        bin = rb_specbin(fx, 3)
     
     
      REVISION HISTORY:
       Written by RB. June 2015
     -
     ------------------------------------------------------------------------------

### lens_sep() (`lens_sep`)

Example code to plot to plot distances for differetn lens separations

    Parameters

        zlist= list of redshifts

    Returns

        Plot of physical separation vs redshift

### __init__() (`LLSFitter`)

Initialize the LLSFitter object.
        
        Parameters:
        -----------
        spectrum_file : str, optional
            Path to the FITS file containing spectrum
        zabs : float, optional
            Absorption redshift

### load_spectrum() (`LLSFitter`)

Load spectrum from FITS file.
        
        Parameters:
        -----------
        spectrum_file : str
            Path to the FITS file

### set_redshift() (`LLSFitter`)

Set the absorption redshift and compute rest-frame wavelengths.
        
        Parameters:
        -----------
        zabs : float
            Absorption redshift

### set_sigma_clip() (`LLSFitter`)

Set the sigma threshold for outlier rejection when creating continuum mask.
        
        Parameters:
        -----------
        sigma : float
            Sigma threshold for outlier rejection

### set_domain_range() (`LLSFitter`)

Set the wavelength domain range for analysis.
        
        Parameters:
        -----------
        wmin : float, optional
            Minimum rest-frame wavelength in Angstroms
        wmax : float, optional
            Maximum rest-frame wavelength in Angstroms

### set_continuum_regions() (`LLSFitter`)

Set the regions used for continuum determination and fitting.
        
        Parameters:
        -----------
        regions : list of tuples, optional
            List of (min, max) wavelength ranges to use for continuum fitting
            If None, use default regions

### get_continuum_mask() (`LLSFitter`)

Create mask for continuum regions and apply sigma clipping to avoid absorption lines.
        
        Parameters:
        -----------
        sigma : float, optional
            Sigma threshold for outlier rejection
        
        Returns:
        --------
        mask : numpy array
            Boolean mask indicating continuum points after sigma clipping

### model_flx() (`LLSFitter`)

Physical model with a continuum and a Lyman limit absorption.
        
        Parameters:
        -----------
        theta : array-like
            Model parameters [C0, C1, logNHI]
        wave : array-like
            Wavelength array in Angstroms
            
        Returns:
        --------
        model : array-like
            Model flux

### model_test() (`LLSFitter`)

Wrapper for model_flx used in curve fitting (like scipy.optimize.curve_fit).

### lnprior() (`LLSFitter`)

Flat prior within given bounds.

### lnlike() (`LLSFitter`)

Log-likelihood function.

### lnprob() (`LLSFitter`)

Log-probability function.

### fit_curve_fit() (`LLSFitter`)

Fit the LLS model using scipy's curve_fit.
        
        Parameters:
        -----------
        theta_init : array-like, optional
            Initial parameter guess [C0, C1, logNHI]
            
        Returns:
        --------
        popt : array-like
            Best-fit parameters
        pcov : array-like
            Covariance matrix

### fit_emcee() (`LLSFitter`)

Fit the LLS model using MCMC with emcee.
        
        Parameters:
        -----------
        nwalkers : int, optional
            Number of walkers
        nsteps : int, optional
            Number of steps per walker
        burnin_frac : float, optional
            Fraction of steps to discard as burn-in
        theta_init : array-like, optional
            Initial parameter guess [C0, C1, logNHI]
            If None, use the results from curve_fit if available
            
        Returns:
        --------
        sampler : emcee.EnsembleSampler
            MCMC sampler object
        samples : array-like
            Flattened chain of samples

### plot_fit() (`LLSFitter`)

Plot the spectrum and the fit.
        
        Parameters:
        -----------
        method : str, optional
            'curve_fit' or 'mcmc'
        show_continuum_regions : bool, optional
            Whether to highlight the continuum regions with gray boxes
            and plot the points used in fitting
        wmin, wmax : float, optional
            Wavelength range to plot (if None, use domain_range)
        ymin, ymax : float, optional
            Flux range to plot (if None, auto-determine from data)
        figsize : tuple, optional
            Figure size
        show_realizations : bool, optional
            Whether to show random realizations (only for MCMC)
        n_realizations : int, optional
            Number of random realizations to show
            
        Returns:
        --------
        fig, ax : matplotlib figure and axis

### plot_corner() (`LLSFitter`)

Plot corner plot of MCMC samples.
        
        Returns:
        --------
        fig : matplotlib figure

### get_results_summary() (`LLSFitter`)

Get a summary of the fitting results.
        
        Returns:
        --------
        dict : Dictionary with fitting results

### compute_EW() (`compute_EW_old`)

Function to compute the equivalent width (EW) within a given velocity window.

    Parameters:
        lam (array): Observed wavelength vector (in Angstroms).
        flx (array): Flux vector (same length as lam, preferably continuum normalized).
        wrest (float): Rest frame wavelength of the absorption line (in Angstroms).
        lmts (list): Velocity limits [vmin, vmax] in km/s.
        flx_err (array): Error spectrum (same length as flx).
        plot (bool): If True, will plot the spectrum and equivalent width. Default is False.

    Optional Parameters:
        f0 (float): Oscillator strength of the transition.
        zabs (float): Absorber redshift. Default is 0.
        sat_limit (float): Limit for saturation. Default is 0.10 (COS specific).
        verbose (bool): If True, will print detailed output. Default is False.
        SNR (bool): If True, computes the Signal-to-Noise Ratio. Default is False.
        _binsize (int): Binning size for SNR calculation. Default is 1.

    Returns:
        dict: A dictionary containing the following keys:
            - 'ew_tot': Total rest frame equivalent width (in Angstroms).
            - 'err_ew_tot': Error on the total equivalent width.
            - 'vel_disp': 1-sigma velocity dispersion.
            - 'vel50_err': Error on the velocity centroid.
            - 'col': AOD column density.
            - 'colerr': Error on the AOD column density.
            - 'Tau_a': Apparent optical depth.
            - 'med_vel': Velocity centroid (EW-weighted velocity within velocity limits).

    Written:
        - Rongmon Bordoloi, 2nd November 2016
        - Translated from Matlab code `compute_EW.m`, which in turn is based on Chris Thom's `eqwrange.pro`.
        - Tested with COS-Halos/Dwarfs data.
        
    Edits:
        - RB, July 5, 2017: Output is a dictionary; edited minor dictionary arrangement.
        - RB, July 25, 2019: Added `med_vel` to the output.
        - RB, April 28, 2021: Modified `med_vel` to be weighted by `EW` and `vel_disp`.
        - RB, February 21, 2025: rewritten for clarity, added SNR keyword to compute signal-to-noise of the spectrum

### rb_iter_contfit() (`rb_iter_contfit`)

Iterative continuum fitter using Legendre polynomials with sigma clipping
    
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
    return_model : bool
        Whether to return the astropy model and fitter objects (default: False)
    
    Returns
    ---------
    fit_final : ndarray
        Final fitted continuum array evaluated at input wavelength points
    resid_final : ndarray
        Residual array (flux/continuum) - useful for normalization
    fit_error : float
        Standard deviation of the residuals (measure of fit quality)
    fit_model : astropy.modeling.Model, optional
        The fitted Legendre polynomial model object (returned if return_model=True)        
    fitter : astropy.modeling.fitting.LevMarLSQFitter, optional
        The fitter object with fit information including covariance matrix
        (returned if return_model=True)
    
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
    
    >>> cont, resid, error = rb_iter_contfit(wavelength, flux, error)
    >>> normalized_flux = flux / cont
    
    With model return and custom parameters:
    
    >>> cont, resid, error, model, fitter = rb_iter_contfit(
    ...     wavelength, flux, error, 
    ...     order=5, 
    ...     sigma=2.5, 
    ...     return_model=True
    ... )
    
    Without error array (will be estimated):
    
    >>> cont, resid, error = rb_iter_contfit(wavelength, flux)
    
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
    Highly modified for optimization April, 2025. 
    --------------------------

### calculate_bic() (`rb_iter_contfit`)

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

### fit_optimal_polynomial() (`rb_iter_contfit`)

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
        Whether to use error spectrum as weights (default: True)
    return_model : bool, optional
        Whether to return the model objects in the result dictionary (default: True)
        Note: During testing of polynomial orders, return_model is always set to True internally.
    plot : bool, optional
        Whether to plot the results, including all tested polynomial fits and BIC values (default: True)
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
        - 'best_order': Optimal polynomial order determined by BIC
        - 'fit_error': Standard deviation of residuals for best fit
        - 'bic_results': List of (order, BIC) tuples for all tested orders
        - 'fit_model': The astropy model for the best fit (if return_model=True)
        - 'fitter': The astropy fitter object (if return_model=True)
    
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
    
    >>> result = fit_optimal_polynomial(wavelength, flux, error, return_model=True)
    >>> new_wave = np.linspace(wavelength.min(), wavelength.max(), 1000)
    >>> new_cont = result['fit_model'](new_wave)

### plot_polynomial_fit_results() (`rb_iter_contfit`)

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

### fit_spectral_region() (`rb_iter_contfit`)

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
    
    Returns
    -------
    dict
        Results from fit_optimal_polynomial

### fit_polynomial_ransac() (`cont_fit_poly_ransac`)

Fits a polynomial using RANSAC to remove outliers and estimates uncertainty either by bootstrap resampling
    or by computing standard error if bootstrap is disabled.
    
    Parameters:
        wave (array): Wavelength values.
        flux (array): Observed flux values.
        error (array): Uncertainty in flux values.
        degree (int): Degree of the polynomial to fit.
        residual_threshold (float): Threshold for identifying outliers based on residuals.
        n_bootstrap (int or bool): Number of bootstrap resampling iterations, or False to skip bootstrap.
        verbose (bool): If True, prints warnings about poor fit or insufficient inliers.
        return_model (bool): If True, returns a function to evaluate the fitted model at new wavelengths.
    
    Returns:
        flux_fit (array): Best-fit polynomial evaluated at wave points.
        model_error (array): Estimated uncertainty in the fit.
        model (callable, optional): Function to compute the continuum model at new wavelength values.

