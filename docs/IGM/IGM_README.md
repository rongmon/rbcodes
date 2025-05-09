# Project Documentation

*Auto-generated documentation from docstrings*

## Modules

### rb_setline

Read in atomic line information for a given or approximate rest frame wavelength.

This module provides functionality to find atomic line information based on different
matching methods including exact wavelength matching, closest wavelength matching,
or matching by species name.

### ransac_contfit

A continuum fitter class. This reads in a 1D spectrum and allows continuum fitting.

### rb_specbin

Rebin 1D spectrum to new pixel scale.

### lens_sep

Example code to plot to plot distances for differetn lens separations

### fit_continuum_full_spec

Full spectrum continuum fitting module that divides the spectrum into chunks, 
applies optimal polynomial fitting to each chunk, and then blends the results.

This module implements an advanced continuum fitting approach for quasar or galaxy spectra
by breaking the spectrum into overlapping chunks, finding the optimal polynomial order
for each chunk using BIC, and then blending the chunks together to create a smooth continuum.

Author: Rongmon Bordoloi
Contributors: Various
Last updated: April 2025

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

Compute physical separation between sightlines in a lensed quasar system.
    
    Parameters
    ----------
    delta_arcsec : float
        Angular separation between two quasars in arcsecond
    zabs_list : array-like
        List or array of absorber redshifts for which to compute physical separation
    z_lens : float
        Lens galaxy redshift
    z_source : float
        Background Quasar redshift
    custom_cosmo : astropy.cosmology object, optional
        Custom cosmology to use. Default is Planck18.
    return_with_units : bool, optional
        If True, returns result with astropy units attached. Default is False.
    
    Returns
    -------
    distlist : numpy.ndarray or astropy.units.Quantity
        Physical Separation for each absorber redshift in kpc
    
    Notes
    -----
    Equation used is equation (5) from Cooke et al. 2010.
    [Cooke, R., Pettini, M., Steidel, C. C., et al. 2010, MNRAS, 409, 679]
    Uses Planck 2018 LambdaCDM cosmology by default.
    
    Examples
    --------
    >>> from rbcodes.IGM import lens_sep_to_kpc as l
    >>> delta_arcsec = 1.  # 1 arcsecond separation
    >>> zabs_list = [.2, .5, 1.2, 2.]
    >>> z_lens = 0.55
    >>> z_source = 3.5
    >>> out = l.lens_sep_to_kpc(delta_arcsec, zabs_list, z_lens, z_source)
    ------------
    Written :- 
        Rongmon Bordoloi                               March 2 2021
        Updated for better error handling and units    April 2025

### rb_setline() (`rb_setline`)

Function to read in atomic line information for a given rest frame wavelength,
    for the line matching the closest wavelength, or by name.

    Parameters
    ----------
    lambda_rest : float
        Rest Frame wavelength (in Å) of the line to match
    method : str
        'closest' -> If set will match the closest line.
        'Exact'   -> If set will match the exact wavelength.
        'Name'    -> Match by name, USE WITH CARE. MUST INPUT OPTIONAL target_name
    linelist : str, optional
        The line list to use. Default is 'atom'.
        Available options: 'atom', 'LLS', 'LLS Small', 'DLA', 'LBG', 'Gal',
        'Eiger_Strong', 'Gal_Em', 'Gal_Abs', 'Gal_long', 'AGN', 'HI_recomb',
        'HI_recomb_light'
    target_name : str, optional
        Required when method='Name'. The name of the target species to match.

    Returns
    -------
    dict
        Dictionary with the following keys:
        - 'wave': Rest frame wavelength(s) of the matched line(s)
        - 'fval': Oscillator strength value(s)
        - 'name': Species name(s)
        - 'gamma': Radiation damping constant (if available in the line list)

    Examples
    --------
    Match the closest line to a given wavelength:
    
    >>> result = rb_setline(2796.3, 'closest')
    >>> print(f"Found: {result['name']} at {result['wave']} Å")
    
    Match a line by exact wavelength:
    
    >>> result = rb_setline(1215.67, 'Exact')
    
    Match a line by name:
    
    >>> result = rb_setline(0, 'Name', target_name='HI 1215')

    Raises
    ------
    ValueError
        If method is not one of 'closest', 'Exact', or 'Name'
        If method='Name' but target_name is not provided
        If method='Exact' but no lines match the wavelength
    FileNotFoundError
        If the requested line list file cannot be found

    Notes
    -----
    - For 'Exact' method, a match is considered if the wavelength is within 0.001 Å
    - The 'Name' method requires the target_name parameter to be set
    - The full list of available line lists can be found in the read_line_list function
    - Line lists are cached to improve performance when called multiple times

    History
    -------
    Written By: Rongmon Bordoloi                Jan 2018, Python 2.7
    Edit:       Rongmon Bordoloi                Sep 2018, Depreciated kwargs to be compatible with python 3
    Edit:       Improved version                Apr 2025, Various improvements while maintaining compatibility

### read_line_list() (`rb_setline`)

Read a line list defined by the label.

    Parameters
    ----------
    label : str
        Label string identifying which line list to read
        Available options: 'atom', 'LLS', 'LLS Small', 'DLA', 'LBG', 'Gal',
        'Eiger_Strong', 'Gal_Em', 'Gal_Abs', 'Gal_long', 'AGN', 'HI_recomb',
        'HI_recomb_light'

    Returns
    -------
    list of dict
        A list of dictionaries containing line data with keys:
        - 'wrest': Rest wavelength
        - 'ion': Ion name
        - 'fval': f-value (oscillator strength)
        - 'gamma': Radiation damping constant (for 'atom' line list only)

    Raises
    ------
    FileNotFoundError
        If the line list file cannot be found
    ValueError
        If an invalid line list label is provided

### from_file() (`ransac_contfit`)

Read spectrum from filename given.

### from_data() (`ransac_contfit`)

read spectrum from input wave,flux,error array.

### rb_specbin() (`rb_specbin`)

This function bins up 1D spectra in integer pixels. The routine returns a
       structure of flux and wavelength and variance that has been rebinned.
    
    Parameters
    -----------
     
        flux      - Input flux array
        nbin      - Number of pixels to bin on
        var=      - Input variance array [Optional]
        wave=     - Input wavelength array [Optional]
    
    Returns
    --------
        output    - Dictionary containing rebinned data with the following keys:
                    'flux': Rebinned flux array
                    'error': Rebinned error array (if var was provided)
                    'wave': Rebinned wavelength array (if wave was provided)
      
    Example
    --------
        output = rb_specbin(fx, 3)
        output = rb_specbin(fx, 3, var=var_array, wave=wave_array)
     
     
    REVISION HISTORY:
       Written by RB. June 2015
       Updated by RB. April 2025 - Improved error handling, documentation, and efficiency

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

### generate_full_coverage_chunks() (`fit_continuum_full_spec`)

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

### create_cosine_taper_weights() (`fit_continuum_full_spec`)

Create a weight array with a cosine taper at both ends.
    
    Parameters
    ----------
    length : int
        Length of the array
    
    Returns
    -------
    array
        Weight array with values between 0 and 1

### plot_continuum_results() (`fit_continuum_full_spec`)

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
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure object

### fit_quasar_continuum() (`fit_continuum_full_spec`)

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

    save_plot : bool, optional
        Whether to save the plot to a file (default: False)
    
    plot_filename : str, optional
        Filename for the plot (default: input filename with _continuum_fit.png suffix)

    
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
        - 'figure': matplotlib figure object (if plot=True)
        - 'plot_file': path to saved plot file (if save_plot=True)

    
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
    include_model : bool
        Whether to include the astropy model and fitter objects in the return dictionary (default: False)
        This replaces the previous return_model parameter for consistency.
    
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

