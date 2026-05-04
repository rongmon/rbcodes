# Project Documentation
[Back to Main Page](../main_readme.md)

*Auto-generated documentation from docstrings*

## Modules

### rb_spectrum

Simple spectrum reader/writer for rbcodes
A lightweight replacement for linetools.XSpectrum1D focused on multispec GUI needs

### rb_utility

Several utility functions. (1) show for loop progress (2) load color list

### rb_x1d_id

Read in HST/COS x1d header files and give some info

### readmultispec

readmultispec.py
Read IRAF (echelle) spectrum in multispec format from a FITS file.
Can read most multispec formats including linear, log, cubic spline,
Chebyshev or Legendre dispersion spectra.
Usage: retdict = readmultispec(fitsfile, reform=True)
Inputs:
fitfile     Name of the FITS file
reform      If true (the default), a single spectrum dimensioned
            [4,1,NWAVE] is returned as flux[4,NWAVE].  If false,
            it is returned as a 3-D array flux[4,1,NWAVE].
Returns a dictionary with these entries:
flux        Array dimensioned [NCOMPONENTS,NORDERS,NWAVE] with the spectra.
            If NORDERS=1, array is [NCOMPONENTS,NWAVE]; if NCOMPONENTS is also
            unity, array is [NWAVE].  (This can be changed
            using the reform keyword.)  Commonly the first dimension
            is 4 and indexes the spectrum, an alternate version of
            the spectrum, the sky, and the error array.  I have also
            seen examples where NCOMPONENTS=2 (probably spectrum and
            error).  Generally I think you can rely on the first element
            flux[0] to be the extracted spectrum.  I don't know of
            any foolproof way to figure out from the IRAF header what the
            various components are.
wavelen     Array dimensioned [NORDERS,NWAVE] with the wavelengths for
            each order.
header      The full FITS header from pyfits.
wavefields  [NORDERS] List with the analytical wavelength
            description (polynomial coefficients, etc.) extracted from
            the header.  This is probably not very useful but is
            included just in case.
History:
Created by Rick White based on my IDL readechelle.pro, 2012 August 15
Apologies for any IDL-isms that remain!

## Classes

### rb_spectrum (`rb_spectrum`)

A simple 1D spectrum class for reading/writing spectral data
    
    Attributes
    ----------
    wavelength : Quantity
        Wavelength array with units
    flux : Quantity  
        Flux array with units
    sig : Quantity or None
        Error array with same units as flux (None if not available)
    co : Quantity or None
        Continuum array with same units as flux (None if not available)
    filename : str
        Source filename
    meta : dict
        Metadata dictionary
    _read_failed : bool
        Flag indicating if reading failed (for fallback to linetools)

### rb_spectrum_collection (`rb_spectrum`)

Container for multiple rb_spectrum objects that behaves like a list
    but maintains compatibility with multispec GUI expectations

## Functions

### __init__() (`rb_spectrum`)

Initialize collection
        
        Parameters
        ----------
        spectrum_list : list of rb_spectrum
            List of spectrum objects

### _ensure_quantity() (`rb_spectrum`)

Convert array to Quantity with given unit

### sig_is_set() (`rb_spectrum`)

Check if error array is available

### co_is_set() (`rb_spectrum`)

Check if continuum array is available

### npix() (`rb_spectrum`)

Number of pixels

### wvmin() (`rb_spectrum`)

Minimum wavelength

### wvmax() (`rb_spectrum`)

Maximum wavelength

### from_file() (`rb_spectrum`)

Load spectrum from file
        
        Parameters
        ----------
        filename : str
            Path to spectrum file
        **kwargs
            Additional arguments passed to readers
            
        Returns
        -------
        rb_spectrum

### from_tuple() (`rb_spectrum`)

Create spectrum from tuple of arrays
        
        Parameters
        ----------
        ituple : tuple
            (wave, flux), (wave, flux, error), or (wave, flux, error, continuum)
        sort : bool, optional
            Sort by wavelength
        **kwargs
            Additional arguments
            
        Returns
        -------
        rb_spectrum

### copy() (`rb_spectrum`)

Create a copy of the spectrum

### from_arrays() (`rb_spectrum`)

Create multiple spectra from 2D arrays
        
        Parameters
        ----------
        wave_2d : array-like, shape (nspec, npix)
            Wavelength arrays for multiple spectra
        flux_2d : array-like, shape (nspec, npix)  
            Flux arrays for multiple spectra
        error_2d : array-like, optional
            Error arrays for multiple spectra (same shape as flux_2d)
        continuum_2d : array-like, optional
            Continuum arrays for multiple spectra (same shape as flux_2d)
        wave_unit : str or astropy.units.Unit, optional
            Unit for wavelength (default: Angstrom)
        flux_unit : str or astropy.units.Unit, optional
            Unit for flux (default: dimensionless)
        filenames : list of str, optional
            Filenames for each spectrum
        **kwargs
            Additional arguments
            
        Returns
        -------
        list of rb_spectrum
            List of spectrum objects

### append() (`rb_spectrum`)

Add a spectrum to the collection

### write() (`rb_spectrum`)

Write spectrum to file (format determined by extension)
        
        Parameters
        ----------
        filename : str
            Output filename
        **kwargs
            Additional arguments passed to writer

### rb_write_fits() (`rb_spectrum`)

Write to multi-extension FITS file
        
        Parameters
        ----------
        filename : str
            Output filename
        clobber : bool
            Overwrite existing file

### rb_write_hdf5() (`rb_spectrum`)

Write to HDF5 file

### rb_write_json() (`rb_spectrum`)

Write to JSON file (rb_spectrum native format)

### rb_write_ascii() (`rb_spectrum`)

Write to ASCII file

### air2vac() (`rb_spectrum`)

Convert wavelengths from air to vacuum

### vac2air() (`rb_spectrum`)

Convert wavelengths from vacuum to air

### extend() (`rb_spectrum`)

Add multiple spectra to the collection

### sort_by_wavelength() (`rb_spectrum`)

Sort spectra by minimum wavelength

### get_wavelength_range() (`rb_spectrum`)

Get overall wavelength range of all spectra

### rb_read_spectrum() (`rb_spectrum`)

Read a spectrum from file with automatic format detection
    
    Parameters
    ----------
    filename : str
        Path to spectrum file
    meta : dict, optional
        User-provided metadata (will override file metadata)
    **kwargs
        Additional arguments
        
    Returns
    -------
    rb_spectrum
        Returns rb_spectrum object with _read_failed=True if reading fails

### _rb_read_fits() (`rb_spectrum`)

Read FITS file with format detection

### _rb_read_json() (`rb_spectrum`)

Read JSON file with auto-detection of rb_spec vs rb_spectrum format

### _rb_read_hdf5() (`rb_spectrum`)

Read HDF5 file

### _rb_read_ascii() (`rb_spectrum`)

Read ASCII file

### _rb_parse_rbspec_json() (`rb_spectrum`)

Parse rb_spec analysis JSON format

### _rb_parse_rbspectrum_json() (`rb_spectrum`)

Parse rb_spectrum native JSON format

### _rb_parse_fits_binary_table() (`rb_spectrum`)

Parse binary FITS table with flexible column detection

### _rb_parse_multi_extension_fits() (`rb_spectrum`)

Parse multi-extension FITS

### _rb_parse_single_extension_fits() (`rb_spectrum`)

Parse single extension FITS and look for separate error file

### _rb_parse_sdss_fits() (`rb_spectrum`)

Parse SDSS-style FITS (2D array with flux/error in different rows)

### _rb_parse_desi_brick() (`rb_spectrum`)

Parse DESI brick format

### _rb_get_table_column() (`rb_spectrum`)

Find column in binary table by trying multiple tag names

### _rb_setwave() (`rb_spectrum`)

Generate wavelength array from FITS header keywords

### rb_read_spec() (`rb_spectrum`)

Convenience function for reading spectra

### rb_write_spec() (`rb_spectrum`)

Convenience function for writing spectra

### filter_2D() (`filter_2d_spec`)

Function to filter a 2D spectrum while masking regions around emission lines with a 2D Gaussian kernel

### estimate_snr() (`compute_SNR_1d`)

Estimate the signal-to-noise ratio (SNR) per pixel for a given 1D spectrum.
    Optionally, the spectrum can be rebinned to a desired resolution before computing the SNR.

    Parameters:
    -----------
    wave : numpy.ndarray
        Wavelength array of the spectrum.
    flux : numpy.ndarray
        Flux array of the spectrum.
    error : numpy.ndarray
        Error (uncertainty) array of the spectrum.
    binsize : int, optional, default=3
        The number of pixels to bin together. If binsize > 1, the spectrum is rebinned.
    snr_range : list, optional, default=[-1, -1]
        The wavelength range [min, max] over which to compute the SNR. 
        If [-1, -1], the entire spectrum is used.
    verbose : bool, optional, default=True
        If True, prints status messages.
    plot : bool, optional, default=False
        If True, plots the SNR as a function of wavelength.
    robust_median : bool, optional, default=False
        If True, computes a sigma-clipped median of the SNR.
    sigma_clip_threshold : float, optional, default=3
        The sigma threshold to use for sigma clipping when computing the robust median SNR.

    Returns:
    --------
    snr_array : dict
        Dictionary containing the wavelength, SNR, flux, error, mean SNR, median SNR, and optionally robust SNR:
        - 'wave': Wavelength array in the selected range.
        - 'snr': Computed SNR per pixel.
        - 'flux': Flux values in the selected range.
        - 'error': Error values in the selected range.
        - 'mean_snr': Mean SNR over the selected range.
        - 'median_snr': Median SNR over the selected range.
        - 'robust_snr': Sigma-clipped median SNR (if robust_median=True).

    WARNING:
    --------
    The user should ensure that `wave`, `flux`, and `error` are proper numpy arrays without bad pixels.

    Example Usage:
    --------------
    >>> from utils import compute_SNR_1d as c
    >>> import numpy as np
    >>> wave = np.linspace(4000, 7000, 3000)  # Example wavelength array
    >>> flux = np.random.random(3000) * 100   # Example flux array
    >>> error = np.random.random(3000) * 10   # Example error array
    >>> snr_result = estimate_snr(wave, flux, error, binsize=5, snr_range=[4500, 6000], verbose=True, plot=True, robust_median=True, sigma_clip_threshold=2.5)

### rb_perccount() (`rb_utility`)

This code reports the percentage of a for loop done to the screen.
    **** NOTE *** This code currently only works on command line. Need to figure out how to make it work in notebook environment.

  
    Parameters
    ----------
       
       I = Current iteration between 0 and Imax
       Imax = Maximum number of iterations
  
    Returns
    ------- 
            Percentage of the for loop done.
            Shows the amount of time elapsed during the job.

    Example
    -------
    Calling Sequence: rb_perccount(I,Imax)
  
    WARNING: Do not print anything to the screen between calls of this function!
  
       Written by Rongmon Bordoloi Nov 15 2016
  --------------------------------------------------------------------

### format_interval() (`rb_utility`)

Formats a number of seconds as a clock time, [H:]MM:SS
        
        Parameters
        ----------
        t  : int
        Number of seconds.
        Returns
        -------
        out  : str
        [H:]MM:SS

### rb_set_color() (`rb_utility`)

Defines an expanded set of colors optimized for dark theme GUIs.
    Returns a dictionary with RGB values (0 to 1).

### rb_colorbar() (`rb_utility`)

Create a custom colorbar for plotting

### compute_galaxy_qso_pa() (`galaxy_qso_pa`)

Compute the position angle between a galaxy's major axis and a background QSO.
    
    Parameters:
    -----------
    fits_file : str
        Path to the FITS image file
    ra_gal, dec_gal : float
        Galaxy position in degrees (RA, Dec)
    ra_qso, dec_qso : float
        QSO position in degrees (RA, Dec)
    cutout_size : tuple, optional
        Size of cutout around galaxy in pixels (default: (100, 100))
    initial_sma : float, optional
        Initial semi-major axis for ellipse fitting in pixels (default: 10)
    initial_eps : float, optional
        Initial ellipticity for ellipse fitting (default: 0.3)
    moments_sigma : float, optional
        Gaussian smoothing sigma for moments method (default: 2.0)
    plot : bool, optional
        Create diagnostic plot (default: False)
    save_plot : str, optional
        Filename to save plot (default: None)
    bootstrap_n : int, optional
        Number of bootstrap samples for moments uncertainty (default: 100)
    
    Returns:
    --------
    dict : PA measurements and diagnostics

### _fit_isophote() (`galaxy_qso_pa`)

Try to fit isophote with multiple parameter attempts.

### _fit_moments() (`galaxy_qso_pa`)

Fit using second moments method with bootstrap uncertainty estimation.

### _calculate_moments() (`galaxy_qso_pa`)

Calculate second moments from image data.

### _unwrap_angles() (`galaxy_qso_pa`)

Handle angle wrapping for uncertainty calculations.

### _pixel_to_celestial_pa() (`galaxy_qso_pa`)

Convert pixel PA to celestial PA using WCS.

### _compute_alignment_angle() (`galaxy_qso_pa`)

Compute alignment angle (0-90°).

### _create_plot() (`galaxy_qso_pa`)

Create single-panel diagnostic plot with correct compass.

### _add_celestial_compass() (`galaxy_qso_pa`)

Add compass with correct celestial North/East directions.

### stiff_rgb() (`rgb_images`)

Emulate STIFF's RGB image generation using gamma + percentile stretch.
    Parameters
    ----------
    r, g, b : 2D numpy arrays
        Red, green, and blue channel images (same shape).
    min_percent : float
        Lower percentile for contrast stretch (e.g., 0.01).
    max_percent : float
        Upper percentile for contrast stretch (e.g., 99.99).
    gamma : float
        Gamma correction value (e.g., 2.2).
    color_balance : bool
        Whether to normalize each channel by its median.
    savefile : str
        If provided, save the RGB image to this path.
    dpi : int
        Resolution for saved image.
    Returns
    -------
    rgb : 3D numpy array
        RGB image array in shape (height, width, 3).

### lupton_rgb() (`rgb_images`)

Create RGB image using Lupton et al. (2004) asinh stretch algorithm.
    
    This algorithm is particularly good for astronomical images with wide
    dynamic range, preserving both faint and bright features.
    
    Parameters
    ----------
    r, g, b : 2D numpy arrays
        Red, green, and blue channel images (same shape).
    Q : float
        Softening parameter. Higher values preserve more faint detail.
        Typical range: 1-20. Default: 8.
    stretch : float
        Linear stretch parameter. Controls the overall brightness.
        Typical range: 0.1-2. Default: 0.5.
    minimum : float
        Minimum value to subtract from all images (background level).
        Default: 0.
    color_balance : bool
        Whether to normalize each channel to have equal contribution
        in bright regions. Default: True.
    savefile : str
        If provided, save the RGB image to this path.
    dpi : int
        Resolution for saved image.
    
    Returns
    -------
    rgb : 3D numpy array
        RGB image array in shape (height, width, 3).
        
    References
    ----------
    Lupton, R., Blanton, M. R., Fekete, G., Hogg, D. W., O'Mullane, W., 
    Szalay, A., & Wherry, N. (2004). Preparing red-green-blue images 
    from CCD data. PASP, 116, 133-137.

### load_and_create_rgb() (`rgb_images`)

Convenience function to load FITS files and create RGB image.
    
    Parameters
    ----------
    r_file, g_file, b_file : str
        Paths to FITS files for red, green, blue channels.
    method : str
        Either 'stiff' or 'lupton' for RGB creation method.
    **kwargs : dict
        Additional parameters passed to the RGB function.
    
    Returns
    -------
    rgb : 3D numpy array
        RGB image array.

### manual_rgb() (`rgb_images`)

Create RGB image with full manual control over each channel's scaling.
    
    Parameters
    ----------
    r, g, b : 2D numpy arrays
        Red, green, and blue channel images (same shape).
    r_min, r_max : float, optional
        Manual min/max values for red channel. If None, uses percentiles.
    g_min, g_max : float, optional  
        Manual min/max values for green channel. If None, uses percentiles.
    b_min, b_max : float, optional
        Manual min/max values for blue channel. If None, uses percentiles.
    r_scale, g_scale, b_scale : float
        Scaling factors to boost/diminish each channel (default: 1.0).
    r_percentiles, g_percentiles, b_percentiles : tuple
        (min_percentile, max_percentile) for each channel if using percentiles.
    use_percentiles : bool
        Whether to use percentiles or manual min/max values.
    stretch : str
        Stretch function: 'linear', 'sqrt', 'log', 'asinh', or 'power'.
    savefile : str
        If provided, save the RGB image to this path.
    dpi : int
        Resolution for saved image.
    
    Returns
    -------
    rgb : 3D numpy array
        RGB image array in shape (height, width, 3).
    scaling_info : dict
        Dictionary containing the scaling parameters used for each channel.

### get_limits() (`rgb_images`)

Get min/max values either from percentiles or manual input.

### apply_stretch() (`rgb_images`)

Apply different stretch functions.

### compare_methods() (`rgb_images`)

Create RGB images using both methods for comparison.
    
    Parameters
    ----------
    r, g, b : 2D numpy arrays
        Red, green, and blue channel images.
    savefile_prefix : str, optional
        If provided, save images as '{prefix}_stiff.png' and '{prefix}_lupton.png'
    
    Returns
    -------
    dict : Dictionary containing both RGB arrays

### write_ds9_regions() (`write_ds9_regions`)

Write RA, DEC coordinates to a DS9 region file with customizable circle size, color, and labels.

    Parameters:
        ra (array-like): Array of Right Ascension values in degrees.
        dec (array-like): Array of Declination values in degrees.
        output_file (str): Name of the output DS9 region file.
        radius (float): Circle radius in arcseconds.
        color (str): Color of the regions (e.g., 'red', 'blue', 'green').
        labels (array-like, optional): Text labels for each region (same length as ra & dec). Default is None.

### print_header() (`rb_x1d_id`)

This function will read all fits files in a folder and print out the header info for raw spectra.
    It will print filename, exposure type, and object type entries from the header.

    Written By: Rongmon Bordoloi April 2018

### nonlinearwave() (`readmultispec`)

Compute non-linear wavelengths from multispec string
    
    Returns wavelength array and dispersion fields.
    Raises a ValueError if it can't understand the dispersion string.

### readmultispec() (`readmultispec`)

Read IRAF echelle spectrum in multispec format from a FITS file
    
    Can read most multispec formats including linear, log, cubic spline,
    Chebyshev or Legendre dispersion spectra
    
    If reform is true, a single spectrum dimensioned 4,1,NWAVE is returned
    as 4,NWAVE (this is the default.)  If reform is false, it is returned as
    a 3-D array.

### compute_offset() (`compute_telescope_offset`)

Compute and print the telescope offset (East, North) to move from (RA1, DEC1) to (RA2, DEC2).

    Parameters:
        ra1, dec1 : float or str
            Right Ascension and Declination of starting point (degrees or sexagesimal).
        ra2, dec2 : float or str
            Right Ascension and Declination of target point (degrees or sexagesimal).
        unit : str, optional
            Unit for output offsets ("arcsec" or "arcmin"). Default is "arcsec".

    Returns:
        tuple : (east_offset, north_offset)
            Offset in chosen unit (arcsec or arcmin).

    -----
    # Example Usage:
    # Using degrees
    ra1, dec1 = 150.1, 2.3  # Start coordinates in degrees
    ra2, dec2 = 150.2, 2.35  # Target coordinates in degrees
    
    # Using sexagesimal (HMS/DMS)
    ra1_s, dec1_s = "10h00m00s", "+02d18m00s"
    ra2_s, dec2_s = "10h00m48s", "+02d21m00s"
    
    # Using colon-separated format
    ra1_colon, dec1_colon = "11:48:15.5526", "+52:51:56.180"
    ra2_colon, dec2_colon = "11:48:17.9932", "+52:52:01.502"
    
    # Compute and print offsets
    compute_offset(ra1_colon, dec1_colon, ra2_colon, dec2_colon, unit="arcsec")

