import numpy as np 
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats, sigma_clip

try:
    from rbcodes.IGM import rb_specbin as b
except:
    from IGM import rb_specbin as b

def estimate_snr(wave, flux, error, binsize=3, snr_range=[-1, -1], verbose=True, plot=False, robust_median=False, sigma_clip_threshold=3):
    """
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
    """
    snr_array = {}  # Initialize output dictionary

    # Bin spectrum if binsize > 1
    if binsize > 1:
        if verbose:
            print(f'Binning the spectrum by {binsize} pixels')
        sp = b.rb_specbin(flux, binsize, wave=wave, var=error**2)
    else:
        sp = {'wave': wave, 'flux': flux, 'error': error}

    # Define SNR computation range
    if snr_range == [-1, -1]:
        snr_range = [min(sp['wave']), max(sp['wave'])]
        if verbose:
            print(f'No wavelength range specified - estimating SNR over: [{snr_range[0]}, {snr_range[1]}]')
    else:
        if verbose:
            print(f'Estimating SNR over [{snr_range[0]}, {snr_range[1]}]')

    # Select data in the specified wavelength range
    mask = (sp['wave'] > snr_range[0]) & (sp['wave'] < snr_range[1])
    
    # Compute SNR
    snr = sp['flux'][mask] / sp['error'][mask]
    
    # Store results
    snr_array['wave'] = sp['wave'][mask]
    snr_array['snr'] = snr
    snr_array['flux'] = sp['flux'][mask]
    snr_array['error'] = sp['error'][mask]
    snr_array['mean_snr'] = np.nanmean(snr)
    snr_array['median_snr'] = np.nanmedian(snr)
    
    if robust_median:
        clipped_snr = sigma_clip(snr, sigma=sigma_clip_threshold, maxiters=5)
        robust_snr = np.nanmedian(clipped_snr.data[~clipped_snr.mask])
        snr_array['robust_snr'] = robust_snr
        clipped_wave = snr_array['wave'][~clipped_snr.mask]
        clipped_snr_values = snr[~clipped_snr.mask]
        if verbose:
            print(f'Robust SNR (sigma-clipped median): {robust_snr}')

    if verbose:
        print(f'Mean SNR within [{snr_range[0]}, {snr_range[1]}] is: {snr_array["mean_snr"]}')

    # Plot SNR if requested
    if plot:
        plt.figure(figsize=(12, 8))
        plt.plot(snr_array['wave'], snr, c='b', marker='.', markersize=8, markerfacecolor='k', 
                 label=f'SNR Binned by {binsize} pixels')
        plt.axhline(snr_array['mean_snr'], c='r', linestyle='-.', linewidth=4, label="Mean SNR")
        
        if robust_median:
            plt.axhline(snr_array['robust_snr'], c='g', linestyle='--', linewidth=3, label="Robust SNR")
            plt.scatter(clipped_wave, clipped_snr_values, c='orange', marker='o', label='Used for Sigma Clipping', alpha=0.6)
        
        plt.legend(fontsize=20)
        plt.xlabel('Wavelength [$\AA$]', size=20)
        plt.ylabel('SNR', size=20)
        plt.tight_layout()
        plt.show()

    return snr_array
