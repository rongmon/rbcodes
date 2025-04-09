import numpy as np
import matplotlib.pyplot as plt

def compute_EW(lam, flx, wrest, lmts, flx_err, plot=False, **kwargs):
    """
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
    """
    # Input Validation
    if not (len(lam) == len(flx) == len(flx_err)):
        raise ValueError("Input arrays (lam, flx, flx_err) must have equal lengths")
    
    # Validate wavelength array is monotonically increasing
    if not np.all(np.diff(lam) > 0):
        raise ValueError("Wavelength array must be monotonically increasing")
    
    # Validate velocity limits
    if len(lmts) != 2 or lmts[0] >= lmts[1]:
        raise ValueError("Invalid velocity limits. Must be [vmin, vmax] with vmin < vmax")

    # Extract optional parameters with default values
    verbose = kwargs.get('verbose', False)
    SNR = kwargs.get('SNR', False)
    zabs = kwargs.get('zabs', 0.0)
    sat_limit = kwargs.get('sat_limit', 'auto')
    f0 = kwargs.get('f0', None)
    _binsize = kwargs.get('_binsize', 1)
    
    # Normalization handling
    normalization = kwargs.get('normalization', 'none')
    if normalization == 'median':
        norm_factor = np.nanmedian(flx)
    elif normalization == 'mean':
        norm_factor = np.nanmean(flx)
    else:
        norm_factor = 1.0
    
    # Robust handling of input arrays
    with np.errstate(invalid='ignore', divide='ignore'):
        # Replace NaN and inf with robust estimates
        flx = np.nan_to_num(flx, 
                             nan=np.nanmedian(flx), 
                             posinf=np.nanmedian(flx), 
                             neginf=np.nanmedian(flx))
        flx_err = np.nan_to_num(flx_err, 
                                 nan=np.nanmedian(flx_err), 
                                 posinf=np.nanmedian(flx_err), 
                                 neginf=np.nanmedian(flx_err))
    
    # Normalize flux
    norm_flx = flx / norm_factor
    norm_flx_err = flx_err / norm_factor

    # Compute velocity array with improved numerical stability
    spl = 2.9979e5  # Speed of light in km/s
    vel = ((lam - wrest * (1.0 + zabs)) * spl / 
           (wrest * (1.0 + zabs) + np.finfo(float).eps))
    lambda_r = lam / (1. + zabs)

    # Select pixels within velocity limits
    pix = np.where((vel >= lmts[0]) & (vel <= lmts[1]))
    
    if len(pix[0]) == 0:
        # If no pixels in velocity range, return default dictionary
        default_output = {
            "ew_tot": np.nan,
            "err_ew_tot": np.nan,
            "vel_disp": np.nan,
            "vel50_err": np.nan,
            "line_saturation": False,
            "saturation_fraction": 0.0
        }
        if f0 is not None:
            default_output.update({
                "col": np.nan,
                "colerr": np.nan,
                "Tau_a": np.nan,
                "med_vel": np.nan
            })
        return default_output

    # New saturation handling based on median error
    if sat_limit == 'auto':
        # Use median error within the integration window as saturation threshold
        med_err = np.nanmedian(norm_flx_err[pix])
        if verbose:
            print(f"Using auto saturation limit: {med_err:.5f}")
        sat_threshold = med_err
    elif sat_limit is None:
        # Disable saturation handling
        sat_threshold = 0.0
    else:
        # Use user-provided saturation threshold
        sat_threshold = float(sat_limit)
    
    # Create a mask for saturated pixels within the integration window
    saturated_mask = norm_flx[pix] <= sat_threshold
    is_saturated = np.any(saturated_mask)
    saturation_fraction = np.sum(saturated_mask) / len(pix[0]) if len(pix[0]) > 0 else 0.0
    
    if verbose and is_saturated:
        print(f"WARNING: Line appears saturated! {saturation_fraction:.1%} of pixels " +
              f"within integration window are below the threshold ({sat_threshold:.5f})")

    # Prevent division by zero and handle low flux regions
    with np.errstate(invalid='ignore', divide='ignore'):
        # Avoid log(0) and division by zero
        norm_flx = np.clip(norm_flx, np.finfo(float).eps, None)
        
        # Handle saturated regions by replacing flux with error values
        if sat_threshold > 0:
            q_saturated = norm_flx <= sat_threshold
            norm_flx[q_saturated] = norm_flx_err[q_saturated] + np.finfo(float).eps

    # Compute wavelength differences with improved handling
    del_lam_j = np.diff(lambda_r)
    del_lam_j = np.pad(del_lam_j, (1, 0), mode='edge')

    # Normalized flux difference for EW computation
    Dj = 1. - norm_flx

    # Equivalent Width Per Pixel (in Angstroms)
    ew = del_lam_j[pix] * Dj[pix]

    # Error on Equivalent Width with robust error propagation
    sig_dj_sq = (norm_flx_err) ** 2.
    err_ew = del_lam_j[pix] * np.sqrt(sig_dj_sq[pix])
    err_ew_tot = np.sqrt(np.sum(err_ew ** 2.))
    ew_tot = np.sum(ew)

    # Compute velocity centroid weighted by equivalent width
    # Use cumulative sum with added numerical stability
    ew_safe = np.maximum(ew, np.finfo(float).eps)
    ew50 = np.cumsum(ew_safe) / np.sum(ew_safe)
    vel50 = np.interp(0.5, ew50, vel[pix])
    vel16 = np.interp(0.16, ew50, vel[pix])
    
    # Robust velocity dispersion calculation
    vel_disp = np.abs(vel50 - vel16)
    vel50_err = vel_disp / np.sqrt(len(ew))

    if verbose:
        print(f"W_lambda = {ew_tot:.3f} +/- {err_ew_tot:.3f} Å over [{lmts[0]:.1f} to {lmts[1]:.1f}] km/s")

    # Prepare output dictionary with saturation information
    output = {
        "ew_tot": ew_tot,
        "err_ew_tot": err_ew_tot,
        "vel_disp": vel_disp,
        "vel50_err": vel50_err,
        "line_saturation": is_saturated,
        "saturation_fraction": saturation_fraction
    }

    # Optional column density and AOD calculations
    if f0 is not None:
        # Compute Apparent Optical Depth (AOD) with improved numerical stability
        Tau_a = -np.log(np.clip(norm_flx, np.finfo(float).eps, None))

        # Compute column density as a function of velocity
        del_vel_j = np.diff(vel)
        del_vel_j = np.pad(del_vel_j, (1, 0), mode='edge')
        
        # Prevent division by zero
        nv = Tau_a / ((2.654e-15) * f0 * lambda_r + np.finfo(float).eps)
        n = nv * del_vel_j

        # Error propagation for column density
        tauerr = norm_flx_err / (norm_flx + np.finfo(float).eps)
        nerr = (tauerr / ((2.654e-15) * f0 * lambda_r + np.finfo(float).eps)) * del_vel_j

        # Compute total column density within selected pixels
        col = np.sum(n[pix])
        colerr = np.sqrt(np.sum(nerr[pix] ** 2.))

        if verbose:
            print(f"Direct N = {np.log10(col):.3f} +/- {np.log10(col + colerr) - np.log10(col):.3f} cm^-2")
            if is_saturated:
                print("WARNING: Column density estimate may be affected by saturation")

        # Update output with column density information
        output.update({
            "col": col,
            "colerr": colerr,
            "Tau_a": Tau_a,
            "med_vel": vel50
        })

    # Optional Signal-to-Noise Ratio calculation
    if SNR:
        from utils import compute_SNR_1d as c
        snr_result = c.estimate_snr(
            lambda_r, 
            norm_flx, 
            norm_flx_err, 
            binsize=_binsize, 
            snr_range=[-1, -1], 
            verbose=verbose, 
            robust_median=True, 
            sigma_clip_threshold=2.
        )
        output["SNR"] = snr_result['robust_snr']

    # Optional plotting
    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

        # Plot normalized flux and error
        ax1.step(vel, norm_flx, color='k', label="Normalized Flux")
        ax1.step(vel, norm_flx_err, color='r', label="Flux Error")
        
        # Highlight saturated regions
        if is_saturated:
            sat_indices = np.where((vel >= lmts[0]) & (vel <= lmts[1]) & (norm_flx <= sat_threshold))
            if len(sat_indices[0]) > 0:
                ax1.scatter(vel[sat_indices], norm_flx[sat_indices], color='orange', 
                         marker='x', s=50, label='Saturated Pixels')
        
        ax1.set_xlim([-600, 600])
        ax1.set_ylim([-0.02, 1.8])
        ax1.plot([-2500, 2500], [0, 0], 'k:')
        ax1.plot([-2500, 2500], [1, 1], 'k:')
        ax1.plot([lmts[0], lmts[0]], [1.5, 1.5], 'r+', markersize=15)
        ax1.plot([lmts[1], lmts[1]], [1.5, 1.5], 'r+', markersize=15)
        
        # Add information about saturation to the plot
        if is_saturated:
            ax1.axhline(y=sat_threshold, color='orange', linestyle='--', alpha=0.7, 
                      label=f'Saturation Threshold ({sat_threshold:.3f})')
            title_text = f'$W_{{rest}} = {ew_tot:.3f} \pm {err_ew_tot:.3f}$ Å [SATURATED: {saturation_fraction:.1%}]'
        else:
            title_text = f'$W_{{rest}} = {ew_tot:.3f} \pm {err_ew_tot:.3f}$ Å'
        
        if SNR:
            ax1.text(-600, 0.2, f'Median SNR: {output.get("SNR", "N/A"):.1f}')
        
        ax1.set_title(title_text)
        ax1.set_xlabel('Velocity [km/s]')
        ax1.legend(loc='upper right')

        # Plot column density as a function of velocity
        if f0 is not None:
            ax2.step(vel[pix], n[pix], color='b', label="Column Density")
            ax2.step(vel, n, color='k', lw=0.5, alpha=0.5)
            
            # Highlight saturated regions in column density plot
            if is_saturated:
                sat_indices = np.where((vel >= lmts[0]) & (vel <= lmts[1]) & (norm_flx <= sat_threshold))
                if len(sat_indices[0]) > 0:
                    ax2.scatter(vel[sat_indices], n[sat_indices], color='orange', 
                             marker='x', s=50, label='Saturated Pixels')
            
            ax2.set_xlim([-600, 600])
            ax2.set_xlabel('Velocity [km/s]')
            ax2.plot([-2500, 2500], [0, 0], 'k:')
            ax2.legend(loc='upper right')

        plt.tight_layout()
        plt.show()

    return output