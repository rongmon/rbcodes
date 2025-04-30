import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Union, Optional, Any, Tuple

def compute_EW(
    lam: np.ndarray,
    flx: np.ndarray,
    wrest: float,
    lmts: List[float],
    flx_err: np.ndarray,
    plot: bool = False,
    zabs: float = 0.0,
    f0: Optional[float] = None,
    sat_limit: Union[float, str, None] = 'auto',
    normalization: str = 'none',
    verbose: bool = False,
    SNR: bool = False,
    _binsize: int = 1,
    **kwargs
) -> Dict[str, Any]:
    """
    Function to compute the equivalent width (EW) within a given velocity window.

    Enhanced version with improved error handling, robustness, and dynamic saturation detection.

    Parameters:
        lam (np.ndarray): Observed wavelength vector (in Angstroms). Must be monotonically increasing.
        flx (np.ndarray): Flux vector (same length as lam, preferably continuum normalized).
        wrest (float): Rest frame wavelength of the absorption line (in Angstroms).
        lmts (List[float]): Velocity limits [vmin, vmax] in km/s. Must have vmin < vmax.
        flx_err (np.ndarray): Error spectrum (same length as flx).
        plot (bool): If True, will plot the spectrum and equivalent width. Default is False.
        zabs (float): Absorber redshift. Default is 0.
        f0 (Optional[float]): Oscillator strength of the transition. Required for column density calculation.
        sat_limit (Union[float, str, None]): Limit for saturation detection:
                                             - 'auto': Uses median error within integration window (default)
                                             - float value: Custom threshold
                                             - None: Disables saturation handling
        normalization (str): Method of flux normalization:
                             - 'none': No normalization (default)
                             - 'median': Normalize by median flux
                             - 'mean': Normalize by mean flux
        verbose (bool): If True, will print detailed output. Default is False.
        SNR (bool): If True, computes the Signal-to-Noise Ratio. Default is False.
        _binsize (int): Binning size for SNR calculation. Default is 1.

    Returns:
        Dict[str, Any]: A dictionary containing equivalent width measurements and related information.
            - 'ew_tot': Total rest frame equivalent width (in Angstroms).
            - 'err_ew_tot': Error on the total equivalent width.
            - 'vel_disp': 1-sigma velocity dispersion.
            - 'vel50_err': Error on the velocity centroid.
            - 'line_saturation': Boolean flag indicating if line is saturated.
            - 'saturation_fraction': Fraction of integration window that is saturated.
            - 'col': AOD column density (only if f0 provided).
            - 'colerr': Error on the AOD column density (only if f0 provided).
            - 'Tau_a': Apparent optical depth (only if f0 provided).
            - 'med_vel': Velocity centroid (only if f0 provided).
            - 'SNR': Signal-to-noise ratio (only if SNR=True).
    
    Examples:
        # Basic usage with minimal parameters
        from rbcodes.IGM.compute_EW import compute_EW
        result = compute_EW(wavelength, flux, 1215.67, [-100, 100], flux_error)
        
        # Full analysis with column density calculation and plotting
        result = compute_EW(
            wavelength, flux, 1215.67, [-150, 150], flux_error, 
            plot=True, zabs=0.1, f0=0.4164, verbose=True
        )
        
        # Check for saturation with custom threshold
        result = compute_EW(
            wavelength, flux, 1215.67, [-100, 100], flux_error,
            sat_limit=0.1, normalization='median'
        )
    
    Written:
        - Rongmon Bordoloi, 2nd November 2016
        - Tested with COS-Halos/Dwarfs data.
        
    Edits:
        - RB, July 5, 2017: Output is a dictionary; edited minor dictionary arrangement.
        - RB, July 25, 2019: Added `med_vel` to the output.
        - RB, April 28, 2021: Modified `med_vel` to be weighted by `EW` and `vel_disp`.
        - RB, February 21, 2025: rewritten for clarity, added SNR keyword to compute signal-to-noise of the spectrum.
        - RB, April 8, 2025: Major improvements.
        - RB, April 9, 2025: Added dynamic saturation detection based on median error.
        - RB, April 26, 2025: plotting updates+type annotations added 

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

    if len(lam) == 0:
        raise ValueError("Input arrays cannot be empty")

    
    # Validate wavelength array is monotonically increasing
    if not np.all(np.diff(lam) > 0):
        raise ValueError("Wavelength array must be monotonically increasing")
    
    # Validate velocity limits
    if len(lmts) != 2 or lmts[0] >= lmts[1]:
        raise ValueError("Invalid velocity limits. Must be [vmin, vmax] with vmin < vmax")


    
    # Normalization handling
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
        if f0 is not None:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        else:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))


        # Plot normalized flux and error
        ax1.step(vel, norm_flx, color='k', label="Normalized Flux",lw=0.5,where='mid')
        ax1.step(vel, norm_flx_err, color='r', label="Flux Error",lw=0.5,alpha=0.5,where='mid')

        # Add velocity centroid marker if available
        if not np.isnan(vel50):
            ax1.axvline(x=vel50, color='green', linestyle='-', alpha=0.7,
                      label=f'Velocity Centroid ({vel50:.1f} km/s)')
        
        

        # Highlight saturated regions
        if is_saturated:
            sat_indices = np.where((vel >= lmts[0]) & (vel <= lmts[1]) & (norm_flx <= sat_threshold))
            if len(sat_indices[0]) > 0:
                ax1.scatter(vel[sat_indices], norm_flx[sat_indices], color='orange', 
                         marker='x', s=50, label='Saturated Pixels')
        
        # Calculate appropriate y limits based on data
        y_max = min(2.0, np.nanpercentile(norm_flx[pix], 95) * 1.5)
        y_min = -0.01
        ax1.set_ylim([y_min, y_max])
        
        # Set x limits to show context around the integration window
        ax1.set_xlim([min(vel),max(vel)])
        
        # Add reference lines
        ax1.axhline(y=0, color='k', linestyle=':', alpha=0.5,lw=0.5)
        ax1.axhline(y=1, color='k', linestyle=':', alpha=0.5,lw=0.5)
        
        # Add vertical lines at integration limits
        ax1.axvline(x=lmts[0], color='blue', linestyle='--', alpha=0.5)
        ax1.axvline(x=lmts[1], color='blue', linestyle='--', alpha=0.5)
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
            ax1.text(min(vel)+200, 0.2, f'Median SNR: {output.get("SNR", "N/A"):.1f}')
        
        ax1.set_title(title_text)
        ax1.set_xlabel('Velocity [km/s]')
        ax1.legend(loc='upper right')

        # Plot column density as a function of velocity
        if f0 is not None:
            ax2.step(vel[pix], n[pix], color='b', label="Column Density",where='mid')
            ax2.step(vel, n, color='k', lw=0.5, alpha=0.5,where='mid')

            # Highlight velocity centroid
            if not np.isnan(vel50):
                idx_v50 = np.argmin(np.abs(vel[pix] - vel50))
                if idx_v50 < len(pix[0]):
                    ax2.axvline(x=vel50, color='green', linestyle='-', alpha=0.7)
            
            # Highlight saturated regions in column density plot
            if is_saturated:
                sat_indices = np.where((vel >= lmts[0]) & (vel <= lmts[1]) & (norm_flx <= sat_threshold))
                if len(sat_indices[0]) > 0:
                    ax2.scatter(vel[sat_indices], n[sat_indices], color='orange', 
                             marker='x', s=50, label='Saturated Pixels')
            # Add total column density annotation
            if not np.isnan(col):
                ax2.annotate(f'log N = {np.log10(col):.3f} ± {np.log10(col + colerr) - np.log10(col):.3f}', 
                           xy=(0.02, 0.85), xycoords='axes fraction',
                           bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
            ax2.set_ylabel(' AOD Column Density (cm$^{-2}$ km$^{-1}$ s)')
        # Plot cumulative EW
        else:
            # Plot cumulative equivalent width if f0 not provided
            vel_pix = vel[pix]
            sort_idx = np.argsort(vel_pix)
            cumulative_ew = np.cumsum(ew[sort_idx]) / np.sum(ew)
            ax2.plot(vel_pix[sort_idx], cumulative_ew, 'b-', label="Cumulative EW")
            
            # Mark 16, 50, and 84 percentiles
            try:
                v16 = np.interp(0.16, cumulative_ew, vel_pix[sort_idx])
                v50 = np.interp(0.50, cumulative_ew, vel_pix[sort_idx])
                v84 = np.interp(0.84, cumulative_ew, vel_pix[sort_idx])
                
                ax2.axvline(x=v16, color='green', linestyle='--', alpha=0.7, label=f'16% ({v16:.1f} km/s)')
                ax2.axvline(x=v50, color='green', linestyle='-', alpha=0.7, label=f'50% ({v50:.1f} km/s)')
                ax2.axvline(x=v84, color='green', linestyle='--', alpha=0.7, label=f'84% ({v84:.1f} km/s)')
                
                # Add dispersion annotation
                disp = (v84 - v16) / 2
                ax2.annotate(f'σ = {disp:.1f} km/s', 
                           xy=(0.02, 0.85), xycoords='axes fraction',
                           bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
            except:
                pass
                
            ax2.set_ylabel('Cumulative EW')
            ax2.set_ylim([0, 1.05])     
        


        ax2.set_xlabel('Velocity [km/s]')
        ax2.axhline(y=0, color='k', linestyle=':', alpha=0.5)
        ax2.legend(loc='upper right')
        

        plt.tight_layout()
        plt.show()

    return output