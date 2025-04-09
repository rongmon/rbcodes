import numpy as np
import matplotlib.pyplot as plt

def compute_EW(lam, flx, wrest, lmts, flx_err, plot=False, **kwargs):
    """
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
    """
    # Set default values for optional parameters
    verbose = kwargs.get('verbose', False)
    SNR = kwargs.get('SNR', False)
    zabs = kwargs.get('zabs', 0.0)  # Default redshift is 0
    sat_limit = kwargs.get('sat_limit', 0.10)  # Default saturation limit
    f0 = kwargs.get('f0', None)  # Default: None, optional for AOD calculation
    _binsize = kwargs.get('_binsize', 1)  # Default bin size for SNR calculation

    spl = 2.9979e5  # Speed of light in km/s
    vel = (lam - wrest * (1.0 + zabs)) * spl / (wrest * (1.0 + zabs))  # Velocity in km/s
    lambda_r = lam / (1. + zabs)

    norm_flx = flx  # Normalize flux (can add further normalization if needed)
    flx_err = flx_err  # Error on flux

    # Handle NaNs and values below saturation limit
    sq = np.isnan(norm_flx)
    norm_flx[sq] = flx_err[sq]  # Replace NaN flux with error values
    q = np.where(norm_flx <= sat_limit)
    norm_flx[q] = flx_err[q]  # Saturated flux set to error values
    q = np.where(norm_flx <= 0.0)
    norm_flx[q] = flx_err[q] + 0.01  # Avoid infinite optical depth

    # Compute the wavelength differences
    del_lam_j = np.diff(lambda_r)
    del_lam_j = np.append([del_lam_j[0]], del_lam_j)

    # Select pixels within velocity limits
    pix = np.where((vel >= lmts[0]) & (vel <= lmts[1]))
    Dj = 1. - norm_flx  # Normalized flux difference for EW computation

    # Equivalent Width Per Pixel (in Angstroms)
    ew = del_lam_j[pix] * Dj[pix]

    # Error on Equivalent Width
    sig_dj_sq = (flx_err) ** 2.
    err_ew = del_lam_j[pix] * np.sqrt(sig_dj_sq[pix])
    err_ew_tot = np.sqrt(np.sum(err_ew ** 2.))
    ew_tot = np.sum(ew)

    # Compute velocity centroid weighted by equivalent width
    ew50 = np.cumsum(ew) / np.max(np.cumsum(ew))
    vel50 = np.interp(0.5, ew50, vel[pix])
    vel16 = np.interp(0.16, ew50, vel[pix])
    vel_disp = np.abs(vel50 - vel16)
    vel50_err = vel_disp / np.sqrt(len(ew))

    if verbose:
        print(f"W_lambda = {ew_tot:.3f} +/- {err_ew_tot:.3f} Å over [{lmts[0]:.1f} to {lmts[1]:.1f}] km/s")

    # Output dictionary
    output = {
        "ew_tot": ew_tot,
        "err_ew_tot": err_ew_tot,
        "vel_disp": vel_disp,
        "vel50_err": vel50_err
    }

    if f0 is not None:
        # Compute Apparent Optical Depth (AOD) and column density
        Tau_a = np.log(1. / norm_flx)

        # Compute the column density as a function of velocity
        del_vel_j = np.diff(vel)
        del_vel_j = np.append([del_vel_j[0]], del_vel_j)
        nv = Tau_a / ((2.654e-15) * f0 * lambda_r)  # Column density per pixel
        n = nv * del_vel_j  # Column density per velocity bin

        tauerr = flx_err / norm_flx
        nerr = (tauerr / ((2.654e-15) * f0 * lambda_r)) * del_vel_j  # Error on column density

        col = np.sum(n[pix])  # Total column density
        colerr = np.sum((nerr[pix]) ** 2.) ** 0.5  # Error on column density

        if verbose:
            print(f"Direct N = {np.log10(col):.3f} +/- {np.log10(col + colerr) - np.log10(col):.3f} cm^-2")

        # Update output dictionary with column density and AOD
        output.update({
            "col": col,
            "colerr": colerr,
            "Tau_a": Tau_a,
            "med_vel": vel50
        })

    if SNR:
        # Compute SNR within the velocity range for the given slice
        from rbcodes.utils import compute_SNR_1d as c
        snr_result = c.estimate_snr(lambda_r, norm_flx, flx_err, binsize=_binsize, snr_range=[-1, -1], verbose=verbose, robust_median=True, sigma_clip_threshold=2.)
        output["SNR"] = snr_result['robust_snr']

    # If 'plot' is set, generate the plot
    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

        # Plot the normalized flux and error
        ax1.step(vel, norm_flx, label="Normalized Flux")
        ax1.step(vel, flx_err, color='r', label="Flux Error")
        ax1.set_xlim([-600, 600])
        ax1.set_ylim([-0.02, 1.8])
        ax1.plot([-2500, 2500], [0, 0], 'k:')
        ax1.plot([-2500, 2500], [1, 1], 'k:')
        ax1.plot([lmts[0], lmts[0]], [1.5, 1.5], 'r+', markersize=15)
        ax1.plot([lmts[1], lmts[1]], [1.5, 1.5], 'r+', markersize=15)
        if SNR:
            ax1.text(-600, 0.2, 'Median SNR: {:.1f}'.format(output["SNR"]))
        ax1.set_title(f'$W_{{rest}} = {ew_tot:.3f} \pm {err_ew_tot:.3f}$ Å')
        ax1.set_xlabel('Velocity [km/s]')
        ax1.legend()

        # Plot the column density as a function of velocity
        ax2.step(vel, n, label="Column Density")
        ax2.set_xlim([-600, 600])
        ax2.set_xlabel('Velocity [km/s]')
        ax2.plot([-2500, 2500], [0, 0], 'k:')
        ax2.legend()

        plt.tight_layout()
        plt.show()

    return output
