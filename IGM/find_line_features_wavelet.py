import numpy as np
from scipy.signal import find_peaks, cwt, ricker
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

import numpy as np
from scipy.signal import find_peaks, cwt, ricker
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit


def gaussian_model(x, amp, center, width, offset=0):
    return amp * np.exp(-(x - center)**2 / (2 * width**2)) + offset


def find_spectral_lines_wavelet(
    wavelength, flux, line_type='both', min_snr=3.0,
    min_scale=None, max_scale=None, num_scales=20,
    fit_lines=True, cluster_tolerance=None, edge_buffer=10
):
    wavelength = np.array(wavelength)
    flux = np.array(flux)
    if len(wavelength) != len(flux):
        raise ValueError("wavelength and flux must be the same length")

    mean_step = np.mean(np.diff(wavelength))
    spectral_length = len(wavelength)

    if min_scale is None:
        min_scale = 1
    if max_scale is None:
        max_scale = max(30, int(spectral_length / 20))

    scales = np.logspace(np.log10(min_scale), np.log10(max_scale), num_scales)

    if cluster_tolerance is None:
        cluster_tolerance = max(3, int(min_scale))

    smoothed_flux = gaussian_filter1d(flux, 7)
    noise = flux - smoothed_flux
    noise_std = np.std(noise)

    cwt_matrix = cwt(flux, ricker, scales)
    abs_cwt = np.abs(cwt_matrix)
    best_scale_idx = np.argmax(abs_cwt, axis=0)
    best_scale_response = np.array([
        cwt_matrix[best_scale_idx[i], i] for i in range(len(flux))
    ])

    edge_mask = np.ones_like(wavelength, dtype=bool)
    if edge_buffer > 0:
        edge_mask[:edge_buffer] = False
        edge_mask[-edge_buffer:] = False

    detected_lines = []

    # Emission lines
    if line_type in ['both', 'emission']:
        peaks, _ = find_peaks(
            best_scale_response,
            height=noise_std * min_snr,
            distance=cluster_tolerance
        )
        peaks = peaks[edge_mask[peaks]]
        for p in peaks:
            detected_lines.append({
                'center': wavelength[p],
                'type': 'emission',
                'snr': abs(best_scale_response[p]) / noise_std,
                'index': p,
                'width': scales[best_scale_idx[p]],
                'strength': best_scale_response[p]
            })

    # Absorption lines
    if line_type in ['both', 'absorption']:
        troughs, _ = find_peaks(
            -best_scale_response,
            height=noise_std * min_snr,
            distance=cluster_tolerance
        )
        troughs = troughs[edge_mask[troughs]]
        for t in troughs:
            detected_lines.append({
                'center': wavelength[t],
                'type': 'absorption',
                'snr': abs(best_scale_response[t]) / noise_std,
                'index': t,
                'width': scales[best_scale_idx[t]],
                'strength': best_scale_response[t]
            })

    detected_lines.sort(key=lambda x: x['center'])

    if not detected_lines:
        return {
            'line_centers': np.array([]),
            'line_types': np.array([]),
            'snr': np.array([]),
            'line_indices': np.array([]),
            'line_widths': np.array([]),
            'fit_params': []
        }


    fit_params = []
    
    for i, line in enumerate(detected_lines):
        idx = line['index']
        width_points = int(max(line['width'] * 3, 10))
        left_idx = max(0, idx - width_points)
        right_idx = min(len(wavelength), idx + width_points + 1)
    
        x_fit = wavelength[left_idx:right_idx]
        y_fit = flux[left_idx:right_idx]
    
        if len(x_fit) < 5:
            fit_params.append(None)
            continue
    
        fit_success = False  # Flag to track fitting outcome
    
        if fit_lines:
            # Initial guess
            offset_init = np.median(y_fit)
            amp_init = (max(y_fit) - offset_init) if line['type'] == 'emission' else (min(y_fit) - offset_init)
            center_init = line['center']
            width_init = abs(mean_step * line['width'] / 2)
            p0 = [amp_init, center_init, width_init, offset_init]
    
            # Parameter bounds
            half_range = (x_fit[-1] - x_fit[0]) / 2
            if line['type'] == 'emission':
                bounds = ([0, x_fit[0], 0, 0],
                          [np.inf, x_fit[-1], half_range, np.inf])
            else:
                bounds = ([-np.inf, x_fit[0], 0, 0],
                          [0, x_fit[-1], half_range, np.inf])
    
            try:
                popt, _ = curve_fit(gaussian_model, x_fit, y_fit, p0=p0, bounds=bounds)
                fit_params.append({
                    'amplitude': popt[0],
                    'center': popt[1],
                    'width': abs(popt[2]),
                    'offset': popt[3],
                    'fit_success': True
                })
                detected_lines[i]['width'] = abs(popt[2]) / mean_step  # Update width estimate
                fit_success = True
            except Exception as e:
                print(f"Fitting failed for line at {line['center']:.2f}: {e}")
    
        if not fit_success:
            # Fallback: estimate width from percentiles
            width_est = estimate_width_percentile(wavelength, flux, idx, window=30)
            if width_est is not None:
                line['width'] = width_est / mean_step  # Update width
                fit_params.append({
                    'amplitude': None,
                    'center': line['center'],
                    'width': width_est,
                    'offset': None,
                    'fit_success': False
                })
            else:
                fit_params.append(None)
    
    
    

    result = {
        'line_centers': np.array([d['center'] for d in detected_lines]),
        'line_types': np.array([d['type'] for d in detected_lines]),
        'snr': np.array([d['snr'] for d in detected_lines]),
        'line_indices': np.array([d['index'] for d in detected_lines]),
        'line_widths': np.array([d['width'] for d in detected_lines]),
        'fit_params': fit_params
    }

    return result

def estimate_width_percentile(wave, flux, center_idx, window=30, lower=1, upper=99):
    """
    Estimate line width using percentile method around a detected line.
    Returns width in wavelength units.
    """
    left = max(0, center_idx - window)
    right = min(len(wave), center_idx + window + 1)
    x = wave[left:right]
    y = flux[left:right]
    
    if len(x) < 5:
        return None

    # Estimate baseline from local median
    baseline = np.median(y)
    
    # Determine whether line is absorption or emission
    if flux[center_idx] < baseline:
        profile = baseline - y  # absorption
    else:
        profile = y - baseline  # emission

    # Mask out negative (non-line) regions
    profile = np.clip(profile, 0, None)
    
    # Abort if no structure
    total = np.sum(profile)
    if total == 0:
        return None

    cumsum = np.cumsum(profile)
    cumsum /= cumsum[-1]  # normalize to 1

    # Function to find closest percentile index
    def closest_percentile_index(target):
        matches = np.where(np.abs(cumsum - target / 100) < 0.01)[0]
        if len(matches) == 0:
            return np.searchsorted(cumsum, target / 100)
        else:
            return matches[np.argmin(np.abs(matches - (center_idx - left)))]

    low_idx = closest_percentile_index(lower)
    high_idx = closest_percentile_index(upper)

    if 0 <= low_idx < len(x) and 0 <= high_idx < len(x) and high_idx > low_idx:
        return abs(x[high_idx] - x[low_idx])
    else:
        return None


def create_feature_mask(wavelength, line_results, width_scale_factor=3.0):
    """
    Create a mask for continuum fitting based on detected spectral lines.
    
    Parameters
    ----------
    wavelength : array-like
        Wavelength array.
    line_results : dict
        Results from find_spectral_lines_wavelet.
    width_scale_factor : float, optional
        Factor to scale the line width for masking. Default is 3.0.
        
    Returns
    -------
    mask : ndarray
        Boolean mask where False indicates regions to mask out (spectral lines).
    """
    # Create mask: True = keep for continuum, False = mask out lines
    mask = np.ones_like(wavelength, dtype=bool)
    
    if len(line_results['line_centers']) == 0:
        return mask
    
    mean_step = np.mean(np.diff(wavelength))
    
    for i, (center, width, snr) in enumerate(zip(
            line_results['line_centers'], 
            line_results['line_widths'], 
            line_results['snr'])):
        
        # Get fitted width if available, otherwise use wavelet width
        fit_width = None
        if line_results['fit_params'][i] is not None and line_results['fit_params'][i]['fit_success']:
            fit_width = line_results['fit_params'][i]['width']
        
        if fit_width is not None:
            # Convert wavelength width to index units
            width_idx = fit_width / mean_step
        else:
            width_idx = width
        
        # Scale mask width based on SNR - stronger lines get wider masks
        snr_factor = 1.0 + 0.2 * np.log10(max(1.0, snr / 3.0))
        scaled_width = width_idx * width_scale_factor * snr_factor
        
        # Find nearest index to center
        idx = np.argmin(np.abs(wavelength - center))
        
        # Calculate mask boundaries
        half_width = int(np.ceil(scaled_width / 2.0))
        left_idx = max(0, idx - half_width)
        right_idx = min(len(wavelength), idx + half_width + 1)
        
        # Apply mask
        mask[left_idx:right_idx] = False
    
    return mask

