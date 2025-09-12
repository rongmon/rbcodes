import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def stiff_rgb(r, g, b,
              min_percent=0.01, max_percent=99.99,
              gamma=2.2,
              color_balance=True,
              savefile=None, dpi=300):
    """
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
    """
    # Stack to handle percentile clipping together
    all_data = np.stack([r, g, b], axis=0)
    vmin = np.percentile(all_data, min_percent)
    vmax = np.percentile(all_data, max_percent)
    
    def stretch(img):
        img = np.clip(img, vmin, vmax)
        img = (img - vmin) / (vmax - vmin)
        img = img ** (1 / gamma)
        return img
    
    r_s = stretch(r)
    g_s = stretch(g)
    b_s = stretch(b)
    
    # Optional color balance (normalize each channel)
    if color_balance:
        r_s /= np.median(r_s)
        g_s /= np.median(g_s)
        b_s /= np.median(b_s)
    
    # Re-normalize to [0, 1] after color balance
    rgb_stack = np.stack([r_s, g_s, b_s], axis=-1)
    rgb_stack /= np.max(rgb_stack)
    
    if savefile:
        if figsize is not None:
            # Use matplotlib figure for custom size control
            fig, ax = plt.subplots(figsize=figsize)
            ax.imshow(rgb_stack, origin='lower')
            ax.axis('off')
            plt.tight_layout()
            plt.savefig(savefile, dpi=dpi, bbox_inches='tight', pad_inches=0)
            plt.close()
        else:
            # Use simple imsave (default behavior)
            plt.imsave(savefile, rgb_stack, dpi=dpi)
    
    return rgb_stack


def lupton_rgb(r, g, b,
               Q=8, stretch=0.5,
               minimum=0,
               color_balance=True,
               savefile=None, dpi=300):
    """
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
    """
    # Subtract minimum (background) from all channels
    r = np.array(r, dtype=float) - minimum
    g = np.array(g, dtype=float) - minimum  
    b = np.array(b, dtype=float) - minimum
    
    # Set negative values to small positive number
    r = np.maximum(r, 1e-10)
    g = np.maximum(g, 1e-10)
    b = np.maximum(b, 1e-10)
    
    # Color balance: normalize channels so they contribute equally in bright regions
    if color_balance:
        # Use 99th percentile for normalization to avoid outliers
        r_norm = np.percentile(r, 99)
        g_norm = np.percentile(g, 99)
        b_norm = np.percentile(b, 99)
        
        # Avoid division by zero
        if r_norm <= 0: r_norm = 1
        if g_norm <= 0: g_norm = 1  
        if b_norm <= 0: b_norm = 1
        
        r = r / r_norm
        g = g / g_norm
        b = b / b_norm
    
    # Calculate intensity (luminance)
    I = (r + g + b) / 3.0
    
    # Apply asinh stretch
    # The key insight: stretch factor alpha depends on intensity
    alpha = stretch * I / Q
    
    # Avoid division by zero in the stretch
    I_safe = np.maximum(I, 1e-10)
    
    # Lupton asinh stretch formula
    f = np.arcsinh(alpha) / alpha
    
    # Handle alpha=0 case (set f=1 for very low intensity regions)
    f = np.where(alpha <= 0, 1, f)
    
    # Apply stretch to each channel
    r_stretched = f * r
    g_stretched = f * g
    b_stretched = f * b
    
    # Normalize to [0, 1] range
    rgb_stack = np.stack([r_stretched, g_stretched, b_stretched], axis=-1)
    rgb_max = np.max(rgb_stack)
    if rgb_max > 0:
        rgb_stack = rgb_stack / rgb_max
    
    # Ensure values are in [0, 1]
    rgb_stack = np.clip(rgb_stack, 0, 1)
    
    if savefile:
        plt.imsave(savefile, rgb_stack, dpi=dpi)
    
    return rgb_stack


# Example usage functions
def load_and_create_rgb(r_file, g_file, b_file, method='stiff', **kwargs):
    """
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
    """
    # Load FITS files
    r_data = fits.getdata(r_file)
    g_data = fits.getdata(g_file) 
    b_data = fits.getdata(b_file)
    
    if method.lower() == 'stiff':
        return stiff_rgb(r_data, g_data, b_data, **kwargs)
    elif method.lower() == 'lupton':
        return lupton_rgb(r_data, g_data, b_data, **kwargs)
    else:
        raise ValueError("Method must be 'stiff' or 'lupton'")


def manual_rgb(r, g, b,
               r_min=None, r_max=None, r_scale=1.0,
               g_min=None, g_max=None, g_scale=1.0,
               b_min=None, b_max=None, b_scale=1.0,
               r_percentiles=(1, 95), g_percentiles=(1, 95), b_percentiles=(1, 95),
               use_percentiles=True,
               stretch='linear',
               savefile=None, dpi=300):
    """
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
    """
    
    def get_limits(data, percentiles, manual_min, manual_max):
        """Get min/max values either from percentiles or manual input."""
        if use_percentiles or (manual_min is None or manual_max is None):
            vmin = np.percentile(data, percentiles[0])
            vmax = np.percentile(data, percentiles[1])
        else:
            vmin = manual_min
            vmax = manual_max
        return vmin, vmax
    
    def apply_stretch(data, stretch_type):
        """Apply different stretch functions."""
        if stretch_type == 'linear':
            return data
        elif stretch_type == 'sqrt':
            return np.sqrt(np.maximum(data, 0))
        elif stretch_type == 'log':
            return np.log10(np.maximum(data, 1e-10))
        elif stretch_type == 'asinh':
            return np.arcsinh(data)
        elif stretch_type == 'power':
            return np.power(np.maximum(data, 0), 0.5)
        else:
            return data
    
    # Get scaling limits for each channel
    r_vmin, r_vmax = get_limits(r, r_percentiles, r_min, r_max)
    g_vmin, g_vmax = get_limits(g, g_percentiles, g_min, g_max)
    b_vmin, b_vmax = get_limits(b, b_percentiles, b_min, b_max)
    
    # Apply scaling to each channel independently
    r_scaled = r * r_scale
    g_scaled = g * g_scale  
    b_scaled = b * b_scale
    
    # Clip and normalize each channel
    r_norm = np.clip((r_scaled - r_vmin) / (r_vmax - r_vmin), 0, None)
    g_norm = np.clip((g_scaled - g_vmin) / (g_vmax - g_vmin), 0, None)
    b_norm = np.clip((b_scaled - b_vmin) / (b_vmax - b_vmin), 0, None)
    
    # Apply stretch function
    r_stretched = apply_stretch(r_norm, stretch)
    g_stretched = apply_stretch(g_norm, stretch)
    b_stretched = apply_stretch(b_norm, stretch)
    
    # Final normalization to [0, 1]
    def final_norm(channel):
        if np.max(channel) > 0:
            return np.clip(channel / np.max(channel), 0, 1)
        else:
            return channel
    
    r_final = final_norm(r_stretched)
    g_final = final_norm(g_stretched)
    b_final = final_norm(b_stretched)
    
    # Create RGB stack
    rgb_stack = np.stack([r_final, g_final, b_final], axis=-1)
    
    # Store scaling info for reference
    scaling_info = {
        'red': {'min': r_vmin, 'max': r_vmax, 'scale': r_scale, 'percentiles': r_percentiles},
        'green': {'min': g_vmin, 'max': g_vmax, 'scale': g_scale, 'percentiles': g_percentiles}, 
        'blue': {'min': b_vmin, 'max': b_vmax, 'scale': b_scale, 'percentiles': b_percentiles},
        'stretch': stretch
    }
    
    if savefile:
        plt.imsave(savefile, rgb_stack, dpi=dpi)
    
    return rgb_stack, scaling_info


def compare_methods(r, g, b, savefile_prefix=None):
    """
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
    """
    results = {}
    
    # Create with STIFF method
    stiff_save = f"{savefile_prefix}_stiff.png" if savefile_prefix else None
    results['stiff'] = stiff_rgb(r, g, b, savefile=stiff_save)
    
    # Create with Lupton method  
    lupton_save = f"{savefile_prefix}_lupton.png" if savefile_prefix else None
    results['lupton'] = lupton_rgb(r, g, b, savefile=lupton_save)
    
    return results