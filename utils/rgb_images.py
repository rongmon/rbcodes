import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

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
def reproject_and_crop(r_file, g_file, b_file,
                      target_grid='r', crop_to='r', custom_crop=None,
                      reproject_method='interp'):
    """
    Reproject three images to common grid and crop to specified coverage.
    
    Parameters
    ----------
    r_file, g_file, b_file : str
        Paths to FITS files for red, green, blue channels.
    target_grid : str
        Which image's pixel grid to use as target: 'r', 'g', 'b', or 'optimal'.
    crop_to : str  
        Which image's coverage to crop to: 'r', 'g', 'b', 'common', or 'union'.
    custom_crop : tuple, optional
        Custom crop region as (ra_min, ra_max, dec_min, dec_max) in degrees.
    reproject_method : str
        'interp' (fast) or 'exact' (precise) reprojection method.
    
    Returns
    -------
    r_reproj, g_reproj, b_reproj : 2D arrays
        Reprojected and cropped image data arrays.
    target_wcs : WCS
        WCS of the final reprojected grid.
    """
    from reproject import reproject_interp, reproject_exact
    from reproject.mosaicking import find_optimal_celestial_wcs
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    # Load all images
    files = {'r': r_file, 'g': g_file, 'b': b_file}
    data = {}
    wcs_dict = {}
    
    for channel, filename in files.items():
        hdu = fits.open(filename)
        raw_data = hdu[0].data
        raw_wcs = WCS(hdu[0].header)
        
        # Handle 3D data/WCS - extract 2D celestial part
        if raw_data.ndim == 3:
            # Take middle slice for 3D cube (or sum/average if preferred)
            print(f"Warning: {channel}_file has 3D data, taking slice {raw_data.shape[0]//2}")
            data[channel] = raw_data[raw_data.shape[0]//2, :, :]
            # Extract 2D celestial WCS
            wcs_dict[channel] = raw_wcs.celestial
        elif raw_data.ndim == 2:
            data[channel] = raw_data
            # Ensure we have 2D celestial WCS even if header has extra axes
            wcs_dict[channel] = raw_wcs.celestial
        else:
            raise ValueError(f"Unsupported data dimensions: {raw_data.ndim}D in {filename}")
        
        hdu.close()
    
    # Determine target WCS
    if target_grid == 'optimal':
        # Create optimal WCS covering all images
        data_wcs_pairs = [(data[ch], wcs_dict[ch]) for ch in ['r', 'g', 'b']]
        target_wcs, target_shape = find_optimal_celestial_wcs(data_wcs_pairs)
        print(f"Using optimal WCS with shape {target_shape}")
    else:
        # Use specified image's WCS
        target_wcs = wcs_dict[target_grid]
        target_shape = data[target_grid].shape
        print(f"Using {target_grid}_file WCS as target grid")
    
    # Reproject all images to target grid
    reproject_func = reproject_exact if reproject_method == 'exact' else reproject_interp
    reprojected = {}
    footprints = {}
    
    for channel in ['r', 'g', 'b']:
        if channel == target_grid and target_grid != 'optimal':
            # No need to reproject if it's already the target
            reprojected[channel] = data[channel]
            footprints[channel] = np.ones_like(data[channel], dtype=bool)
        else:
            print(f"Reprojecting {channel}_file to target grid...")
            reprojected[channel], footprints[channel] = reproject_func(
                (data[channel], wcs_dict[channel]), 
                target_wcs, 
                shape_out=target_shape
            )
    
    # Handle cropping
    if custom_crop is not None:
        # Custom RA/Dec crop box
        ra_min, ra_max, dec_min, dec_max = custom_crop
        print(f"Applying custom crop: RA({ra_min:.3f}, {ra_max:.3f}), Dec({dec_min:.3f}, {dec_max:.3f})")
        
        # Create coordinate arrays
        ny, nx = target_shape
        x_pix, y_pix = np.meshgrid(np.arange(nx), np.arange(ny))
        coords = target_wcs.pixel_to_world(x_pix, y_pix)
        
        # Create crop mask
        crop_mask = ((coords.ra.deg >= ra_min) & (coords.ra.deg <= ra_max) & 
                    (coords.dec.deg >= dec_min) & (coords.dec.deg <= dec_max))
        
        # Find bounding box
        y_indices, x_indices = np.where(crop_mask)
        if len(y_indices) == 0:
            raise ValueError("Custom crop region contains no pixels!")
        
        y_min, y_max = y_indices.min(), y_indices.max() + 1
        x_min, x_max = x_indices.min(), x_indices.max() + 1
        
    elif crop_to in ['r', 'g', 'b']:
        # Crop to specific image's footprint
        print(f"Cropping to {crop_to}_file coverage...")
        footprint = footprints[crop_to]
        
        # Ensure boolean types for logical operations
        footprint_bool = footprint.astype(bool)
        nan_mask = ~np.isnan(reprojected[crop_to])
        
        # Find bounding box of valid data
        valid_pixels = np.where(footprint_bool & nan_mask)
        if len(valid_pixels[0]) == 0:
            raise ValueError(f"No valid data found in {crop_to}_file!")
        
        y_min, y_max = valid_pixels[0].min(), valid_pixels[0].max() + 1
        x_min, x_max = valid_pixels[1].min(), valid_pixels[1].max() + 1
        
    elif crop_to == 'common':
        # Crop to intersection of all footprints
        print("Cropping to common coverage area...")
        
        # Ensure all boolean types
        footprint_r = footprints['r'].astype(bool)
        footprint_g = footprints['g'].astype(bool) 
        footprint_b = footprints['b'].astype(bool)
        
        nan_mask_r = ~np.isnan(reprojected['r'])
        nan_mask_g = ~np.isnan(reprojected['g'])
        nan_mask_b = ~np.isnan(reprojected['b'])
        
        common_footprint = (footprint_r & footprint_g & footprint_b &
                           nan_mask_r & nan_mask_g & nan_mask_b)
        
        valid_pixels = np.where(common_footprint)
        if len(valid_pixels[0]) == 0:
            raise ValueError("No common coverage area found!")
        
        y_min, y_max = valid_pixels[0].min(), valid_pixels[0].max() + 1
        x_min, x_max = valid_pixels[1].min(), valid_pixels[1].max() + 1
        
    elif crop_to == 'union':
        # Keep all coverage (no cropping)
        print("Keeping full union coverage...")
        y_min, y_max = 0, target_shape[0]
        x_min, x_max = 0, target_shape[1]
        
    else:
        raise ValueError("crop_to must be 'r', 'g', 'b', 'common', 'union', or use custom_crop")
    
    # Apply cropping
    r_cropped = reprojected['r'][y_min:y_max, x_min:x_max]
    g_cropped = reprojected['g'][y_min:y_max, x_min:x_max] 
    b_cropped = reprojected['b'][y_min:y_max, x_min:x_max]
    
    # Update WCS for cropped region
    target_wcs_cropped = target_wcs.deepcopy()
    target_wcs_cropped.wcs.crpix[0] -= x_min
    target_wcs_cropped.wcs.crpix[1] -= y_min
    
    print(f"Final cropped shape: {r_cropped.shape}")
    print(f"Pixel scale: {target_wcs.proj_plane_pixel_scales()[0]*3600:.3f} arcsec/pixel")
    
    return r_cropped, g_cropped, b_cropped, target_wcs_cropped


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