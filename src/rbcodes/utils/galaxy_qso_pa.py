from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from photutils.isophote import Ellipse, EllipseGeometry
from photutils.aperture import EllipticalAperture
from photutils.centroids import centroid_com
from astropy.nddata import Cutout2D
from scipy.ndimage import gaussian_filter

def compute_galaxy_qso_pa(fits_file, ra_gal, dec_gal, ra_qso, dec_qso,
                          cutout_size=(100, 100), initial_sma=10, initial_eps=0.3,
                          moments_sigma=2.0, plot=False, save_plot=None, 
                          bootstrap_n=100):
    """
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
    """
    
    # Load data and create cutout
    hdu = fits.open(fits_file)[0]
    data = hdu.data
    wcs = WCS(hdu.header)
    
    gal_coord = SkyCoord(ra_gal, dec_gal, unit='deg')
    qso_coord = SkyCoord(ra_qso, dec_qso, unit='deg')
    gal_pix = wcs.world_to_pixel(gal_coord)
    qso_pix = wcs.world_to_pixel(qso_coord)
    
    cutout = Cutout2D(data, position=gal_pix, size=cutout_size, wcs=wcs)
    cut_data = cutout.data
    
    # Find galaxy centroid
    try:
        x_centroid, y_centroid = centroid_com(cut_data)
        print(f"Galaxy centroid: ({x_centroid:.2f}, {y_centroid:.2f})")
    except:
        x_centroid, y_centroid = cutout_size[1]/2, cutout_size[0]/2
        print("Using geometric center")
    
    # Fit galaxy shape: try isophote first, fallback to moments
    try:
        method_used, shape_results = _fit_isophote(cut_data, initial_sma, initial_eps, 
                                                  x_centroid, y_centroid)
        print(f"Isophote fitting successful")
    except Exception as e:
        print(f"Isophote fitting failed: {e}")
        print("Switching to moments method")
        method_used, shape_results = _fit_moments(cut_data, moments_sigma, bootstrap_n)
    
    # Convert pixel PA to celestial PA and propagate errors
    pa_major_celestial = _pixel_to_celestial_pa(shape_results['pa_pixel'], cutout, cutout_size)
    
    # Estimate uncertainty in celestial PA from pixel PA uncertainty
    if 'pa_pixel_err' in shape_results and shape_results['pa_pixel_err'] > 0:
        # For small errors, celestial error ≈ pixel error (WCS projection is approximately linear)
        pa_major_err = shape_results['pa_pixel_err']
    else:
        pa_major_err = 0.0
    
    # Calculate QSO direction in celestial coordinates (no uncertainty - coordinates are exact)
    pa_to_qso_celestial = gal_coord.position_angle(qso_coord).to(u.deg).value
    
    # Compute alignment angle (0-90°) and propagate uncertainty
    delta_pa = _compute_alignment_angle(pa_major_celestial, pa_to_qso_celestial)
    delta_pa_err = pa_major_err  # Uncertainty propagates directly for alignment angle
    alignment_type = "major axis aligned" if delta_pa <= 45 else "minor axis aligned"
    
    # Build results
    results = {
        'pa_major': pa_major_celestial,
        'pa_major_err': pa_major_err,
        'pa_to_qso': pa_to_qso_celestial,
        'delta_pa': delta_pa,
        'delta_pa_err': delta_pa_err,
        'alignment_type': alignment_type,
        'method_used': method_used,
        'galaxy_pix': gal_pix,
        'qso_pix': qso_pix,
        'cutout': cutout,
        'centroid': (x_centroid, y_centroid),
        'shape_results': shape_results
    }
    
    # Print results with uncertainties
    print(f"Method: {method_used}")
    if pa_major_err > 0:
        print(f"Galaxy PA: {pa_major_celestial:.1f} ± {pa_major_err:.1f}°, QSO PA: {pa_to_qso_celestial:.1f}°")
        print(f"Alignment: {delta_pa:.1f} ± {delta_pa_err:.1f}° ({alignment_type})")
    else:
        print(f"Galaxy PA: {pa_major_celestial:.1f}°, QSO PA: {pa_to_qso_celestial:.1f}°")
        print(f"Alignment: {delta_pa:.1f}° ({alignment_type})")
        print("No uncertainty estimate available")
    
    if plot:
        _create_plot(results, cut_data, cutout_size, save_plot)
    
    return results

def _fit_isophote(data, initial_sma, initial_eps, x_centroid, y_centroid):
    """Try to fit isophote with multiple parameter attempts."""
    
    param_sets = [
        (initial_sma, initial_eps),
        (initial_sma * 0.5, 0.1),
        (initial_sma * 1.5, 0.5),
        (initial_sma * 2.0, 0.3)
    ]
    
    for i, (sma, eps) in enumerate(param_sets):
        try:
            geometry = EllipseGeometry(x0=x_centroid, y0=y_centroid, sma=sma, eps=eps, pa=0)
            ellipse = Ellipse(data, geometry)
            isolist = ellipse.fit_image()
            
            if len(isolist) > 0:
                print(f"Isophote successful on attempt {i+1}")
                
                # Get best isophote and extract uncertainties
                valid_isos = [iso for iso in isolist if hasattr(iso, 'valid') and iso.valid]
                best_iso = valid_isos[-1] if valid_isos else isolist[-1]
                
                # Extract PA uncertainty from isophote fitting
                pa_pixel_err = 0.0
                if hasattr(best_iso, 'pa_err') and best_iso.pa_err is not None:
                    pa_pixel_err = np.rad2deg(abs(best_iso.pa_err))  # Convert to degrees
                
                return 'isophote', {
                    'pa_pixel': np.rad2deg(best_iso.pa),
                    'pa_pixel_err': pa_pixel_err,
                    'ellipticity': best_iso.eps,
                    'center': (best_iso.x0, best_iso.y0),
                    'sma': best_iso.sma,
                    'isolist': isolist
                }
        except:
            continue
    
    raise ValueError("All isophote attempts failed")

def _fit_moments(data, sigma, bootstrap_n=100):
    """Fit using second moments method with bootstrap uncertainty estimation."""
    
    # Apply Gaussian smoothing
    if sigma > 0:
        data_smooth = gaussian_filter(data, sigma=sigma)
    else:
        data_smooth = data
    
    # Main moments calculation
    pa_degrees, ellipticity, x_center, y_center, r_eff = _calculate_moments(data_smooth)
    
    # Bootstrap uncertainty estimation
    if bootstrap_n > 0:
        pa_bootstrap = []
        ny, nx = data_smooth.shape
        
        # Create pixel coordinates for resampling
        y_coords, x_coords = np.indices(data_smooth.shape)
        pixel_coords = np.column_stack([x_coords.flatten(), y_coords.flatten()])
        pixel_values = data_smooth.flatten()
        
        # Only use positive values for bootstrap
        positive_mask = pixel_values > 0
        if np.sum(positive_mask) < 10:  # Need minimum pixels
            print("Warning: Too few positive pixels for bootstrap")
            pa_pixel_err = 0.0
        else:
            pixel_coords_pos = pixel_coords[positive_mask]
            pixel_values_pos = pixel_values[positive_mask]
            
            print(f"Bootstrap uncertainty estimation with {bootstrap_n} samples...")
            
            for i in range(bootstrap_n):
                # Resample pixels with replacement, weighted by intensity
                weights = pixel_values_pos / np.sum(pixel_values_pos)
                n_pixels = len(pixel_values_pos)
                
                # Bootstrap sample
                sample_indices = np.random.choice(n_pixels, size=n_pixels, replace=True, p=weights)
                
                # Create bootstrap image
                bootstrap_image = np.zeros_like(data_smooth)
                for idx in sample_indices:
                    x_coord, y_coord = pixel_coords_pos[idx].astype(int)
                    if 0 <= x_coord < nx and 0 <= y_coord < ny:
                        bootstrap_image[y_coord, x_coord] += pixel_values_pos[idx]
                
                # Calculate PA for bootstrap sample
                try:
                    pa_boot, _, _, _, _ = _calculate_moments(bootstrap_image)
                    pa_bootstrap.append(pa_boot)
                except:
                    continue
            
            # Calculate uncertainty from bootstrap distribution
            if len(pa_bootstrap) > 10:  # Need minimum successful samples
                pa_bootstrap = np.array(pa_bootstrap)
                # Handle angle wrapping for uncertainty calculation
                pa_bootstrap = _unwrap_angles(pa_bootstrap, pa_degrees)
                pa_pixel_err = np.std(pa_bootstrap)
                print(f"Bootstrap PA uncertainty: ±{pa_pixel_err:.1f}°")
            else:
                pa_pixel_err = 0.0
                print("Warning: Bootstrap failed, no uncertainty estimate")
    else:
        pa_pixel_err = 0.0
    
    print(f"Moments: PA={pa_degrees:.1f} ± {pa_pixel_err:.1f}°, ellipticity={ellipticity:.3f}")
    
    return 'moments', {
        'pa_pixel': pa_degrees,
        'pa_pixel_err': pa_pixel_err,
        'ellipticity': ellipticity,
        'center': (x_center, y_center),
        'r_eff': r_eff
    }

def _calculate_moments(data_smooth):
    """Calculate second moments from image data."""
    
    # Calculate center of mass
    total_flux = np.sum(data_smooth)
    if total_flux <= 0:
        raise ValueError("No positive flux found")
    
    y_indices, x_indices = np.indices(data_smooth.shape)
    x_center = np.sum(x_indices * data_smooth) / total_flux
    y_center = np.sum(y_indices * data_smooth) / total_flux
    
    # Calculate second moments
    dx = x_indices - x_center
    dy = y_indices - y_center
    
    Mxx = np.sum(dx * dx * data_smooth) / total_flux
    Myy = np.sum(dy * dy * data_smooth) / total_flux
    Mxy = np.sum(dx * dy * data_smooth) / total_flux
    
    # Eigenvalue analysis
    M = np.array([[Mxx, Mxy], [Mxy, Myy]])
    eigenvals, eigenvecs = np.linalg.eigh(M)
    
    # Sort by eigenvalue (largest first)
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]
    
    # Calculate PA and ellipticity
    major_axis_vec = eigenvecs[:, 0]
    pa_radians = np.arctan2(major_axis_vec[0], major_axis_vec[1])
    pa_degrees = np.degrees(pa_radians) % 360
    
    ellipticity = 1 - np.sqrt(eigenvals[1] / eigenvals[0]) if eigenvals[0] > 0 else 0
    r_eff = np.sqrt(np.sqrt(eigenvals[0] * eigenvals[1]))
    
    return pa_degrees, ellipticity, x_center, y_center, r_eff

def _unwrap_angles(angles, reference_angle):
    """Handle angle wrapping for uncertainty calculations."""
    
    # Unwrap angles relative to reference
    unwrapped = []
    for angle in angles:
        diff = angle - reference_angle
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        unwrapped.append(reference_angle + diff)
    
    return np.array(unwrapped)

def _pixel_to_celestial_pa(pa_pixel, cutout, cutout_size):
    """Convert pixel PA to celestial PA using WCS."""
    
    center_x, center_y = cutout_size[1]/2, cutout_size[0]/2
    pa_rad = np.deg2rad(pa_pixel)
    
    # Create two points along major axis
    dx = 10 * np.cos(pa_rad)
    dy = 10 * np.sin(pa_rad)
    
    point1 = (center_x, center_y)
    point2 = (center_x + dx, center_y + dy)
    
    # Convert to world coordinates
    cutout_wcs = cutout.wcs
    point1_world = cutout_wcs.pixel_to_world(point1[0], point1[1])
    point2_world = cutout_wcs.pixel_to_world(point2[0], point2[1])
    
    return point1_world.position_angle(point2_world).to(u.deg).value

def _compute_alignment_angle(pa_major, pa_to_qso):
    """Compute alignment angle (0-90°)."""
    
    delta_pa = abs(pa_to_qso - pa_major)
    if delta_pa > 180:
        delta_pa = 360 - delta_pa
    if delta_pa > 90:
        delta_pa = 180 - delta_pa
    return delta_pa

def _create_plot(results, cut_data, cutout_size, save_plot):
    """Create single-panel diagnostic plot with correct compass."""
    
    fig = plt.figure(figsize=(10, 10))
    
    # Single Panel: Sky view with celestial coordinates and shape analysis
    cutout = results['cutout']
    ax = plt.subplot(111, projection=cutout.wcs)
    im = ax.imshow(cutout.data, origin='lower', cmap='viridis')
    plt.colorbar(im, ax=ax, shrink=0.8)
    
    shape = results['shape_results']
    method = results['method_used']
    
    # Manual ellipse drawing
    pa_pixel = shape['pa_pixel']
    ellipticity = shape['ellipticity']
    center = shape['center']
    
    # Estimate semi-major axis size
    if method == 'moments' and 'r_eff' in shape:
        sma = shape['r_eff'] * 2
    elif method == 'isophote' and 'sma' in shape:
        sma = shape['sma']
    else:
        sma = 25  # Default size
    
    # Draw ellipse
    pa_rad = np.deg2rad(pa_pixel)
    semi_major = sma
    semi_minor = sma * (1 - ellipticity)
    
    ellipse = EllipticalAperture(center, semi_major, semi_minor, theta=pa_rad)
    ellipse.plot(ax=ax, color='red', linewidth=2, label=f'{method.title()} ellipse', 
                 transform=ax.get_transform('pixel'))
    
    # Draw major axis line in yellow
    dx = semi_major * np.cos(pa_rad)
    dy = semi_major * np.sin(pa_rad)
    ax.plot([center[0] - dx, center[0] + dx], 
             [center[1] - dy, center[1] + dy],
             'y--', linewidth=3, alpha=0.9, label='Major axis', 
             transform=ax.get_transform('pixel'))
    
    # Mark galaxy center
    ax.plot(center[0], center[1], 'r+', markersize=12, markeredgewidth=3, label='Galaxy center',
            transform=ax.get_transform('pixel'))
    
    # Galaxy and QSO positions
    gal_pix_orig = results['galaxy_pix']
    qso_pix_orig = results['qso_pix']
    
    # Use actual fitted galaxy center (not geometric center of cutout)
    fitted_center = shape['center']  # This is the real galaxy center from fitting
    
    dx_qso = qso_pix_orig[0] - gal_pix_orig[0]
    dy_qso = qso_pix_orig[1] - gal_pix_orig[1]
    qso_pix_cutout = (fitted_center[0] + dx_qso, fitted_center[1] + dy_qso)
    
    # Plot galaxy position at fitted center
    ax.plot(fitted_center[0], fitted_center[1], 'r+', markersize=15, 
             markeredgewidth=3, label='Galaxy', transform=ax.get_transform('pixel'))
    
    # Plot QSO or direction arrow
    if (0 <= qso_pix_cutout[0] < cutout_size[1] and 0 <= qso_pix_cutout[1] < cutout_size[0]):
        ax.plot(qso_pix_cutout[0], qso_pix_cutout[1], 'bo', markersize=8, 
                 markerfacecolor='cyan', markeredgewidth=2, label='QSO', 
                 transform=ax.get_transform('pixel'))
        ax.plot([fitted_center[0], qso_pix_cutout[0]], 
                 [fitted_center[1], qso_pix_cutout[1]], 
                 'b--', linewidth=2, alpha=0.7, label='Galaxy→QSO', 
                 transform=ax.get_transform('pixel'))
    else:
        # QSO outside cutout - draw direction arrow from fitted galaxy center
        qso_direction = np.array([dx_qso, dy_qso])
        qso_direction = qso_direction / np.linalg.norm(qso_direction) * 45
        
        ax.arrow(fitted_center[0], fitted_center[1], 
                  qso_direction[0], qso_direction[1],
                  head_width=4, head_length=4, fc='yellow', ec='yellow',
                  linewidth=3, alpha=0.9, label='QSO direction',
                  transform=ax.get_transform('pixel'))
    
    # Add correct compass with celestial directions
    _add_celestial_compass(ax, results, cutout)
    
    # Celestial coordinate setup
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('Dec')
    ax.coords[0].set_major_formatter('hh:mm:ss.s')
    ax.coords[1].set_major_formatter('dd:mm:ss')
    ax.coords.grid(color='white', alpha=0.3, linestyle='-', linewidth=0.5)
    
    ax.set_title(f'Galaxy-QSO Analysis ({method}) - PA: {results["pa_major"]:.1f}°, Alignment: {results["delta_pa"]:.1f}°')
    ax.legend(loc='upper right', framealpha=0.8)
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(save_plot, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {save_plot}")
    
    plt.show()

def _add_celestial_compass(ax, results, cutout):
    """Add compass with correct celestial North/East directions."""
    
    # Position in upper left corner using relative coordinates
    compass_x = 0.15  # 15% from left edge
    compass_y = 0.85  # 85% from bottom (near top)
    compass_size = 0.08  # 8% of plot size
    
    # Convert relative coordinates to data coordinates
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    compass_x_data = xlim[0] + compass_x * (xlim[1] - xlim[0])
    compass_y_data = ylim[0] + compass_y * (ylim[1] - ylim[0])
    compass_size_data = compass_size * min((xlim[1] - xlim[0]), (ylim[1] - ylim[0]))
    
    # Calculate true celestial North and East directions in pixel coordinates
    cutout_wcs = cutout.wcs
    
    # Center point in pixel coordinates
    center_pix = (compass_x_data, compass_y_data)
    
    # Convert center to world coordinates
    center_world = cutout_wcs.pixel_to_world(center_pix[0], center_pix[1])
    
    # Calculate points offset in celestial North and East directions
    offset_deg = 0.01  # Small offset in degrees
    
    # North direction: increase Dec, keep RA same
    north_world = SkyCoord(center_world.ra, center_world.dec + offset_deg * u.deg)
    north_pix = cutout_wcs.world_to_pixel(north_world)
    north_dx = north_pix[0] - center_pix[0]
    north_dy = north_pix[1] - center_pix[1]
    north_norm = np.sqrt(north_dx**2 + north_dy**2)
    if north_norm > 0:
        north_dx = (north_dx / north_norm) * compass_size_data
        north_dy = (north_dy / north_norm) * compass_size_data
    
    # East direction: increase RA, keep Dec same  
    east_world = SkyCoord(center_world.ra + offset_deg * u.deg / np.cos(center_world.dec.radian), 
                         center_world.dec)
    east_pix = cutout_wcs.world_to_pixel(east_world)
    east_dx = east_pix[0] - center_pix[0]
    east_dy = east_pix[1] - center_pix[1]
    east_norm = np.sqrt(east_dx**2 + east_dy**2)
    if east_norm > 0:
        east_dx = (east_dx / east_norm) * compass_size_data
        east_dy = (east_dy / east_norm) * compass_size_data
    
    # Draw North arrow
    ax.arrow(compass_x_data, compass_y_data, north_dx, north_dy,
             head_width=compass_size_data*0.15, head_length=compass_size_data*0.15, 
             fc='white', ec='white', alpha=1.0, linewidth=2)
    
    # Draw East arrow
    ax.arrow(compass_x_data, compass_y_data, east_dx, east_dy,
             head_width=compass_size_data*0.15, head_length=compass_size_data*0.15, 
             fc='white', ec='white', alpha=1.0, linewidth=2)
    
    # Compass labels in white
    label_offset = compass_size_data * 1.3
    
    # North label
    ax.text(compass_x_data + north_dx * 1.5, compass_y_data + north_dy * 1.5, 'N', 
            ha='center', va='center', fontsize=12, fontweight='bold', color='white')
    
    # East label  
    ax.text(compass_x_data + east_dx * 1.5, compass_y_data + east_dy * 1.5, 'E', 
            ha='center', va='center', fontsize=12, fontweight='bold', color='white')
    
    # Angle measurement
    angle_text = f'{results["delta_pa"]:.1f}°\n{results["alignment_type"]}'
    ax.text(compass_x_data, compass_y_data - compass_size_data * 1.8, angle_text,
            ha='center', va='center', fontsize=10, color='white', fontweight='bold')

# Example usage
if __name__ == "__main__":
    fits_file = 'your_image.fits'
    ra_gal, dec_gal = 150.114, 2.205
    ra_qso, dec_qso = 150.120, 2.208
    
    # Simple usage
    results = compute_galaxy_qso_pa(
        fits_file, ra_gal, dec_gal, ra_qso, dec_qso,
        cutout_size=(160, 160),
        initial_sma=10,
        initial_eps=0.3,
        plot=True,
        save_plot='galaxy_qso_analysis.png'
    )
    
    print(f"Final alignment: {results['delta_pa']:.1f}° ({results['alignment_type']})")