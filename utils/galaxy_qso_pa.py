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
                          moments_sigma=2.0, plot=False, save_plot=None):
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
        method_used, shape_results = _fit_moments(cut_data, moments_sigma)
    
    # Convert pixel PA to celestial PA
    pa_major_celestial = _pixel_to_celestial_pa(shape_results['pa_pixel'], cutout, cutout_size)
    
    # Calculate QSO direction in celestial coordinates
    pa_to_qso_celestial = gal_coord.position_angle(qso_coord).to(u.deg).value
    
    # Compute alignment angle (0-90°)
    delta_pa = _compute_alignment_angle(pa_major_celestial, pa_to_qso_celestial)
    alignment_type = "major axis aligned" if delta_pa <= 45 else "minor axis aligned"
    
    # Build results
    results = {
        'pa_major': pa_major_celestial,
        'pa_to_qso': pa_to_qso_celestial,
        'delta_pa': delta_pa,
        'alignment_type': alignment_type,
        'method_used': method_used,
        'galaxy_pix': gal_pix,
        'qso_pix': qso_pix,
        'cutout': cutout,
        'centroid': (x_centroid, y_centroid),
        'shape_results': shape_results
    }
    
    # Print results
    print(f"Method: {method_used}")
    print(f"Galaxy PA: {pa_major_celestial:.1f}°, QSO PA: {pa_to_qso_celestial:.1f}°")
    print(f"Alignment: {delta_pa:.1f}° ({alignment_type})")
    
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
                
                # Get best isophote
                valid_isos = [iso for iso in isolist if hasattr(iso, 'valid') and iso.valid]
                best_iso = valid_isos[-1] if valid_isos else isolist[-1]
                
                return 'isophote', {
                    'pa_pixel': np.rad2deg(best_iso.pa),
                    'ellipticity': best_iso.eps,
                    'center': (best_iso.x0, best_iso.y0),
                    'sma': best_iso.sma,
                    'isolist': isolist
                }
        except:
            continue
    
    raise ValueError("All isophote attempts failed")

def _fit_moments(data, sigma):
    """Fit using second moments method."""
    
    # Apply Gaussian smoothing
    if sigma > 0:
        data_smooth = gaussian_filter(data, sigma=sigma)
    else:
        data_smooth = data
    
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
    
    print(f"Moments: PA={pa_degrees:.1f}°, ellipticity={ellipticity:.3f}")
    
    return 'moments', {
        'pa_pixel': pa_degrees,
        'ellipticity': ellipticity,
        'center': (x_center, y_center),
        'r_eff': r_eff
    }

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
    """Create two-panel diagnostic plot."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Left Panel: Galaxy shape analysis
    im1 = ax1.imshow(cut_data, origin='lower', cmap='viridis')
    plt.colorbar(im1, ax=ax1, shrink=0.8)
    
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
    ellipse.plot(ax=ax1, color='red', linewidth=2, label=f'{method.title()} ellipse')
    
    # Draw major axis line
    dx = semi_major * np.cos(pa_rad)
    dy = semi_major * np.sin(pa_rad)
    ax1.plot([center[0] - dx, center[0] + dx], 
             [center[1] - dy, center[1] + dy],
             'r--', linewidth=2, alpha=0.8, label='Major axis')
    
    # Mark center
    ax1.plot(center[0], center[1], 'r+', markersize=12, markeredgewidth=3, label='Center')
    
    ax1.set_title(f'Galaxy Shape ({method})')
    ax1.set_xlabel('Pixels')
    ax1.set_ylabel('Pixels')
    ax1.legend()
    
    # Right Panel: Sky view with celestial coordinates
    cutout = results['cutout']
    ax2 = plt.subplot(122, projection=cutout.wcs)
    im2 = ax2.imshow(cutout.data, origin='lower', cmap='viridis')
    plt.colorbar(im2, ax=ax2, shrink=0.8)
    
    # Galaxy and QSO positions
    gal_pix_orig = results['galaxy_pix']
    qso_pix_orig = results['qso_pix']
    
    cutout_center = (cutout_size[1]//2, cutout_size[0]//2)
    gal_pix_cutout = cutout_center
    
    dx_qso = qso_pix_orig[0] - gal_pix_orig[0]
    dy_qso = qso_pix_orig[1] - gal_pix_orig[1]
    qso_pix_cutout = (cutout_center[0] + dx_qso, cutout_center[1] + dy_qso)
    
    # Plot galaxy
    ax2.plot(gal_pix_cutout[0], gal_pix_cutout[1], 'r+', markersize=15, 
             markeredgewidth=3, label='Galaxy', transform=ax2.get_transform('pixel'))
    
    # Plot QSO or direction arrow
    if (0 <= qso_pix_cutout[0] < cutout_size[1] and 0 <= qso_pix_cutout[1] < cutout_size[0]):
        ax2.plot(qso_pix_cutout[0], qso_pix_cutout[1], 'bo', markersize=8, 
                 markerfacecolor='cyan', markeredgewidth=2, label='QSO', 
                 transform=ax2.get_transform('pixel'))
        ax2.plot([gal_pix_cutout[0], qso_pix_cutout[0]], 
                 [gal_pix_cutout[1], qso_pix_cutout[1]], 
                 'b--', linewidth=2, alpha=0.7, label='Galaxy→QSO', 
                 transform=ax2.get_transform('pixel'))
    else:
        # QSO outside cutout - draw direction arrow
        qso_direction = np.array([dx_qso, dy_qso])
        qso_direction = qso_direction / np.linalg.norm(qso_direction) * 20
        ax2.arrow(gal_pix_cutout[0], gal_pix_cutout[1], 
                  qso_direction[0], qso_direction[1],
                  head_width=3, head_length=3, fc='blue', ec='blue',
                  linewidth=2, alpha=0.8, label='QSO direction', 
                  transform=ax2.get_transform('pixel'))
    
    # Draw galaxy major axis
    center_sky = shape['center']
    pa_rad_sky = np.deg2rad(shape['pa_pixel'])
    axis_length = sma
    
    dx_maj = axis_length * np.cos(pa_rad_sky)
    dy_maj = axis_length * np.sin(pa_rad_sky)
    ax2.plot([center_sky[0] - dx_maj, center_sky[0] + dx_maj], 
             [center_sky[1] - dy_maj, center_sky[1] + dy_maj],
             'r--', linewidth=3, alpha=0.8, label='Major axis', 
             transform=ax2.get_transform('pixel'))
    
    # Add PA compass
    _add_pa_compass(ax2, results, cutout_size)
    
    # Celestial coordinate setup
    ax2.coords[0].set_axislabel('RA')
    ax2.coords[1].set_axislabel('Dec')
    ax2.coords[0].set_major_formatter('hh:mm:ss.s')
    ax2.coords[1].set_major_formatter('dd:mm:ss')
    ax2.coords.grid(color='white', alpha=0.3, linestyle='-', linewidth=0.5)
    
    ax2.set_title('Sky View with Position Angles')
    ax2.legend(loc='upper right', framealpha=0.8)
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(save_plot, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {save_plot}")
    
    plt.show()

def _add_pa_compass(ax, results, cutout_size):
    """Add PA compass to corner of sky plot."""
    
    # Position in upper left corner
    compass_x = cutout_size[1] * 0.15
    compass_y = cutout_size[0] * 0.85
    compass_size = 15
    
    # PA vectors
    pa_maj_rad = np.deg2rad(results['pa_major'])
    pa_qso_rad = np.deg2rad(results['pa_to_qso'])
    
    # Major axis vector (red)
    ax.arrow(compass_x, compass_y, 
             compass_size * np.sin(pa_maj_rad), compass_size * np.cos(pa_maj_rad),
             head_width=2, head_length=2, fc='red', ec='red', alpha=0.8,
             linewidth=2, transform=ax.get_transform('pixel'))
    
    # QSO direction vector (blue)
    ax.arrow(compass_x, compass_y,
             compass_size * np.sin(pa_qso_rad), compass_size * np.cos(pa_qso_rad),
             head_width=2, head_length=2, fc='blue', ec='blue', alpha=0.8,
             linewidth=2, transform=ax.get_transform('pixel'))
    
    # Compass labels
    ax.text(compass_x, compass_y + compass_size + 5, 'N', 
            ha='center', va='center', fontsize=10, fontweight='bold', color='white',
            transform=ax.get_transform('pixel'))
    ax.text(compass_x + compass_size + 5, compass_y, 'E', 
            ha='center', va='center', fontsize=10, fontweight='bold', color='white',
            transform=ax.get_transform('pixel'))
    
    # Angle measurement
    angle_text = f'{results["delta_pa"]:.1f}°\n({results["alignment_type"]})'
    ax.text(compass_x, compass_y - compass_size - 10, angle_text,
            ha='center', va='center', fontsize=9, color='white',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.6),
            transform=ax.get_transform('pixel'))

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