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
                          method='auto', moments_sigma=2.0, plot=False, save_plot=None,
                          background_method='sigma_clip', background_value=None, 
                          background_percentile=25, annulus_r_in=None, annulus_r_out=None,
                          two_component=False, transition_radius=20):
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
    method : str, optional
        Method: 'isophote', 'moments', or 'auto' (default: 'auto')
    two_component : bool, optional
        Fit inner and outer components separately (default: False)
    transition_radius : float, optional
        Radius separating inner and outer components in pixels (default: 20)
    background_method : str, optional
        Background: 'sigma_clip', 'percentile', 'annulus', 'manual' (default: 'sigma_clip')
    plot : bool, optional
        Create diagnostic plot (default: False)
    
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
    
    # Background subtraction
    background = _estimate_background(cut_data, background_method, background_value, 
                                    background_percentile, annulus_r_in, annulus_r_out)
    cut_data_bkg_sub = cut_data - background
    
    # Handle negative values by adding positive offset
    min_value = np.min(cut_data_bkg_sub)
    if min_value < 0:
        offset = abs(min_value) + 0.001  # Small positive buffer
        cut_data_clean = cut_data_bkg_sub + offset
        print(f"Background subtracted: {background:.6f}")
        print(f"Added offset {offset:.6f} to handle negative values (min was {min_value:.6f})")
    else:
        offset = 0.0
        cut_data_clean = np.maximum(cut_data_bkg_sub, 0.001)  # Ensure minimum positive value
        print(f"Background subtracted: {background:.6f}")
        if np.min(cut_data_clean) < 0.001:
            print("Applied minimum threshold of 0.001 to avoid zero values")
    
    # Find galaxy centroid
    try:
        x_centroid, y_centroid = centroid_com(cut_data_clean)
        print(f"Galaxy centroid: ({x_centroid:.2f}, {y_centroid:.2f})")
    except:
        x_centroid, y_centroid = cutout_size[1]/2, cutout_size[0]/2
        print("Using geometric center")
    
    # Fit galaxy shape
    method_used, shape_results = _fit_galaxy_shape(
        cut_data_clean, method, initial_sma, initial_eps, moments_sigma,
        x_centroid, y_centroid, two_component, transition_radius
    )
    
    # Convert to celestial coordinates and compute alignments
    pa_to_qso_celestial = gal_coord.position_angle(qso_coord).to(u.deg).value
    
    # Build results
    results = {
        'pa_to_qso': pa_to_qso_celestial,
        'method_used': method_used,
        'two_component': two_component,
        'galaxy_pix': gal_pix,
        'qso_pix': qso_pix,
        'cutout': cutout,
        'background_subtracted': background,
        'offset_applied': offset,
        'centroid': (x_centroid, y_centroid)
    }
    
    if two_component and 'pa_inner' in shape_results:
        # Two component results
        pa_inner_cel = _pixel_to_celestial_pa(shape_results['pa_inner'], cutout, cutout_size)
        pa_outer_cel = _pixel_to_celestial_pa(shape_results['pa_outer'], cutout, cutout_size)
        
        delta_inner = _compute_alignment_angle(pa_inner_cel, pa_to_qso_celestial)
        delta_outer = _compute_alignment_angle(pa_outer_cel, pa_to_qso_celestial)
        
        results.update({
            'pa_major': pa_outer_cel,      # Use outer for main result
            'pa_inner': pa_inner_cel,
            'pa_outer': pa_outer_cel,
            'delta_pa': delta_outer,       # Use outer for main result
            'delta_pa_inner': delta_inner,
            'delta_pa_outer': delta_outer,
            'transition_radius': transition_radius
        })
        
        # Print two-component results
        print(f"Inner PA: {pa_inner_cel:.1f}°, Outer PA: {pa_outer_cel:.1f}°")
        print(f"QSO alignment - Inner: {delta_inner:.1f}°, Outer: {delta_outer:.1f}°")
        
    else:
        # Single component results
        pa_major_cel = _pixel_to_celestial_pa(shape_results['pa_pixel'], cutout, cutout_size)
        delta_pa = _compute_alignment_angle(pa_major_cel, pa_to_qso_celestial)
        
        results.update({
            'pa_major': pa_major_cel,
            'delta_pa': delta_pa
        })
        
        print(f"Galaxy PA: {pa_major_cel:.1f}°, QSO alignment: {delta_pa:.1f}°")
    
    # Add alignment type
    alignment_type = "major axis aligned" if results['delta_pa'] <= 45 else "minor axis aligned"
    results['alignment_type'] = alignment_type
    results['shape_results'] = shape_results
    
    print(f"Method: {method_used}, Alignment: {results['delta_pa']:.1f}° ({alignment_type})")
    
    if plot:
        _create_plot(results, cut_data_clean, cutout_size, save_plot)
    
    return results

def _estimate_background(data, method, value, percentile, r_in, r_out):
    """Estimate background using specified method."""
    
    if method == 'manual':
        return value
    elif method == 'percentile':
        return np.percentile(data, percentile)
    elif method == 'sigma_clip':
        try:
            from astropy.stats import sigma_clipped_stats
            _, median, _ = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
            return median
        except ImportError:
            return np.percentile(data, 25)
    elif method == 'annulus':
        return _annulus_background(data, r_in, r_out)
    else:
        return np.percentile(data, 25)

def _annulus_background(data, r_in, r_out):
    """Background from annulus around center."""
    
    ny, nx = data.shape
    center_x, center_y = nx // 2, ny // 2
    max_radius = min(center_x, center_y)
    
    if r_in is None:
        r_in = max_radius * 0.6
    if r_out is None:
        r_out = max_radius * 0.9
    if r_out >= max_radius:
        r_out = max_radius - 1
    if r_in >= r_out:
        r_in = r_out * 0.7
    
    y_indices, x_indices = np.indices(data.shape)
    distances = np.sqrt((x_indices - center_x)**2 + (y_indices - center_y)**2)
    annulus_mask = (distances >= r_in) & (distances <= r_out)
    
    if np.sum(annulus_mask) == 0:
        return np.percentile(data, 25)
    
    annulus_pixels = data[annulus_mask]
    return np.median(annulus_pixels)

def _fit_galaxy_shape(data, method, initial_sma, initial_eps, moments_sigma,
                      x_centroid, y_centroid, two_component, transition_radius):
    """Fit galaxy shape using specified method."""
    
    if method == 'moments':
        return _use_moments_method(data, moments_sigma, x_centroid, y_centroid)
    elif method == 'isophote':
        return _use_isophote_method(data, initial_sma, initial_eps, x_centroid, y_centroid,
                                  two_component, transition_radius)
    elif method == 'auto':
        try:
            return _use_isophote_method(data, initial_sma, initial_eps, x_centroid, y_centroid,
                                      two_component, transition_radius)
        except Exception as e:
            print(f"Isophote failed: {e}, switching to moments")
            return _use_moments_method(data, moments_sigma, x_centroid, y_centroid)
    else:
        raise ValueError("Method must be 'isophote', 'moments', or 'auto'")

def _use_isophote_method(data, initial_sma, initial_eps, x_centroid, y_centroid,
                        two_component, transition_radius):
    """Fit isophotes."""
    
    if not two_component:
        # Single component
        isolist = _fit_single_isophote(data, initial_sma, initial_eps, x_centroid, y_centroid)
        best_iso = isolist[-1]
        return 'isophote', {
            'pa_pixel': np.rad2deg(best_iso.pa),
            'ellipticity': best_iso.eps,
            'center': (best_iso.x0, best_iso.y0)
        }
    
    else:
        # Two component
        results = {}
        
        # Create masks
        y_indices, x_indices = np.indices(data.shape)
        distances = np.sqrt((x_indices - x_centroid)**2 + (y_indices - y_centroid)**2)
        inner_mask = distances <= transition_radius
        outer_mask = distances > transition_radius
        
        # Fit inner
        try:
            inner_data = np.where(inner_mask, data, 0)
            isolist_inner = _fit_single_isophote(inner_data, min(initial_sma, transition_radius*0.7), 
                                                initial_eps, x_centroid, y_centroid)
            best_inner = isolist_inner[-1]
            results['pa_inner'] = np.rad2deg(best_inner.pa)
            print(f"Inner fitted: PA={results['pa_inner']:.1f}°")
        except:
            results['pa_inner'] = 0.0
            print("Inner fitting failed")
        
        # Fit outer
        try:
            outer_data = np.where(outer_mask, data, 0)
            isolist_outer = _fit_single_isophote(outer_data, max(initial_sma, transition_radius*1.2),
                                                initial_eps, x_centroid, y_centroid)
            best_outer = isolist_outer[-1]
            results['pa_outer'] = np.rad2deg(best_outer.pa)
            results['pa_pixel'] = results['pa_outer']  # Use outer as main
            print(f"Outer fitted: PA={results['pa_outer']:.1f}°")
        except:
            # Fallback to single component
            isolist = _fit_single_isophote(data, initial_sma, initial_eps, x_centroid, y_centroid)
            best_iso = isolist[-1]
            return 'isophote', {
                'pa_pixel': np.rad2deg(best_iso.pa),
                'ellipticity': best_iso.eps,
                'center': (best_iso.x0, best_iso.y0)
            }
        
        results['center'] = (x_centroid, y_centroid)
        return 'isophote', results

def _fit_single_isophote(data, initial_sma, initial_eps, x_centroid, y_centroid):
    """Fit single isophote with multiple attempts."""
    
    param_sets = [
        (initial_sma, initial_eps),
        (initial_sma * 0.5, 0.1),
        (initial_sma * 1.5, 0.5),
        (initial_sma * 2.0, 0.3)
    ]
    
    for sma, eps in param_sets:
        try:
            geometry = EllipseGeometry(x0=x_centroid, y0=y_centroid, sma=sma, eps=eps, pa=0)
            ellipse = Ellipse(data, geometry)
            isolist = ellipse.fit_image()
            if len(isolist) > 0:
                return isolist
        except:
            continue
    
    raise ValueError("All isophote attempts failed")

def _use_moments_method(data, sigma, x_centroid, y_centroid):
    """Use second moments method."""
    
    if sigma > 0:
        data = gaussian_filter(data, sigma=sigma)
    
    total_flux = np.sum(data)
    if total_flux <= 0:
        raise ValueError("No positive flux")
    
    # Calculate moments
    y_indices, x_indices = np.indices(data.shape)
    x_center = np.sum(x_indices * data) / total_flux
    y_center = np.sum(y_indices * data) / total_flux
    
    dx = x_indices - x_center
    dy = y_indices - y_center
    
    Mxx = np.sum(dx * dx * data) / total_flux
    Myy = np.sum(dy * dy * data) / total_flux
    Mxy = np.sum(dx * dy * data) / total_flux
    
    # Eigenanalysis
    M = np.array([[Mxx, Mxy], [Mxy, Myy]])
    eigenvals, eigenvecs = np.linalg.eigh(M)
    
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]
    
    major_axis_vec = eigenvecs[:, 0]
    pa_radians = np.arctan2(major_axis_vec[0], major_axis_vec[1])
    pa_degrees = np.degrees(pa_radians) % 360
    
    ellipticity = 1 - np.sqrt(eigenvals[1] / eigenvals[0]) if eigenvals[0] > 0 else 0
    
    print(f"Moments: PA={pa_degrees:.1f}°, ellipticity={ellipticity:.3f}")
    
    return 'moments', {
        'pa_pixel': pa_degrees,
        'ellipticity': ellipticity,
        'center': (x_center, y_center)
    }

def _pixel_to_celestial_pa(pa_pixel, cutout, cutout_size):
    """Convert pixel PA to celestial PA."""
    
    center_x, center_y = cutout_size[1]/2, cutout_size[0]/2
    pa_rad = np.deg2rad(pa_pixel)
    
    dx = 10 * np.cos(pa_rad)
    dy = 10 * np.sin(pa_rad)
    
    point1 = (center_x, center_y)
    point2 = (center_x + dx, center_y + dy)
    
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
    """Create three-panel diagnostic plot."""
    
    fig = plt.figure(figsize=(18, 6))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133, projection=results['cutout'].wcs)
    
    # Panel 1: Cutout with shape
    im1 = ax1.imshow(cut_data, origin='lower', cmap='viridis')
    plt.colorbar(im1, ax=ax1, shrink=0.8)
    
    shape = results['shape_results']
    center = shape['center']
    
    # Draw ellipse and major axis
    if results['method_used'] == 'isophote' and 'pa_outer' in shape:
        # Two component
        ax1.plot(center[0], center[1], 'r+', markersize=12, markeredgewidth=3)
        
        # Draw transition radius circle
        circle = plt.Circle(center, results['transition_radius'], 
                          fill=False, color='yellow', linestyle='--', linewidth=2)
        ax1.add_patch(circle)
        
        # Draw major axes
        for component, color, label in [('inner', 'orange', 'Inner'), ('outer', 'red', 'Outer')]:
            pa_key = f'pa_{component}'
            if pa_key in shape:
                pa_rad = np.deg2rad(shape[pa_key])
                r = results['transition_radius'] if component == 'inner' else results['transition_radius'] * 1.5
                dx = r * np.cos(pa_rad)
                dy = r * np.sin(pa_rad)
                ax1.plot([center[0] - dx, center[0] + dx], 
                        [center[1] - dy, center[1] + dy],
                        color=color, linewidth=3, label=f'{label} axis')
        
        ax1.legend()
        ax1.set_title(f'Two-Component Shape ({results["method_used"]})')
        
    else:
        # Single component
        ax1.plot(center[0], center[1], 'r+', markersize=12, markeredgewidth=3)
        
        pa_rad = np.deg2rad(shape['pa_pixel'])
        r = 30
        dx = r * np.cos(pa_rad)
        dy = r * np.sin(pa_rad)
        ax1.plot([center[0] - dx, center[0] + dx], 
                [center[1] - dy, center[1] + dy],
                'r--', linewidth=3, label='Major axis')
        
        ax1.legend()
        ax1.set_title(f'Galaxy Shape ({results["method_used"]})')
    
    ax1.set_xlabel('Pixels')
    ax1.set_ylabel('Pixels')
    
    # Panel 2: PA vectors
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1.5, 1.5)
    ax2.set_aspect('equal')
    
    pa_maj_rad = np.deg2rad(results['pa_major'])
    pa_qso_rad = np.deg2rad(results['pa_to_qso'])
    
    ax2.arrow(0, 0, np.sin(pa_maj_rad), np.cos(pa_maj_rad),
              head_width=0.1, head_length=0.1, fc='red', ec='red',
              linewidth=2, label=f'Major ({results["pa_major"]:.1f}°)')
    
    ax2.arrow(0, 0, np.sin(pa_qso_rad), np.cos(pa_qso_rad),
              head_width=0.1, head_length=0.1, fc='blue', ec='blue',
              linewidth=2, label=f'To QSO ({results["pa_to_qso"]:.1f}°)')
    
    # Angle arc
    angle_range = np.linspace(pa_maj_rad, pa_qso_rad, 50)
    if abs(results['pa_to_qso'] - results['pa_major']) > 180:
        if results['pa_major'] > results['pa_to_qso']:
            angle_range = np.linspace(pa_maj_rad, pa_qso_rad + 2*np.pi, 50)
        else:
            angle_range = np.linspace(pa_maj_rad - 2*np.pi, pa_qso_rad, 50)
    
    arc_r = 0.5
    ax2.plot(arc_r * np.sin(angle_range), arc_r * np.cos(angle_range), 'g-', linewidth=2)
    
    # Labels
    mid_angle = (pa_maj_rad + pa_qso_rad) / 2
    if abs(results['pa_to_qso'] - results['pa_major']) > 180:
        mid_angle += np.pi
    
    ax2.text(0.7 * np.sin(mid_angle), 0.7 * np.cos(mid_angle),
             f'{results["delta_pa"]:.1f}°\n({results["alignment_type"]})',
             fontsize=10, ha='center', va='center',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax2.text(0, 1.3, 'N', ha='center', va='center', fontsize=12, fontweight='bold')
    ax2.text(1.3, 0, 'E', ha='center', va='center', fontsize=12, fontweight='bold')
    
    ax2.set_title('Position Angles')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Sky view
    cutout = results['cutout']
    im3 = ax3.imshow(cutout.data, origin='lower', cmap='viridis')
    plt.colorbar(im3, ax=ax3, shrink=0.8)
    
    # Galaxy and QSO positions
    gal_pix_orig = results['galaxy_pix']
    qso_pix_orig = results['qso_pix']
    
    cutout_center = (cutout_size[1]//2, cutout_size[0]//2)
    gal_pix_cutout = cutout_center
    
    dx_qso = qso_pix_orig[0] - gal_pix_orig[0]
    dy_qso = qso_pix_orig[1] - gal_pix_orig[1]
    qso_pix_cutout = (cutout_center[0] + dx_qso, cutout_center[1] + dy_qso)
    
    ax3.plot(gal_pix_cutout[0], gal_pix_cutout[1], 'r+', markersize=15, 
             markeredgewidth=3, label='Galaxy', transform=ax3.get_transform('pixel'))
    
    # QSO position or direction
    if (0 <= qso_pix_cutout[0] < cutout_size[1] and 0 <= qso_pix_cutout[1] < cutout_size[0]):
        ax3.plot(qso_pix_cutout[0], qso_pix_cutout[1], 'bo', markersize=8, 
                 markerfacecolor='cyan', label='QSO', transform=ax3.get_transform('pixel'))
        ax3.plot([gal_pix_cutout[0], qso_pix_cutout[0]], 
                 [gal_pix_cutout[1], qso_pix_cutout[1]], 
                 'b--', linewidth=2, alpha=0.7, transform=ax3.get_transform('pixel'))
    else:
        qso_direction = np.array([dx_qso, dy_qso])
        qso_direction = qso_direction / np.linalg.norm(qso_direction) * 20
        ax3.arrow(gal_pix_cutout[0], gal_pix_cutout[1], 
                  qso_direction[0], qso_direction[1],
                  head_width=3, head_length=3, fc='blue', ec='blue',
                  linewidth=2, label='QSO direction', transform=ax3.get_transform('pixel'))
    
    # Coordinate labels
    ax3.coords[0].set_axislabel('RA')
    ax3.coords[1].set_axislabel('Dec')
    ax3.coords.grid(color='white', alpha=0.5, linestyle='-', linewidth=0.5)
    
    ax3.set_title('Sky View')
    ax3.legend()
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(save_plot, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {save_plot}")
    
    plt.show()

# Example usage
if __name__ == "__main__":
    fits_file = 'your_image.fits'
    ra_gal, dec_gal = 150.114, 2.205
    ra_qso, dec_qso = 150.120, 2.208
    
    # Two-component analysis for noisy galaxies
    results = compute_galaxy_qso_pa(
        fits_file, ra_gal, dec_gal, ra_qso, dec_qso,
        cutout_size=(150, 150),
        two_component=True,
        transition_radius=25,
        background_method='sigma_clip',
        plot=True
    )
    
    print(f"Final alignment: {results['delta_pa']:.1f}°")
    if 'pa_inner' in results:
        print(f"Inner vs Outer PA: {results['pa_inner']:.1f}° vs {results['pa_outer']:.1f}°")