from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from photutils.isophote import Ellipse, EllipseGeometry
from photutils.aperture import EllipticalAperture
from astropy.nddata import Cutout2D

def compute_galaxy_qso_pa(fits_file, ra_gal, dec_gal, ra_qso, dec_qso, 
                          cutout_size=(100, 100), initial_sma=10, 
                          initial_eps=0.3, plot=False, save_plot=None,
                          try_multiple_params=True):
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
    plot : bool, optional
        Whether to create a diagnostic plot (default: False)
    save_plot : str, optional
        Filename to save the plot (default: None, display only)
    try_multiple_params : bool, optional
        Try different initial parameters if first attempt fails (default: True)
    
    Returns:
    --------
    dict : Dictionary containing:
        - 'pa_major': Galaxy major axis PA in degrees (east of north)
        - 'pa_to_qso': Galaxy to QSO PA in degrees (east of north)
        - 'delta_pa': Alignment angle (0°=major axis, 90°=minor axis)
        - 'alignment_type': String describing alignment type
        - 'galaxy_pix': Galaxy pixel coordinates in original image
        - 'qso_pix': QSO pixel coordinates in original image
        - 'cutout': Cutout2D object for further analysis
        - 'isolist': Isophote fitting results
    """
    
    # Load data
    hdu = fits.open(fits_file)[0]
    data = hdu.data
    wcs = WCS(hdu.header)
    
    # Get pixel positions
    gal_coord = SkyCoord(ra_gal, dec_gal, unit='deg')
    qso_coord = SkyCoord(ra_qso, dec_qso, unit='deg')
    gal_pix = wcs.world_to_pixel(gal_coord)
    qso_pix = wcs.world_to_pixel(qso_coord)
    
    # Make cutout around galaxy
    cutout = Cutout2D(data, position=gal_pix, size=cutout_size, wcs=wcs)
    cut_data = cutout.data
    
    # Initial geometry for ellipse fitting
    geometry = EllipseGeometry(
        x0=cutout_size[1]/2, 
        y0=cutout_size[0]/2, 
        sma=initial_sma, 
        eps=initial_eps, 
        pa=0
    )
    
    # Try fitting isophotes with error handling
    isolist = None
    attempt_params = [
        (initial_sma, initial_eps),
        (initial_sma * 0.5, 0.1),     # Smaller, more circular
        (initial_sma * 1.5, 0.5),     # Larger, more elliptical  
        (initial_sma * 0.3, 0.2),     # Very small start
        (initial_sma * 2.0, 0.3),     # Large start
    ]
    
    for i, (sma, eps) in enumerate(attempt_params):
        try:
            geometry = EllipseGeometry(
                x0=cutout_size[1]/2, 
                y0=cutout_size[0]/2, 
                sma=sma, 
                eps=eps, 
                pa=0
            )
            ellipse = Ellipse(cut_data, geometry)
            isolist = ellipse.fit_image()
            
            if len(isolist) > 0:
                print(f"Isophote fitting successful on attempt {i+1}")
                break
                
        except Exception as e:
            if i == 0:
                print(f"Initial fitting attempt failed: {e}")
            continue
    
    if isolist is None or len(isolist) == 0:
        raise ValueError("All isophote fitting attempts failed. Check your data and parameters.")
    
    # Check if fitting was successful and get best isophote
    valid_isos = [iso for iso in isolist if hasattr(iso, 'valid') and iso.valid]
    if len(valid_isos) == 0:
        print("Warning: No explicitly valid isophotes found, using last fitted isophote")
        best_iso = isolist[-1]
    else:
        best_iso = valid_isos[-1]  # Use outermost valid isophote
        
    pa_major = np.rad2deg(best_iso.pa) % 360  # Convert to degrees and normalize
    
    # Compute vector PA to QSO in celestial coordinates (proper way)
    # Get the celestial coordinates and compute proper PA
    gal_skycoord = SkyCoord(ra_gal, dec_gal, unit='deg')
    qso_skycoord = SkyCoord(ra_qso, dec_qso, unit='deg')
    
    # Compute position angle from galaxy to QSO in celestial coordinates
    pa_to_qso_celestial = gal_skycoord.position_angle(qso_skycoord).to(u.deg).value
    
    # Convert galaxy major axis PA from pixel coordinates to celestial coordinates
    # We need to account for image rotation relative to celestial coordinates
    
    # Method 1: Use WCS transformation
    # Create two points along the major axis in cutout coordinates
    center_x, center_y = cutout_size[1]/2, cutout_size[0]/2
    pa_rad = best_iso.pa
    
    # Points along major axis in cutout coordinates
    dx = 10 * np.cos(pa_rad)  # Small offset along major axis
    dy = 10 * np.sin(pa_rad)
    
    point1_cut = (center_x, center_y)
    point2_cut = (center_x + dx, center_y + dy)
    
    # Convert cutout coordinates back to original image coordinates
    # The cutout has its own WCS
    cutout_wcs = cutout.wcs
    point1_world = cutout_wcs.pixel_to_world(point1_cut[0], point1_cut[1])
    point2_world = cutout_wcs.pixel_to_world(point2_cut[0], point2_cut[1])
    
    # Calculate celestial PA of major axis
    pa_major_celestial = point1_world.position_angle(point2_world).to(u.deg).value
    
    # Now we have both PAs in proper celestial coordinates
    pa_major = pa_major_celestial
    pa_to_qso = pa_to_qso_celestial
    
    # Compute relative angle (0° = along major axis, 90° = along minor axis)
    delta_pa = abs(pa_to_qso - pa_major)
    
    # Handle wraparound (e.g., 350° - 10° = 340°, but should be 20°)
    if delta_pa > 180:
        delta_pa = 360 - delta_pa
    
    # Force to 0-90° range (major/minor axis alignment is symmetric)
    if delta_pa > 90:
        delta_pa = 180 - delta_pa
    
    # Determine alignment type
    if delta_pa <= 45:
        alignment_type = "major axis aligned"
    else:
        alignment_type = "minor axis aligned" 
    
    # Results dictionary
    results = {
        'pa_major': pa_major,
        'pa_to_qso': pa_to_qso,
        'delta_pa': delta_pa,
        'alignment_type': alignment_type,
        'galaxy_pix': gal_pix,
        'qso_pix': qso_pix,
        'cutout': cutout,
        'isolist': isolist
    }
    
    # Print results
    print(f"Galaxy major axis PA = {pa_major:.2f}° (celestial, east of north)")
    print(f"Galaxy → QSO PA = {pa_to_qso:.2f}° (celestial)")
    print(f"Alignment angle = {delta_pa:.2f}° ({alignment_type})")
    print(f"  • 0° = perfect major axis alignment")
    print(f"  • 90° = perfect minor axis alignment") 
    print(f"Image rotation from celestial north = {_get_image_rotation(wcs):.2f}°")
    
    # Create diagnostic plot if requested
    if plot:
        _create_diagnostic_plot(results, cut_data, cutout_size, save_plot)
    
    return results

def _get_image_rotation(wcs):
    """Calculate how much the image is rotated from celestial north."""
    try:
        # Get CD matrix or PC matrix
        if hasattr(wcs.wcs, 'cd') and wcs.wcs.cd is not None:
            cd = wcs.wcs.cd
        elif hasattr(wcs.wcs, 'pc') and wcs.wcs.pc is not None:
            cd = wcs.wcs.pc * np.array([[wcs.wcs.cdelt[0], 0], [0, wcs.wcs.cdelt[1]]])
        else:
            return 0.0
        
        # Calculate rotation angle from CD matrix
        rotation = np.degrees(np.arctan2(-cd[0,1], cd[1,1]))
        return rotation % 360
    except:
        return 0.0

def _create_diagnostic_plot(results, cut_data, cutout_size, save_plot=None):
    """Create diagnostic plot showing galaxy, fitted ellipse, and QSO direction."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Cutout with fitted ellipse
    im1 = ax1.imshow(cut_data, origin='lower', cmap='viridis')
    plt.colorbar(im1, ax=ax1, shrink=0.8)
    
    # Plot fitted ellipse from best isophote
    best_iso = results['isolist'][-1]
    ellipse_aperture = EllipticalAperture(
        (best_iso.x0, best_iso.y0),
        best_iso.sma,
        best_iso.sma * (1 - best_iso.eps),
        theta=best_iso.pa
    )
    ellipse_aperture.plot(ax=ax1, color='red', linewidth=2, label='Fitted ellipse')
    
    # Mark galaxy center
    ax1.plot(cutout_size[1]/2, cutout_size[0]/2, 'r+', markersize=10, markeredgewidth=2)
    
    # Draw major axis line
    center_x, center_y = cutout_size[1]/2, cutout_size[0]/2
    sma = best_iso.sma
    pa_rad = best_iso.pa
    
    # Major axis endpoints
    dx_maj = sma * np.cos(pa_rad)
    dy_maj = sma * np.sin(pa_rad)
    ax1.plot([center_x - dx_maj, center_x + dx_maj], 
             [center_y - dy_maj, center_y + dy_maj], 
             'r--', linewidth=2, label=f'Major axis (PA={results["pa_major"]:.1f}°)')
    
    ax1.set_title('Galaxy Cutout with Fitted Ellipse')
    ax1.set_xlabel('Pixels')
    ax1.set_ylabel('Pixels')
    ax1.legend()
    
    # Plot 2: Vector diagram
    ax2.set_xlim(-1.5, 1.5)
    ax2.set_ylim(-1.5, 1.5)
    ax2.set_aspect('equal')
    
    # Draw major axis vector
    pa_maj_rad = np.deg2rad(results['pa_major'])
    ax2.arrow(0, 0, np.sin(pa_maj_rad), np.cos(pa_maj_rad), 
             head_width=0.1, head_length=0.1, fc='red', ec='red', 
             linewidth=2, label=f'Major axis ({results["pa_major"]:.1f}°)')
    
    # Draw QSO direction vector  
    pa_qso_rad = np.deg2rad(results['pa_to_qso'])
    ax2.arrow(0, 0, np.sin(pa_qso_rad), np.cos(pa_qso_rad), 
             head_width=0.1, head_length=0.1, fc='blue', ec='blue', 
             linewidth=2, label=f'To QSO ({results["pa_to_qso"]:.1f}°)')
    
    # Draw angle arc
    angle_range = np.linspace(pa_maj_rad, pa_qso_rad, 50)
    if abs(results['pa_to_qso'] - results['pa_major']) > 180:
        # Handle wrap-around case
        if results['pa_major'] > results['pa_to_qso']:
            angle_range = np.linspace(pa_maj_rad, pa_qso_rad + 2*np.pi, 50)
        else:
            angle_range = np.linspace(pa_maj_rad - 2*np.pi, pa_qso_rad, 50)
    
    arc_r = 0.5
    ax2.plot(arc_r * np.sin(angle_range), arc_r * np.cos(angle_range), 
             'g-', linewidth=2)
    
    # Add angle label
    mid_angle = (pa_maj_rad + pa_qso_rad) / 2
    if abs(results['pa_to_qso'] - results['pa_major']) > 180:
        mid_angle += np.pi
    
    # Label with alignment info
    angle_text = f'{results["delta_pa"]:.1f}°\n({results["alignment_type"]})'
    ax2.text(0.7 * np.sin(mid_angle), 0.7 * np.cos(mid_angle), 
             angle_text, fontsize=10, ha='center', va='center',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add compass
    ax2.text(0, 1.3, 'N', ha='center', va='center', fontsize=12, fontweight='bold')
    ax2.text(1.3, 0, 'E', ha='center', va='center', fontsize=12, fontweight='bold')
    
    ax2.set_title('Position Angle Comparison')
    ax2.set_xlabel('East →')
    ax2.set_ylabel('North →')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plot:
        plt.savefig(save_plot, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {save_plot}")
    
    plt.show()

# Example usage:
if __name__ == "__main__":
    # Example parameters
    fits_file = 'your_image.fits'
    ra_gal, dec_gal = 150.114, 2.205    # Galaxy position (deg)
    ra_qso, dec_qso = 150.120, 2.208    # QSO position (deg)
    
    # Run analysis with plotting
    results = compute_galaxy_qso_pa(
        fits_file, ra_gal, dec_gal, ra_qso, dec_qso,
        cutout_size=(100, 100),
        initial_sma=10,
        initial_eps=0.3,
        plot=True,
        save_plot='galaxy_qso_analysis.png'
    )
    
    # Access individual results
    print(f"\nDetailed results:")
    print(f"Galaxy pixel position: {results['galaxy_pix']}")
    print(f"QSO pixel position: {results['qso_pix']}")
    print(f"Number of fitted isophotes: {len(results['isolist'])}")