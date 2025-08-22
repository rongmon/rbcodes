from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np  # Import numpy for cosine calculation

def compute_offset(ra1, dec1, ra2, dec2, unit="arcsec"):
    """
    Compute and print the telescope offset (East, North) to move from (RA1, DEC1) to (RA2, DEC2).

    Parameters:
        ra1, dec1 : float or str
            Right Ascension and Declination of starting point (degrees or sexagesimal).
        ra2, dec2 : float or str
            Right Ascension and Declination of target point (degrees or sexagesimal).
        unit : str, optional
            Unit for output offsets ("arcsec" or "arcmin"). Default is "arcsec".

    Returns:
        tuple : (east_offset, north_offset)
            Offset in chosen unit (arcsec or arcmin).

    -----
    # Example Usage:
    # Using degrees
    ra1, dec1 = 150.1, 2.3  # Start coordinates in degrees
    ra2, dec2 = 150.2, 2.35  # Target coordinates in degrees
    
    # Using sexagesimal (HMS/DMS)
    ra1_s, dec1_s = "10h00m00s", "+02d18m00s"
    ra2_s, dec2_s = "10h00m48s", "+02d21m00s"
    
    # Using colon-separated format
    ra1_colon, dec1_colon = "11:48:15.5526", "+52:51:56.180"
    ra2_colon, dec2_colon = "11:48:17.9932", "+52:52:01.502"
    
    # Compute and print offsets
    compute_offset(ra1_colon, dec1_colon, ra2_colon, dec2_colon, unit="arcsec")
    """
    # Validate unit
    if unit not in ["arcsec", "arcmin"]:
        raise ValueError("Invalid unit. Choose 'arcsec' or 'arcmin'.")

    # Detect if input coordinates are in sexagesimal format (contain ':')
    def is_sexagesimal(coord):
        return isinstance(coord, str) and ":" in coord

    # Determine unit based on input format
    unit_type = (u.hourangle, u.deg) if is_sexagesimal(ra1) else (u.deg, u.deg)

    # Create SkyCoord objects (auto-detects format)
    coord1 = SkyCoord(ra1, dec1, unit=unit_type, frame="icrs")
    coord2 = SkyCoord(ra2, dec2, unit=unit_type, frame="icrs")

    # Convert declination to radians for cosine calculation
    cos_dec = np.cos(coord1.dec.radian)

    # Compute offsets in arcseconds
    delta_ra = ((coord2.ra - coord1.ra).to(u.arcsec) * cos_dec).value  # East offset
    delta_dec = (coord2.dec - coord1.dec).to(u.arcsec).value  # North offset

    # Convert to arcminutes if requested
    if unit == "arcmin":
        delta_ra /= 60
        delta_dec /= 60

    # Print the result
    print(f"Move telescope: East Offset = {delta_ra:.3f} {unit}, North Offset = {delta_dec:.3f} {unit}")

    return delta_ra, delta_dec  # Return values

