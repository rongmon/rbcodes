def write_ds9_regions(ra, dec, output_file="regions.reg", radius=1.0, color="green", labels=None):
    """
    Write RA, DEC coordinates to a DS9 region file with customizable circle size, color, and labels.

    Parameters:
        ra (array-like): Array of Right Ascension values in degrees.
        dec (array-like): Array of Declination values in degrees.
        output_file (str): Name of the output DS9 region file.
        radius (float): Circle radius in arcseconds.
        color (str): Color of the regions (e.g., 'red', 'blue', 'green').
        labels (array-like, optional): Text labels for each region (same length as ra & dec). Default is None.
    """
    if labels is None:
        labels = ["" for _ in ra]  # Empty labels if none provided

    with open(output_file, "w") as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color={} dashlist=8 3 width=1 font='helvetica 10 normal' "
                "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n".format(color))
        f.write("fk5\n")  # Use the FK5 coordinate system (J2000)

        for r, d, label in zip(ra, dec, labels):
            region_line = f"circle({r},{d},{radius}\") # color={color}"
            if label:
                region_line += f" text={{ {label} }}"
            f.write(region_line + "\n")

    print(f"DS9 region file '{output_file}' created successfully.")

