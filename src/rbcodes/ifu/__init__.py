"""
IFU cube tools: reading, images, moment maps, extraction, regions and line fitting.

    from rbcodes.ifu import load_fits
    cube = load_fits('cube.fits')          # cube.flux (n_wave, ny, nx), .var, .wave, .wcs (2D)

Modules
-------
io               load_fits (auto-detects KCWI/MUSE/generic cubes and 2D images), load_image
cube_collapse    white-light, narrow-band and continuum-subtracted images
moment_maps      moment 0/1/2 maps, S/N maps, linear continuum subtraction
aperture_extract aperture masks and 1D spectrum extraction
spatial_mask     ds9 region parsing, region -> pixel mask, overlays, region export
noise            noise cube from the variance, or estimated from sky spaxels
linefit          emission-line fits -> velocity / dispersion / flux maps, quality cuts
"""
from .io import load_fits
from .io.image2d import load_image

__all__ = ['load_fits', 'load_image']
