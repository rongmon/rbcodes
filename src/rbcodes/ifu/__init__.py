"""
IFU cube tools shared by rb_ifuview, rb_align and rbcodes.turbulence.

    from rbcodes.ifu import load_fits
    cube = load_fits('cube.fits')          # cube.flux (n_wave, ny, nx), .var, .wave, .wcs (2D)

Modules
-------
io               load_fits (auto-detects KCWI/MUSE/generic cubes and 2D images), load_image
cube_collapse    white-light, narrow-band and continuum-subtracted images
moment_maps      moment 0/1/2 maps, S/N maps, linear continuum subtraction
aperture_extract aperture masks and 1D spectrum extraction
spatial_mask     ds9 region parsing, region -> pixel mask, overlays, region export
"""
from .io import load_fits
from .io.image2d import load_image

__all__ = ['load_fits', 'load_image']
