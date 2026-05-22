"""
IFUCube base class — defines the interface all loaders must implement.
"""
import numpy as np
from pathlib import Path


class IFUCube:
    """
    Base class for IFU data cubes loaded from FITS files.

    Attributes
    ----------
    flux : np.ndarray, shape (n_wave, ny, nx)
        Flux data cube.
    var : np.ndarray, same shape as flux, or None
        Variance cube (or None if not present).
    wave : np.ndarray, shape (n_wave,)
        Wavelength array in Angstroms.
    header : astropy.io.fits.Header
        Primary FITS header (full 3D).
    wcs : astropy.wcs.WCS or None
        2D spatial WCS (3rd axis stripped).
    name : str
        Filename stem (no directory, no extension).
    path : str
        Full file path.
    """

    def __init__(self):
        self.flux = None
        self.var = None
        self.wave = None
        self.header = None
        self.wcs = None
        self.name = ''
        self.path = ''

    def load(self, path):
        """Load data from *path*. Must set flux, var, wave, header, wcs, name."""
        raise NotImplementedError

    def spatial_header(self):
        """Return a 2D FITS header suitable for ds9 / region WCS transforms."""
        if self.header is None:
            return None
        from astropy.io import fits
        h = self.header.copy()
        # Remove 3rd-axis keywords so the header describes a 2D image
        for key in ('NAXIS3', 'CRVAL3', 'CRPIX3', 'CDELT3', 'CD3_3',
                    'CTYPE3', 'CUNIT3', 'CROTA3'):
            h.remove(key, ignore_missing=True)
        h['NAXIS'] = 2
        return h

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def n_wave(self):
        return self.flux.shape[0] if self.flux is not None else 0

    @property
    def ny(self):
        return self.flux.shape[1] if self.flux is not None else 0

    @property
    def nx(self):
        return self.flux.shape[2] if self.flux is not None else 0

    @property
    def wave_range(self):
        """(wave_min, wave_max) tuple in Angstroms."""
        if self.wave is None or len(self.wave) == 0:
            return (None, None)
        return (float(self.wave[0]), float(self.wave[-1]))

    def __repr__(self):
        shape = self.flux.shape if self.flux is not None else 'not loaded'
        wr = self.wave_range
        return (f"<{self.__class__.__name__} '{self.name}' "
                f"shape={shape} wave={wr[0]:.1f}–{wr[1]:.1f} Å>")
