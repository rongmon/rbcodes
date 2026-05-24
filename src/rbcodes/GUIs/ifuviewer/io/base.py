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

    def crop(self, x1, y1, x2, y2):
        """
        Return a new instance spatially cropped to pixel region x1:x2, y1:y2.

        Works for both 3-D cubes (flux shape n_wave, ny, nx) and 2-D images
        (flux shape ny, nx).  The wavelength array is unchanged.  WCS CRPIX
        values are updated automatically.

        Parameters
        ----------
        x1, y1 : int  lower-left corner (inclusive, 0-indexed)
        x2, y2 : int  upper-right corner (exclusive, numpy convention)

        Returns
        -------
        New instance of the same class with cropped data and updated header/WCS.
        """
        import re
        cropped = type(self)()

        # Flux + variance
        if self.flux.ndim == 3:
            cropped.flux = self.flux[:, y1:y2, x1:x2].copy()
            cropped.var  = (self.var[:, y1:y2, x1:x2].copy()
                            if self.var is not None else None)
        else:  # 2-D image
            cropped.flux = self.flux[y1:y2, x1:x2].copy()
            cropped.var  = None

        cropped.wave = self.wave   # shared ref — wavelength axis unchanged

        # WCS — update CRPIX manually (pixel offset is 0-indexed; CRPIX is 1-indexed)
        if self.wcs is not None:
            try:
                import copy as _copy
                wcs_new = _copy.deepcopy(self.wcs)
                wcs_new.wcs.crpix[0] -= x1
                wcs_new.wcs.crpix[1] -= y1
                wcs_new.wcs.set()
                cropped.wcs = wcs_new
            except Exception:
                cropped.wcs = None
        else:
            cropped.wcs = None

        # Header — update NAXIS1/2 and CRPIX1/2
        if self.header is not None:
            h = self.header.copy()
            h['NAXIS1'] = x2 - x1
            h['NAXIS2'] = y2 - y1
            if 'CRPIX1' in h:
                h['CRPIX1'] = float(h['CRPIX1']) - x1
            if 'CRPIX2' in h:
                h['CRPIX2'] = float(h['CRPIX2']) - y1
            if self.flux.ndim == 3:
                h['NAXIS3'] = self.flux.shape[0]
            cropped.header = h
        else:
            cropped.header = None

        # Name — strip any existing _crop suffix then re-append
        base = re.sub(r'_crop\d*$', '', self.name)
        cropped.name = base + '_crop'
        cropped.path = ''

        return cropped

    def spatial_header(self):
        """Return a 2D FITS header suitable for ds9 / region WCS transforms."""
        if self.header is None:
            return None
        from astropy.io import fits
        h = self.header.copy()
        # Remove all 3rd-axis WCS keywords (FITS standard + CD/PC matrix terms)
        scalar_keys = ('NAXIS3', 'CRVAL3', 'CRPIX3', 'CDELT3', 'CD3_3',
                       'CTYPE3', 'CUNIT3', 'CROTA3',
                       'PC3_1', 'PC3_2', 'PC3_3',
                       'PC1_3', 'PC2_3',
                       'CD1_3', 'CD2_3', 'CD3_1', 'CD3_2',
                       'LTVL3', 'LTV3', 'LTM3_3',
                       'WCSAXES')
        for key in scalar_keys:
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
