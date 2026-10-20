"""
FITSImage — wraps a plain 2D FITS image (HST, ground-based, etc.).

Same interface as IFUCube but flux is 2D and wave/var are None.
"""
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS

from .base import IFUCube


class FITSImage(IFUCube):
    """
    A 2D FITS image with the same interface as IFUCube.

    Attributes
    ----------
    flux : np.ndarray, shape (ny, nx)
        Image data.  (Not 3D — callers must check.)
    var : None
        Not applicable for 2D images.
    wave : None
        Not applicable for 2D images.
    """

    def load(self, path):
        path = Path(path)
        self.path = str(path)
        self.name = path.stem
        self.var = None
        self.wave = None

        with fits.open(path) as hdul:
            hdu = _find_2d_hdu(hdul)
            self.header = hdu.header.copy()
            self.flux = hdu.data.astype(np.float64)
            try:
                self.wcs = WCS(self.header, naxis=2)
            except Exception:
                self.wcs = None

        return self

    # Override base properties — not meaningful for 2D
    @property
    def n_wave(self):
        return 0

    @property
    def ny(self):
        return self.flux.shape[0] if self.flux is not None else 0

    @property
    def nx(self):
        return self.flux.shape[1] if self.flux is not None else 0

    @property
    def wave_range(self):
        return (None, None)

    def spatial_header(self):
        return self.header.copy() if self.header is not None else None

    def __repr__(self):
        shape = self.flux.shape if self.flux is not None else 'not loaded'
        return f"<FITSImage '{self.name}' shape={shape}>"


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _find_2d_hdu(hdul):
    """Return the first HDU that has 2D data."""
    # Prefer primary if it has data
    if hdul[0].data is not None and hdul[0].data.ndim == 2:
        return hdul[0]
    # Some files store image in extension 1 with empty primary
    for hdu in hdul[1:]:
        if hdu.data is not None and hdu.data.ndim == 2:
            return hdu
    # Last resort: squeeze extra dimensions
    for hdu in hdul:
        if hdu.data is not None and hdu.data.squeeze().ndim == 2:
            hdu.data = hdu.data.squeeze()
            return hdu
    raise ValueError("No 2D image extension found in FITS file.")


def load_image(path):
    """Convenience loader — returns a FITSImage."""
    img = FITSImage()
    img.load(path)
    return img
