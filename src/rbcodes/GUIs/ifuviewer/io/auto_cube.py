"""
GenericCube — auto-detecting FITS cube loader.

Entry point: load_fits(path) → IFUCube subclass instance.
"""
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS

from .base import IFUCube


# Map INSTRUME header values to subclass names (populated at bottom of file)
_INSTRUMENT_MAP = {}   # filled after subclasses are defined


class GenericCube(IFUCube):
    """
    Loads any 3-axis FITS cube using generic WCS keywords.

    Strategy
    --------
    - Look for flux in extension 0 (NAXIS=3); fall back to extension 1.
    - Build wavelength array from CRVAL3/CDELT3 (or CD3_3)/CRPIX3/NAXIS3.
    - Variance: search extensions named VAR, VARIANCE, STAT, ERR.
    - Strip 3rd WCS axis to get 2D spatial WCS.
    """

    def load(self, path, var_path=None):
        path = Path(path)
        self.path = str(path)
        self.name = path.stem

        with fits.open(path) as hdul:
            flux_hdu, flux_data = _find_flux_hdu(hdul)
            self.header = flux_hdu.header.copy()
            self.flux = flux_data.astype(np.float64)
            self.var = _find_variance(hdul, self.flux.shape)
            self.wave = _build_wave(self.header)
            self.wcs = _make_2d_wcs(self.header)

        if var_path is not None:
            self.var = _load_external_var(var_path, self.flux.shape)

        return self


# ---------------------------------------------------------------------------
# KCWI subclass
# ---------------------------------------------------------------------------

class KCWICube(GenericCube):
    """KCWI-specific loader. Currently uses GenericCube logic."""

    def load(self, path, var_path=None):
        super().load(path, var_path=var_path)
        # KCWI stores flux in e-/s; variance in (e-/s)^2 — no unit conversion needed
        return self


# ---------------------------------------------------------------------------
# MUSE subclass
# ---------------------------------------------------------------------------

class MUSECube(GenericCube):
    """MUSE-specific loader. MUSE stores flux in DATA ext, variance in STAT."""

    def load(self, path, var_path=None):
        path = Path(path)
        self.path = str(path)
        self.name = path.stem

        with fits.open(path) as hdul:
            # MUSE primary is empty; data is in 'DATA'
            if 'DATA' in [h.name for h in hdul]:
                flux_hdu = hdul['DATA']
            else:
                flux_hdu, _ = _find_flux_hdu(hdul)
            self.header = flux_hdu.header.copy()
            self.flux = flux_hdu.data.astype(np.float64)
            # MUSE variance extension
            var_hdu = _get_extension(hdul, ('STAT', 'VAR', 'VARIANCE', 'ERR'))
            if var_hdu is not None and var_hdu.data is not None:
                if var_hdu.data.shape == self.flux.shape:
                    self.var = var_hdu.data.astype(np.float64)
            self.wave = _build_wave(self.header)
            self.wcs = _make_2d_wcs(self.header)

        if var_path is not None:
            self.var = _load_external_var(var_path, self.flux.shape)

        return self


# ---------------------------------------------------------------------------
# Instrument map (must come after subclasses)
# ---------------------------------------------------------------------------

_INSTRUMENT_MAP = {
    'KCWI': KCWICube,
    'MUSE': MUSECube,
    'MUSEWFM': MUSECube,
    'MUSE-WFM': MUSECube,
    'MUSE-NFM': MUSECube,
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def load_fits(path, var=None):
    """
    Load a FITS file and return the appropriate IFUCube subclass.

    Reads the INSTRUME keyword and dispatches to the matching subclass.
    Falls back to GenericCube if the instrument is unknown or absent.

    Parameters
    ----------
    path : str or Path
        Path to the flux FITS cube.
    var : str or Path, optional
        Path to a separate variance FITS file.  When provided, overrides any
        variance extension found inside *path*.

    Returns
    -------
    IFUCube subclass instance (GenericCube, KCWICube, MUSECube, …)
    """
    from .image2d import FITSImage
    path = Path(path)

    # Peek at the header to decide cube vs. 2D image
    with fits.open(path) as hdul:
        instrume = _read_instrume(hdul)
        has_3d = any(
            hdu.data is not None and hdu.data.ndim == 3 for hdu in hdul
        )

    if not has_3d:
        img = FITSImage()
        img.load(path)
        return img

    cls = _INSTRUMENT_MAP.get(instrume.upper() if instrume else '', GenericCube)
    cube = cls()
    cube.load(path, var_path=var)
    return cube


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _read_instrume(hdul):
    """Return INSTRUME keyword from whichever header has it, or ''."""
    for hdu in hdul:
        val = hdu.header.get('INSTRUME', None)
        if val:
            return str(val).strip()
    return ''


def _find_flux_hdu(hdul):
    """
    Return (hdu, data) for the first 3D extension found.
    Tries extension 0 first, then extension 1, then all remaining.
    """
    candidates = [hdul[0], hdul[1]] if len(hdul) > 1 else [hdul[0]]
    for hdu in candidates:
        if hdu.data is not None and hdu.data.ndim == 3:
            return hdu, hdu.data
    # Fallback: scan all extensions
    for hdu in hdul:
        if hdu.data is not None and hdu.data.ndim == 3:
            return hdu, hdu.data
    raise ValueError("No 3D extension found in FITS file.")


def _find_variance(hdul, flux_shape):
    """Search extensions for a variance array matching flux_shape."""
    var_hdu = _get_extension(hdul, ('VAR', 'VARIANCE', 'STAT', 'ERR'))
    if var_hdu is not None and var_hdu.data is not None:
        if var_hdu.data.shape == flux_shape:
            return var_hdu.data.astype(np.float64)
    return None


def _get_extension(hdul, names):
    """Return the first HDU whose .name is in *names*, or None."""
    name_set = {n.upper() for n in names}
    for hdu in hdul:
        if hdu.name.upper() in name_set:
            return hdu
    return None


def _build_wave(header):
    """
    Build a 1-D wavelength array (Angstroms) from 3-axis WCS keywords.

    Supports CRVAL3/CDELT3/CRPIX3 and CD3_3 variants.
    Converts from nm or µm to Å if CUNIT3 says so.
    """
    naxis3 = header.get('NAXIS3')
    if naxis3 is None:
        raise ValueError("Header has no NAXIS3 — cannot build wavelength array.")

    crval3 = float(header.get('CRVAL3', 0.0))
    crpix3 = float(header.get('CRPIX3', 1.0))
    cdelt3 = float(header.get('CDELT3', header.get('CD3_3', 1.0)))
    cunit3 = str(header.get('CUNIT3', 'Angstrom')).strip().lower()

    pixels = np.arange(1, naxis3 + 1)  # 1-indexed, FITS convention
    wave = crval3 + cdelt3 * (pixels - crpix3)

    # Unit conversion to Angstroms
    if cunit3 in ('nm', 'nanometer', 'nanometers'):
        wave *= 10.0
    elif cunit3 in ('um', 'micron', 'microns', 'micrometer', 'micrometers'):
        wave *= 1e4

    return wave


def _load_external_var(var_path, flux_shape):
    """
    Load variance from a separate FITS file.

    Searches for the first 3D extension whose shape matches *flux_shape*.
    Raises ValueError if none is found.
    """
    var_path = Path(var_path)
    with fits.open(var_path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.shape == flux_shape:
                return hdu.data.astype(np.float64)
    raise ValueError(
        f"No array with shape {flux_shape} found in variance file '{var_path}'."
    )


def _make_2d_wcs(header):
    """
    Build a 2D spatial astropy WCS from a 3D cube header.
    Uses naxis=2 to read only the first two (spatial) axes directly,
    avoiding the WCS-axis-count mismatch warning.
    """
    try:
        wcs2d = WCS(header, naxis=2)
        return wcs2d
    except Exception:
        return None
