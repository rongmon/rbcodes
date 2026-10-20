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

    def load(self, path, var_path=None, wave_range=None, air_to_vacuum=False):
        path = Path(path)
        self.path = str(path)
        self.name = path.stem

        with fits.open(path) as hdul:
            flux_hdu, _ = _find_flux_hdu(hdul)
            self._read(hdul, flux_hdu, ('VAR', 'VARIANCE', 'STAT', 'ERR'), wave_range)

        if var_path is not None:
            self.var = _load_external_var(var_path, self._full_shape, self._channels)
        _finish(self, air_to_vacuum)
        return self

    def _read(self, hdul, flux_hdu, var_names, wave_range):
        """Read flux (and variance) — only the channels in wave_range if given."""
        header = flux_hdu.header.copy()
        full_wave = _build_wave(header)
        sl = _channel_slice(full_wave, wave_range)
        self._channels = sl
        self._full_shape = flux_hdu.shape
        self.flux = np.asarray(flux_hdu.section[sl], dtype=np.float64)
        self.var = None
        var_hdu = _get_extension(hdul, var_names)
        if var_hdu is not None and var_hdu.header.get('NAXIS', 0) == 3 and var_hdu.shape == flux_hdu.shape:
            self.var = np.asarray(var_hdu.section[sl], dtype=np.float64)
        self.wave = full_wave[sl]
        if sl.start or sl.stop != len(full_wave):
            header['NAXIS3'] = self.flux.shape[0]
            header['CRPIX3'] = float(header.get('CRPIX3', 1.0)) - sl.start
            header.add_history(f'rbcodes: channels {sl.start}-{sl.stop - 1} of {len(full_wave)} read')
        self.header = header
        self.wcs = _make_2d_wcs(header)


# ---------------------------------------------------------------------------
# KCWI subclass
# ---------------------------------------------------------------------------

class KCWICube(GenericCube):
    """KCWI-specific loader. Currently uses GenericCube logic."""

    def load(self, path, var_path=None, wave_range=None, air_to_vacuum=False):
        super().load(path, var_path=var_path, wave_range=wave_range, air_to_vacuum=air_to_vacuum)
        # KCWI stores flux in e-/s; variance in (e-/s)^2 — no unit conversion needed
        return self


# ---------------------------------------------------------------------------
# MUSE subclass
# ---------------------------------------------------------------------------

class MUSECube(GenericCube):
    """MUSE-specific loader. MUSE stores flux in DATA ext, variance in STAT."""

    def load(self, path, var_path=None, wave_range=None, air_to_vacuum=False):
        path = Path(path)
        self.path = str(path)
        self.name = path.stem

        with fits.open(path) as hdul:
            # MUSE primary is empty; data is in 'DATA'
            if 'DATA' in [h.name for h in hdul]:
                flux_hdu = hdul['DATA']
            else:
                flux_hdu, _ = _find_flux_hdu(hdul)
            self._read(hdul, flux_hdu, ('STAT', 'VAR', 'VARIANCE', 'ERR'), wave_range)

        if var_path is not None:
            self.var = _load_external_var(var_path, self._full_shape, self._channels)
        _finish(self, air_to_vacuum)
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

def load_fits(path, var=None, wave_range=None, air_to_vacuum=False):
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
    wave_range : (lo, hi) in Angstrom, optional
        Read only the channels with lo <= wavelength <= hi (in the file's own air/vacuum
        convention). Default: the whole cube. Useful for large cubes: only those channels
        are read from disk.
    air_to_vacuum : bool
        Convert the wavelength axis from air to vacuum (``rbcodes.utils.rb_spectrum.air_to_vacuum``)
        when the cube is in air (CTYPE3 = 'AWAV'). A cube already in vacuum (CTYPE3 = 'WAVE') is
        left unchanged with a warning. ``cube.airvac`` records the result ('air', 'vac' or 'unknown').

    Returns
    -------
    IFUCube subclass instance (GenericCube, KCWICube, MUSECube, …)
    """
    from .image2d import FITSImage
    path = Path(path)

    # Peek at the headers to decide cube vs. 2D image (without reading the data)
    with fits.open(path) as hdul:
        instrume = _read_instrume(hdul)
        has_3d = any(hdu.header.get('NAXIS', 0) == 3 for hdu in hdul)

    if not has_3d:
        img = FITSImage()
        img.load(path)
        return img

    cls = _INSTRUMENT_MAP.get(instrume.upper() if instrume else '', GenericCube)
    cube = cls()
    cube.load(path, var_path=var, wave_range=wave_range, air_to_vacuum=air_to_vacuum)
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
    Return (hdu, None) for the first 3D extension found (the data are not read here).
    Tries extension 0 first, then extension 1, then all remaining.
    """
    candidates = [hdul[0], hdul[1]] if len(hdul) > 1 else [hdul[0]]
    for hdu in candidates:
        if hdu.header.get('NAXIS', 0) == 3:
            return hdu, None
    # Fallback: scan all extensions
    for hdu in hdul:
        if hdu.header.get('NAXIS', 0) == 3:
            return hdu, None
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


def _load_external_var(var_path, flux_shape, channels=slice(None)):
    """
    Load variance from a separate FITS file (channels ``channels`` only).

    Searches for the first 3D extension whose shape matches *flux_shape* (the full cube).
    Raises ValueError if none is found.
    """
    var_path = Path(var_path)
    with fits.open(var_path) as hdul:
        for hdu in hdul:
            if hdu.header.get('NAXIS', 0) == 3 and hdu.shape == tuple(flux_shape):
                return np.asarray(hdu.section[channels], dtype=np.float64)
    raise ValueError(
        f"No array with shape {flux_shape} found in variance file '{var_path}'."
    )


def _channel_slice(wave, wave_range):
    """Slice of the channels with lo <= wave <= hi (all channels if wave_range is None)."""
    if wave_range is None:
        return slice(0, len(wave))
    lo, hi = wave_range
    idx = np.flatnonzero((wave >= lo) & (wave <= hi))
    if len(idx) == 0:
        raise ValueError(f'No channels in wave_range {wave_range}; the cube covers '
                         f'{wave.min():.1f}-{wave.max():.1f} A.')
    return slice(int(idx[0]), int(idx[-1]) + 1)


def _finish(cube, air_to_vacuum):
    """Record the air/vacuum convention and convert to vacuum if asked."""
    ctype = str(cube.header.get('CTYPE3', '')).strip().upper()
    cube.airvac = {'AWAV': 'air', 'WAVE': 'vac'}.get(ctype, 'unknown')
    if not air_to_vacuum:
        return
    if cube.airvac == 'vac':
        import warnings
        warnings.warn('air_to_vacuum: the cube is already in vacuum (CTYPE3 = WAVE); not converted.')
        return
    if cube.airvac == 'unknown':
        import warnings
        warnings.warn('air_to_vacuum: CTYPE3 does not say air or vacuum; converting as asked.')
    from rbcodes.utils.rb_spectrum import air_to_vacuum as _a2v
    cube.wave = _a2v(cube.wave)
    cube.airvac = 'vac'
    cube.header['HISTORY'] = 'rbcodes: wavelengths converted air -> vacuum (cube.wave); header WCS is still air'



def _make_2d_wcs(header):
    """
    Build a 2D spatial astropy WCS from a 3D cube header.
    Uses naxis=2 to read only the first two (spatial) axes directly,
    avoiding the WCS-axis-count mismatch warning.

    Falls back to manually constructing a WCS from scalar CRPIX/CRVAL/CDELT
    keywords when the header uses a mixed CD/CDELT convention that confuses
    wcslib (e.g. CDELT1/2 for spatial + CD3_3 for wavelength, no CD1_1 etc.).
    """
    try:
        wcs2d = WCS(header, naxis=2)
        return wcs2d
    except Exception:
        pass
    # Fallback: build a clean 2D WCS from individual header keywords
    try:
        w = WCS(naxis=2)
        w.wcs.crpix = [header['CRPIX1'], header['CRPIX2']]
        w.wcs.crval = [header['CRVAL1'], header['CRVAL2']]
        w.wcs.cdelt = [header['CDELT1'], header['CDELT2']]
        w.wcs.ctype = [header.get('CTYPE1', 'RA---TAN'),
                       header.get('CTYPE2', 'DEC--TAN')]
        if 'CROTA2' in header:
            w.wcs.crota = [0.0, header['CROTA2']]
        w.wcs.set()
        return w
    except Exception:
        return None
