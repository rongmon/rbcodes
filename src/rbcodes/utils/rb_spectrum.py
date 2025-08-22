"""
Simple spectrum reader/writer for rbcodes
A lightweight replacement for linetools.XSpectrum1D focused on multispec GUI needs
"""

import numpy as np
import warnings
import os
import json
import traceback
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.units import Quantity, UnitBase

# Optional dependencies
try:
    import h5py
    HDF5_AVAILABLE = True
except ImportError:
    HDF5_AVAILABLE = False


class rb_spectrum:
    """
    A simple 1D spectrum class for reading/writing spectral data
    
    Attributes
    ----------
    wavelength : Quantity
        Wavelength array with units
    flux : Quantity  
        Flux array with units
    sig : Quantity or None
        Error array with same units as flux (None if not available)
    co : Quantity or None
        Continuum array with same units as flux (None if not available)
    filename : str
        Source filename
    meta : dict
        Metadata dictionary
    _read_failed : bool
        Flag indicating if reading failed (for fallback to linetools)
    """
    
    def __init__(self, wavelength, flux, error=None, continuum=None, 
                 filename=None, meta=None, units=None, _read_failed=False):
        """
        Initialize rb_spectrum
        
        Parameters
        ----------
        wavelength : array or Quantity
            Wavelength array
        flux : array or Quantity
            Flux array  
        error : array or Quantity, optional
            Error array
        continuum : array or Quantity, optional
            Continuum array
        filename : str, optional
            Source filename
        meta : dict, optional
            Metadata dictionary
        units : dict, optional
            Units dictionary with 'wave' and 'flux' keys
        _read_failed : bool, optional
            Flag indicating if reading failed
        """
        self._read_failed = _read_failed
        
        # If read failed, create minimal dummy object
        if _read_failed:
            self.wavelength = Quantity([0, 1], unit=u.AA)
            self.flux = Quantity([0, 0], unit=u.dimensionless_unscaled)
            self.sig = None
            self.co = None
            self.filename = filename or 'failed'
            self.meta = {'read_error': True}
            self.units = {'wave': u.AA, 'flux': u.dimensionless_unscaled}
            return
        
        # Handle units
        if units is None:
            units = {'wave': u.AA, 'flux': u.dimensionless_unscaled}
            
        # Convert to Quantity arrays
        self.wavelength = self._ensure_quantity(wavelength, units['wave'])
        self.flux = self._ensure_quantity(flux, units['flux'])
        
        if error is not None:
            self.sig = self._ensure_quantity(error, units['flux'])
        else:
            self.sig = None
            
        if continuum is not None:
            self.co = self._ensure_quantity(continuum, units['flux'])
        else:
            self.co = None
            
        self.filename = filename or 'none'
        self.meta = meta or {'airvac': 'vac'}
        self.units = units
        
    def _ensure_quantity(self, array, unit):
        """Convert array to Quantity with given unit"""
        if hasattr(array, 'unit'):
            return array
        else:
            return Quantity(array, unit=unit)
    
    @property
    def sig_is_set(self):
        """Check if error array is available"""
        return self.sig is not None
        
    @property
    def co_is_set(self):
        """Check if continuum array is available"""
        return self.co is not None
        
    @property
    def npix(self):
        """Number of pixels"""
        return len(self.wavelength)
        
    @property
    def wvmin(self):
        """Minimum wavelength"""
        return np.min(self.wavelength)
        
    @property
    def wvmax(self):
        """Maximum wavelength"""
        return np.max(self.wavelength)
    
    @classmethod
    def from_file(cls, filename, **kwargs):
        """
        Load spectrum from file
        
        Parameters
        ----------
        filename : str
            Path to spectrum file
        **kwargs
            Additional arguments passed to readers
            
        Returns
        -------
        rb_spectrum
        """
        return rb_read_spectrum(filename, **kwargs)
    
    @classmethod
    def from_tuple(cls, ituple, sort=True, **kwargs):
        """
        Create spectrum from tuple of arrays
        
        Parameters
        ----------
        ituple : tuple
            (wave, flux), (wave, flux, error), or (wave, flux, error, continuum)
        sort : bool, optional
            Sort by wavelength
        **kwargs
            Additional arguments
            
        Returns
        -------
        rb_spectrum
        """
        try:
            if not isinstance(ituple, tuple):
                raise ValueError("Input must be a tuple")
                
            if len(ituple) < 2:
                raise ValueError("Need at least wavelength and flux")
                
            # Parse wavelength
            wave = ituple[0]
            if hasattr(wave, 'unit'):
                wv_unit = wave.unit
                wave_vals = wave.value
            else:
                wv_unit = u.AA
                wave_vals = wave
                warnings.warn("Assuming wavelength unit is Angstroms")
                
            # Parse flux
            flux = ituple[1]
            if hasattr(flux, 'unit'):
                fx_unit = flux.unit
                flux_vals = flux.value
            else:
                fx_unit = u.dimensionless_unscaled
                flux_vals = flux
                
            # Parse error
            error = ituple[2] if len(ituple) > 2 else None
            error_vals = error.value if hasattr(error, 'value') else error
            
            # Parse continuum
            continuum = ituple[3] if len(ituple) > 3 else None
            continuum_vals = continuum.value if hasattr(continuum, 'value') else continuum
            
            # Sort if requested
            if sort:
                srt = np.argsort(wave_vals)
                wave_vals = wave_vals[srt]
                flux_vals = flux_vals[srt]
                if error_vals is not None:
                    error_vals = error_vals[srt]
                if continuum_vals is not None:
                    continuum_vals = continuum_vals[srt]
                    
            # Create spectrum
            units = {'wave': wv_unit, 'flux': fx_unit}
            return cls(wave_vals, flux_vals, error_vals, continuum_vals, 
                      units=units, **kwargs)
        
        except Exception as e:
            print(f"rb_spectrum.from_tuple failed: {str(e)}")
            return cls(None, None, _read_failed=True, **kwargs)
    
    def copy(self):
        """Create a copy of the spectrum"""
        if self._read_failed:
            return rb_spectrum(None, None, _read_failed=True, filename=self.filename)
            
        return rb_spectrum(
            self.wavelength.copy(),
            self.flux.copy(), 
            self.sig.copy() if self.sig_is_set else None,
            self.co.copy() if self.co_is_set else None,
            filename=self.filename,
            meta=self.meta.copy(),
            units=self.units.copy()
        )


    @classmethod
    def from_arrays(cls, wave_2d, flux_2d, error_2d=None, continuum_2d=None,
                   wave_unit=None, flux_unit=None, filenames=None, **kwargs):
        """
        Create multiple spectra from 2D arrays
        
        Parameters
        ----------
        wave_2d : array-like, shape (nspec, npix)
            Wavelength arrays for multiple spectra
        flux_2d : array-like, shape (nspec, npix)  
            Flux arrays for multiple spectra
        error_2d : array-like, optional
            Error arrays for multiple spectra (same shape as flux_2d)
        continuum_2d : array-like, optional
            Continuum arrays for multiple spectra (same shape as flux_2d)
        wave_unit : str or astropy.units.Unit, optional
            Unit for wavelength (default: Angstrom)
        flux_unit : str or astropy.units.Unit, optional
            Unit for flux (default: dimensionless)
        filenames : list of str, optional
            Filenames for each spectrum
        **kwargs
            Additional arguments
            
        Returns
        -------
        list of rb_spectrum
            List of spectrum objects
        """
        try:
            # Convert to numpy arrays
            wave_2d = np.asarray(wave_2d)
            flux_2d = np.asarray(flux_2d)
            
            if error_2d is not None:
                error_2d = np.asarray(error_2d)
            if continuum_2d is not None:
                continuum_2d = np.asarray(continuum_2d)
                
            # Handle units
            if wave_unit is None:
                wave_unit = u.AA
                warnings.warn("No wavelength unit specified, assuming Angstroms")
            elif isinstance(wave_unit, str):
                wave_unit = getattr(u, wave_unit)
                
            if flux_unit is None:
                flux_unit = u.dimensionless_unscaled
            elif isinstance(flux_unit, str):
                flux_unit = getattr(u, flux_unit)
                
            # Validate arrays
            if wave_2d.ndim != 2 or flux_2d.ndim != 2:
                raise ValueError("Input arrays must be 2D")
                
            if wave_2d.shape != flux_2d.shape:
                raise ValueError("Wave and flux arrays must have same shape")
                
            nspec, npix = wave_2d.shape
            
            # Validate error and continuum shapes if provided
            if error_2d is not None and error_2d.shape != (nspec, npix):
                raise ValueError("Error array shape doesn't match wave/flux arrays")
            if continuum_2d is not None and continuum_2d.shape != (nspec, npix):
                raise ValueError("Continuum array shape doesn't match wave/flux arrays")
                
            # Handle filenames
            if filenames is None:
                filenames = [f'spectrum_{i}' for i in range(nspec)]
            elif len(filenames) != nspec:
                raise ValueError(f"Number of filenames ({len(filenames)}) doesn't match number of spectra ({nspec})")
                
            # Create spectrum objects
            spectra = []
            for i in range(nspec):
                wave_i = wave_2d[i, :]
                flux_i = flux_2d[i, :]
                error_i = error_2d[i, :] if error_2d is not None else None
                continuum_i = continuum_2d[i, :] if continuum_2d is not None else None
                
                # Convert wavelength to Angstroms for consistency
                if wave_unit != u.AA:
                    wave_quantity = wave_i * wave_unit
                    wave_i = wave_quantity.to(u.AA).value
                    internal_wave_unit = u.AA
                else:
                    internal_wave_unit = wave_unit
                    
                internal_units = {'wave': internal_wave_unit, 'flux': flux_unit}
                
                spec = cls(wave_i, flux_i, error_i, continuum_i,
                          filename=filenames[i], units=internal_units, **kwargs)
                spectra.append(spec)
                
            return spectra
            
        except Exception as e:
            print(f"rb_spectrum.from_arrays failed: {str(e)}")
            traceback.print_exc()
            # Return list with failed spectrum objects
            return [cls(None, None, _read_failed=True, filename=f'failed_spectrum_{i}') 
                   for i in range(len(filenames) if filenames else 1)]
    

    @classmethod  
    def append(cls, spectrum_list, sort_by_wavelength=True):
        """
        Combine multiple rb_spectrum objects into a single multi-spectrum container
        
        Parameters
        ----------
        spectrum_list : list of rb_spectrum
            List of spectrum objects to combine
        sort_by_wavelength : bool, optional
            Sort spectra by their minimum wavelength (default: True)
            
        Returns
        -------
        rb_spectrum_collection
            Container object that behaves like a list of spectra
        """
        if not spectrum_list:
            raise ValueError("Cannot append empty spectrum list")
            
        if not all(isinstance(spec, rb_spectrum) for spec in spectrum_list):
            raise ValueError("All items must be rb_spectrum objects")
            
        # Sort by wavelength if requested
        if sort_by_wavelength:
            spectrum_list = sorted(spectrum_list, key=lambda s: s.wvmin.value)
            
        return rb_spectrum_collection(spectrum_list)

    
    def write(self, filename, **kwargs):
        """
        Write spectrum to file (format determined by extension)
        
        Parameters
        ----------
        filename : str
            Output filename
        **kwargs
            Additional arguments passed to writer
        """
        if self._read_failed:
            raise RuntimeError("Cannot write spectrum that failed to read")
            
        ext = filename.split('.')[-1].lower()
        
        if ext in ['fits', 'fit']:
            self.rb_write_fits(filename, **kwargs)
        elif ext == 'hdf5':
            if HDF5_AVAILABLE:
                self.rb_write_hdf5(filename, **kwargs)
            else:
                raise ImportError("h5py not available for HDF5 writing")
        elif ext == 'json':
            self.rb_write_json(filename, **kwargs)
        elif ext in ['txt', 'ascii', 'dat']:
            self.rb_write_ascii(filename, **kwargs)
        else:
            # Default to FITS
            self.rb_write_fits(filename, **kwargs)
    
    def rb_write_fits(self, filename, clobber=True):
        """
        Write to multi-extension FITS file
        
        Parameters
        ----------
        filename : str
            Output filename
        clobber : bool
            Overwrite existing file
        """
        # Primary HDU with flux
        hdu_list = [fits.PrimaryHDU(self.flux.value)]
        hdu_list[0].name = 'FLUX'
        
        # Error extension
        if self.sig_is_set:
            hdu_list.append(fits.ImageHDU(self.sig.value, name='ERROR'))
            
        # Wavelength extension  
        hdu_list.append(fits.ImageHDU(self.wavelength.value, name='WAVELENGTH'))
        
        # Continuum extension
        if self.co_is_set:
            hdu_list.append(fits.ImageHDU(self.co.value, name='CONTINUUM'))
            
        # Add metadata
        hdu_list[0].header['NPIX'] = self.npix
        if self.meta:
            hdu_list[0].header['METADATA'] = json.dumps(self.meta)
        hdu_list[0].header['UNITS'] = json.dumps({
            'wave': str(self.units['wave']),
            'flux': str(self.units['flux'])
        })
        
        # Write
        hdulist = fits.HDUList(hdu_list)
        hdulist.writeto(filename, overwrite=clobber)
        print(f'Wrote spectrum to {filename}')
    
    def rb_write_hdf5(self, filename, clobber=True):
        """Write to HDF5 file"""
        if not HDF5_AVAILABLE:
            raise ImportError("h5py not available")
            
        if not clobber and os.path.exists(filename):
            raise IOError("File exists and clobber=False")
            
        with h5py.File(filename, 'w') as f:
            f.create_dataset('wavelength', data=self.wavelength.value)
            f.create_dataset('flux', data=self.flux.value)
            
            if self.sig_is_set:
                f.create_dataset('error', data=self.sig.value)
            if self.co_is_set:
                f.create_dataset('continuum', data=self.co.value)
                
            # Metadata
            f.attrs['units'] = json.dumps({
                'wave': str(self.units['wave']),
                'flux': str(self.units['flux'])
            })
            if self.meta:
                f.attrs['metadata'] = json.dumps(self.meta)
                
        print(f'Wrote spectrum to {filename}')
    
    def rb_write_json(self, filename, clobber=True):
        """Write to JSON file (rb_spectrum native format)"""
        if not clobber and os.path.exists(filename):
            raise IOError("File exists and clobber=False")
            
        # Prepare data structure
        data = {
            'wavelength': self.wavelength.value.tolist(),
            'flux': self.flux.value.tolist(),
            'units': {
                'wave': str(self.units['wave']),
                'flux': str(self.units['flux'])
            },
            'metadata': self.meta
        }
        
        if self.sig_is_set:
            data['error'] = self.sig.value.tolist()
        if self.co_is_set:
            data['continuum'] = self.co.value.tolist()
            
        # Write to file
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
            
        print(f'Wrote spectrum to {filename}')
    
    def rb_write_ascii(self, filename, format='ascii.ecsv'):
        """Write to ASCII file"""
        from astropy.table import QTable, Column
        
        table = QTable([self.wavelength, self.flux], names=('WAVE', 'FLUX'))
        
        if self.sig_is_set:
            table.add_column(Column(self.sig, name='ERROR'))
        if self.co_is_set:
            table.add_column(Column(self.co, name='CONTINUUM'))
            
        table.write(filename, format=format, overwrite=True)
        print(f'Wrote spectrum to {filename}')
    
    def air2vac(self):
        """Convert wavelengths from air to vacuum"""
        if self._read_failed:
            return
            
        if self.meta.get('airvac') == 'vac':
            warnings.warn("Already in vacuum. Not applying correction")
            return
            
        # Convert to Angstroms
        wave_aa = self.wavelength.to(u.AA).value
        
        # Standard conversion
        sigma_sq = (1e4 / wave_aa)**2
        factor = (1 + (5.792105e-2/(238.0185-sigma_sq)) + 
                 (1.67918e-3/(57.362-sigma_sq)))
        factor = factor * (wave_aa >= 2000.) + 1. * (wave_aa < 2000.)
        
        # Apply conversion
        new_wave = wave_aa * factor * u.AA
        self.wavelength = new_wave.to(self.wavelength.unit)
        self.meta['airvac'] = 'vac'
        
    def vac2air(self):
        """Convert wavelengths from vacuum to air"""
        if self._read_failed:
            return
            
        if self.meta.get('airvac') == 'air':
            warnings.warn("Already in air. Not applying correction")
            return
            
        # Convert to Angstroms
        wave_aa = self.wavelength.to(u.AA).value
        
        # Standard conversion (inverse of air to vac)
        sigma_sq = (1e4 / wave_aa)**2
        factor = (1 + (5.792105e-2/(238.0185-sigma_sq)) + 
                 (1.67918e-3/(57.362-sigma_sq)))
        factor = factor * (wave_aa >= 2000.) + 1. * (wave_aa < 2000.)
        
        # Apply inverse conversion
        new_wave = wave_aa / factor * u.AA
        self.wavelength = new_wave.to(self.wavelength.unit)
        self.meta['airvac'] = 'air'
    
    def __repr__(self):
        if self._read_failed:
            return f'<rb_spectrum: READ FAILED - file={self.filename}>'
        return (f'<rb_spectrum: file={self.filename}, '
                f'npix={self.npix}, '
                f'wvmin={self.wvmin:.1f}, wvmax={self.wvmax:.1f}>')


class rb_spectrum_collection:
    """
    Container for multiple rb_spectrum objects that behaves like a list
    but maintains compatibility with multispec GUI expectations
    """
    
    def __init__(self, spectrum_list):
        """
        Initialize collection
        
        Parameters
        ----------
        spectrum_list : list of rb_spectrum
            List of spectrum objects
        """
        self.spectra = spectrum_list
        self.nspec = len(spectrum_list)
        
    def __len__(self):
        return len(self.spectra)
        
    def __getitem__(self, index):
        return self.spectra[index]
        
    def __iter__(self):
        return iter(self.spectra)
        
    def append(self, spectrum):
        """Add a spectrum to the collection"""
        if not isinstance(spectrum, rb_spectrum):
            raise ValueError("Can only append rb_spectrum objects")
        self.spectra.append(spectrum)
        self.nspec += 1
        
    def extend(self, spectrum_list):
        """Add multiple spectra to the collection"""
        for spec in spectrum_list:
            self.append(spec)
            
    def sort_by_wavelength(self):
        """Sort spectra by minimum wavelength"""
        self.spectra.sort(key=lambda s: s.wvmin.value)
        
    def get_wavelength_range(self):
        """Get overall wavelength range of all spectra"""
        if not self.spectra:
            return None, None
        min_wave = min(spec.wvmin for spec in self.spectra)
        max_wave = max(spec.wvmax for spec in self.spectra)
        return min_wave, max_wave
        
    def __repr__(self):
        return f'<rb_spectrum_collection: {self.nspec} spectra>'


# =============================================================================
# File reading functions (adapted from linetools io.py)
# =============================================================================

def rb_read_spectrum(filename, meta=None, **kwargs):
    """
    Read a spectrum from file with automatic format detection
    
    Parameters
    ----------
    filename : str
        Path to spectrum file
    meta : dict, optional
        User-provided metadata (will override file metadata)
    **kwargs
        Additional arguments
        
    Returns
    -------
    rb_spectrum
        Returns rb_spectrum object with _read_failed=True if reading fails
    """
    try:
        # Check file exists
        if not os.path.exists(filename):
            print(f"rb_spectrum: File not found: {filename}")
            return rb_spectrum(None, None, _read_failed=True, filename=filename)
        
        # Determine file type by extension
        ext = filename.split('.')[-1].lower()
        
        if ext == 'hdf5':
            spectrum = _rb_read_hdf5(filename, **kwargs)
        elif ext == 'json':
            spectrum = _rb_read_json(filename, **kwargs)
        elif ext in ['fits', 'fit']:
            spectrum = _rb_read_fits(filename, **kwargs)
        else:
            # Default to ASCII for unknown extensions
            spectrum = _rb_read_ascii(filename, **kwargs)
        
        # Apply user-provided metadata if specified
        if meta is not None and not spectrum._read_failed:
            # Update spectrum metadata with user values (user takes precedence)
            spectrum.meta.update(meta)
        
        return spectrum
    
    except Exception as e:
        print(f"rb_spectrum failed to read {filename}: {str(e)}")
        print(f"Full traceback: {traceback.format_exc()}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_read_fits(filename, **kwargs):
    """Read FITS file with format detection"""
    try:
        hdulist = fits.open(filename)
        header = hdulist[0].header
        
        # First check for binary table (regardless of NAXIS)
        if (header.get('XTENSION') == 'BINTABLE' or 
            (header.get('NAXIS', 0) == 0 and len(hdulist) > 1 and 
             hasattr(hdulist[1], 'data') and hasattr(hdulist[1].data, 'dtype'))):
            return _rb_parse_fits_binary_table(hdulist, filename=filename, **kwargs)
        
        # Then proceed with NAXIS-based logic for image data
        naxis = header.get('NAXIS', 0)
        
        if naxis == 0:
            # True binary table case
            return _rb_parse_fits_binary_table(hdulist, filename=filename, **kwargs)
        elif naxis == 1:
            if len(hdulist) == 1:
                # Single extension - look for separate error file
                return _rb_parse_single_extension_fits(filename, hdulist, **kwargs)
            else:
                # Multi-extension
                return _rb_parse_multi_extension_fits(hdulist, filename=filename, **kwargs)
        elif naxis == 2:
            # Check for DESI brick format
            if (len(hdulist) >= 3 and hdulist[0].name == 'FLUX' and 
                hdulist[2].name == 'WAVELENGTH'):
                return _rb_parse_desi_brick(hdulist, filename=filename, **kwargs)
            else:
                # SDSS-style format (true 2D image)
                return _rb_parse_sdss_fits(hdulist, filename=filename, **kwargs)
        else:
            raise ValueError(f"Unsupported FITS format: NAXIS={naxis}")
    
    except Exception as e:
        print(f"rb_spectrum FITS parsing failed for {filename}: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_read_json(filename, **kwargs):
    """Read JSON file with auto-detection of rb_spec vs rb_spectrum format"""
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        
        # Auto-detect format based on content
        if 'wave_slice' in data and 'flux_slice' in data:
            # rb_spec analysis format
            return _rb_parse_rbspec_json(data, filename=filename, **kwargs)
        elif 'wavelength' in data and 'flux' in data:
            # rb_spectrum native format  
            return _rb_parse_rbspectrum_json(data, filename=filename, **kwargs)
        else:
            raise ValueError("Unrecognized JSON format. Expected rb_spec or rb_spectrum format.")
    
    except Exception as e:
        print(f"rb_spectrum JSON parsing failed for {filename}: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_read_hdf5(filename, **kwargs):
    """Read HDF5 file"""
    try:
        if not HDF5_AVAILABLE:
            raise ImportError("h5py not available for HDF5 reading")
            
        with h5py.File(filename, 'r') as f:
            wave = f['wavelength'][:]
            flux = f['flux'][:]
            sig = f['error'][:] if 'error' in f else None
            co = f['continuum'][:] if 'continuum' in f else None
            
            # Extract metadata
            meta = {'airvac': 'vac'}
            if 'metadata' in f.attrs:
                meta.update(json.loads(f.attrs['metadata']))
            
            # Extract units
            units = None
            if 'units' in f.attrs:
                units_data = json.loads(f.attrs['units'])
                units = {
                    'wave': getattr(u, units_data.get('wave', 'AA')),
                    'flux': getattr(u, units_data.get('flux', 'dimensionless_unscaled'))
                }
            
            return rb_spectrum(wave, flux, sig, co, filename=filename, 
                              meta=meta, units=units)
    
    except Exception as e:
        print(f"rb_spectrum HDF5 parsing failed for {filename}: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_read_ascii(filename, **kwargs):
    """Read ASCII file"""
    try:
        table = Table.read(filename)
        
        # Handle column names
        if 'WAVE' in table.colnames:
            wave = table['WAVE']
        elif 'WAVELENGTH' in table.colnames:
            wave = table['WAVELENGTH']
        else:
            wave = table[table.colnames[0]]  # First column
            
        if 'FLUX' in table.colnames:
            flux = table['FLUX']
        else:
            flux = table[table.colnames[1]]  # Second column
            
        sig = table['ERROR'] if 'ERROR' in table.colnames else None
        co = table['CONTINUUM'] if 'CONTINUUM' in table.colnames else None
        
        return rb_spectrum(wave, flux, sig, co, filename=filename)
    
    except Exception as e:
        print(f"rb_spectrum ASCII parsing failed for {filename}: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


# =============================================================================
# JSON format parsers
# =============================================================================

def _rb_parse_rbspec_json(data, filename=None, **kwargs):
    """Parse rb_spec analysis JSON format"""
    # Extract core spectral arrays
    wave = np.array(data['wave_slice'])*(1+data['zabs'])
    flux = np.array(data['flux_slice']) 
    error = np.array(data['error_slice']) if 'error_slice' in data else None
    
    # Handle continuum
    if 'cont' in data:
        continuum = np.array(data['cont'])
    elif 'fnorm' in data and 'flux_slice' in data:
        # Reconstruct continuum from normalized flux: cont = flux / fnorm
        fnorm = np.array(data['fnorm'])
        # Avoid division by zero
        continuum = np.zeros_like(flux)
        good = fnorm != 0
        continuum[good] = flux[good] / fnorm[good]
    else:
        continuum = None
    
    # Create metadata from analysis results
    meta = {'airvac': 'vac'}  # Default
    
    # Preserve key analysis results in metadata
    analysis_data = {}
    preserve_keys = ['zabs', 'linelist', 'trans', 'fval', 'trans_wave', 'vmin', 'vmax',
                     'W', 'W_e', 'N', 'N_e', 'logN', 'logN_e', 'vel_centroid', 
                     'vel_disp', 'SNR', 'continuum_fit_params']
    
    for key in preserve_keys:
        if key in data:
            analysis_data[key] = data[key]
    
    if analysis_data:
        meta['rb_spec_analysis'] = analysis_data
    
    # Preserve original metadata if available
    if 'metadata' in data:
        meta['original_metadata'] = data['metadata']
    
    return rb_spectrum(wave, flux, error, continuum, filename=filename, meta=meta)


def _rb_parse_rbspectrum_json(data, filename=None, **kwargs):
    """Parse rb_spectrum native JSON format"""
    wave = np.array(data['wavelength'])
    flux = np.array(data['flux'])
    error = np.array(data['error']) if 'error' in data else None
    continuum = np.array(data['continuum']) if 'continuum' in data else None
    
    # Extract metadata
    meta = data.get('metadata', {'airvac': 'vac'})
    
    # Extract units if available
    units = None
    if 'units' in data:
        units_data = data['units']
        if isinstance(units_data, str):
            units_data = json.loads(units_data)
        units = {
            'wave': getattr(u, units_data.get('wave', 'AA')),
            'flux': getattr(u, units_data.get('flux', 'dimensionless_unscaled'))
        }
    
    return rb_spectrum(wave, flux, error, continuum, filename=filename, 
                      meta=meta, units=units)


# =============================================================================
# FITS format parsers
# =============================================================================

def _rb_parse_fits_binary_table(hdulist, filename=None, **kwargs):
    """Parse binary FITS table with flexible column detection"""
    try:
        table = Table(hdulist[1].data)
        
        # Check for required columns with case insensitivity
        if not any(col in table.colnames for col in ['WAVE', 'wave', 'WAVELENGTH', 'wavelength', 'loglam', 'LOGLAM']):
            raise ValueError(f"Wavelength column not found in FITS file: {filename}")
        if not any(col in table.colnames for col in ['FLUX', 'flux', 'SPEC', 'spec']):
            raise ValueError(f"Flux column not found in FITS file: {filename}")
        
        # Get column names with case insensitivity
        wave_col = None
        for col in ['WAVE', 'wave', 'WAVELENGTH', 'wavelength', 'loglam', 'LOGLAM']:
            if col in table.colnames:
                wave_col = col
                break
                
        flux_col = None  
        for col in ['FLUX', 'flux', 'SPEC', 'spec']:
            if col in table.colnames:
                flux_col = col
                break
                
        error_col = None
        for col in ['ERROR', 'error', 'ERR', 'err', 'SIGMA', 'sigma', 'ivar', 'IVAR']:
            if col in table.colnames:
                error_col = col
                break
        
        # Extract arrays (handles both single-row arrays and multi-row scalars)
        wave = np.array(table[wave_col]).flatten()
        flux = np.array(table[flux_col]).flatten()
        
        # Handle log wavelength
        if wave_col in ['loglam', 'LOGLAM']:
            wave = 10**wave
            
        if error_col:
            error = np.array(table[error_col]).flatten()
        else:
            warnings.warn(f"⚠️ Warning: No error column found. Assuming 10% error.", 
                         category=UserWarning)
            error = 0.1 * np.abs(flux)
        
        # Handle inverse variance if needed
        if error_col in ['ivar', 'IVAR']:
            sig = np.zeros_like(error)
            good = error > 0
            sig[good] = 1.0 / np.sqrt(error[good])
            error = sig
        
        return rb_spectrum(wave, flux, error, filename=filename)
        
    except Exception as e:
        print(f"rb_spectrum binary table parsing failed: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_parse_multi_extension_fits(hdulist, filename=None, **kwargs):
    """Parse multi-extension FITS"""
    try:
        flux = hdulist[0].data
        
        # Look for error in extension 1
        sig = None
        if len(hdulist) > 1 and hdulist[1].data is not None:
            sig = hdulist[1].data
            
        # Look for wavelength in extension 2
        if len(hdulist) > 2 and hdulist[2].data is not None:
            wave = hdulist[2].data
        else:
            # Generate wavelength from header keywords
            wave = _rb_setwave(hdulist[0].header)
            
        # Look for continuum in extension 3
        co = None
        if len(hdulist) > 3 and hdulist[3].data is not None:
            co = hdulist[3].data
            
        return rb_spectrum(wave, flux, sig, co, filename=filename)
    
    except Exception as e:
        print(f"rb_spectrum multi-extension parsing failed: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_parse_single_extension_fits(filename, hdulist, **kwargs):
    """Parse single extension FITS and look for separate error file"""
    try:
        flux = hdulist[0].data
        wave = _rb_setwave(hdulist[0].header)
        
        # Look for error file using common naming conventions
        sig = None
        base = filename.replace('.fits', '').replace('.fit', '')
        error_files = [f"{base}e.fits", f"{base}E.fits", f"{base}_err.fits", 
                       f"{base}.err.fits"]
        
        for err_file in error_files:
            if os.path.exists(err_file):
                try:
                    sig = fits.getdata(err_file)
                    break
                except:
                    continue
                
        return rb_spectrum(wave, flux, sig, filename=filename)
    
    except Exception as e:
        print(f"rb_spectrum single extension parsing failed: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_parse_sdss_fits(hdulist, filename=None, **kwargs):
    """Parse SDSS-style FITS (2D array with flux/error in different rows)"""
    try:
        flux = hdulist[0].data[0, :].flatten()
        sig = hdulist[0].data[2, :].flatten()
        wave = _rb_setwave(hdulist[0].header)
        
        return rb_spectrum(wave, flux, sig, filename=filename)
    
    except Exception as e:
        print(f"rb_spectrum SDSS parsing failed: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


def _rb_parse_desi_brick(hdulist, select=0, filename=None, **kwargs):
    """Parse DESI brick format"""
    try:
        flux = hdulist[0].data[select, :]
        
        # Handle error or inverse variance
        if hdulist[1].name == 'ERROR':
            sig = hdulist[1].data[select, :]
        else:
            # Assume inverse variance
            ivar = hdulist[1].data[select, :]
            sig = np.zeros_like(ivar)
            good = ivar > 0
            sig[good] = 1.0 / np.sqrt(ivar[good])
            
        wave = hdulist[2].data
        
        return rb_spectrum(wave, flux, sig, filename=filename)
    
    except Exception as e:
        print(f"rb_spectrum DESI brick parsing failed: {str(e)}")
        return rb_spectrum(None, None, _read_failed=True, filename=filename)


# =============================================================================
# Utility functions
# =============================================================================

def _rb_get_table_column(tags, hdulist, idx=1):
    """Find column in binary table by trying multiple tag names"""
    try:
        if idx >= len(hdulist):
            return None, None
            
        table = Table(hdulist[idx].data)
        names = set(table.colnames)
        
        for tag in tags:
            if tag in names:
                return table[tag], tag
                
        return None, None
    
    except Exception:
        return None, None


def _rb_setwave(header):
    """Generate wavelength array from FITS header keywords"""
    try:
        # Required keywords
        if 'NAXIS1' not in header:
            raise ValueError("NAXIS1 keyword not found in header")
        if 'CRVAL1' not in header:
            raise ValueError("CRVAL1 keyword not found in header")
            
        npix = header['NAXIS1']
        crpix1 = header.get('CRPIX1', 1.0)  # Default to 1.0 if missing
        crval1 = header['CRVAL1']
        
        # Get wavelength step size - try multiple keywords
        cdelt1 = None
        if 'CDELT1' in header:
            cdelt1 = header['CDELT1']
        elif 'CD1_1' in header:
            cdelt1 = header['CD1_1']  # SDSS style
        else:
            raise ValueError("No wavelength step size found in header (CDELT1 or CD1_1)")
            
        # Check for log-linear wavelength scale
        dc_flag = header.get('DC-FLAG', 0)
        
        # Auto-detect log scale if CDELT1 is very small
        if cdelt1 < 1e-4:
            dc_flag = 1
            warnings.warn("Small CDELT1 detected, assuming log wavelength scale")
            
        # Generate wavelength array
        wave = crval1 + cdelt1 * (np.arange(npix) + 1.0 - crpix1)
        
        if dc_flag == 1:
            wave = 10**wave  # Convert from log wavelength
            
        return wave
    
    except Exception as e:
        print(f"Error generating wavelength array from header: {str(e)}")
        raise


# =============================================================================
# Convenience functions
# =============================================================================

def rb_read_spec(filename, **kwargs):
    """Convenience function for reading spectra"""
    return rb_read_spectrum(filename, **kwargs)


def rb_write_spec(spectrum, filename, **kwargs):
    """Convenience function for writing spectra"""
    spectrum.write(filename, **kwargs)