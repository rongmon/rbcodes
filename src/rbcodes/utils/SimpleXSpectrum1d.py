import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.table import Table
import warnings
import os


class SimpleXSpectrum1D:
    """A simplified replacement for linetools.XSpectrum1D.
    
    This class handles 1D spectra with wavelength, flux, and error arrays.
    It provides basic functionality for loading, manipulating, and analyzing spectra.
    
    Parameters
    ----------
    wave : array-like
        Wavelength array
    flux : array-like
        Flux array
    sig : array-like, optional
        Error array
    mask : array-like, optional
        Mask array (1=masked)
    meta : dict, optional
        Metadata dictionary
    units : dict, optional
        Dictionary with units (e.g., {'wave': 'Angstrom', 'flux': 'erg/s/cm^2/Angstrom'})
    """
    
    def __init__(self, wave, flux, sig=None, mask=None, meta=None, units=None, continuum=None, co_is_set=False):
        # Convert to numpy arrays
        self.wavelength = np.array(wave)
        self.flux = np.array(flux)
        
        # Error array - default to zeros if not provided
        if sig is not None:
            self.sig = np.array(sig)
        else:
            self.sig = np.zeros_like(self.flux)
        
        # Mask array - default to zeros (no masking) if not provided
        if mask is not None:
            self.mask = np.array(mask, dtype=bool)
        else:
            self.mask = np.zeros_like(self.flux, dtype=bool)
        
        # Metadata
        if meta is None:
            self.meta = {}
        else:
            self.meta = meta
            
        # Units
        default_units = {'wave': 'Angstrom', 'flux': 'erg/s/cm^2/Angstrom'}
        if units is None:
            self.units = default_units
        else:
            self.units = units
        
        # Continuum
        self.co_is_set = co_is_set
        if continuum is not None:
            self.continuum = np.array(continuum)
            self.co_is_set = True
            
            # Check continuum array size
            if len(self.continuum) != len(self.wavelength):
                raise ValueError("Continuum array dimensions do not match wavelength array")
        elif co_is_set:
            # If co_is_set is True but no continuum is provided, initialize with ones
            self.continuum = np.ones_like(self.flux)
        else:
            self.continuum = None
            
        # Check array sizes
        if not (len(self.wavelength) == len(self.flux) == len(self.sig) == len(self.mask)):
            raise ValueError("Array dimensions do not match")
    
    @property
    def npix(self):
        """Number of pixels in the spectrum"""
        return len(self.wavelength)
    
    @classmethod
    def from_file(cls, filename, **kwargs):
        """Read spectrum from file based on file extension
        
        Parameters
        ----------
        filename : str
            File to read
        
        Returns
        -------
        spec : SimpleXSpectrum1D
            Spectrum object
        
        Notes
        -----
        - FITS files should have extensions with WAVE, FLUX, and optionally ERR, MASK
        - ASCII files should have columns for wavelength, flux, and optionally error
        """
        # Determine file type from extension
        if filename.endswith(('.fits', '.fit', '.FITS', '.FIT')):
            return cls.from_fits(filename, **kwargs)
        elif filename.endswith(('.txt', '.dat', '.ascii')):
            return cls.from_ascii(filename, **kwargs)
        else:
            raise ValueError(f"Unrecognized file format: {filename}")
    
    @classmethod
    def from_fits(cls, filename, ext=1, wave_key='WAVE', flux_key='FLUX', 
                  sig_key='ERR', mask_key=None, cont_key='CONTINUUM'):
        """Read spectrum from a FITS file
        
        Parameters
        ----------
        filename : str
            FITS file name
        ext : int, optional
            FITS extension number
        wave_key, flux_key, sig_key, mask_key, cont_key : str, optional
            Keys for wavelength, flux, error, mask, and continuum arrays
        
        Returns
        -------
        spec : SimpleXSpectrum1D
            Spectrum object
        """
        with fits.open(filename) as hdul:
            header = hdul[ext].header
            
            # Try different ways to extract wavelength
            try:
                wave = hdul[ext].data[wave_key]
            except (KeyError, TypeError):
                try:
                    # Try as table column
                    wave = hdul[ext].data.field(wave_key)
                except:
                    # Try to reconstruct from header WCS
                    naxis1 = header.get('NAXIS1')
                    crval1 = header.get('CRVAL1')
                    cdelt1 = header.get('CDELT1', header.get('CD1_1'))
                    if None in (naxis1, crval1, cdelt1):
                        raise ValueError("Could not extract wavelength from header")
                    wave = crval1 + cdelt1 * np.arange(naxis1)
            
            # Get flux
            try:
                flux = hdul[ext].data[flux_key]
            except (KeyError, TypeError):
                try:
                    flux = hdul[ext].data.field(flux_key)
                except:
                    raise ValueError(f"Could not find flux array with key {flux_key}")
            
            # Get error (optional)
            sig = None
            if sig_key:
                try:
                    sig = hdul[ext].data[sig_key]
                except (KeyError, TypeError):
                    try:
                        sig = hdul[ext].data.field(sig_key)
                    except:
                        warnings.warn(f"Could not find error array with key {sig_key}")
            
            # Get mask (optional)
            mask = None
            if mask_key:
                try:
                    mask = hdul[ext].data[mask_key]
                except (KeyError, TypeError):
                    try:
                        mask = hdul[ext].data.field(mask_key)
                    except:
                        warnings.warn(f"Could not find mask array with key {mask_key}")
            
            # Get continuum (optional)
            continuum = None
            co_is_set = header.get('CO_IS_SET', False)
            if cont_key:
                try:
                    continuum = hdul[ext].data[cont_key]
                    co_is_set = True
                except (KeyError, TypeError):
                    try:
                        continuum = hdul[ext].data.field(cont_key)
                        co_is_set = True
                    except:
                        if co_is_set:
                            warnings.warn(f"CO_IS_SET flag is True but could not find continuum array with key {cont_key}")
            
            # Create metadata from header
            meta = {'header': header}
            
            # Try to determine units
            units = {'wave': header.get('CUNIT1', 'Angstrom'),
                     'flux': header.get('CUNIT2', 'erg/s/cm^2/Angstrom')}
            
            return cls(wave, flux, sig=sig, mask=mask, meta=meta, units=units, 
                      continuum=continuum, co_is_set=co_is_set)
    
    @classmethod
    def from_ascii(cls, filename, wave_col=0, flux_col=1, sig_col=2, 
                   comment='#', **kwargs):
        """Read spectrum from an ASCII file
        
        Parameters
        ----------
        filename : str
            ASCII file name
        wave_col, flux_col, sig_col : int, optional
            Column indices for wavelength, flux, and error
        comment : str, optional
            Comment character
        **kwargs : 
            Additional arguments passed to np.loadtxt
        
        Returns
        -------
        spec : SimpleXSpectrum1D
            Spectrum object
        """
        # Read data
        try:
            data = np.loadtxt(filename, comments=comment, **kwargs)
        except Exception as e:
            raise IOError(f"Error reading ASCII file: {str(e)}")
        
        # Extract columns
        try:
            wave = data[:, wave_col]
            flux = data[:, flux_col]
        except IndexError:
            raise ValueError("Could not extract wavelength or flux columns")
        
        # Error array (optional)
        try:
            sig = data[:, sig_col]
        except IndexError:
            sig = None
        
        # Create metadata
        meta = {'filename': os.path.basename(filename)}
        
        return cls(wave, flux, sig=sig, meta=meta)
    
    def rebin(self, new_wave, do_sig=True):
        """Rebin the spectrum to a new wavelength grid
        
        Parameters
        ----------
        new_wave : array-like
            New wavelength array
        do_sig : bool, optional
            Whether to rebin the error array
        
        Returns
        -------
        spec : SimpleXSpectrum1D
            Rebinned spectrum
        """
        from scipy.interpolate import interp1d
        
        # Skip if the wavelength grids are identical
        if np.array_equal(self.wavelength, new_wave):
            return self
        
        # Create interpolation function
        f = interp1d(self.wavelength, self.flux, bounds_error=False, fill_value=np.nan)
        new_flux = f(new_wave)
        
        # Handle error array
        if do_sig and self.sig is not None:
            f_sig = interp1d(self.wavelength, self.sig, bounds_error=False, fill_value=np.nan)
            new_sig = f_sig(new_wave)
        else:
            new_sig = None
            
        # Handle continuum array if available
        new_continuum = None
        if self.co_is_set and self.continuum is not None:
            f_co = interp1d(self.wavelength, self.continuum, bounds_error=False, fill_value=np.nan)
            new_continuum = f_co(new_wave)
        
        # Create new spectrum
        return self.__class__(new_wave, new_flux, sig=new_sig, meta=self.meta.copy(), 
                             units=self.units, continuum=new_continuum, co_is_set=self.co_is_set)
    
    def write(self, filename, format=None, overwrite=False):
        """Write spectrum to file
        
        Parameters
        ----------
        filename : str
            Output file name
        format : str, optional
            Format ('fits' or 'ascii'). If None, determined from file extension.
        overwrite : bool, optional
            Whether to overwrite existing file
        """
        # Determine format from extension if not specified
        if format is None:
            if filename.endswith(('.fits', '.fit', '.FITS', '.FIT')):
                format = 'fits'
            elif filename.endswith(('.txt', '.dat', '.ascii')):
                format = 'ascii'
            else:
                format = 'fits'  # Default to FITS
        
        # Write according to format
        if format.lower() == 'fits':
            self._write_fits(filename, overwrite=overwrite)
        elif format.lower() == 'ascii':
            self._write_ascii(filename, overwrite=overwrite)
        else:
            raise ValueError(f"Unrecognized format: {format}")
    
    def _write_fits(self, filename, overwrite=False):
        """Write spectrum to FITS file"""
        # Create a table for the data
        t = Table()
        t['WAVE'] = self.wavelength
        t['FLUX'] = self.flux
        if self.sig is not None:
            t['ERR'] = self.sig
        if self.mask is not None:
            t['MASK'] = self.mask.astype(int)
        if self.co_is_set and self.continuum is not None:
            t['CONTINUUM'] = self.continuum
        
        # Convert to HDU
        hdu = fits.BinTableHDU(t)
        
        # Add metadata to header
        for key, value in self.meta.items():
            if key != 'header':  # Avoid header recursion
                if isinstance(value, (str, int, float, bool)):
                    hdu.header[key] = value
        
        # Add continuum flag to header
        hdu.header['CO_IS_SET'] = self.co_is_set
        
        # Add unit information
        if 'wave' in self.units:
            hdu.header['CUNIT1'] = self.units['wave']
        if 'flux' in self.units:
            hdu.header['CUNIT2'] = self.units['flux']
        
        # Create HDU list and write
        hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
        hdul.writeto(filename, overwrite=overwrite)
    
    def _write_ascii(self, filename, overwrite=False):
        """Write spectrum to ASCII file"""
        # Check if file exists and overwrite flag
        if os.path.exists(filename) and not overwrite:
            raise IOError(f"File exists: {filename}. Use overwrite=True to overwrite.")
        
        # Create data array
        if self.sig is not None:
            data = np.column_stack((self.wavelength, self.flux, self.sig))
            header = f"# WAVE({self.units.get('wave', 'Angstrom')}) FLUX({self.units.get('flux', 'erg/s/cm^2/Angstrom')}) ERROR"
        else:
            data = np.column_stack((self.wavelength, self.flux))
            header = f"# WAVE({self.units.get('wave', 'Angstrom')}) FLUX({self.units.get('flux', 'erg/s/cm^2/Angstrom')})"
        
        # Write to file
        np.savetxt(filename, data, header=header)
    
    def normalize(self, window_size=5, iterations=3, sigma_lower=3.0, sigma_upper=3.0):
        """Normalize the spectrum using sigma-clipping
        
        Parameters
        ----------
        window_size : int, optional
            Window size for the median filter (in pixels)
        iterations : int, optional
            Number of sigma-clipping iterations
        sigma_lower, sigma_upper : float, optional
            Lower and upper sigma clipping thresholds
        
        Returns
        -------
        spec : SimpleXSpectrum1D
            Normalized spectrum
        """
        # If continuum is already set, use it
        if self.co_is_set and self.continuum is not None:
            continuum = self.continuum
        else:
            from scipy.signal import medfilt
            from astropy.stats import sigma_clip
            
            # Copy flux for manipulation
            flux_work = self.flux.copy()
            
            # Iterative sigma-clipping to find continuum
            for _ in range(iterations):
                # Apply median filter
                filtered = medfilt(flux_work, window_size)
                
                # Compute residuals
                residuals = flux_work - filtered
                
                # Sigma clip
                mask = sigma_clip(residuals, sigma_lower=sigma_lower, sigma_upper=sigma_upper).mask
                
                # Replace masked values with filtered values
                flux_work[mask] = filtered[mask]
            
            # Final median filter is our continuum
            continuum = medfilt(flux_work, window_size)
        
        # Normalize
        norm_flux = self.flux / continuum
        
        # Also normalize errors if present
        if self.sig is not None:
            norm_sig = self.sig / continuum
        else:
            norm_sig = None
        
        # Create new metadata with normalization info
        new_meta = self.meta.copy()
        new_meta['normalized'] = True
        
        # Update units
        new_units = self.units.copy()
        new_units['flux'] = 'normalized'
        if 'sig' in new_units:
            new_units['sig'] = 'normalized'
        
        return self.__class__(self.wavelength, norm_flux, sig=norm_sig, 
                              mask=self.mask, meta=new_meta, units=new_units, 
                              continuum=continuum, co_is_set=True)
    
    def select_region(self, wvmin, wvmax):
        """Select a wavelength region of the spectrum
        
        Parameters
        ----------
        wvmin, wvmax : float
            Minimum and maximum wavelength
        
        Returns
        -------
        spec : SimpleXSpectrum1D
            Spectrum in the selected region
        """
        # Create mask for region
        region_mask = (self.wavelength >= wvmin) & (self.wavelength <= wvmax)
        
        # Select continuum if available
        continuum = None
        if self.co_is_set and self.continuum is not None:
            continuum = self.continuum[region_mask]
        
        # Return subset
        return self.__class__(self.wavelength[region_mask], self.flux[region_mask],
                             sig=None if self.sig is None else self.sig[region_mask],
                             mask=None if self.mask is None else self.mask[region_mask],
                             meta=self.meta.copy(), units=self.units,
                             continuum=continuum, co_is_set=self.co_is_set)
    
    def __add__(self, other):
        """Add two spectra or add a constant to flux"""
        if isinstance(other, (int, float)):
            new_flux = self.flux + other
            new_sig = self.sig  # Error doesn't change when adding constant
            return self.__class__(self.wavelength, new_flux, sig=new_sig, 
                                  mask=self.mask, meta=self.meta.copy(), units=self.units)
        elif isinstance(other, SimpleXSpectrum1D):
            # Check if wavelength grids are identical
            if not np.array_equal(self.wavelength, other.wavelength):
                # Rebin the other spectrum to match this one
                other = other.rebin(self.wavelength)
            
            # Add fluxes
            new_flux = self.flux + other.flux
            
            # Combine errors in quadrature if available
            if self.sig is not None and other.sig is not None:
                new_sig = np.sqrt(self.sig**2 + other.sig**2)
            else:
                new_sig = None
            
            # Combine masks (regions masked in either spectrum)
            new_mask = self.mask | other.mask
            
            return self.__class__(self.wavelength, new_flux, sig=new_sig,
                                  mask=new_mask, meta=self.meta.copy(), units=self.units)
        else:
            raise TypeError(f"Unsupported operand type: {type(other)}")
    
    def __mul__(self, other):
        """Multiply spectrum by a constant or another spectrum"""
        if isinstance(other, (int, float)):
            new_flux = self.flux * other
            new_sig = None if self.sig is None else self.sig * abs(other)
            return self.__class__(self.wavelength, new_flux, sig=new_sig,
                                  mask=self.mask, meta=self.meta.copy(), units=self.units)
        elif isinstance(other, SimpleXSpectrum1D):
            # Check if wavelength grids are identical
            if not np.array_equal(self.wavelength, other.wavelength):
                # Rebin the other spectrum to match this one
                other = other.rebin(self.wavelength)
            
            # Multiply fluxes
            new_flux = self.flux * other.flux
            
            # Propagate errors if available
            if self.sig is not None and other.sig is not None:
                # Error propagation for multiplication: (f*g)² * ((df/f)² + (dg/g)²)
                rel_err1 = self.sig / self.flux
                rel_err2 = other.sig / other.flux
                # Handle division by zero
                rel_err1[~np.isfinite(rel_err1)] = 0
                rel_err2[~np.isfinite(rel_err2)] = 0
                
                new_sig = new_flux * np.sqrt(rel_err1**2 + rel_err2**2)
            else:
                new_sig = None
            
            # Combine masks
            new_mask = self.mask | other.mask
            
            return self.__class__(self.wavelength, new_flux, sig=new_sig,
                                  mask=new_mask, meta=self.meta.copy(), units=self.units)
        else:
            raise TypeError(f"Unsupported operand type: {type(other)}")
    
    def __rmul__(self, other):
        """Right multiplication"""
        return self.__mul__(other)
    
    def compute_equivalent_width(self, wvmin, wvmax):
        """Compute the equivalent width of a line
        
        Parameters
        ----------
        wvmin, wvmax : float
            Wavelength range to integrate over
        
        Returns
        -------
        ew : float
            Equivalent width
        ew_err : float
            Error in equivalent width
        """
        # Select region
        region = self.select_region(wvmin, wvmax)
        
        # Check if this is a normalized spectrum or has a continuum set
        is_normalized = self.meta.get('normalized', False)
        
        # Wavelength step sizes
        dwave = np.diff(region.wavelength)
        # Add the last step to make lengths match
        dwave = np.append(dwave, dwave[-1])
        
        # Use proper continuum if available or spectrum is normalized
        if region.co_is_set and region.continuum is not None:
            # Use the stored continuum
            ew = np.sum((1.0 - region.flux/region.continuum) * dwave)
            
            # Error calculation
            if region.sig is not None:
                ew_err = np.sqrt(np.sum((region.sig/region.continuum * dwave)**2))
            else:
                ew_err = np.nan
        elif is_normalized:
            # For normalized spectra, EW = sum((1-flux)*dwave)
            ew = np.sum((1.0 - region.flux) * dwave)
            
            # Error calculation if available
            if region.sig is not None:
                ew_err = np.sqrt(np.sum((region.sig * dwave)**2))
            else:
                ew_err = np.nan
        else:
            # Warn if not normalized and no continuum
            warnings.warn("Computing EW on non-normalized spectrum without continuum may give incorrect results")
            
            # For non-normalized spectra without continuum, need to estimate continuum
            # This is a very simple linear continuum estimate
            left_idx, right_idx = 0, -1
            cont_level = (region.flux[left_idx] + region.flux[right_idx]) / 2.0
            
            # Calculate EW
            ew = np.sum((1.0 - region.flux/cont_level) * dwave)
            
            # Error calculation (simplified)
            if region.sig is not None:
                ew_err = np.sqrt(np.sum((region.sig/cont_level * dwave)**2))
            else:
                ew_err = np.nan
        
        return ew, ew_err
    
    def smooth(self, window_size=5, method='boxcar'):
        """Smooth the spectrum
        
        Parameters
        ----------
        window_size : int, optional
            Size of the smoothing window
        method : str, optional
            Smoothing method ('boxcar', 'gaussian', or 'median')
        
        Returns
        -------
        spec : SimpleXSpectrum1D
            Smoothed spectrum
        """
        if method == 'boxcar':
            kernel = np.ones(window_size) / window_size
            new_flux = np.convolve(self.flux, kernel, mode='same')
        elif method == 'gaussian':
            from scipy.ndimage import gaussian_filter
            new_flux = gaussian_filter(self.flux, sigma=window_size/2.355)
        elif method == 'median':
            from scipy.signal import medfilt
            new_flux = medfilt(self.flux, kernel_size=window_size)
        else:
            raise ValueError(f"Unknown smoothing method: {method}")
        
        # Create new metadata with smoothing info
        new_meta = self.meta.copy()
        new_meta['smoothed'] = True
        new_meta['smooth_method'] = method
        new_meta['smooth_window'] = window_size
        
        return self.__class__(self.wavelength, new_flux, sig=self.sig,
                              mask=self.mask, meta=new_meta, units=self.units)