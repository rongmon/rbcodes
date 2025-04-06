#!/usr/bin/env python
"""
Spectrum data model - Core component for representing spectrum data.
"""
import numpy as np
from astropy.convolution import convolve, Box1DKernel

class Spectrum:
    """
    A class to represent a spectrum with wavelength, flux, and error arrays.
    
    This class provides an interface for manipulating spectrum data, including
    smoothing, normalization, and other common operations.
    
    Parameters
    ----------
    wavelength : array-like
        Wavelength array
    flux : array-like
        Flux array
    error : array-like
        Error array
    redshift : float, optional
        Redshift of the spectrum, by default 0.0
    """
    
    def __init__(self, wavelength, flux, error, redshift=0.0):
        """Initialize a Spectrum object with data arrays and redshift."""
        self.wavelength = np.asarray(wavelength)
        self.flux = np.asarray(flux)
        self.error = np.asarray(error)
        self.redshift = redshift
        
        # Create a copy for smooth operations
        self.smoothed_flux = self.flux.copy()
        self.smoothed_error = self.error.copy()
        
        # Track smooth level
        self.smooth_level = 1
        
        # Store initial view limits for reset operations
        self._initial_xlim = (np.min(self.wavelength), np.max(self.wavelength))
        
        # Calculate a reasonable initial y-limit
        median_flux = np.nanmedian(self.flux)
        self._initial_ylim = (0.0, median_flux*2.5 if median_flux > 0 else 1.0)
    
    @property
    def rest_wavelength(self):
        """Get wavelength array in rest frame."""
        return self.wavelength / (1.0 + self.redshift)
    
    @property
    def initial_xlim(self):
        """Get initial x limits for plotting."""
        return self._initial_xlim
    
    @property
    def initial_ylim(self):
        """Get initial y limits for plotting."""
        return self._initial_ylim
    
    def set_redshift(self, redshift):
        """
        Set a new redshift value.
        
        Parameters
        ----------
        redshift : float
            New redshift value
        """
        self.redshift = redshift
    
    def smooth(self, increment=True):
        """
        Apply smoothing to the spectrum.
        
        Parameters
        ----------
        increment : bool, optional
            If True, increase smoothing level; if False, decrease
            
        Returns
        -------
        tuple
            Smoothed flux and error arrays
        """
        if increment:
            self.smooth_level += 2
        else:
            self.smooth_level = max(1, self.smooth_level - 2)
        
        # Apply box filter smoothing
        kernel_size = int(self.smooth_level)
        self.smoothed_flux = convolve(self.flux, Box1DKernel(kernel_size))
        self.smoothed_error = convolve(self.error, Box1DKernel(kernel_size))
        
        return self.smoothed_flux, self.smoothed_error
    
    def reset_smooth(self):
        """Reset smoothing to original data."""
        self.smooth_level = 1
        self.smoothed_flux = self.flux.copy()
        self.smoothed_error = self.error.copy()
        
        return self.smoothed_flux, self.smoothed_error
    
    def compute_ew(self, wave_limits, flux_limits, use_rest_frame=True):
        """
        Compute equivalent width between two wavelength points.
        
        Parameters
        ----------
        wave_limits : array-like
            Two wavelength points defining the region
        flux_limits : array-like
            Two flux points defining the continuum level
        use_rest_frame : bool, optional
            If True, treat wavelengths as rest frame, by default True
            
        Returns
        -------
        dict
            Results including EW, error, continuum, and wavelength slice
        """
        # Sort limits if they're not in ascending order
        sorted_indices = np.argsort(wave_limits)
        wave_limits = np.array(wave_limits)[sorted_indices]
        flux_limits = np.array(flux_limits)[sorted_indices]
        
        # Determine which wavelength array to use
        if use_rest_frame:
            wave = self.rest_wavelength
            flux = self.flux
            error = self.error
        else:
            wave = self.wavelength
            flux = self.flux
            error = self.error
        
        # Select data within the specified wavelength range
        mask = (wave >= wave_limits[0]) & (wave <= wave_limits[1])
        wave_slice = wave[mask]
        flux_slice = flux[mask]
        error_slice = error[mask]
        
        # Interpolate continuum
        from scipy.interpolate import splrep, splev
        spline = splrep(wave_limits, flux_limits, k=1)
        continuum = splev(wave_slice, spline)
        
        # Normalize flux and error by continuum
        norm_flux = flux_slice / continuum
        norm_error = error_slice / continuum
        
        # Handle NaN values
        nan_mask = np.isnan(norm_flux)
        if np.any(nan_mask):
            norm_flux[nan_mask] = norm_error[nan_mask]
        
        # Clip flux to avoid infinite optical depth
        sat_lim = -0.01
        sat_mask = norm_flux <= sat_lim
        if np.any(sat_mask):
            norm_flux[sat_mask] = norm_error[sat_mask]
            
        zero_mask = norm_flux <= 0.0
        if np.any(zero_mask):
            norm_flux[zero_mask] = norm_error[zero_mask] + 0.01
        
        # Calculate equivalent width
        ew = np.trapz(1.0 - norm_flux, x=wave_slice)
        
        # Calculate error on equivalent width
        delta_wave = np.median(np.diff(wave_slice))
        sig_w = delta_wave * norm_error
        sig_ew = np.sqrt(np.sum(sig_w**2.0))
        
        # Prepare output data
        if not use_rest_frame:
            wave_slice = wave_slice * (1.0 + self.redshift)
            
        return {
            'ew': ew,
            'error': sig_ew,
            'continuum': continuum,
            'wave_slice': wave_slice * (1 + self.redshift) if use_rest_frame else wave_slice,
            'norm_flux': norm_flux
        }
    
    def compute_column_density(self, wave_limits, flux_limits, f_value, rest_wavelength):
        """
        Compute column density for an absorption line.
        
        Parameters
        ----------
        wave_limits : array-like
            Two wavelength points defining the region
        flux_limits : array-like
            Two flux points defining the continuum level
        f_value : float
            Oscillator strength of the transition
        rest_wavelength : float
            Rest wavelength of the transition
            
        Returns
        -------
        dict
            Results including column density, error, EW, and other parameters
        """
        # First calculate equivalent width
        ew_result = self.compute_ew(wave_limits, flux_limits)
        
        # Convert wavelength to velocity space
        c = 2.9979e5  # Speed of light in km/s
        observed_wavelength = rest_wavelength * (1.0 + self.redshift)
        vel = (ew_result['wave_slice'] - observed_wavelength) * c / observed_wavelength
        
        # Calculate optical depth
        tau_a = np.log(1.0 / ew_result['norm_flux'])
        
        # Calculate median optical depth weighted velocity
        tau_50 = np.cumsum(tau_a) / np.max(tau_a) if np.max(tau_a) > 0 else np.zeros_like(tau_a)
        vel_50 = np.interp(0.5, tau_50, vel) if np.any(tau_a > 0) else 0.0
        
        # Calculate velocity differences
        delta_vel = np.diff(vel)
        delta_vel = np.append([delta_vel[0]], delta_vel)
        
        # Column density calculation
        lambda_r = rest_wavelength
        nv = tau_a / ((2.654e-15) * f_value * lambda_r)  # cm^-2 / (km s^-1)
        n = nv * delta_vel  # column density per bin
        
        # Error calculation
        tau_err = ew_result['error'] / ew_result['norm_flux']
        n_err = (tau_err / ((2.654e-15) * f_value * lambda_r)) * delta_vel
        
        # Total column density and error
        col_density = np.sum(n)
        col_density_err = np.sqrt(np.sum(n_err**2))
        
        return {
            'column_density': col_density,
            'error': col_density_err,
            'ew': ew_result['ew'] * 1000,  # Convert to mÅ
            'ew_error': ew_result['error'] * 1000,  # Convert to mÅ
            'continuum': ew_result['continuum'],
            'wave_slice': ew_result['wave_slice'],
            'vel_50': vel_50
        }