import numpy as np
import matplotlib.pyplot as plt
from linetools.spectra.xspectrum1d import XSpectrum1D
import emcee
import corner
from scipy.optimize import curve_fit
from tqdm import tqdm
from astropy.stats import sigma_clip

class LLSFitter:
    """
    Class for measuring the column density of Lyman Limit Systems (LLS)
    using both curve_fit and MCMC methods.
    """
    
    def __init__(self, spectrum_file=None, zabs=None):
        """
        Initialize the LLSFitter object.
        
        Parameters:
        -----------
        spectrum_file : str, optional
            Path to the FITS file containing spectrum
        zabs : float, optional
            Absorption redshift
        """
        self.spectrum = None
        self.zabs = None
        self.wave = None
        self.flux = None
        self.error = None
        self.rest_wave = None
        self.continuum_regions = None
        self.domain_range = [800, 1100]  # Default domain range in Angstroms
        self.sigma_clip = 3.0  # Default sigma value for outlier rejection
        
        # Default fitting parameters
        self.theta_init = np.array([1.0, 0.0, 17.0])  # [C0, C1, logNHI]
        self.bounds_lower = np.array([-10.0, -10.0, 14.0])
        self.bounds_upper = np.array([10.0, 10.0, 19.0])
        
        # Results storage
        self.curve_fit_results = None
        self.mcmc_results = None
        self.samples = None
        
        # If spectrum_file and zabs are provided, load and process immediately
        if spectrum_file is not None and zabs is not None:
            self.load_spectrum(spectrum_file)
            self.set_redshift(zabs)
    
    def load_spectrum(self, spectrum_file):
        """
        Load spectrum from FITS file.
        
        Parameters:
        -----------
        spectrum_file : str
            Path to the FITS file
        """
        self.spectrum = XSpectrum1D.from_file(spectrum_file)
        
        # Extract data and normalize
        self.wave = self.spectrum.wavelength.value
        c = np.nanmedian(self.spectrum.flux.value)
        self.flux = self.spectrum.flux.value / c
        self.error = self.spectrum.sig.value / c
        
        return self
    
    def set_redshift(self, zabs):
        """
        Set the absorption redshift and compute rest-frame wavelengths.
        
        Parameters:
        -----------
        zabs : float
            Absorption redshift
        """
        self.zabs = zabs
        if self.wave is not None:
            self.rest_wave = self.wave / (1.0 + self.zabs)
        
        return self
    def set_sigma_clip(self, sigma):
        """
        Set the sigma threshold for outlier rejection when creating continuum mask.
        
        Parameters:
        -----------
        sigma : float
            Sigma threshold for outlier rejection
        """
        self.sigma_clip = sigma
        return self
    
    def set_domain_range(self, wmin=800, wmax=1100):
        """
        Set the wavelength domain range for analysis.
        
        Parameters:
        -----------
        wmin : float, optional
            Minimum rest-frame wavelength in Angstroms
        wmax : float, optional
            Maximum rest-frame wavelength in Angstroms
        """
        self.domain_range = [wmin, wmax]
        return self
    
    def set_continuum_regions(self, regions=None):
        """
        Set the regions used for continuum determination and fitting.
        
        Parameters:
        -----------
        regions : list of tuples, optional
            List of (min, max) wavelength ranges to use for continuum fitting
            If None, use default regions
        """
        if regions is None:
            # Default regions from original code
            self.continuum_regions = [
                (860, 872), (893, 899), (897, 899), (905, 910),
                (940, 944), (946, 948), (927, 928), (918.5, 919),
                (931, 933), (934, 936), (951, 970)
            ]
        else:
            self.continuum_regions = regions
        
        return self
    
    def get_continuum_mask(self, sigma=None):
        """
        Create mask for continuum regions and apply sigma clipping to avoid absorption lines.
        
        Parameters:
        -----------
        sigma : float, optional
            Sigma threshold for outlier rejection
        
        Returns:
        --------
        mask : numpy array
            Boolean mask indicating continuum points after sigma clipping
        """
        
        if sigma is None:
            sigma = self.sigma_clip

        if self.continuum_regions is None:
            self.set_continuum_regions()

        
        # Initialize final mask
        final_mask = np.zeros_like(self.rest_wave, dtype=bool)
        
        # Process each continuum region separately
        for rmin, rmax in self.continuum_regions:
            # Create mask for this region
            region_mask = (self.rest_wave >= rmin) & (self.rest_wave <= rmax)
            
            # Skip if no points in this region
            if not np.any(region_mask):
                continue
                
            # Apply sigma clipping to this region only
            flux_in_region = self.flux[region_mask]
            clipped_data = sigma_clip(flux_in_region, sigma=sigma, maxiters=5)
            
            # Create mask for non-clipped points in this region
            region_final_mask = np.zeros_like(region_mask)
            region_final_mask[region_mask] = ~clipped_data.mask
            
            # Add to the final mask
            final_mask = final_mask | region_final_mask
        
        return final_mask
    
    # --- Models ---
    def model_flx(self, theta, wave):
        """
        Physical model with a continuum and a Lyman limit absorption.
        
        Parameters:
        -----------
        theta : array-like
            Model parameters [C0, C1, logNHI]
        wave : array-like
            Wavelength array in Angstroms
            
        Returns:
        --------
        model : array-like
            Model flux
        """
        c = theta[0] + theta[1] * (wave - 911.0)
        NHI = 10.0 ** theta[2]
        Tau = NHI * 6.3e-18 * (wave / 912.0) ** 2.75

        model = np.copy(c)
        q = np.where(wave <= 912)
        model[q] = c[q] * np.exp(-Tau[q])
        return model
    
    def model_test(self, wave, *params):
        """
        Wrapper for model_flx used in curve fitting (like scipy.optimize.curve_fit).
        """
        theta = np.array(params)
        return self.model_flx(theta, wave)
    
    # --- Likelihood Components for MCMC ---
    def lnprior(self, theta, lb, ub):
        """Flat prior within given bounds."""
        if np.any(theta < lb) or np.any(theta > ub):
            return -np.inf
        return 0.0

    def lnlike(self, theta, x, y, yerr):
        """Log-likelihood function."""
        model = self.model_flx(theta, x)
        inv_sigma2 = 1.0 / (yerr ** 2)
        return -0.5 * np.sum((y - model) ** 2 * inv_sigma2 + np.log(2 * np.pi * yerr**2))

    def lnprob(self, theta, lb, ub, x, y, yerr):
        """Log-probability function."""
        lp = self.lnprior(theta, lb, ub)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.lnlike(theta, x, y, yerr)
    
    def fit_curve_fit(self, theta_init=None):
        """
        Fit the LLS model using scipy's curve_fit.
        
        Parameters:
        -----------
        theta_init : array-like, optional
            Initial parameter guess [C0, C1, logNHI]
            
        Returns:
        --------
        popt : array-like
            Best-fit parameters
        pcov : array-like
            Covariance matrix
        """
        if theta_init is None:
            theta_init = self.theta_init
        
        # Get continuum mask
        mask = self.get_continuum_mask()
        
        # Set up bounds
        bounds = (self.bounds_lower, self.bounds_upper)
        
        # Perform fit
        popt, pcov = curve_fit(
            self.model_test, 
            self.rest_wave[mask], 
            self.flux[mask], 
            p0=theta_init,
            sigma=self.error[mask],
            bounds=bounds
        )
        
        # Store results
        self.curve_fit_results = {
            'parameters': popt,
            'covariance': pcov,
            'errors': np.sqrt(np.diag(pcov))
        }
        
        return popt, pcov
    
    def fit_emcee(self, nwalkers=50, nsteps=500, burnin_frac=0.2, theta_init=None):
        """
        Fit the LLS model using MCMC with emcee.
        
        Parameters:
        -----------
        nwalkers : int, optional
            Number of walkers
        nsteps : int, optional
            Number of steps per walker
        burnin_frac : float, optional
            Fraction of steps to discard as burn-in
        theta_init : array-like, optional
            Initial parameter guess [C0, C1, logNHI]
            If None, use the results from curve_fit if available
            
        Returns:
        --------
        sampler : emcee.EnsembleSampler
            MCMC sampler object
        samples : array-like
            Flattened chain of samples
        """
        # If no initial parameters provided, use curve_fit results if available
        if theta_init is None:
            if self.curve_fit_results is not None:
                theta_init = self.curve_fit_results['parameters']
            else:
                theta_init = self.theta_init
                
        # Get continuum mask
        mask = self.get_continuum_mask()
        
        # Setup for MCMC
        ndim = len(theta_init)
        
        # Generate initial positions for walkers
        pos = [theta_init + 1.e-6 * np.random.randn(ndim) for _ in range(nwalkers)]
        
        # Calculate burn-in steps
        burnin = int(round(nsteps * burnin_frac))
        
        # Initialize sampler
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, self.lnprob,
            args=(self.bounds_lower, self.bounds_upper, 
                  self.rest_wave[mask], self.flux[mask], self.error[mask])
        )
        
        # Run burn-in
        print(f"🔥 Burn-in Phase ({burnin} steps)...")
        pos = sampler.run_mcmc(pos, burnin, progress=True)
        sampler.reset()
        print("✅ Burn-in complete.\n")
        
        # Run production
        print(f"🚀 Main Sampling Phase ({nsteps} steps)...")
        for _ in tqdm(sampler.sample(pos, iterations=nsteps), 
                     total=nsteps, desc="Sampling", ncols=80):
            pass
        print("✅ Sampling complete.")
        
        # Extract and store samples (discarding first 20% as burn-in)
        discard = int(nsteps * 0.2)
        self.samples = sampler.chain[:, discard:, :].reshape((-1, ndim))
        
        # Calculate median and standard deviation of parameters
        params_median = np.median(self.samples, axis=0)
        params_std = np.std(self.samples, axis=0)
        
        # Store results
        self.mcmc_results = {
            'parameters': params_median,
            'errors': params_std,
            'sampler': sampler,
            'samples': self.samples
        }
        
        return sampler, self.samples
    
    def plot_fit(self, method='curve_fit', show_continuum_regions=True, 
               wmin=None, wmax=None, ymin=None, ymax=None, figsize=(12, 6), 
               show_realizations=False, n_realizations=100):
        """
        Plot the spectrum and the fit.
        
        Parameters:
        -----------
        method : str, optional
            'curve_fit' or 'mcmc'
        show_continuum_regions : bool, optional
            Whether to highlight the continuum regions with gray boxes
            and plot the points used in fitting
        wmin, wmax : float, optional
            Wavelength range to plot (if None, use domain_range)
        ymin, ymax : float, optional
            Flux range to plot (if None, auto-determine from data)
        figsize : tuple, optional
            Figure size
        show_realizations : bool, optional
            Whether to show random realizations (only for MCMC)
        n_realizations : int, optional
            Number of random realizations to show
            
        Returns:
        --------
        fig, ax : matplotlib figure and axis
        """
        if wmin is None:
            wmin = self.domain_range[0]
        if wmax is None:
            wmax = self.domain_range[1]
        
        # Get parameters based on method
        if method == 'curve_fit':
            if self.curve_fit_results is None:
                raise ValueError("Run fit_curve_fit first")
            params = self.curve_fit_results['parameters']
            title = f"Curve Fit: logNHI = {params[2]:.2f} ± {self.curve_fit_results['errors'][2]:.2f}"
        elif method == 'mcmc':
            if self.mcmc_results is None:
                raise ValueError("Run fit_emcee first")
            params = self.mcmc_results['parameters']
            title = f"MCMC Fit: logNHI = {params[2]:.2f} ± {self.mcmc_results['errors'][2]:.2f}"
        else:
            raise ValueError("Method must be 'curve_fit' or 'mcmc'")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot continuum regions if requested - do this first so they're in the background
        if show_continuum_regions:
            # Get continuum mask - these are the points actually used in fitting
            continuum_mask = self.get_continuum_mask()
            
            # Plot each continuum region as a gray box
            y_box_min = 0
            y_box_max = 2  # Default range, will be adjusted later
                
            for rmin, rmax in self.continuum_regions:
                # Only show regions within the plot range
                if (rmin <= wmax and rmax >= wmin):
                    rect = plt.Rectangle((rmin, y_box_min), rmax-rmin, y_box_max-y_box_min, 
                                        color='gray', alpha=0.15, zorder=0)
                    ax.add_patch(rect)
        
        # Plot MCMC realizations if requested
        if method == 'mcmc' and show_realizations and self.samples is not None:
            # Generate random indices for the realizations
            n_samples = len(self.samples)
            indices = np.random.randint(0, n_samples, n_realizations)
            
            # Plot each realization
            for idx in indices:
                params_sample = self.samples[idx]
                model_realization = self.model_flx(params_sample, self.rest_wave)
                ax.plot(self.rest_wave, model_realization, color='lightgreen', alpha=0.1, zorder=1)
        
        # Plot data
        ax.step(self.rest_wave, self.flux, 'b-', alpha=0.7, where='mid', label='Spectrum', zorder=2)
        
        # Generate and plot model
        model_flux = self.model_flx(params, self.rest_wave)
        ax.plot(self.rest_wave, model_flux, 'g-', linewidth=2, label='Model Fit', zorder=4)
        
        # Add vertical line at Lyman limit
        ax.axvline(x=912, color='k', linestyle='--', alpha=0.7, label='Lyman Limit (912 Å)', zorder=3)
        
        # Plot continuum points used (do this last so they're on top)
        if show_continuum_regions:
            ax.plot(self.rest_wave[continuum_mask], self.flux[continuum_mask], 'r.', 
                    alpha=0.7, markersize=3, label='Continuum Points Used', zorder=5)
        
        # Set axes limits for x-axis
        ax.set_xlim(wmin, wmax)
        
        # Set axes limits for y-axis - either use provided values or auto-calculate
        if ymin is None or ymax is None:
            # Auto-calculate y limits from visible data
            visible_flux = self.flux[(self.rest_wave >= wmin) & (self.rest_wave <= wmax)]
            auto_ymin = np.percentile(visible_flux, 1)
            auto_ymax = np.percentile(visible_flux, 99)
            margin = 0.2 * (auto_ymax - auto_ymin)
            
            # Use provided values if available, otherwise use auto-calculated ones
            y_lower = ymin if ymin is not None else max(0, auto_ymin - margin)
            y_upper = ymax if ymax is not None else auto_ymax + margin
        else:
            y_lower = ymin
            y_upper = ymax
        
        # Update the rectangle heights for continuum regions
        if show_continuum_regions:
            for patch in ax.patches:
                if isinstance(patch, plt.Rectangle):
                    patch.set_height(y_upper - y_lower)
                    patch.set_y(y_lower)
        
        ax.set_ylim(y_lower, y_upper)
        
        # Add labels and legend
        ax.set_xlabel('Rest Wavelength (Å)')
        ax.set_ylabel('Normalized Flux')
        ax.set_title(title)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig, ax    
    def plot_corner(self, figsize=(10, 10)):
        """
        Plot corner plot of MCMC samples.
        
        Returns:
        --------
        fig : matplotlib figure
        """
        if self.mcmc_results is None:
            raise ValueError("Run fit_emcee first")
        
        # Create corner plot
        fig = corner.corner(
            self.mcmc_results['samples'], 
            labels=["C0", "C1", "logNHI"],
            truths=self.mcmc_results['parameters'],
            show_titles=True,
            title_kwargs={"fontsize": 14},
            figsize=figsize
        )
        
        return fig
    
    def get_results_summary(self):
        """
        Get a summary of the fitting results.
        
        Returns:
        --------
        dict : Dictionary with fitting results
        """
        results = {}
        
        if self.curve_fit_results is not None:
            params = self.curve_fit_results['parameters']
            errors = self.curve_fit_results['errors']
            results['curve_fit'] = {
                'C0': params[0],
                'C0_err': errors[0],
                'C1': params[1],
                'C1_err': errors[1],
                'logNHI': params[2],
                'logNHI_err': errors[2]
            }
        
        if self.mcmc_results is not None:
            params = self.mcmc_results['parameters']
            errors = self.mcmc_results['errors']
            results['mcmc'] = {
                'C0': params[0],
                'C0_err': errors[0],
                'C1': params[1],
                'C1_err': errors[1],
                'logNHI': params[2],
                'logNHI_err': errors[2]
            }
        
        return results


# Example usage
if __name__ == "__main__":
    # Create LLS fitter
    lls = LLSFitter('J1154+4635_nbin3_coadd.fits', zabs=0.52859638)
    #set sigma clip to 5 sigma
    lls.set_sigma_clip(25)
    
    # Set continuum regions (using default regions)
    lls.set_continuum_regions()
    
    # Fit with curve_fit
    popt, pcov = lls.fit_curve_fit()
    print("Curve fit results:")
    print(f"logNHI = {popt[2]:.2f} ± {np.sqrt(pcov[2,2]):.2f}")
    
    # Plot the fit
    fig, ax = lls.plot_fit(method='curve_fit', wmin=880, wmax=975)
    plt.savefig('lls_curvefit.png')
    plt.close(fig)
    
    # Fit with MCMC
    sampler, samples = lls.fit_emcee(nwalkers=50, nsteps=500)
    
    # Plot the MCMC fit with realizations
    fig, ax = lls.plot_fit(method='mcmc', wmin=880, wmax=975, 
                          show_realizations=True, n_realizations=100)
    plt.savefig('lls_mcmc_with_realizations.png')
    plt.close(fig)
    
    # Plot corner plot
    fig = lls.plot_corner()
    plt.savefig('lls_corner_plot.png')
    plt.close(fig)
    
    # Print summary
    results = lls.get_results_summary()
    print("\nSummary of results:")
    print(f"Curve fit: logNHI = {results['curve_fit']['logNHI']:.2f} ± {results['curve_fit']['logNHI_err']:.2f}")
    print(f"MCMC: logNHI = {results['mcmc']['logNHI']:.2f} ± {results['mcmc']['logNHI_err']:.2f}")