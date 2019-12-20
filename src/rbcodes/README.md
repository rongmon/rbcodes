# rbcodes
This is a public release of my python codes
Add this folder to your pythonpath and you are good to go.

Dependencies:  astropy, lmfit, scipy, numpy, matplotlib

The GUIs folder contains several graphical user interfaces for astrophysics. 

	1) rb_cont.py :  A simple interactive continuum fitter. 

	2) rb_spec.py :  An absorption line analysis pipeline. Allowing to create an object from an 1D spectrum do continuum fitting, measure equivalenth widths/column densities and save the final data and measurements. Also allows the fitting of a simple Voigt profile. 

	3) rb_interactive_vpfit_singlet.py : Interactive simple Voigt profile fitter [used with rb_spec.py]
	[ Look for rbvfit github repository for an interactive Voigt profile fitter, suited for any complex profiles.]

	4) rb_plot_spec.py:  To plot 1D spectrum, smooth, pan, zoom and explore spectrum. Allows to plot absorption lines at any redshifts, and do simple equivalent width and Gaussian fitting. My version of x_specplot routine. 


The IGM folder contains several python routines for intergalactic medium and circumgalactic medium calculations.

	1) compute_EW.py :  Compute equivalent width and column densities of any absorption/emission lines.
	2) rb_setline.py :  For any given rest frame wavelength find corresponding atomic transition and the fvalues. [Can find any line nearest to the input wavelength].
	3) rb_iter_contfit.py : Iteratively fit [with sigma clipping] continuum to a small slice of 1D spectrum using Legendre polynomials.
	4) rb_specbin.py  :  Rebin 1D spectrum to increase S/N

The halo folder contains modules to compute NFW halo profiles

       1) rb_nfw.py      : Compute NFW halo profile
       2) mstar2mhalo.py : Convert stellar mass to Halo mass using Moster et al. 2010.
The stat folder contains several statistical codes:

       1) rb_wilsonscore.py :  Wilsonscore confidence intervals.
       2) rb_boot.py .      :  Bootstram function

The utils folder contains several utility modules:

      1) rb_utility.py .  :  Several utility functions
	



