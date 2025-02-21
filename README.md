rbcodes
-------

This is a public release of python codes commonly used for astrophysical reserach by Rongmon Bordoloi.
This package is constantly under development and will be periodically updated. 



Installation:
-------------
Add this folder to your pythonpath and you are good to go.
To make sure the GUIs function properly we recommend that you create a new environemnt as shown in requirements_simple.txt

	Dependencies:  astropy, lmfit, scipy, numpy, matplotlib, linetools
	Partial Dependencies: PysimpleGUI for some GUIs. 


Contents:
---------

The GUIs folder contains several graphical user interfaces for astrophysics. 

	1) rb_cont.py :  A simple interactive continuum fitter. 

	2) rb_spec.py :  An absorption line analysis pipeline. Allowing to create an object from an 1D 
	                 spectrum do continuum fitting, measure equivalenth widths/column densities and
	                 save the final data and measurements. Also allows the fitting of a simple Voigt profile. 
	                 NOTE: Major overhaul to calling sequence since April 2022. See example in the file.

	3) rb_interactive_vpfit_singlet.py : Interactive simple Voigt profile fitter [used with rb_spec.py]
	                                    [ Look for rbvfit github repository for an interactive Voigt 
	                                    profile fitter, suited for any complex profiles.]

	4) rb_plot_spec.py:  To plot 1D spectrum, smooth, pan, zoom and explore spectrum. Allows to plot
	                     absorption lines at any redshifts, and do simple equivalent width and Gaussian
	                     fitting. My version of x_specplot routine. Now updated to plot multiple absorbers
	                     and load previously identified absorbers. 
			 
	5) abstools: A complex absorption line analysis GUI to simultaneously continuum fit and 
		     measure EW/logN of several absorption lines. This is under GUIs/abstools. 
		     Look up readme file in that folder for a working example.
		     
	6) rb_specgui: Wrapper for a complex 1d spectrum viwer and line identifier [Plotspec_integrated.py]. 
	               This is the upgraded pyqt5 version of rb_plot_spec package.
        
	7) zgui/main.py : This is a new GUI to measure redshifts of galaxies using 1d or 2d spectrum. 
			     Can be used in batch mode and is optimized for JWST NIRCam/Grism spectroscopy.



The IGM folder contains several python routines for intergalactic medium and circumgalactic medium calculations.

	1) compute_EW.py :  Compute equivalent width and column densities of any absorption/emission lines.

	2) rb_setline.py :  For any given rest frame wavelength find corresponding atomic transition and 
	                    the fvalues. [Can find any line nearest to the input wavelength].

	3) rb_iter_contfit.py : Iteratively fit [with sigma clipping] continuum to a small slice of 1D
	                        spectrum using Legendre polynomials.

	4) rb_specbin.py  :  Rebin 1D spectrum to increase S/N

	5) ransac_contfit.py : Powerful hybrid continuum fitter. Initial fit is done by an automated ransac 
	                       algorithm, then allows user to tweak the fit using linetools.fit_continuum 
	                       routine. Saves the fitted continuum and the 1D spectrum as a xspectrum1D object
	                       Very helpful to fit a full quasar continuum.
			       
	6) lens_sep_to_kpc.py :Function to compute physical separation between sightlines in a lensed quasar system under the thin lens approximation.

The halo folder contains modules to compute NFW halo profiles

       1) rb_nfw.py      : Compute NFW halo profile

       2) mstar2mhalo.py : Convert stellar mass to Halo mass using Moster et al. 2010.
       
The stat folder contains several statistical codes:

       1) rb_wilsonscore.py :  Wilsonscore confidence intervals.

       2) rb_boot.py        :  Bootstrap function
 
The lensing folder contains code to ray trace positions to source plane for a given deflection matrix.

       1) lens_ang_sep.py  :        Inititate the delensing object to load the required deflection matrices

The utils folder contains several utility modules:

      1) rb_utility.py     :  Several utility functions

      2) rb_x1d_id.py      : This function will read all HST/COS x1d fits files in a folder and print out
                             the header info for raw spectra. It will print filename, exposure type, 
                             and object type entries from the header.

      3) readmultispec.py  : Read IRAF (echelle) spectrum in multispec format from a FITS file. Can read 
                             most multispec formats including linear, log, cubic spline, Chebyshev or 
                             Legendre dispersion spectra. I got this code from https://github.com/kgullikson88/General.

     4) cos_inspec.py      : Routine to quickly inspect HST/COS x1d files.

     5) filter_2d_spec.py  : Custom routine to mask and filter JWST/NIRCam 2D grism spectra from the EIGER survey. [Optimized for a custom data format]

     6) compute_SNR_1d.py  : Compute signal-to-noise ratio of an 1D spectrum and plot the result.
			    
The catalog folder contains codes for catalog search:

	1)rb_search.py.    : Function to do a cone search around any (ra,dec) pointing with respect to a list of ra,dec entries.
	2) convert_FIRE_coordinates.py:  Custom code only to take Magellan 1D FIRE spectrum and use the header information
                                         to transform the co-ordinates to J2000 epoch.
					 
We are constantly developing...

	
DOI
---
[![DOI](https://zenodo.org/badge/192408573.svg)](https://zenodo.org/badge/latestdoi/192408573)




