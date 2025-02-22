""" Spectrum class to read in, analyze and measure absorption lines."""
import numpy as np
from scipy.interpolate import splrep,splev
from numpy.polynomial.legendre import Legendre
import sys
import os
import pdb

import json


def load_rb_spec_object(filename, verbose=True):
    """
    Load an rb_spec object from a JSON file and populate its attributes with precomputed values.

    This function reads a JSON file containing precomputed spectral data, initializes an 
    `rb_spec` object using key spectral arrays, and then dynamically sets all remaining 
    attributes from the JSON data. Lists in the JSON file are converted to `numpy` arrays 
    for consistency and efficient numerical operations.

    Parameters
    ----------
    filename : str
        Path to the JSON file containing the spectral data.
    verbose : bool, optional
        If True, prints a message when loading is complete. Default is True.

    Returns
    -------
    rb_spec
        An instance of the `rb_spec` class with all attributes loaded from the JSON file.
        Returns None if there is an error in loading the JSON file.

    """
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error loading JSON file {filename}: {e}")
        return None

    spec_object = rb_spec.from_data(np.array(data['wave_slice'])*(1.+data['zabs']), np.array(data['flux_slice']), np.array(data['error_slice']))

    # Set each value in json file as an attribute of the object
    for key, value in data.items():
        if isinstance(value, list):
            value = np.array(value)  # Convert lists to numpy arrays
        setattr(spec_object, key, value)

    if verbose:
        print('---Finished loading saved rb_spec object----')

    return spec_object

# Calculate the confidence bounds
def calculate_confidence_bounds(x, model, cov_matrix):
    # Evaluate the Legendre basis functions at the given x values
    P = np.polynomial.legendre.legvander(x, model.degree)

    # Propagate the parameter uncertainties to the fitted values
    y_err = np.sqrt(np.sum((P @ cov_matrix) * P, axis=1))

    return y_err

class rb_spec(object):
    """A spectrum read into a class, spectrum will have following properties.

    Attributes 
    ----------
        wave: wavelength.
        flux: flux.
        error: error
        filename=filename and location
        filetype = False [default] : other options 
                 ascii, fits, HSLA, xfits, p [pickle], temp, and linetools [uses linetools.io routine for this]

        Optional: 
            All only valid for filetype=linetools option
            efil= errorfile [Default None]

    Returns
    -------
        This gives a rb_spec object with following attributes:

        self.zabs= Absorber redshift

        self.wave_slice= sliced observed wavelength vector
        self.flux_slice= sliced observed flux vector
        self.error_slice= sliced velocity spectra 
        self.linelist=. LineList used
        self.velo=  sliced velocity vector
        self.cont = Fitted continuum
        self.fnorm= Normalized flux
        self.enorm= Normalized error
        self.trans=  Name of the Transition
        self.fval= fvalue of transition
        self.trans_wave= rest frame wavelength of transition
        self.vmin=     velocity minimum used for equivalent width calculation
        self.vmax=    velocity maximum used for equivalent width calculation
        self.W=    Rest Frame equivalenth width
        self.W_e=  uncertainty on rest frame equivalent width
        self.logN=  AOD column density
        self.logN_e= AOD column density uncertainty
        self.Tau= Apparant optical depth as a function of velocity
        self.vel_centroid= EW weighted velocity centroid of the absorption line
        self.vel_disp=    1sigma velocity dispersion
        self.vel50_err = error on velocity centroid



        Written : Rongmon Bordoloi      April 2018
        Edit    : Rongmon Bordoloi      September 2018 Changed kwargs to be compatible to python 3   
        Edit    : Rongmon Bordoloi      Aug 2020: added linetools.io.readspec file
        Edit    : Rongmon Bordoloi      April 2021: added velocity centroid estimates
        Edit    : Rongmon Bordoloi      March 2022: Added more continuum fitting methods
        Edit    : Rongmon Bordoloi      April 2022: Added velocity centroid error
        Edit    : Rongmon Bordoloi      April 2022: Added different calling sequence to ingest numpy arrays directly. 
        Edit    : Rongmon Bordoloi      April 2022: Small updates to have all continuum fitting routines working
        Edit    : Rongmon Bordoloi      April 2022: Added plotting sliced spectrum option
        Edit    : Rongmon Bordoloi      September 2024: Added saving as json option, and loading the json dictionary as a rb_spec object.
        
        # WARNING: CALLING SEQUENCE HAS CHANGED SINCE APRIL 2022.
        CAREFULLY LOOK AT THE EXAMPLE BELOW

    Example
    -------
        import numpy as np
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt
        from GUIs.rb_spec import rb_spec as r 
        #List of absorber redshifts
        zabs=[0.511020,1.026311,1.564481]
        transition= 2796.3
        #Which absorber to analyze
        index=0
        filename='Quasar_Spectrum.fits'
        #Read in file
        s=r.from_file(filename,filetype='linetools')

        #-------------------------------------------------------------------------------
        # DETOUR --->
        #ALTERNATIVE 
        #IF YOU WANT TO DIRECTLY INJEST NUMPY ARRAYS DO THE FOLLOWING
        s=r.from_data(wave,flux,error)
        # HERE wave,flux,error are numpy arrays of wavelength,flux and error respectively.
        #-------------------------------------------------------------------------------
                
        #Shift spectra to rest frame
        s.shift_spec(zabs[index]);
        #Velocity window around transition
        xlim=[-1500,1500]
        
        #Slice Spectrum within that window
        s.slice_spec(transition,xlim[0],xlim[1],use_vel=True);
        
        #Fit continuum Mask the regions defined by velocity
        s.fit_continuum(mask=[-200,300,500,1100],domain=xlim,Legendre=3)
        
        #-------------------------------------------------------------------------------
        # DETOUR 1--->
        #Alternative Fit continuum methods.
        #s.fit_continuum_ransac(window=149,mednorm=False)  
        
        #-------------------------------------------------------------------------------
        # DETOUR 2--->
        #Aternate continuum fitting method [interactive]
        s.fit_continuum(Interactive=True)
        #Aternate continuum fitting method [input prefit continuum]
        # Length of prefit continuum array = length of sliced spectrum
        s.fit_continuum(Legendre=False,prefit_cont=cont_arrary)
        #-------------------------------------------------------------------------------
        # DETOUR 3--->
        #Aternate continuum fitting method, using sigma clipping
        s.fit_continuum(domain=xlim,Legendre=3,sigma_clip=True)
        #-------------------------------------------------------------------------------
        
       
        
        #Compute EW
        #Compute equivalent width within a velocity window
        s.compute_EW(transition,vmin=-200.,vmax=360.);
        
        # Saving the analysis
        #--------------------
        # There are two options:
        # First method: save everything as a pickle file [default]
        s.save_slice('outfile.p')




        #---------------------
        # Second method: Saving information as a json file
        s.save_slice('outfile.json',file_format='json')
        


        #-----------------------------------------------------------
        #Loading the saved rb_spec object from the above two methods

        # Loading the rb_spec object back from the saved pickle file
        import pickle
        with open('outfile.p', 'rb') as f:
            # Load the pickled data
            sp_test = pickle.load(f)



        ----------------------------------------------------------
        # Loading the rb_spec object back from the saved json file
        from GUIs.rb_spec import load_rb_spec_object as r_load
        
        f='outfile.json'

        sp_test=r_load(f)






        #-------------------------------------------------------------------------------
        # Additional inspection routines        
        #plot the Full spectrum
        s.plot_spec()
        
        #Plot the sliced spectrum with the fitted continuum
        s.plot_slice()
    
        #plot the sliced and fitted continuum
        #Plot stuff
        plt.subplot(2,1,1)
        plt.step(s.velo,s.flux_slice)
        plt.step(s.velo,s.flux_slice/s.fnorm)
        plt.xlim(xlim)
        plt.subplot(2,1,2)
        plt.step(s.velo,s.fnorm)
        plt.plot([-1500,1500],[1,1],'--')
        plt.xlim(xlim)
        plt.show()
    """
    def __init__(self,wave,flux,error,filename=False):#,filetype=False, efil=None,**kwargs):
        """ creates the spectrum object """
        #print('Initializing rb_spec object for absorption line analysis!')
        self.wave=wave
        self.flux=flux/np.nanmedian(flux)
        self.error=error/np.nanmedian(flux)
        #self.wrest=wave*(1.+0.)
        #self.zabs=0.
        self.filename=filename


    @classmethod
    def from_file(cls,filename,filetype=False,efil=None,**kwargs):

        if filetype==False:
            #Take File Extention and try
            if filename==False:
                if 'wave' in kwargs:
                    wave=kwargs['wave']
                else:
                    raise IOError("Input wavelength array")
                if 'flux' in kwargs:
                    flux=kwargs['flux']
                else:
                    raise IOError("Input flux array")
                if 'error' in kwargs:
                    error=kwargs['error']
                else:
                    raise IOError("Input error array")

            else:
                tt=os.path.splitext(filename)[1]
                if (tt=='txt')| (tt=='dat'):
                    filetype='ascii'
                else:
                    filetype=tt[1:len(tt)] 


        # Read in Files in differet formats
        if filetype=='ascii':
            from astropy.io import ascii
            dat=ascii.read(filename)
            tab=dat.keys()
            wave=np.array(dat[tab[0]])
            flux=np.array(dat[tab[1]])
            if (len(dat.keys())>=3):
                error=dat[tab[2]]
            else:
                error=0.*flux

        elif filetype=='fits':
            from astropy.io import fits
            file=fits.open(filename)#(cwd+'/'+filename)
            dat=file[1].data
            tab=dat.names
            wave=np.array(dat['wave'][0])
            flux=np.array(dat['flux'][0])
            if (len(tab)>=3):
                error=np.array(dat['error'][0])
            else:
                error=0.*flux
        elif filetype=='HSLA':
            from astropy.io import fits
            file=fits.open(filename)#(cwd+'/'+filename)
            dat=file[1].data
            tab=dat.names
            wave=np.array(dat['WAVE'])
            flux=np.array(dat['FLUX'])#/np.median(np.array(dat['FLUX']))
            if (len(tab)>=3):
                error=np.array(dat['ERROR'])#/np.median(np.array(dat['FLUX']))
            else:
                error=0.*flux


        if filetype=='xfits':
            from linetools.spectra.xspectrum1d import XSpectrum1D  
            sp=XSpectrum1D.from_file(filename)

            wave=sp.wavelength.value
            flux=sp.flux.value
            error=sp.sig.value


            if sp.co_is_set == True:
                print('Normalizing spectrum using given continuum...')
                flux=sp.flux.value/sp.co.value
                error=sp.sig.value/sp.co.value
 

        elif filetype=='p':
            import pickle
            dat=pickle.load( open(filename, "rb" ))
            #tab=dat.keys()
            wave=np.array(dat['wave'])
            flux=np.array(dat['flux'])
            if (len(tab)>=3):
                error=np.array(dat['error'])
            else:
                error=0.*flux
        elif filetype=='temp':
            from astropy.io import fits
            file=fits.open(filename)#(cwd+'/'+filename)
            #a=fits.open(path+'spec_knotA.fits')

            wave=file[2].data
            flux=file[0].data
            error=file[1].data
        #Use linetools.io.readspec to read file
        elif filetype =='linetools':
            from linetools.spectra import io as tio
            sp=tio.readspec(filename,inflg=None, efil=efil,**kwargs)
            wave=sp.wavelength.value
            flux=sp.flux.value

            if sp.sig_is_set == False:
                print('Assuiming arbiarbitrary 10% error on flux')
                error=0.1*flux
            else:
                error=sp.sig.value

        return cls(wave,flux,error,filename=filename)


    @classmethod
    def from_data(cls,wave,flux,error):

        return cls(wave,flux,error,filename=None)


        

        #self.wave=wave
        #self.flux=flux
        #self.error=error
        #self.wrest=wave*(1.+0.)
        #self.zabs=0.
    
    def shift_spec(self,zabs):
        """ Shifts wavelength to absorber rest frame"""
        self.wrest=self.wave/(1.+zabs)
        self.zabs=zabs
        return self.wrest, self.zabs

    def slice_spec(self,lam_rest,lam_min,lam_max,method='closest',linelist='LLS',use_vel=False):
        """
        Slice the spectrum around a central wavelength and convert it to velocity space
        lam_rest : approximate rest wavelength of a transition
        lam_min  : minimum wavelength/velocity to slice 
        lam_max  : maximum wavelength/velocity to slice 

        Keywords:   method = 'closest' [default] -> sets lam_rest to closest atomic transition
                    method = 'Exact' -> uses given lam_rest value to look for transition
                    linelist= Default LLS line linelist, otherwise uses the specified line list

                    use_vel = True -> uses velocity space to slice.
                                   here inputs are lam_min = vel_min [in km/sec]
                                                   lam_max =vel_max [km/s]


        """
        from IGM import rb_setline as s

        str=s.rb_setline(lam_rest,method,linelist=linelist)

        spl=2.9979e5;  #speed of light
        vel = (self.wrest-str['wave']*(1.0 + 0.))*spl/(str['wave']*(1.0 + 0.))

        # Now slice spectrum either velocity or wave space
        if use_vel==False:
            q=np.where((self.wrest >= lam_min) & (self.wrest <= lam_max))
        else:
            vel_min=lam_min
            vel_max=lam_max
            q=np.where((vel >= vel_min) & (vel <= vel_max))

        self.wave_slice=self.wrest[q]
        self.flux_slice=self.flux[q]
        self.error_slice=self.error[q]
        self.linelist=linelist
        self.slice_spec_lam_min=lam_min
        self.slice_spec_lam_max=lam_max
        self.slice_spec_method=use_vel

        
        self.velo=vel[q]
        self.transition=str['wave']
        self.transition_name=str['name']
        self.line_sel_flag=method


    def fit_continuum(self,mask=False,domain=False,Legendre=False,**kwargs):
        """ By default calls an interactive continuum fitter to the sliced spectrum.
            Or an automated Legendre polynomial fitter if keyword set Legendre.
            Order is given by Legendre=order
        """

        verbose = kwargs.get('verbose', False)  # Default is False if not provided

        n_sigma=kwargs.get('n_sigma',3) # sigma clipping level 

        if Legendre==False:
            #pdb.set_trace()
            if 'Interactive' in kwargs:
                Interactive=kwargs['Interactive']
            else:
                Interactive=False

            if 'prefit_cont' in kwargs:
                prefit_cont=kwargs['prefit_cont']
                if verbose:
                    print('Using prefitted continuum...')
                if len(prefit_cont)==1:
                    prefit_cont=prefit_cont*np.ones(len(self.velo),)
                cont=prefit_cont
            else:
                if verbose:
                    print('Initializing interactive continuum fitter...')
                from GUIs import rb_fit_interactive_continuum as f
                s=f.rb_fit_interactive_continuum(self.wave_slice,self.flux_slice,self.error_slice)
                cont=s.cont



        else:
            from astropy.modeling import models, fitting

            order=Legendre
            weight= 1./(self.error_slice**2.)
            
            if domain==False:
                domain=[-600.,600.]

            if mask==False:
                #No Mask
                q=0.*self.wave_slice +1. 

            else:
                #mask=kwargs['mask']            
                #Number of Masks 
                nmsk=int(len(mask)/2)
                vmin=np.zeros(nmsk,)
                vmax=np.zeros(nmsk,)

                for i in range(0,nmsk):
                    vmin[i]=mask[2*i]
                    vmax[i]=mask[2*i+1]

                q=0.*self.wave_slice +1. 

                for i in range(0,nmsk):
                    sq=np.where((self.velo >= vmin[i]) & (self.velo <= vmax[i]))
                    q[sq]=0.
            #Select Unmasked part of the spectrum
            qtq=np.where((q ==1))

            # Fitting the masked Data
            #e=Legendre.fit(self.velo[qtq],self.flux_slice[qtq],order,w=weight[qtq],domain=domain);
            #cont=e(self.velo)

            legendre_init = models.Legendre1D(degree=order)
            # Fit the model using Levenberg-Marquardt minimization
            fitter = fitting.LevMarLSQFitter()
            #fitter=fitting.LinearLSQFitter()
            # Sigma clip the date for a better continuum fit. [Optional]
            if 'sigma_clip' in kwargs:
                from astropy.stats import sigma_clip
                # Apply sigma clipping on the flux array
                sel_flux=self.flux_slice[qtq]
                sel_weight=weight[qtq]
                sel_velo=self.velo[qtq]
                clipped_flux = sigma_clip(sel_flux, sigma=n_sigma, maxiters=5,cenfunc=np.nanmedian)
                # Mask for unclipped (valid) values
                unclipped_mask = ~clipped_flux.mask
                # Select unclipped flux values and corresponding wave values
                unclipped_flux = sel_flux[unclipped_mask]
                unclipped_vel = sel_velo[unclipped_mask]
                unclipped_weigths=sel_weight[unclipped_mask]
                legendre_fit = fitter(legendre_init, unclipped_vel,unclipped_flux,weights=unclipped_weigths)

            else:
                legendre_fit = fitter(legendre_init, self.velo[qtq],self.flux_slice[qtq],weights=weight[qtq])

            cont=legendre_fit(self.velo)
            # Calculate uncertainties (standard deviations) of the fitted parameters
            cov_matrix = fitter.fit_info['param_cov']
            if cov_matrix is not None:
                param_uncertainties = np.sqrt(np.diag(cov_matrix))
                print("Both statistical and continuum fitting error included.")
                # Calculate the 1-sigma confidence bounds
                self.cont_err = calculate_confidence_bounds(self.velo, legendre_fit, cov_matrix)
                self.error_slice=np.sqrt((self.error_slice**2)+(self.cont_err**2))

            else:
                print("Covariance matrix is not available. The fit might be poorly constrained.")
                print("Using only statistical error.")





            # Now mask the part of spectrum that we don't want to fit. 
            # Mask is created to have multiple low vel, high vel ranges.
            # e.g. mask = [-300.,-250.,100.,120.] will exclude -300,-250 and 100,120 km/s parts of the spectrum in the fit


        self.cont=cont
        self.fnorm=self.flux_slice/self.cont
        self.enorm=self.error_slice/self.cont



        #return self.cont,self.fnorm,self.enorm

    def fit_polynomial_ransac(self,degree=3,residual_threshold=0.1,**kwargs):
        """
          Alternate continuum fitting, using ransac to fit polynomial
          degree: polynomial order

        residual_threshold : important RANSAC paramerter to identify inmask points to fit

        _n_bootstrap= number of bootstrap sampling to estimate fitting uncertainty [default = 100]
        """
        _n_bootstrap = kwargs.get('_n_bootstrap', 100)
 
        from IGM import cont_fit_poly_ransac as cf 
        # Fit Legendre polynomial with RANSACdegree 
        # Fit polynomial using RANSAC
        self.cont, self.cont_err = cf.fit_polynomial_ransac(self.wave_slice, self.flux_slice, self.error_slice, degree,residual_threshold=residual_threshold,n_bootstrap=_n_bootstrap)

        self.fnorm=self.flux_slice/self.cont
        self.error_slice=np.sqrt((self.error_slice**2)+(self.cont_err**2))
        self.enorm=self.error_slice/self.cont




    def fit_continuum_ransac(self,window=149,mednorm=False):
        """Alternate continuum fitting method. Does iterative ransac continumm fitting.

        """
        from IGM import ransac_contfit as cf 
        sp=cf.cont_fitter()
        sp=cf.cont_fitter.from_data(self.wave_slice,self.flux_slice,error=self.error_slice,mednorm=mednorm)
        sp.fit_continuum(window=window)        


        self.cont=sp.cont
        self.fnorm=self.flux_slice/self.cont
        self.enorm=self.error_slice/self.cont



    def compute_EW(self, lam_cen, vmin=-50., vmax=50., method='closest', plot=False, **kwargs):
        """Computes rest frame equivalent width and column density for a desired atomic line.
        Around the species lam_cen and given vmin and vmax keyword values. 

        """
        verbose = kwargs.get('verbose', False)  # Default is False if not provided
        SNR=kwargs.get('SNR', False)  # Default is False if not provided
        _binsize = kwargs.get('_binsize', 1)


        from IGM import rb_setline as s
        str=s.rb_setline(lam_cen,method,linelist=self.linelist)

        from IGM import compute_EW as EW
        out = EW.compute_EW(self.wave_slice,self.fnorm,str['wave'],[vmin,vmax],self.enorm,f0=str['fval'],zabs=0.,plot=plot, verbose=verbose,SNR=SNR,_binsize=_binsize)


        self.trans=str['name']
        self.fval=str['fval']
        self.trans_wave=str['wave']
        self.vmin=vmin
        self.vmax=vmax


        self.W= out['ew_tot']
        self.W_e=out['err_ew_tot']
        self.logN=out['col']
        self.logN_e=out['colerr']

        self.Tau=out['Tau_a']
        self.vel_centroid=out['med_vel']
        self.vel_disp=out['vel_disp']
        self.vel50_err = out['vel50_err'] 
        if SNR:
            self.SNR=out['SNR']
        else:
            self.SNR=-99

        #return self.trans,self.fval,self.vmin,self.vmax,self.trans_wave,self.W,self.W_e,self.logN,self.logN_e,self.Tau


    def plot_spec(self):
        """Quick wrapper to call an interactive plotter for the full spectrum as given in input file.
        """
        #from GUIs import rb_plot_spec as sp
        #tt=sp.rb_plot_spec(self.wave,self.flux,self.error)
        from GUIs import PlotSpec_Integrated as sp
        tt=sp.rb_plotspec(self.wave,self.flux,self.error)

    def plot_slice(self):
        """Quick wrapper to call an interactive plotter for the full spectrum as given in input file.
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex = ax1)
        #fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        ax1.step(self.velo,self.flux_slice,'k-')
        ax1.step(self.velo,self.cont,'b-')

        ax1.step(self.velo,self.error_slice,'r-')

        ax1.set_xlim([min(self.velo),max(self.velo)])
        ax1.set_ylim([min(self.flux_slice)-0.02*min(self.flux_slice),max(self.flux_slice)+.1*max(self.flux_slice)])
        ax1.plot([-2500,2500],[0,0],'k:')
        ax1.plot([-2500,2500],[1,1],'k:')       
        ax1.set_xlabel('vel [km/s]')

        #ax2=fig.add_subplot(212)
        ax2.step(self.velo,self.fnorm,'k')
        ax2.step(self.velo,self.enorm,color='r')

        ax2.set_xlim([min(self.velo),max(self.velo)])
        ax2.set_ylim([-0.02,1.8])
        ax2.plot([-2500,2500],[0,0],'k:')
        ax2.plot([-2500,2500],[1,1],'k:')       
        ax2.set_xlabel('vel [km/s]')
        plt.show()


    def save_slice(self, outfilename, file_format='json', verbose=True):
        """Saves the slice object for future processing.

        Parameters:
        -----------
        outfilename : str
            The file path to save the slice object.
        file_format : str, optional
            Format to save the object. Options: 'pickle' (default) or 'json'.

        Notes:
        ------
        - Pickle saves the entire object for later editing.
        - JSON saves only the output data, not the object, so it cannot be reloaded for editing.
        """
        if file_format == 'pickle':
            import pickle
            with open(outfilename, 'wb') as output:
                pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

        elif file_format == 'json':
            import json
            import numpy as np

            # Helper function to convert non-serializable objects
            def convert_for_json(obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                else:
                    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")

            # Create a dictionary of data to save
            data_out = {
                'zabs': self.zabs,
                'linelist': self.linelist,
                'line_sel_flag': self.line_sel_flag,
                'trans': self.trans,
                'fval': self.fval,
                'trans_wave': self.trans_wave,
                'vmin': self.vmin,
                'vmax': self.vmax,
                'W': self.W,
                'W_e': self.W_e,
                'logN': self.logN,
                'logN_e': self.logN_e,
                'vel_centroid': self.vel_centroid,
                'vel_disp': self.vel_disp,
                'vel50_err': self.vel50_err,
                'SNR':self.SNR,
                'wave_slice': self.wave_slice,
                'flux_slice': self.flux_slice,
                'error_slice': self.error_slice,
                'velo': self.velo,
                'cont': self.cont,
                'fnorm': self.fnorm,
                'enorm': self.enorm,
                'Tau': self.Tau,
                'slice_spec_lam_min': self.slice_spec_lam_min,
                'slice_spec_lam_max': self.slice_spec_lam_max,
                'slice_spec_method': self.slice_spec_method
            }

            # Convert arrays to lists before saving
            for key, value in data_out.items():
                if isinstance(value, np.ndarray):
                    data_out[key] = value.tolist()

            # Write JSON data to a file with error handling
            try:
                with open(outfilename, 'w') as json_file:
                    json.dump(data_out, json_file, indent=4, default=convert_for_json)
                if verbose:
                    print(f"File saved to {outfilename} successfully!")
            except TypeError as e:
                print(f"Error saving to JSON: {e}")




    def plot_doublet(self,lam1,lam2,vmin=-600.,vmax=600.,method='closest'):
        """Plot a given doublet defined by the lam1 and lam2 wavelength centers.
        """

        from IGM import rb_setline as s
        str1=s.rb_setline(lam1,method,linelist=self.linelist)
        str2=s.rb_setline(lam2,method,linelist=self.linelist)

        spl=2.9979e5;  #speed of light
        vel1 = (self.wave_slice-str1['wave']*(1.0 + 0.))*spl/(str1['wave']*(1.0 + 0.))
        vel2 = (self.wave_slice-str2['wave']*(1.0 + 0.))*spl/(str2['wave']*(1.0 + 0.))
        
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1=fig.add_subplot(211)
        ax1.step(vel1,self.fnorm)
        ax1.step(vel1,self.enorm,color='r')

        ax1.set_xlim([vmin,vmax])
        ax1.set_ylim([-0.02,1.8])
        ax1.plot([-2500,2500],[0,0],'k:')
        ax1.plot([-2500,2500],[1,1],'k:')       
        ax1.set_xlabel('vel [km/s]')
   
        ax2=fig.add_subplot(212)
        ax2.step(vel2,self.fnorm)
        ax2.step(vel2,self.enorm,color='r')

        ax2.set_xlim([vmin,vmax])
        ax2.set_ylim([-0.02,1.8])
        ax2.plot([-2500,2500],[0,0],'k:')
        ax2.plot([-2500,2500],[1,1],'k:')       
        ax2.set_xlabel('vel [km/s]')
        plt.show()

    def vpfit_singlet(self,FWHM=6.5):
        """Test Wrapper to call vpfit GUI
        """
        from GUIs import rb_interactive_vpfit_singlet as vf 
        vt=vf.rb_interactive_vpfit_singlet(self.wave_slice,self.fnorm,self.enorm,self.transition);    






