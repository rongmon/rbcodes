""" Spectrum class to read in, analyze and measure absorption lines."""
import numpy as np
from scipy.interpolate import splrep,splev
from numpy.polynomial.legendre import Legendre
import sys
import os
import pdb
import warnings

import json

#rbcodes imports 

try:
    from rbcodes.IGM import compute_EW as EW
    from rbcodes.IGM import rb_setline as s


except:
    from IGM import compute_EW as EW
    from IGM import rb_setline as s




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
        from rbcodes.GUIs.rb_spec import rb_spec as r 
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
        
        #Fit continuum Mask the regions defined by velocity/ use_weight= True/False toggles if error is used for weighting the fit
        s.fit_continuum(mask=[-200,300,500,1100],domain=xlim,Legendre=3, use_weights=False)
        
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
        from rbcodes.GUIs.rb_spec import load_rb_spec_object as r_load
        
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
    def __init__(self, wave, flux, error, filename=False):
        """
        Initializes a Spectrum object for absorption line analysis.
    
        Parameters:
        wave (array-like): Wavelength values.
        flux (array-like): Flux values.
        error (array-like): Error values.
        filename (str, optional): Filename if reading from a file. Default is False.
        
        Raises:
        ValueError: If input arrays have mismatched lengths or contain invalid values.
        """
        # Verify input arrays have the same length
        if len(wave) != len(flux):
            raise ValueError("Wavelength and flux arrays must have the same length")
        if len(wave) != len(error):
            raise ValueError("Wavelength and error arrays must have the same length")
        
        # Check for NaN or infinity values
        if np.any(np.isnan(wave)) or np.any(np.isinf(wave)):
            raise ValueError("Wavelength array contains NaN or infinity values")
        
        # Check if arrays are empty
        if len(wave) == 0:
            raise ValueError("Input arrays cannot be empty")
        
        # Calculate median flux carefully to handle edge cases
        try:
            median_flux = np.nanmedian(flux)
            if median_flux == 0 or np.isnan(median_flux):
                # If median is zero or NaN, use 1.0 for scaling to avoid division by zero
                median_flux = 1.0
                #print("Warning: Median flux is zero or NaN. Using unscaled flux values.")
                warnings.warn("⚠️ Warning: Median flux is zero or NaN. Using unscaled flux values.", category=UserWarning)
        except Exception as e:
            #print(f"Warning: Could not calculate median flux: {e}. Using unscaled values.")
            warnings.warn(f"⚠️ Warning: Could not calculate median flux: {e}. Using unscaled values.", category=UserWarning)

            median_flux = 1.0
        
        self.wave = np.array(wave)
        self.flux = np.array(flux) / median_flux
        self.error = np.array(error) / median_flux
        self.filename = filename
    
    @classmethod
    def from_file(cls, filename, filetype=False, efil=None, **kwargs):
        """
        Creates a Spectrum object from a file.
    
        Parameters:
        filename (str): Path to the file.
        filetype (str, optional): Type of file (e.g., 'ascii', 'fits', 'HSLA', etc.). If False, it is inferred.
        efil (str, optional): Error file if separate.
        kwargs: Additional arguments for specific file readers.
    
        Returns:
        Spectrum: An instance of the Spectrum class.
    
        Raises:
        FileNotFoundError: If the file does not exist.
        ImportError: If required modules cannot be imported.
        ValueError: If file type is not supported or file content is invalid.
        IOError: If there's an error reading the file.
        """
        # Check if filename is provided
        if filename is None:
            if 'wave' in kwargs and 'flux' in kwargs and 'error' in kwargs:
                return cls(kwargs['wave'], kwargs['flux'], kwargs['error'])
            raise IOError("Input wavelength, flux, and error arrays are required.")
        
        # Check if the file exists
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File not found: {filename}")
        
        # Infer filetype if not specified
        if filetype is False:
            ext = os.path.splitext(filename)[1].lower()
            if ext in ['.txt', '.dat']:
                filetype = 'ascii'
            elif ext == '.fits':
                filetype = 'fits'
            elif ext == '.p' or ext == '.pkl':
                filetype = 'p'
            else:
                filetype = ext[1:] if ext.startswith('.') else ext
            print(f"Inferred file type: {filetype} based on extension")
        
        try:
            if filetype == 'ascii':
                try:
                    from astropy.io import ascii
                except ImportError:
                    raise ImportError("astropy.io.ascii is required to read ASCII files")
                
                try:
                    data = ascii.read(filename)
                    if len(data) == 0:
                        raise ValueError(f"Empty data file: {filename}")
                    
                    keys = data.keys()
                    if len(keys) < 2:
                        raise ValueError(f"Not enough columns in file: {filename}. Need at least wavelength and     flux.")
                    
                    wave = np.array(data[keys[0]])
                    flux = np.array(data[keys[1]])
                    
                    # Handle error column or create error array
                    if len(keys) >= 3:
                        error = np.array(data[keys[2]])
                    else:
                        #print("Warning: No error column found. Assuming 10% error.")
                        warnings.warn(f"⚠️ Warning: No error column found. Assuming 10% error.", category=UserWarning)

                        error = 0.1 * np.abs(flux)
                except Exception as e:
                    raise IOError(f"Error reading ASCII file: {e}")
            
            elif filetype in ['fits', 'HSLA']:
                try:
                    from astropy.io import fits
                except ImportError:
                    raise ImportError("astropy.io.fits is required to read FITS files")
                
                try:
                    with fits.open(filename) as file:
                        if len(file) < 2:
                            raise ValueError(f"FITS file does not have expected HDU structure: {filename}")
                        
                        data = file[1].data
                        if data is None or len(data) == 0:
                            raise ValueError(f"Empty data in FITS file: {filename}")
                        
                        # Check for required columns
                        if 'WAVE' not in data.names and 'wave' not in data.names:
                            raise ValueError(f"Wavelength column not found in FITS file: {filename}")
                        if 'FLUX' not in data.names and 'flux' not in data.names:
                            raise ValueError(f"Flux column not found in FITS file: {filename}")
                        
                        # Get column names with case insensitivity
                        wave_col = 'WAVE' if 'WAVE' in data.names else 'wave'
                        flux_col = 'FLUX' if 'FLUX' in data.names else 'flux'
                        error_col = None
                        for col in ['ERROR', 'error', 'ERR', 'err']:
                            if col in data.names:
                                error_col = col
                                break
                        
                        wave = np.array(data[wave_col])
                        flux = np.array(data[flux_col])
                        
                        if error_col:
                            error = np.array(data[error_col])
                        else:
                            #print("Warning: No error column found in FITS file. Assuming 10% error.")
                            warnings.warn(f"⚠️ Warning: No error column found. Assuming 10% error.", category=UserWarning)

                            error = 0.1 * np.abs(flux)
                except Exception as e:
                    raise IOError(f"Error reading FITS file: {e}")
            
            elif filetype == 'xfits':
                try:
                    from linetools.spectra.xspectrum1d import XSpectrum1D
                except ImportError:
                    raise ImportError("linetools is required to read XFITS files")
                
                try:
                    sp = XSpectrum1D.from_file(filename)
                    wave = sp.wavelength.value
                    flux = sp.flux.value
                    
                    if sp.sig_is_set:
                        error = sp.sig.value
                    else:
                        #print("Warning: No error information in XSpectrum. Assuming 10% error.")
                        warnings.warn(f"⚠️ Warning: No error information in XSpectrum. Assuming 10% error.", category=UserWarning)

                        error = 0.1 * np.abs(flux)
                    
                    if sp.co_is_set:
                        print("Normalizing by continuum from file")
                        flux /= sp.co.value
                        error /= sp.co.value
                except Exception as e:
                    raise IOError(f"Error reading XFITS file: {e}")
            
            elif filetype == 'p':
                try:
                    import pickle
                except ImportError:
                    raise ImportError("pickle module is required to read pickle files")
                
                try:
                    with open(filename, "rb") as f:
                        data = pickle.load(f)
                    
                    if not isinstance(data, dict):
                        raise ValueError(f"Pickle file does not contain a dictionary: {filename}")
                    
                    if 'wave' not in data:
                        raise ValueError(f"'wave' key not found in pickle file: {filename}")
                    if 'flux' not in data:
                        raise ValueError(f"'flux' key not found in pickle file: {filename}")
                    
                    wave = np.array(data['wave'])
                    flux = np.array(data['flux'])
                    
                    if 'error' in data:
                        error = np.array(data['error'])
                    else:
                        #print("Warning: No 'error' key found in pickle file. Assuming 10% error.")
                        warnings.warn(f"⚠️ Warning: No 'error' key found in pickle file. Assuming 10% error.", category=UserWarning)

                        error = 0.1 * np.abs(flux)
                except Exception as e:
                    raise IOError(f"Error reading pickle file: {e}")
            
            elif filetype == 'temp':
                try:
                    from astropy.io import fits
                except ImportError:
                    raise ImportError("astropy.io.fits is required to read TEMP files")
                
                try:
                    with fits.open(filename) as file:
                        if len(file) < 3:
                            raise ValueError(f"TEMP FITS file does not have expected HDU structure: {filename}")
                        
                        wave = file[2].data
                        flux = file[0].data
                        error = file[1].data
                        
                        if wave is None or flux is None or error is None:
                            raise ValueError(f"Missing data in TEMP FITS file: {filename}")
                except Exception as e:
                    raise IOError(f"Error reading TEMP FITS file: {e}")
            
            elif filetype == 'linetools':
                try:
                    from linetools.spectra import io as tio
                except ImportError:
                    raise ImportError("linetools is required to read files with linetools")
                
                try:
                    sp = tio.readspec(filename, inflg=None, efil=efil, **kwargs)
                    wave = sp.wavelength.value
                    flux = sp.flux.value
                    
                    if sp.sig_is_set:
                        error = sp.sig.value
                    else:
                        #print("Warning: No error information. Assuming 10% error.")
                        warnings.warn(f"⚠️ Warning: No error information. Assuming 10% error.", category=UserWarning)
                        error = 0.1 * flux
                except Exception as e:
                    raise IOError(f"Error reading file with linetools: {e}")
            
            else:
                raise ValueError(f"Unsupported file type: {filetype}")
            
            # Verify data integrity
            if len(wave) == 0 or len(flux) == 0:
                raise ValueError(f"Empty arrays read from file: {filename}")
            
            # Check for NaN or Inf values
            if np.any(np.isnan(wave)) or np.any(np.isinf(wave)):
                #print("Warning: Wavelength array contains NaN or infinity values")
                warnings.warn(f"⚠️ Warning: Wavelength array contains NaN or infinity values", category=UserWarning)
            
            if np.any(np.isnan(flux)) or np.any(np.isinf(flux)):
                #print("Warning: Flux array contains NaN or infinity values")
                warnings.warn(f"⚠️ Warning: Flux array contains NaN or infinity values", category=UserWarning)
                # Replace NaN/Inf with interpolated or zero values
                bad_indices = np.isnan(flux) | np.isinf(flux)
                if np.all(bad_indices):
                    #print("All flux values are NaN or Inf. Replacing with zeros.")
                    warnings.warn(f"⚠️ Warning: Either Check original spectrum or proceed as stated below!", category=UserWarning)
                    warnings.warn(f"⚠️ All flux values are NaN or Inf. Replacing with zeros.", category=UserWarning)
                    flux = np.zeros_like(flux)
                else:
                    good_indices = ~bad_indices
                    if np.any(good_indices):
                        try:
                            from scipy.interpolate import interp1d
                            x_good = wave[good_indices]
                            y_good = flux[good_indices]
                            f = interp1d(x_good, y_good, bounds_error=False, fill_value=0)
                            flux[bad_indices] = f(wave[bad_indices])
                            print(f"Replaced {np.sum(bad_indices)} NaN/Inf flux values with interpolated values")
                        except:
                            flux[bad_indices] = 0
                            print(f"Replaced {np.sum(bad_indices)} NaN/Inf flux values with zeros")
            
            if np.any(np.isnan(error)) or np.any(np.isinf(error)):
                warnings.warn(f"⚠️ Warning: Either Check original spectrum or proceed as stated below!", category=UserWarning)

                #print("Warning: Error array contains NaN or infinity values")
                warnings.warn(f"⚠️ Warning: Error array contains NaN or infinity values.", category=UserWarning)
                # Replace NaN/Inf with median or percentage values
                bad_indices = np.isnan(error) | np.isinf(error) | (error <= 0)
                if np.all(bad_indices):
                    #print("All error values are invalid. Using 10% of flux as error.")
                    warnings.warn(f"⚠️ All error values are invalid. Using 10% of flux as error.", category=UserWarning)
                    error = 0.1 * np.abs(flux)
                else:
                    good_indices = ~bad_indices
                    if np.any(good_indices):
                        median_error = np.median(error[good_indices])
                        error[bad_indices] = median_error if median_error > 0 else 0.1 * np.abs(flux[    bad_indices])
                        #print(f"Replaced {np.sum(bad_indices)} invalid error values with {median_error:.2e}")
                        warnings.warn(f"⚠️ Replaced {np.sum(bad_indices)} invalid error values with {median_error:.2e}", UserWarning)

            return cls(wave, flux, error, filename=filename)
        
        except Exception as e:
            # Re-raise the exception with additional context
            raise type(e)(f"{str(e)} while processing file: {filename}").with_traceback(sys.exc_info()[2])
    
    @classmethod
    def from_data(cls, wave, flux, error):
        """
        Creates a Spectrum object from given data arrays.
    
        Parameters:
        wave (array-like): Wavelength values.
        flux (array-like): Flux values.
        error (array-like): Error values.
    
        Returns:
        Spectrum: An instance of the Spectrum class.
    
        Raises:
        ValueError: If input arrays have mismatched lengths or contain invalid values.
        """
        # Convert inputs to numpy arrays if they aren't already
        wave = np.asarray(wave)
        flux = np.asarray(flux)
        error = np.asarray(error)
        
        # Verify shapes
        if wave.shape != flux.shape:
            raise ValueError(f"Wave and flux arrays must have the same shape: {wave.shape} vs {flux.shape}")
        if wave.shape != error.shape:
            raise ValueError(f"Wave and error arrays must have the same shape: {wave.shape} vs {error.shape}")
        
        # Check for empty arrays
        if wave.size == 0:
            raise ValueError("Input arrays cannot be empty")
        
        # Check for all-NaN arrays
        if np.all(np.isnan(wave)):
            raise ValueError("Wavelength array contains only NaN values")
        if np.all(np.isnan(flux)):
            raise ValueError("Flux array contains only NaN values")
        
        # Check for negative wavelengths
        if np.any(wave <= 0):
            warnstr = "Wavelength array contains zero or negative values"
            warnings.warn(f"⚠️ {warnstr}", category=UserWarning)
        
        # Check for negative error values
        if np.any(error < 0):
            warnstr="Warning: Error array contains negative values. Taking absolute values."
            warnings.warn(f"⚠️ {warnstr}", category=UserWarning)
        
            error = np.abs(error)
        
        # If error array is all zeros, set to default
        if np.all(error == 0):
            warnstr="Warning: Error array contains all zeros. Setting to 10% of flux."
            warnings.warn(f"⚠️ {warnstr}", category=UserWarning)
        
            error = 0.1 * np.abs(flux)
        
        return cls(wave, flux, error)    

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
                try:
                    from rbcodes.GUIs import rb_fit_interactive_continuum as f
                except:
                    from GUIs import rb_fit_interactive_continuum as f
                s=f.rb_fit_interactive_continuum(self.wave_slice,self.flux_slice,self.error_slice)
                cont=s.cont



        else:

            #Moving onto using my iterative continuum fitting module. More versatile.
            try:
                from rbcodes.IGM.rb_iter_contfit import rb_iter_contfit
            except:
                from IGM.rb_iter_contfit import rb_iter_contfit
            
            order=Legendre
            # Setup parameters
            if domain == False:
                domain = [-600., 600.]
            
            # Handle masks
            if mask == False:
                # No Mask
                q = 0. * self.wave_slice + 1. 
            else:
                # Process the mask regions
                nmsk = int(len(mask)/2)
                vmin = np.zeros(nmsk,)
                vmax = np.zeros(nmsk,)
                for i in range(0, nmsk):
                    vmin[i] = mask[2*i]
                    vmax[i] = mask[2*i+1]
                q = 0. * self.wave_slice + 1. 
                for i in range(0, nmsk):
                    sq = np.where((self.velo >= vmin[i]) & (self.velo <= vmax[i]))
                    q[sq] = 0.
            
            # Get unmasked indices
            unmasked_indices = np.where(q == 1)
            
            # Extract only unmasked data points
            velo_unmasked = self.velo[unmasked_indices]  # Use velocity for fitting
            wave_unmasked=self.wave_slice[unmasked_indices]
            flux_unmasked = self.flux_slice[unmasked_indices]
            error_unmasked = self.error_slice[unmasked_indices]
            
            # Set up fitting parameters
            use_weights = kwargs.get('use_weights', False)
            use_sigma_clip = 'sigma_clip' in kwargs
            sigma_level = kwargs.get('n_sigma', 3.0) if use_sigma_clip else 3.0
            
            # Call rb_iter_contfit with only unmasked points
            _, _, fit_error, fit_model,fitter = rb_iter_contfit(
                velo_unmasked,
                flux_unmasked,
                error_unmasked,
                order=order,
                sigma=sigma_level,
                use_weights=use_weights,
                return_model=True,
                maxiter=kwargs.get('maxiter', 25)
            )
            
            # Generate continuum over the full range
            cont = fit_model(self.velo)
            
            # Calculate uncertainties from the fit
            if use_weights:

                if hasattr(fitter, 'fit_info') and 'param_cov' in fitter.fit_info:
                    cov_matrix = fitter.fit_info['param_cov']
                    if cov_matrix is not None:
                        param_uncertainties = np.sqrt(np.diag(cov_matrix))
                        print("Both statistical and continuum fitting error included.")
                        # Calculate the 1-sigma confidence bounds
                        self.cont_err = calculate_confidence_bounds(self.velo, fit_model, cov_matrix)
                        self.error_slice = np.sqrt((self.error_slice**2) + (self.cont_err**2))
                    else:
                        print("Covariance matrix is not available. The fit might be poorly constrained.")
                        print("Using fit-error as proxy.")
                        self.cont_err= np.ones(len(self.velo),)*fit_error
                        self.error_slice = np.sqrt((self.error_slice**2) + (self.cont_err**2))                        
                else:
                    self.cont_err=fit_error
                    print("Using only statistical error.")
            else:
                self.cont_err=fit_error
                print("Unweighted fit performed - using only statistical error.")
            
            """
            #This is the older implementation. Now a cleaner implementation is used.
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

            """




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
        try:
            from rbcodes.IGM import cont_fit_poly_ransac as cf 
        except:
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
        
        try:
            from rbcodes.IGM import ransac_contfit as cf 
        except:
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


        str=s.rb_setline(lam_cen,method,linelist=self.linelist)

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
        try:
            from GUIs import PlotSpec_Integrated as sp
        except:
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
        try:
            from rbcodes.GUIs import rb_interactive_vpfit_singlet as vf 
        except:
            from GUIs import rb_interactive_vpfit_singlet as vf 
        vt=vf.rb_interactive_vpfit_singlet(self.wave_slice,self.fnorm,self.enorm,self.transition);    






