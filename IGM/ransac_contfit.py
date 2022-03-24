""" A continuum fitter class. This reads in a 1D spectrum and allows continuum fitting."""
from __future__ import print_function, absolute_import, division, unicode_literals
import os
import time
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg') #python 3 only
import ipdb
from linetools.spectra import io as tio
from scipy.signal import medfilt
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from sklearn.linear_model import RANSACRegressor
from linetools.spectra.xspectrum1d import XSpectrum1D


class cont_fitter(object):
    """
   A continuum fitter class. This reads in a 1D spectrum and allows continuum fitting.
   The initial continuum is fitted using a random sample consensus model.
   Then user has the option to tweak the fitted continuum by hand usling interactive 
   linetools API and/or save the spectrum object.


    Attributes:

        Input:
              filename=input spectrum filename
              efil = input error spectrum name [If it exists, otherwise None]
              window=    default [149], smoothing window
              mednorm = False [default], if set normalizes the spectrum with median value


        output: This gives a cont_fitter object with following attributes:

        self.wave= wavelength

        self.flux= Flux
        self.error= error
        self.cont= Fitted continuum 


    Written : Rongmon Bordoloi      August 2020
    Based on RANSAC continuum fitter written by Bin Liu Summer 2020.
    --------------------------------------------------------------------------------------------
    EXAMPLE: 
               from IGM import ransac_contfit as c 

               sp=c.cont_fitter()

            Two ways to read in spectrum, from file: 
                 #efil = optional if error spectrum is defined in another file
               sp=c.from_file(fluxfilename,efil=errorfilename)

            or from input wave,flux,error array. 
               sp=c.from_data(wave,flux,error=error)
            

            Now fit continuum
               sp=c.fit_continuum()

               #AND YOU ARE DONE.

               #OPTIONAL:

               #Tweak the fitted continuum 
               sp.tweak_continuum()

               #Show new continuum
               sp.plot_spectrum()

               #Save it as a fits file
               sp.save_spectrum(outputfilename)


    --------------------------------------------------------------------------------------------
    """
     

    def from_file(self,filename,efil=None,mednorm=False,**kwargs):
        """
            Read spectrum from filename given.        
        """
        sp=tio.readspec(filename,inflg=None, efil=efil,**kwargs)
        wave=sp.wavelength.value
        if mednorm ==True:
            cnt=np.nanmedian(sp.flux.value)
        else:
            cnt=1.

        flux=sp.flux.value/cnt

        if sp.sig_is_set == False:
            print('Assuiming arbitrary 10% error on flux')
            error=0.1*flux/cnt
        else:
            error=sp.sig.value/cnt

        self.wave=wave
        self.flux=flux
        self.error=error
        self.mednorm=mednorm
        

    def from_data(self, wave,flux,mednorm=False, **kwargs):
        """ read spectrum from input wave,flux,error array. 
        """

        if mednorm ==True:
            cnt=np.nanmedian(flux)
        else:
            cnt=1.

        flux=flux/cnt

        if 'error' in kwargs:
            error=kwargs['error']/cnt
        else:
            print('Assuiming arbitrary 10% error on flux')
            error=0.1*flux/cnt
            

        # Generate
        self.flux=flux
        self.error=error
        self.wave=wave
        self.mednorm=mednorm

    def __init__(self,**kwargs):
        print('Initializing RANSAC continuum fitter')
        np.random.seed(1)



    def fit_continuum(self,window=149):
        self.prepare_data(window=window)
        self.run_ransac(window=window)
        self.spectrum= XSpectrum1D.from_tuple((self.wave, self.flux, self.error,self.cont),verbose=False)


    

    def prepare_data(self,window=149):
        self.flx_md=medfilt(self.flux,window)
        (self.peaks, prop) = find_peaks(self.flux, height=self.flx_md)

        #interpolate the peaks
        peak_interp = np.interp(self.wave, self.wave[self.peaks],self.flux[self.peaks]) 
        peak_interp = np.array(peak_interp)

        ####4. RANSAC RANSAC(spec-evenlope)
        self.sp_diff = self.flux - peak_interp
        self.window=window


    def run_ransac(self,window=149):
        #inlier_masks = []
        #outlier_masks = []
        ransac = RANSACRegressor()
        ransac.fit(self.wave.reshape(-1,1), self.sp_diff)
        #inlier_masks.append(ransac.inlier_mask_)
        inlier_masks=ransac.inlier_mask_

        #outlier_masks.append(np.logical_not(ransac.inlier_mask_))
        outlier_masks=np.logical_not(ransac.inlier_mask_)

        ####5. Use `inlier_masks` to interpolate
        spec_inliers = np.interp(self.wave,self.wave[inlier_masks],self.flux[inlier_masks]) 
        self.cont = medfilt(spec_inliers, window)


    def save_spectrum(self,filename):
        spec= XSpectrum1D.from_tuple((self.wave, self.flux, self.error,self.cont), masking='none')
        spec.write_to_fits(filename)

    def plot_spectrum(self):
        plt.figure(figsize=(10,6))
        plt.plot(self.wave, self.flux, c='grey', label='Original Spectrum')
        plt.plot(self.wave, self.flx_md, c='blue', label='RunningMedFilter, '+str(self.window))
        plt.plot(self.wave, self.cont, c='red', label='RANSAC fit')
        plt.xlabel('Wavelength ($\AA$)', fontsize=15)
        plt.ylabel('Relative Flux ', fontsize=15)
        plt.legend(loc='upper right', fontsize=12)
        plt.ylim([-0.01, self.flux.max()+2])

        plt.show()


    def tweak_continuum(self):
        deltax=self.wave[50]-self.wave[0]
        knot_x=np.arange(min(self.wave),max(self.wave), deltax)
        knot_y =np.interp(knot_x,self.wave,self.cont)

        knots=[[knot_x[i], knot_y[i]] for i in range(len(knot_x))]

        # Get rid of previous knots 
        if os.path.exists("_knots.jsn"):
            os.remove("_knots.jsn")        

        self.spectrum.fit_continuum(knots=knots)
        self.cont=self.spectrum.co







