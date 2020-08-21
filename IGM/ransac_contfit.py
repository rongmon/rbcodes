""" A continuum fitter class. This reads in a 1D spectrum and allows continuum fitting."""
import os
import time
import numpy as np
import matplotlib
matplotlib.use('TkAgg') #python 3 only
import ipdb
np.random.seed(1)
from linetools.spectra import io as tio
from scipy.signal import medfilt
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from sklearn.linear_model import RANSACRegressor
from linetools.spectra.xspectrum1d import XSpectrum1D


class cont_fitter(object):

    def __init__(self,filename, efil=None,window=149,mednorm=False,**kwargs):

        '''
            """A continuum fitter class. This reads in a 1D spectrum and allows continuum fitting.
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


               #efil = optional if error spectrum is defined in another file
               sp=c.c.cont_fitter(fluxfilename,efil=errorfilename)

               #AND YOU ARE DONE.

               #OPTIONAL:

               #Tweak the fitted continuum 
               sp.tweak_continuum()

               #Show new continuum
               sp.plot_spectrum()

               #Save it as a fits file
               sp.save_spectrum(outputfilename)


    --------------------------------------------------------------------------------------------

 
        '''
        self.mednorm=mednorm
        self.read_spectrum(filename,efil=efil)
        self.prepare_data(window=window)
        self.run_ransac()
        self.spectrum= XSpectrum1D.from_tuple((self.wave, self.flux, self.error,self.cont),verbose=False)


    def read_spectrum(self,filename,efil=None,**kwargs):
        sp=tio.readspec(filename,inflg=None, efil=efil,**kwargs)

        wave=sp.wavelength.value
        if self.mednorm ==True:
            cnt=np.nanmedian(sp.flux.value)
        else:
            cnt=1.

        flux=sp.flux.value/cnt

        if sp.sig_is_set == False:
            print('Assuiming arbiarbitrary 10% error on flux')
            error=0.1*flux/cnt
        else:
            error=sp.sig.value/cnt

        self.wave=wave
        self.flux=flux
        self.error=error


    def prepare_data(self,window=149):
        self.flx_md=medfilt(self.flux,window)
        (self.peaks, prop) = find_peaks(self.flux, height=self.flx_md)

        #interpolate the peaks
        peak_interp = np.interp(self.wave, self.wave[self.peaks],self.flux[self.peaks]) 
        peak_interp = np.array(peak_interp)

        ####4. RANSAC RANSAC(spec-evenlope)
        self.sp_diff = self.flux - peak_interp
        self.window=window


    def run_ransac(self):
        inlier_masks = []
        outlier_masks = []
        ransac = RANSACRegressor()
        ransac.fit(self.wave.reshape(-1,1), self.sp_diff)
        inlier_masks.append(ransac.inlier_mask_)
        outlier_masks.append(np.logical_not(ransac.inlier_mask_))
        ####5. Use `inlier_masks` to interpolate
        spec_inliers = np.interp(self.wave,self.wave[inlier_masks],self.flux[inlier_masks]) 
        self.cont = medfilt(spec_inliers, 99)

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







