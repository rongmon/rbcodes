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

    def __init__(self,filename, efil=None,medfilter_window=149,**kwargs):
        self.read_spectrum(filename,efil=efil)
        self.prepare_data(medfilter_window=medfilter_window)
        self.run_ransac()
        self.spectrum= XSpectrum1D.from_tuple((self.wave, self.flux, self.error,self.cont),verbose=False)


    def read_spectrum(self,filename,efil=None,**kwargs):
        sp=tio.readspec(filename,inflg=None, efil=efil,**kwargs)
        wave=sp.wavelength.value
        flux=sp.flux.value
        if sp.sig_is_set == False:
            print('Assuiming arbiarbitrary 10% error on flux')
            error=0.1*flux
        else:
            error=sp.sig.value

        self.wave=wave
        self.flux=flux
        self.error=error


    def prepare_data(self,medfilter_window=149):
        self.flx_md=medfilt(self.flux,medfilter_window)
        (self.peaks, prop) = find_peaks(self.flux, height=self.flx_md)

        #interpolate the peaks
        peak_interp = np.interp(self.wave, self.wave[self.peaks],self.flux[self.peaks]) 
        peak_interp = np.array(peak_interp)

        ####4. RANSAC RANSAC(spec-evenlope)
        self.sp_diff = self.flux - peak_interp
        self.medfilter_window=medfilter_window


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
        spec.write_to_fits(filename)

    def plot_spectrum(self):
        plt.figure(figsize=(10,6))
        plt.plot(self.wave, self.flux, c='grey', label='Original Spectrum')
        plt.plot(self.wave, self.flx_md, c='blue', label='RunningMedFilter, '+str(self.medfilter_window))
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







