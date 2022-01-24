""" Rebin 1D spectrum to new pixel scale."""

import math
import numpy as np
#binning function
#bins up 1D spectra in integer pixels. 
#The routine returns a structure of flux and wavelength and variance that has been rebinned.
def rb_specbin(flux,nbin,**kwargs):
    '''
    %--------------------------------------------------------------------------------
    %   This is a translation of x_specbin.pro from xidl. 
    %  
    %   This function bins up 1D spectra in integer pixels. The routine returns a
    %   structure of flux and wavelength and variance that has been rebinned.
    % 
    %  CALLING SEQUENCE:
    %    bin = x_specbin(fx, nbin WAV, var)
    % 
    %  INPUTS:
    % 
    %    fx       - Flux
    %    nbin     - Number of pixels to bin on
    % 
    %  RETURNS:
    %    bin       - Structure of data
    % 
    % 
    %   Optional Input: 
    %   VAR=  -- Input variance array
    %   WAV=  -- Input wavelength array
    % 
    %  EXAMPLES:
    %    bin = rb_specbin(fx, 3)
    % 
    % 
    %  REVISION HISTORY:
    %   Written by RB. June 2015
    % -
    % ------------------------------------------------------------------------------
    '''
    TF=0
    TFF=0
    if 'var' in kwargs:
        var = kwargs['var']
        TF = 1
    if 'wave' in kwargs:
        wave = kwargs['wave']
        TFF = 1
        

    wavePixels = len(flux)
    if (wavePixels % nbin) != 0:
        numPix = math.floor(wavePixels/nbin) + 1
    else:
        numPix = math.floor(wavePixels/nbin)
    first = wavePixels % nbin
    if first == 0:
        first = nbin
    newFlux = np.zeros(numPix,)
    newVar = np.zeros(numPix,)
    newWave = np.zeros(numPix,) 
    # Binning
    for qq in range(numPix-1): 
        ii = qq*nbin
        index = np.array(range(ii,ii+nbin))
        if qq != numPix: 
            #add them up
            newFlux[qq] = np.mean(flux[index])
            if TF == 1:
                newVar[qq] = np.mean(var[index])
            if TFF == 1:
                newWave[qq] = np.mean(wave[index])
        else: #last pixel
            index = np.array(range(ii,ii+first))
            newFlux[qq] = np.mean(flux[index])
            if TF == 1:
                newVar[qq] = np.mean(var[index])
            if TFF == 1:
                newWave[qq] = np.mean(wave[index])
    output={}
    output['flux']=newFlux

    #returns = newFlux
    if TF == 1:
        #newVar is in second column of returns
        #returns = np.append(returns,newVar,axis = 1)
        output['error']=np.sqrt(newVar/nbin)
    if TFF == 1:
        #returns = np.append(returns,newWave,axis = 1)
        output['wave']=newWave



    return output

