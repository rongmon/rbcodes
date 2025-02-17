""" Rebin 1D spectrum to new pixel scale."""

import math
import numpy as np
def rb_specbin(flux,nbin,**kwargs):
    """This function bins up 1D spectra in integer pixels. The routine returns a
       structure of flux and wavelength and variance that has been rebinned.
    
    Parameters
    -----------
     
        fx       - Flux
        nbin     - Number of pixels to bin on
        VAR=  -- Input variance array [Optional]
        WAV=  -- Input wavelength array [Optional]
    
    Returns
    --------
        bin       - Structure of data
      
    Example
    --------
        bin = rb_specbin(fx, 3)
     
     
      REVISION HISTORY:
       Written by RB. June 2015
     -
     ------------------------------------------------------------------------------
    """
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
    output['flux']=newFlux[:-1]

    #returns = newFlux
    if TF == 1:
        #newVar is in second column of returns
        #returns = np.append(returns,newVar,axis = 1)
        output['error']=np.sqrt(newVar[:-1]/nbin)
    if TFF == 1:
        #returns = np.append(returns,newWave,axis = 1)
        output['wave']=newWave[:-1]



    return output

