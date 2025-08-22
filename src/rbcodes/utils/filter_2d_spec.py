from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, Box1DKernel, interpolate_replace_nans, AiryDisk2DKernel
import copy
import numpy as np


def filter_1D(wav_array,oned,method='convolution'):
    dx=9.
    smoothing_kernel=19

    #sel_mask=(wav_array>5008.24-1.5*dx)*(wav_array< 5008.24+1.5*dx) + (wav_array>4960.-1.2*dx)*(wav_array< 4960.+2*dx)+(wav_array>4864.-2*dx)*(wav_array< 4864.+2*dx) + (wav_array>5877.28-2*dx)*(wav_array< 5877.28+2*dx) + (wav_array>4340.-dx)*(wav_array< 4340.+dx)
    sel_mask=(wav_array>5008.24-1.5*dx)*(wav_array< 5008.24+1.5*dx)+ (wav_array>4960.-1.2*dx)*(wav_array< 4960.+2*dx)+(wav_array>4864.-2*dx)*(wav_array< 4864.+2*dx) + (wav_array>5877.28-2*dx)*(wav_array< 5877.28+2*dx) 
    + (wav_array>4341.-dx)*(wav_array< 4341.+dx)+ (wav_array>5016.5-dx)*(wav_array< 5016.5+dx)

    thisoned=copy.deepcopy(oned)
    thisoned[sel_mask]=np.nan


    if method =='median_filter':
        #Widekernel
        kx=71
        kx_gap=9        
        n_x=len(wav_array)
        idx_x = np.fromfunction(lambda  i, j: i+j, (n_x, kx), dtype=np.int64) - kx // 2
    
        # Exclude the center gap
        if kx_gap>0:
                #print('Gap size: ', kx_gap, flush=True)
                cols_remain=np.append(np.arange((kx - kx_gap)/2., dtype=np.int64),
                               np.flip(kx-1-np.arange((kx - kx_gap)/2., dtype=np.int64)))
                idx_x = idx_x[:, cols_remain]
    
        idx_x[idx_x < 0]=0
        idx_x[idx_x > n_x-1]=n_x-1
        cont = np.nanmedian(thisoned[idx_x], axis=1)
    elif method=='convolution':
        kernel = Gaussian1DKernel(stddev=smoothing_kernel,mode='oversample', factor=10)
        #kernel=Box1DKernel(width=dx)
        cont = convolve(thisoned, kernel)
   

    filtered=cont
    return filtered



def median_filter(wav_array,twod):

    img_back=copy.deepcopy(twod)

    ly,lx=np.shape(twod)
    for i in range(0,ly):
        #img_back[i,:]=cont_fit_1d(wave,flux_sampled_A[i,:],e_sampled_A[i,:])
        img_back[i,:]=filter_1D(wav_array,twod[i,:],method='median_filter')
    return img_back

def Gaussian_filter(wav_array,twod):

    img_back=copy.deepcopy(twod)

    ly,lx=np.shape(twod)
    for i in range(0,ly):
        #img_back[i,:]=cont_fit_1d(wave,flux_sampled_A[i,:],e_sampled_A[i,:])
        img_back[i,:]=filter_1D(wav_array,twod[i,:],method='convolution')
    return img_back



def filter_2D(wav_array,twod,pix_radius=15,smoothing_kernel=19,dx=1,y_stddev=0.2):
    """
       Function to filter a 2D spectrum while masking regions around emission lines with a 2D Gaussian kernel

    """
    #dx=1.
    #pix_radius=15
    #smoothing_kernel=19
    #HeI 5016,5017.079
    #HeI 5876,5877.299
    #,7067.198,7283.356,8752.876,8865.217,9017.385,9071.1,9533.2,9231.547,9548.591,10052.13,10833.22,10941.09,12821.59,18756.13,19450.87,21661.2
    linelist=[4341.,4960.,4864.,5008.24,5017.,5877.28,6548.,6565.,6584.,6618.,6717.,6732.,]

    #ly,lx=np.shape(twod)
    nrows, ncols = twod.shape
    row, col = np.ogrid[:nrows, :ncols]

    thistwod=copy.deepcopy(twod)


    #First mask out regions of interst
    for index in range(0,len(linelist)):
        sel_mask=np.where( (wav_array>linelist[index]-1.5*dx) & (wav_array<linelist[index]+1.5*dx))
        
        #Start Masking
        if len(sel_mask[0])>0:
                this_index=sel_mask[0][0]
                cen_x,cen_y= this_index,int(nrows/2)

                ROI=((row - cen_y)**2 + (col - cen_x)**2 < (pix_radius)**2)
                thistwod[ROI]=np.nan

    kernel =  Gaussian2DKernel(x_stddev=smoothing_kernel,y_stddev=y_stddev)
    #kernel=AiryDisk2DKernel(smoothing_kernel)
    cont = convolve_fft(thistwod, kernel,nan_treatment='interpolate')

    # create a "fixed" image with NaNs replaced by interpolated values
    #cont = interpolate_replace_nans(thistwod, kernel)



    filtered=cont
    return filtered
