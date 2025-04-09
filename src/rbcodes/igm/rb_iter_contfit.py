import numpy as np
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
import warnings
import pdb

def rb_iter_contfit(wave,flux,error,**kwargs):
    """Iterative continuum fitter using Legendre polynomials
    
    Parameters
    ----------
                wavelength array
                flux array 
                error array
        optional input:
                maxiter :-  maximum iteration [25 default]
                order   :-  polynomial order of fit [4 default]

    Returns
    ---------
               fit_final : Final fitted continuum array
               resid_final : residual error array
               fit_error  : error on the fit [standard deviation of the residual]

    Written by:  Rongmon Bordoloi
    Tested on Python 3.7  Sep 4 2019
    --------------------------

    Example
    --------
    from IGM import rb_iter_contfit as r
        out= r.rb_iter_contfit(wave,flux,error,order=5)

        out[0] = fitted continuum
    """

    if 'maxiter' in kwargs:
        maxiter=kwargs['maxiter']
    else:
        maxiter=25
    if 'order' in kwargs:
        order=kwargs['order']
    else:
        order=4

    #Initialize a mask
    mask=np.ones((np.size(wave),))

    # Looking for chip gaps 
    chip_gap= np.where(((error==0) & (flux==0)) | (error-flux==0));
    if (np.size(chip_gap)>0):
        #error[chip_gap]=1;
        #flux[chip_gap]=1;
        mask[chip_gap]=0

    # Now get rid of negative error values
    qq=np.where(error <= 0)
    if (np.size(qq) >0):
        error[qq]=np.median(error)

    #Do a sanity check to avoid bad flux values
    q=np.where(flux <= 0)
    if (np.size(q)>0):
        flux[q]=error[q]

    w=1/error**2. # Weight

    index=0 #Initialize the counter
    outside_chip_gap= np.where(((error!=0) & (flux!=0)) | (error-flux!=0));
    med_err=np.median(error[outside_chip_gap])
    med_flux=np.median(flux[outside_chip_gap])

    # Now take out any possible emission or absorption features
    #%Flux values less than the median error are 1sigma within zero, so should be masked
    #%Try to exclude emission. anything above 2sigma over the median flux should be masked
    bd=np.where((flux < med_err)) #| (flux > (med_flux+3*med_err)))
    nbd = len(bd[0])
    if (nbd>0):
        mask[bd]=0;
    mask=np.array(mask)
    qq=np.where(mask==1)
    #Here I'm cutting out any chip gap or any other absorption feature. These are permanently excluded from the fit
    flux_new=flux[qq]
    wave_new=wave[qq]
    weights=np.array(w[qq])
    # Fit the data using Legedre Polynomials
    g_init = models.Legendre1D(order)
    # initialize fitters
    fit = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore') 
        new_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,niter=maxiter, sigma=3.0)
        filtered_fit,filtered_data = new_fit(g_init, wave_new,flux_new)#,weights=weights)

    ''''
    while (index < maxiter):
        print("Index="+np.str(index))
        index=index+1
        oldmask=mask
        # Fit the data using Legedre Polynomials
        g_init = models.Legendre1D(order)
        # initialize fitters
        fit = fitting.LevMarLSQFitter()
        new_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,niter=maxiter, sigma=3.0)
        filtered_data, new_fit = new_fit(g_init, wave_new,flux_new,weights=weights)


        #fit_g = fitting.LevMarLSQFitter()
        #new_fit = fit_g(g_init, wave_new,flux_new) # Learn how to use the weights....
        cont=new_fit(wave_new)
        resid=flux_new/np.double(cont)
        #mederr=np.sqrt(scipy.stats.moment(resid,2)) #;;Use the median error to set "sigma" for the clipping. If you median is skewed, i think you're pretty much hosed
        mederr=np.std(resid)
        # Do some sigma clipping here
        median_val=np.median(resid)
        kappa=3.;  # 3 sigma clipping
        bd = np.where( (resid > (median_val+kappa*mederr)) | (resid < (median_val - mederr*kappa))) 
        gd= np.where( (resid <= (median_val+kappa*mederr)) & (resid >= (median_val - mederr*kappa))) 
        print np.size(gd)
        print np.size(bd)
        if (np.size(gd) > 0):
            mask[gd] = 1  # Allowing previously rejected points to reenter
        if (np.size(bd) > 0):
            mask[bd] = 0
 
        
        if (index >0):
            diff=oldmask-mask
            qq = np.where(diff != 0.)
            if (np.size(qq)==0):
                print index
                print "No more Clipping"
                index=maxiter # ;;mask and oldmask are identical if all elements match
        
        # Updating everything
        qq=np.where(mask==1)
        flux_new=flux[qq]
        wave_new=wave[qq]
        weights=w[qq]

    '''
    #pdb.set_trace()
    fit_final=filtered_fit(wave)
    resid_final=flux/fit_final
    fit_error=np.std(resid_final)#np.sqrt(scipy.stats.moment(resid_final,2))

    
    return fit_final,resid_final,fit_error

