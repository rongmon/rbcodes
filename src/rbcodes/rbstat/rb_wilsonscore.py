"""
Wilson score confidence interval for binomial proportions.

Standalone module — not imported by other package modules; available for direct use.

Computes the Wilson score interval for a binomial distribution, which is more
accurate than the normal approximation (Wald interval) for small samples or
extreme probabilities.

Example
-------
    from rbcodes.rbstat.rb_wilsonscore import rb_wilsonscore
    center, hi, lo = rb_wilsonscore(10, 20, 0.95)
"""
from scipy.special import ndtri
import numpy as np
def rb_wilsonscore(count,nobs,confint):
	"""This function computes the wilson score confidence intervals Score Interval for a binomial distribution. 
	
	    Paramters
	    ---------
	        count   =    Number of successes
	        nobs    =    Number of total Trials
	        confint =    Confindence interval for which Wilson Score is computed [e.g. confint =0.95 2\sigma]
	
	    Returns
	    -------
			center =  gives the center of the score intervals given the data
	        hi     =   Upper bound for given confint
	        lo     =   Lower bound for given confint
	
	    Example
	    -------

	        import rb_wilsonscore as w
	        XC, hi, lo = w.rb_wilsonscore(10.,20.,.95)
	
	    Written by :   Rongmon Bordoloi  Nov 2017
	    Tested on  : Python 2.7, 3.x
	-----------------------------------------------------------------------------------
	"""
	count=np.double(count)
	nobs=np.double(nobs)
	confint=np.double(confint)
	# Written by RB
	if nobs == 0.0: return (0.0,0.5 ,1.0)
	z = ndtri(1. - 0.5 * (1.-confint))
	p=count/nobs
	# now do it with Wilson score interval
	alpha=(p+ (z*z)/(2.*nobs))
	beta=((p*(1.-p)/ nobs) + ((z**2.)/(4.*(nobs**2.))))**0.5 
	center = (alpha) / (1. + ((z**2.)/ nobs))
	hi = (alpha+ (z*beta))/ (1. + ((z**2.)/ nobs))
	lo = (alpha - z*(beta)) / (1. + ((z**2.)/ nobs))
	return (center,hi,lo)
