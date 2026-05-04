# Project Documentation
[Back to Main Page](../main_readme.md)

*Auto-generated documentation from docstrings*

## Modules

### rb_boot

Bootstrap resampling utility.

Standalone module — not imported by other package modules; available for direct use.

Provides case-bootstrap resampling on numpy arrays, with optional support for
a user-supplied summary statistic function (``bootfunc``). Returns bootstrap
resamples or the distribution of the statistic, which can be used to estimate
confidence intervals.

Example
-------
    from rbcodes.rbstat.rb_boot import bootstrap
    import numpy as np
    data = np.random.normal(0, 1, 100)
    boot_means = bootstrap(data, bootnum=1000, bootfunc=np.mean)

### rb_wilsonscore

Wilson score confidence interval for binomial proportions.

Standalone module — not imported by other package modules; available for direct use.

Computes the Wilson score interval for a binomial distribution, which is more
accurate than the normal approximation (Wald interval) for small samples or
extreme probabilities.

Example
-------
    from rbcodes.rbstat.rb_wilsonscore import rb_wilsonscore
    center, hi, lo = rb_wilsonscore(10, 20, 0.95)

## Functions

### bootstrap() (`rb_boot`)

Performs bootstrap resampling on numpy arrays.

    Bootstrap resampling is used to understand confidence intervals of sample
    estimates. This function returns versions of the dataset resampled with
    replacement ("case bootstrapping"). These can all be run through a function
    or statistic to produce a distribution of values which can then be used to
    find the confidence intervals.

    Parameters
    ----------
    data : numpy.ndarray
        N-D array. The boostrap resampling will be performed on the first
        index, so the first index should access the relevant information
        to be bootstrapped.
    bootnum : int
        Number of bootstrap resamples
    samples : int
        Number of samples in each resample. The default None sets samples to
        the number of datapoints
    bootfunc : function
        Function to reduce the resampled data. Each bootstrap resample will
        be put through this function and the results returned. If None, the
        bootstrapped data will be returned

    Returns
    -------
    boot : numpy.ndarray
        Bootstrapped data. Each row is a bootstrap resample of the data.

### rb_wilsonscore() (`rb_wilsonscore`)

This function computes the wilson score confidence intervals Score Interval for a binomial distribution. 
	
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

