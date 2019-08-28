import numpy as np
def bootstrap(data, bootnum=100, samples=None, bootfunc=None):
    """Performs bootstrap resampling on numpy arrays.

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

    """
    if samples is None:
        samples = data.shape[0]

    #make sure the input is sane
    assert samples > 0, "samples cannot be less than one"
    assert bootnum > 0, "bootnum cannot be less than one"

    if bootfunc is None:
        resultdims = (bootnum,) + (samples,) + data.shape[1:]
        boot = np.empty(resultdims)
    else:
        resultdims = (bootnum,)
        boot = np.empty(resultdims)

    for i in range(bootnum):
        bootarr = np.random.randint(low=0,high=data.shape[0],size=samples)
        if bootfunc is None:
            boot[i] = data[bootarr]
        else:
            boot[i] = bootfunc(data[bootarr])

    return boot