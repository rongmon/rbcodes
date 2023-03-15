from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
import matplotlib
import Spec_Inspect as S




filename='/Users/bordoloi/Dropbox/COS-Pairs/Targets/J2228-0950/Data/lris_longslit/long_radd.fits'
a=fits.open(filename)
data = a[0].data

two_d_spec=np.transpose(deepcopy(data[:,2100:2200]))
S.Spec_Inspect(two_d_spec)
