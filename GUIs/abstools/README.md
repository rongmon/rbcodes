# Absorption toolbox

A new toolbox for absorption line analysis. 
Can take one spectrum and interactively and simultaneously analyze several
ions

Example:

path='/Users/bordoloi/Dropbox/COS-Pairs/Targets/PG0832+251/Data/'
filename=path+'LineList_Identified.log'


```python

from astropy.io import fits
from GUIs.abstools import Absorber as A
from GUIs.abstools import Metal_Plot as M   
filename='/Users/bordoloi/Dropbox/COS-Pairs/Targets/PG0832+251/Data/PG0832+251_nbin3_coadd.fits'


a=fits.open(filename)
wave=a[1].data['wave'][0]
flux=a[1].data['flux'][0]   
error=a[1].data['error'][0]  


z=0.017
lines=[1215.67,1260,1334,1393.9,1402,1548.3,1550,1526,1670]
absys=A.Absorber(z,wave,flux,error,lines)   
Abs=absys.ions

# Interevening absorption line list
path='/Users/bordoloi/Dropbox/COS-Pairs/Targets/PG0832+251/Data/'
filename=path+'LineList_Identified.log'

# Intervening absorption line list is optional to show any intervening absorbers
M.Transitions(Abs,intervening=filename)
```


import pickle
pfile='test.p'
with open(pfile,'rb') as pklfile:
    absys=pickle.load(pklfile)

from GUIs.abstools import Metal_Plot as M   
M.Transitions(absys)