# Absorption toolbox

A GUI toolbox for CGM/IGM absorption line analysis. 
Can take one spectrum and interactively and simultaneously analyze several
ions.

Example:

First Example shows how to start the GUI up from an ipython terminal

```python

#Read in packages 
from astropy.io import fits
from GUIs.abstools import Absorber as A
from GUIs.abstools import Metal_Plot as M   

# Read in the 1D spectrum to be analyzed
filename='/Users/bordoloi/Dropbox/COS-Pairs/Targets/PG0832+251/Data/PG0832+251_nbin3_coadd.fits'


a=fits.open(filename)
wave=a[1].data['wave'][0]
flux=a[1].data['flux'][0]   
error=a[1].data['error'][0]  
#------------------------------


#Specify redshift at which to perform analysis
z=0.017

#Give approximate absorption line rest frame wavelengths to be analyzed
lines=[1215.67,1260,1334,1393.9,1402,1548.3,1550,1526,1670]

#Preprocessing:
# Create an absorber class to feed into the main GUI
absys=A.Absorber(z,wave,flux,error,lines=lines,window_lim=[-2000,2000])   
Abs=absys.ions


#Optional: if you have a pre-identified intervening absorption line list, feed it in.
# Interevening absorption line list
path='/Users/bordoloi/WORK/python/rbcodes/GUIs/abstools/SpecPlot_Projects/'
filename=path+'Identified_Linelist.txt'

# Intervening absorption line list is optional to show any intervening absorbers

#Run the main GUI
M.Transitions(Abs,intervening=filename)
```

Second Example shows how to load up a saved analysis file

```python

#First read in the saved file
import pickle
pfile='Spectrum_Analysis_z_0.017.p'
with open(pfile,'rb') as pklfile:
    absys=pickle.load(pklfile)

#Run the Master GUI
from GUIs.abstools import Metal_Plot as M   
M.Transitions(absys)
```


-----------------------------------------------------------------------------


Dependencies:- PyQt5, astropy, rbcodes, numpy, matpltolib, pickle

May 2021.

Brought to you by Sean Clark and Rongmon Bordoloi.

Developed from a basic original code from Tom Cooper.
