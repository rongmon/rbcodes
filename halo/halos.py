#!/usr/bin/python
import sys
import math
import re
from scipy.interpolate import interp1d
from astropy.cosmology import Planck18
from astropy.constants import G
import astropy.units as u
from numpy import log10

"""NB: Stellar-mass to halo-mass relations from Behroozi et al 
Universe Machine. It is possible to calculate scatter in the
Stellar to halo mass relation but shoulf assume that in the
mass range of interest, the scatter in this correlation is
about 0.1 dex (in reality, 0.04-0.13 dex)

example
from rbcode.halo import halos as h
import numpy as np
import matplotlib.pyplot as plt

#define stellar mass in log
sm=np.arange(8,10,.2)



hm1=np.zeros(len(sm),)
hm2=np.zeros(len(sm),)
hm3=np.zeros(len(sm),)
hm4=np.zeros(len(sm),)

z1=np.ones(len(sm),)*5.3
z2=np.ones(len(sm),)*5.6
z3=np.ones(len(sm),)*6
z4=np.ones(len(sm),)*6.5

for i in range(0,len(sm)):
    hm1[i]=h.stellarToHaloMass(z1[i],sm[i])    
    hm2[i]=h.stellarToHaloMass(z2[i],sm[i])
    hm3[i]=h.stellarToHaloMass(z3[i],sm[i])
    hm4[i]=h.stellarToHaloMass(z4[i],sm[i])




plt.plot(sm,hm1,'.',label='z=5.3')
plt.plot(sm,hm2,'.',label='z=5.6')
plt.plot(sm,hm3,'.',label='z=6')
plt.plot(sm,hm4,'.',label='z=6.5')

plt.xlabel('log Stellar Mass')
plt.ylabel('log halo Mass')

plt.show()


"""
def stellarToHaloMass(z,stellar_mass):
    
    #Load params
    param_file = open('smhm_med_params.txt', "r")
    param_list = []
    allparams = []
    for line in param_file:
        param_list.append(float((line.split(" "))[1]))
        allparams.append(line.split(" "))
        
    if (len(param_list) != 20):
        print("Parameter file not correct length.  (Expected 20 lines, got %d)." % len(param_list))
        quit()

    if (stellar_mass > 20):
        stellar_mass = log10(stellar_mass)
        
    names = "EFF_0 EFF_0_A EFF_0_A2 EFF_0_Z M_1 M_1_A M_1_A2 M_1_Z ALPHA ALPHA_A ALPHA_A2 ALPHA_Z BETA BETA_A BETA_Z DELTA GAMMA GAMMA_A GAMMA_Z CHI2".split(" ");
    params = dict(zip(names, param_list))

    #Print SMHM relation
    a = 1.0/(1.0+z)
    a1 = a - 1.0
    lna = math.log(a)
    zparams = {}
    zparams['m_1'] = params['M_1'] + a1*params['M_1_A'] - lna*params['M_1_A2'] \
        + z*params['M_1_Z']
    zparams['sm_0'] = zparams['m_1'] + params['EFF_0'] + a1*params['EFF_0_A'] \
        - lna*params['EFF_0_A2'] + z*params['EFF_0_Z']
    zparams['alpha'] = params['ALPHA'] + a1*params['ALPHA_A'] \
        - lna*params['ALPHA_A2'] + z*params['ALPHA_Z']
    zparams['beta'] = params['BETA'] + a1*params['BETA_A'] + z*params['BETA_Z']
    zparams['delta'] = params['DELTA']
    zparams['gamma'] = 10**(params['GAMMA'] + a1*params['GAMMA_A'] + z*params['GAMMA_Z'])
    
    smhm_max = 14.5-0.35*z
    if (params['CHI2']>200):
        print('#Warning: chi^2 > 200 implies that not all features are well fit.')
        print('Comparison with the raw data (in data/smhm/median_raw/) is crucial.')

    Mhalo = [x*0.05 for x in range(int(10.5*20),int(smhm_max*20+1),1)]
    Mstar = []
        
    for m in [x*0.05 for x in range(int(10.5*20),int(smhm_max*20+1),1)]:
        dm = m-zparams['m_1'];
        dm2 = dm/zparams['delta'];
        sm = zparams['sm_0'] - math.log10(10**(-zparams['alpha']*dm) + \
                                          10**(-zparams['beta']*dm)) + \
                                          zparams['gamma']*math.exp(-0.5*(dm2*dm2));
        Mstar.append(sm)
        # print("%.2f %.6f %.6f" % (m,sm,sm-m))

    f = interp1d(Mstar, Mhalo, fill_value='extrapolate')
    Mh = f(stellar_mass).flatten()

    if (len(Mh) == 1):
        return(Mh[0])
    else:
        return(Mh)

##########################################
#
#  Takes a linear (not log) halo mass and returns Virial radius

def R_200(Mh, z):

    #if (Mh > 20):
    #    Mh = log10(Mh)
    
    H = Planck18.H(z)
    r200 = (((G * (10**Mh * u.M_sun).to(u.kg) / (100 * H.cgs**2))**(1/3)).to(u.kpc)).flatten()

    if (len(r200) == 1):
        return(r200[0])
    else:
        return(r200)




    
# Osterbrock 2006
# L_Ha = 1.37e-12 * Qion
# SFR_Ha = 5.5e-42 * L_Ha (Calzetti et al 2012)
# L(Ha) = 2.86 * L(Hb) (Case B)

# Qion = 7.3e11 * L(Ha) = 2.04e12 * L(Hb)
# SFR(Hb) = 1.57e-41 * L(Hb)

