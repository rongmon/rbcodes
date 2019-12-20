import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np
import matplotlib.pyplot as plt 
from rbvfit import model as m
from importlib import reload
reload(m)
from scipy.optimize import curve_fit
import pdb
'''
n_clouds = 1
logN=np.array([15.])
b=np.array([20.])
v=np.array([0.])

zabs     = np.zeros((n_clouds))
lambda_rest = 1526.9 * np.ones((n_clouds))

s=m.create_voigt(zabs,lambda_rest)
theta=np.concatenate((logN,b,v))
outflux= s.model_flux(theta,wave)
'''


def vel2shift(vel):
    c = 299792.458;  #% speed of light [km/sec]
    Beta  = Vel/c;
    z = np.sqrt((1.+Beta)/(1.-Beta)) - 1.;
    return z



def quick_nv_estimate(wave,norm_flx,wrest,f0):
    # All done in rest frame
    spl=2.9979e5  #speed of light
    vel = (wave-wrest*(1.0 + 0.))*spl/(wrest*(1.0 + 0.))
    lambda_r=wave/(1+0.)    
    #compute apparent optical depth
    Tau_a =np.log(1./norm_flx);
    # REMEMBER WE ARE SWITCHING TO VELOCITY HERE
    del_vel_j = np.diff(vel);
    del_vel_j = np.append([del_vel_j[0]], del_vel_j)
    # Column density per pixel as a function of velocity
    nv = Tau_a / ((2.654e-15) * f0 * lambda_r) # in units cm^-2 / (km s^-1), SS91 
    n = nv * del_vel_j  # column density per bin obtained by multiplying differential Nv by bin width
    return vel, n 

def rb_interactive_vpfit_singlet(wave,flux,error,wrest,custom_guess=False,FWHM=15.):


    # Fix infinities

    sq=np.isnan(flux);
    flux[sq]=0;
    sqq=flux<=0;
    flux[sqq]=0;
    q=flux<=0;
    flux[q]=error[q];

    # Converting To velocities
    from IGM import rb_setline as s
    str=s.rb_setline(wrest,'closest','atom')
    spl=2.9979e5;  #speed of light
    vel = (wave-str['wave']*(1.0 + 0.))*spl/(str['wave']*(1.0 + 0.))


    vel,nv=quick_nv_estimate(wave,flux,str['wave'],str['fval']);

    # Start guessing initial fit interactively
    if custom_guess == False:
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.step(vel,flux)
        ax.step(vel,error,color='r')
        plt.xlim([np.min(vel),np.max(vel)])
        plt.ylim([-0.02,1.8])
        ax.plot([-2500,2500],[0,0],'k:')
        ax.plot([-2500,2500],[1,1],'k:')       
        ax.set_xlabel('vel [km/s]')
        ax.set_ylabel('Normalized Flux')
        ax.set_title('click to select line centers!')
        xx = plt.ginput(n=-1)
        print('Once you are done, press enter to continue!')
        n_clouds=len(xx)


        # Create the guesses for starting a fit
        bguess=np.zeros(n_clouds,)
        vguess=np.zeros(n_clouds,)
        nguess=np.zeros(n_clouds,)
        plt.close()

        for i in range(0,n_clouds):
            vguess[i]=xx[i][0]
            qq=np.where( (vel < vguess[i]+ 80.) & (vel > vguess[i]-80.))
            nguess[i]=np.log10(sum(nv[qq]))        

            #Now ask interactively for b values     
            prompt='Guess  b  for line ' +np.str(i+1)+ '/'+np.str(n_clouds) +', vel guess = ' + np.str('%.1f' % vguess[i])  +', col guess= '+ np.str('%.1f' % nguess[i])+ ': '
            tmp_b =  input(prompt)
            bguess[i]= tmp_b 


    # Now starting to set up the Voigt Profile model
    zabs  = np.zeros((n_clouds))
    lambda_rest = str['wave'] * np.ones((n_clouds))
    s=m.create_voigt(zabs,lambda_rest)

    theta=np.concatenate((nguess,bguess,vguess))
    outflux= s.model_flux(theta,wave)
    def test_fit(wave,*params):
        return  s.model_flux(params,wave)

    #pdb.set_trace()
    popt, pcov = curve_fit(test_fit, wave, flux, p0=(theta))
    print(popt)


    fig = plt.figure()
    ax=fig.add_subplot(111)
    ax.step(vel,flux)
    ax.step(vel,error,color='r')
    ax.plot(vel,outflux,color='k',lw=4)
    ax.plot(vel,s.model_flux(popt,wave),color='g',lw=6)
    plt.xlim([np.min(vel),np.max(vel)])
    plt.ylim([-0.02,1.8])
    ax.plot([-2500,2500],[0,0],'k:')
    ax.plot([-2500,2500],[1,1],'k:')       
    ax.set_xlabel('vel [km/s]')
    ax.set_ylabel('Normalized Flux')
    plt.show()





    
