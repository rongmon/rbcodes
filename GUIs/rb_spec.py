""" Spectrum class to read in, analyze and measure absorption lines."""

import numpy as np
from scipy.interpolate import splrep,splev
import sys
import os
import pdb
class read_spec(object):
    """A spectrum read into a class. spectrum will have
    following properties:

    Attributes:
        wave: wavelength.
        flux: flux.
        error: error
        filename=filename and location
        filetype = False [default] : other options 
                 ascii, fits, HSLA, xfits, p [pickle], temp, and linetools [uses linetools.io routine for this]

    Optional:  All only valid for filetype=linetools option
         efil= errorfile [Default None]

        output: This gives a read_spec object with following attributes:

        self.zabs= Absorber redshift

        self.wave_slice= sliced observed wavelength vector
        self.flux_slice= sliced observed flux vector
        self.error_slice= sliced velocity spectra 
        self.linelist=. LineList used

        
        self.velo=  sliced velocity vector



        self.cont = Fitted continuum
        self.fnorm= Normalized flux
        self.enorm= Normalized error

        self.trans=  Name of the Transition
        self.fval= fvalue of transition
        self.trans_wave= rest frame wavelength of transition
        self.vmin=     velocity minimum used for equivalent width calculation
        self.vmax=    velocity maximum used for equivalent width calculation
        self.W=    Rest Frame equivalenth width
        self.W_e=  uncertainty on rest frame equivalent width
        self.logN=  AOD column density
        self.logN_e= AOD column density uncertainty
        self.Tau= Apparant optical depth as a function of velocity



    Written : Rongmon Bordoloi      April 2018
    Edit    : Rongmon Bordoloi      September 2018 Changed kwargs to be compatible to python 3   
    Edit    : Rongmon Bordoloi      Aug 2020: added linetools.io.readspec file

    --------------------------------------------------------------------------------------------
    EXAMPLE: import numpy as np
            import matplotlib
            matplotlib.use('TkAgg')
            import matplotlib.pyplot as plt
            from GUIs import rb_spec as r 

            # List of absorber redshifts
            zabs=[0.511020,1.026311,1.564481]
            transition= 2796.3
            #Which absorber to analyze
            index=0
            filename='Quasar_Spectrum.fits'
            #Read in file
            s=r.read_spec(filename,filetype='linetools')

            #Shift spectra to rest frame
            s.shift_spec(zabs[index]);
            #Velocity window around transition
            xlim=[-1500,1500]
            
            #Slice Spectrum within that window
            s.slice_spec(transition,xlim[0],xlim[1],use_vel=True);

            #Fit continuum Mask the regions defined by velocity
            s.fit_continuum(mask=[-200,300,500,1100],domain=xlim,Legendre=3)

            #Compute EW
            #Compute equivalent width within a velocity window
            s.compute_EW(transition,vmin=-200.,vmax=360.);

            #save everything as a pickle file
            s.save_slice('outfile.p')

            #plot the Full spectrum
            s.plot_spec()

            #plot the sliced and fitted continuum
            #Plot stuff
            plt.subplot(2,1,1)
            plt.step(s.velo,s.flux_slice)
            plt.step(s.velo,s.flux_slice/s.fnorm)
            plt.xlim(xlim)
            plt.subplot(2,1,2)
            plt.step(s.velo,s.fnorm)
            plt.plot([-1500,1500],[1,1],'--')
            plt.xlim(xlim)
            plt.show()
    --------------------------------------------------------------------------------------------

    """
    def __init__(self,filename,filetype=False, efil=None,**kwargs):
        """ creates the spectrum object """
        self.filename=filename

        if filetype==False:
            #Take File Extention and try
            tt=os.path.splitext(filename)[1]
            if (tt=='txt')| (tt=='dat'):
                filetype='ascii'
            else:
                filetype=tt[1:len(tt)] 


        # Read in Files in differet formats
        if filetype=='ascii':
            from astropy.io import ascii
            dat=ascii.read(filename)
            tab=dat.keys()
            wave=np.array(dat[tab[0]])
            flux=np.array(dat[tab[1]])
            if (len(dat.keys())>=3):
                error=dat[tab[2]]
            else:
                error=0.*flux

        elif filetype=='fits':
            from astropy.io import fits
            file=fits.open(filename)#(cwd+'/'+filename)
            dat=file[1].data
            tab=dat.names
            wave=np.array(dat['wave'][0])
            flux=np.array(dat['flux'][0])
            if (len(tab)>=3):
                error=np.array(dat['error'][0])
            else:
                error=0.*flux
        elif filetype=='HSLA':
            from astropy.io import fits
            file=fits.open(filename)#(cwd+'/'+filename)
            dat=file[1].data
            tab=dat.names
            wave=np.array(dat['WAVE'])
            flux=np.array(dat['FLUX'])#/np.median(np.array(dat['FLUX']))
            if (len(tab)>=3):
                error=np.array(dat['ERROR'])#/np.median(np.array(dat['FLUX']))
            else:
                error=0.*flux


        if filetype=='xfits':
            from linetools.spectra.xspectrum1d import XSpectrum1D  
            sp=XSpectrum1D.from_file(filename)

            wave=sp.wavelength.value
            flux=sp.flux.value
            error=sp.sig.value

        elif filetype=='p':
            import pickle
            dat=pickle.load( open(filename, "rb" ))
            #tab=dat.keys()
            wave=np.array(dat['wave'])
            flux=np.array(dat['flux'])
            if (len(tab)>=3):
                error=np.array(dat['error'])
            else:
                error=0.*flux
        elif filetype=='temp':
            from astropy.io import fits
            file=fits.open(filename)#(cwd+'/'+filename)
            #a=fits.open(path+'spec_knotA.fits')

            wave=file[2].data
            flux=file[0].data
            error=file[1].data
        #Use linetools.io.readspec to read file
        elif filetype =='linetools':
            from linetools.spectra import io as tio
            sp=tio.readspec(filename,inflg=None, efil=efil,**kwargs)
            wave=sp.wavelength.value
            flux=sp.flux.value

            if sp.sig_is_set == False:
                print('Assuiming arbiarbitrary 10% error on flux')
                error=0.1*flux
            else:
                error=sp.sig.value


        

        self.wave=wave
        self.flux=flux
        self.error=error
        self.wrest=wave*(1.+0.)
        self.zabs=0.
    
    def shift_spec(self,zabs):
        """ Shifts wavelength to absorber rest frame"""
        self.wrest=self.wave/(1.+zabs)
        self.zabs=zabs
        return self.wrest, self.zabs

    def slice_spec(self,lam_rest,lam_min,lam_max,method='closest',linelist='LLS',use_vel=False):
        """
        Slice the spectrum around a central wavelength and convert it to velocity space
        lam_rest : approximate rest wavelength of a transition
        lam_min  : minimum wavelength/velocity to slice 
        lam_max  : maximum wavelength/velocity to slice 

        Keywords:   method = 'closest' [default] -> sets lam_rest to closest atomic transition
                    method = 'Exact' -> uses given lam_rest value to look for transition
                    linelist= Default LLS line linelist, otherwise uses the specified line list

                    use_vel = True -> uses velocity space to slice.
                                   here inputs are lam_min = vel_min [in km/sec]
                                                   lam_max =vel_max [km/s]


        """
        '''
        if kwargs.has_key('method'):
            method=kwargs['method']
        else:
            method='closest'

        if kwargs.has_key('linelist'):
            linelist=kwargs['linelist']
        else:
            linelist='LLS'
        '''

        from IGM import rb_setline as s

        str=s.rb_setline(lam_rest,method,linelist=linelist)

        spl=2.9979e5;  #speed of light
        vel = (self.wrest-str['wave']*(1.0 + 0.))*spl/(str['wave']*(1.0 + 0.))

        # Now slice spectrum either velocity or wave space
        if use_vel==False:
            q=np.where((self.wrest >= lam_min) & (self.wrest <= lam_max))
        else:
            vel_min=lam_min
            vel_max=lam_max
            q=np.where((vel >= vel_min) & (vel <= vel_max))

        self.wave_slice=self.wrest[q]
        self.flux_slice=self.flux[q]
        self.error_slice=self.error[q]
        self.linelist=linelist

        
        self.velo=vel[q]

        return self.wave_slice,self.error_slice,self.flux_slice,self.velo,self.linelist

    def fit_continuum(self,mask=False,domain=False,Legendre=False,prefit_cont=[1.]):
        """
        By default calls an interactive continuum fitter to the sliced spectrum.
        Or an automated Legendre polynomial fitter if keyword set Legendre.
            Order is given by Legendre=order

        """
        if Legendre==False:
            #pdb.set_trace()
            if prefit_cont == 'False':
                from GUIs import rb_fit_interactive_continuum as f
                s=f.rb_fit_interactive_continuum(self.wave_slice,self.flux_slice,self.error_slice)
                cont=s.cont
            else:
                cont=prefit_cont



        else:
            import numpy.polynomial.legendre as L
            order=Legendre
            weight= 1./(self.error_slice**2.)
            
            if domain==False:
                domain=[-600.,600.]

            if mask==False:
                #No Mask
                q=0.*self.wave_slice +1. 

            else:
                #mask=kwargs['mask']            
                #Number of Masks 
                nmsk=int(len(mask)/2)
                vmin=np.zeros(nmsk,)
                vmax=np.zeros(nmsk,)

                for i in range(0,nmsk):
                    vmin[i]=mask[2*i]
                    vmax[i]=mask[2*i+1]

                q=0.*self.wave_slice +1. 

                for i in range(0,nmsk):
                    sq=np.where((self.velo >= vmin[i]) & (self.velo <= vmax[i]))
                    q[sq]=0.
            #Select Unmasked part of the spectrum
            qtq=np.where((q ==1))

            # Fitting the masked Data
            e=L.Legendre.fit(self.velo[qtq],self.flux_slice[qtq],order,w=weight[qtq],domain=domain);
            cont=e(self.velo)
            # Now mask the part of spectrum that we don't want to fit. 
            # Mask is created to have multiple low vel, high vel ranges.
            # e.g. mask = [-300.,-250.,100.,120.] will exclude -300,-250 and 100,120 km/s parts of the spectrum in the fit


        self.cont=cont
        self.fnorm=self.flux_slice/cont
        self.enorm=self.error_slice/cont



        return self.cont,self.fnorm,self.enorm

    def compute_EW(self,lam_cen,vmin=-50.,vmax=50.,method='closest',plot=False):
        """
        Computes rest frame equivalent width and column density for a desired atomic line.
        Around the species lam_cen and given vmin and vmax keyword values. 

        """

        from IGM import rb_setline as s
        str=s.rb_setline(lam_cen,method,linelist=self.linelist)

        from IGM import compute_EW as EW
        out=EW.compute_EW(self.wave_slice,self.fnorm,str['wave'],[vmin,vmax],self.enorm,f0=str['fval'],zabs=0.,plot=plot)

        self.trans=str['name']
        self.fval=str['fval']
        self.trans_wave=str['wave']
        self.vmin=vmin
        self.vmax=vmax


        self.W= out['ew_tot']
        self.W_e=out['err_ew_tot']
        self.logN=out['col']
        self.logN_e=out['colerr']

        self.Tau=out['Tau_a']

        return self.trans,self.fval,self.vmin,self.vmax,self.trans_wave,self.W,self.W_e,self.logN,self.logN_e,self.Tau


    def plot_spec(self):
        """
        Quick wrapper to call an interactive plotter for the full spectrum as given in input file.
        """
        from GUIs import rb_plot_spec as sp
        tt=sp.rb_plot_spec(self.wave,self.flux,self.error)

    def save_slice(self,outfilename):
        """
        Saves the slice object for future processing.
        """
        import pickle
        with open(outfilename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)


    def plot_doublet(self,lam1,lam2,vmin=-600.,vmax=600.,method='closest'):
        """
        Plot a given doublet defined by the lam1 and lam2 wavelength centers.
        """

        from IGM import rb_setline as s
        str1=s.rb_setline(lam1,method,linelist=self.linelist)
        str2=s.rb_setline(lam2,method,linelist=self.linelist)

        spl=2.9979e5;  #speed of light
        vel1 = (self.wave_slice-str1['wave']*(1.0 + 0.))*spl/(str1['wave']*(1.0 + 0.))
        vel2 = (self.wave_slice-str2['wave']*(1.0 + 0.))*spl/(str2['wave']*(1.0 + 0.))
        
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1=fig.add_subplot(211)
        ax1.step(vel1,self.fnorm)
        ax1.step(vel1,self.enorm,color='r')

        ax1.set_xlim([vmin,vmax])
        ax1.set_ylim([-0.02,1.8])
        ax1.plot([-2500,2500],[0,0],'k:')
        ax1.plot([-2500,2500],[1,1],'k:')       
        ax1.set_xlabel('vel [km/s]')
   
        ax2=fig.add_subplot(212)
        ax2.step(vel2,self.fnorm)
        ax2.step(vel2,self.enorm,color='r')

        ax2.set_xlim([vmin,vmax])
        ax2.set_ylim([-0.02,1.8])
        ax2.plot([-2500,2500],[0,0],'k:')
        ax2.plot([-2500,2500],[1,1],'k:')       
        ax2.set_xlabel('vel [km/s]')
        plt.show()

    def vpfit_singlet(self):
        from GUIs import rb_interactive_vpfit_singlet as vf 
        vt=vf.rb_interactive_vpfit_singlet(self.wave_slice,self.fnorm,self.enorm,self.velo);    






