import matplotlib
matplotlib.use('Tkagg')
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import ascii
from scipy.interpolate import splrep,splev
import sys
import os
from scipy.signal import medfilt
from astropy.convolution import convolve, Box1DKernel
from numpy import sqrt, pi, exp, linspace, loadtxt
from lmfit import  Model
from matplotlib.widgets import TextBox, RadioButtons
from pkg_resources import resource_filename


class rb_plot_spec(object):

    def __init__(self,wave,flux,error,zabs=0.):
        self.wave=wave
        self.flux=flux
        self.error=error
        self.zabs=zabs
        self.label='None' # Initializing a label
        
        fig, ax = plt.subplots(1)
        self.fig=fig
        plt.subplots_adjust(bottom=0.2)
        initial_text = str(self.zabs)

        spectrum, = ax.step(self.wave, self.flux, 'k-',lw=1)
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        xr=[min(self.wave),max(self.wave)]
        yr=[0.,np.median(flux)*2.5]
        ax.set_ylim(yr)
        ax.set_xlim(xr)
        self.ax=ax
        self.vel=np.array([1.])
        self.lam_lim=[]
        self.lam_ylim=[]
        self.FXval=[]
        self.FYval=[]

        # Box to input redshift
        axbox = plt.axes([0.1, 0.05, 0.2, 0.075])
        text_box = TextBox(axbox, 'Redshift', initial=initial_text)
        text_box.on_submit(self.set_redshift)


        # Radio Button for Line Selection
        rax = plt.axes([0.7, 0.01, 0.2, 0.1])
        radio = RadioButtons(rax,  ('None','LLS', 'LLS Small', 'DLA'))
        radio.on_clicked(self.DrawLineList)

    
        plt.gcf().canvas.mpl_connect('key_press_event',self.ontype)
        plt.draw()
        plt.show()


    def ontype(self,event):
        zabs=np.double(0.)
        ax=self.ax
            # when the user hits 'r': clear the axes and plot the original spectrum
        if event.key=='r':
            ax.cla()
            ax.step(self.wave/(1.+zabs),self.flux,'k-',linewidth=1)
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')
            xr=[np.min(self.wave),np.max(self.wave)]/(1.+zabs)
            yr=[np.min(self.flux)-0.2,np.max(self.flux)+.2]
            ax.set_ylim([yr[0],yr[1]])
            ax.set_xlim([xr[0], xr[1]])
    
        # Set top y max
        elif event.key=='t':
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.set_ylim([ylim[0],event.ydata])
            ax.set_xlim(xlim)
        elif event.key=='j':
            print('Testing to see if we can get a pop up window to work')
            #fig1, ax1 = plt.subplots()
            #print(event.xdata)
            ax.plot([event.xdata,event.xdata],[-10,10])
            #import ipywidgets as widgets
            #widgets.IntSlider()

        # Set top y min
        elif event.key=='b':
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.set_ylim([event.ydata,ylim[1]])
            ax.set_xlim(xlim)


        # Smooth spectrum
        elif event.key=='S':
            self.vel[0] += 2
            Filter_size=np.int(self.vel[0]) 
            smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.cla()
            ax.step(self.wave/(1.+zabs),smoothed_spectrum,'k-',lw=1,label='smooth')
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        #Unsmooth Spectrum
        elif event.key=='U':
            self.vel[0] -= 2
            if self.vel[0] <= 0:
                self.vel[0]=1;
            Filter_size=np.int(self.vel[0]) 
            smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.cla()
            ax.step(self.wave/(1.+zabs),smoothed_spectrum,'k-',lw=1,label='smooth')
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Flux')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            # Set X max
        elif event.key=='X':
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.set_xlim([xlim[0],event.xdata])
            ax.set_ylim(ylim)
            plt.draw()
        # Set x min
        elif event.key=='x':
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.set_xlim([event.xdata,xlim[1]])
            ax.set_ylim(ylim)

        # Set pan spectrum
        elif event.key==']':
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            delx=(xlim[1]-xlim[0])
            ax.set_xlim([xlim[1],xlim[1]+delx])
            ax.set_ylim(ylim)

        # Set pan spectrum
        elif event.key=='[':
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            delx=(xlim[1]-xlim[0])
            ax.set_xlim([xlim[0]-delx,xlim[0]])
            ax.set_ylim(ylim)

        # Compute Equivalent Width between two points
        elif event.key=='E':
            #ekeycounts +=1
            self.lam_lim=np.append(self.lam_lim,event.xdata)
            self.lam_ylim=np.append(self.lam_ylim,event.ydata)

    
            # Keep running tab of all E clicks
            eclick=len(self.lam_lim);

            ax.plot(np.array([event.xdata]),np.array([event.ydata]),'ro',ms=5,markeredgecolor='k')
            plt.draw()




            if eclick==2:
                # Check if the wave entries are monotonously increasing
                tab=self.lam_lim.argsort()


                EW,sig_EW=compute_EW(self.wave/(1.+zabs),self.flux,self.lam_lim[tab],self.lam_ylim[tab],self.error,self.ax)
                EW=np.array(EW)*1000.
                sig_EW=np.array(sig_EW)*1000.
                print('---------------------- Equivalent Width -------------------------------------')
                Wval='EW [mAA]: '+ '%.1f' % EW + ' +/- ' + '%.1f' % sig_EW
                ax.text(np.mean([self.lam_lim]),np.max(self.lam_ylim)+0.2,Wval, rotation=90,verticalalignment='bottom')
                print(Wval)
                print('---------------------------------------------------------------------------')
                plt.draw()


    
                self.lam_lim=[]
                self.lam_ylim=[]

        
        # Fit a Gaussian
        elif event.key=='F':
            self.FXval=np.append(self.FXval, event.xdata)
            self.FYval=np.append(self.FYval, event.ydata)
    
            fclick=len(self.FXval)    
            ax.plot(event.xdata,event.ydata,'rs',ms=5,picker=5,label='EW_pt',markeredgecolor='k')
            self.fig.canvas.draw()
            plt.show()

            #Start Fitting
            if fclick==3:
                # First fit a quick continuum
                qtq=np.where( ( (self.wave/(1.+zabs)) >= self.FXval[0] ) & ( (self.wave/(1.+zabs)) <= self.FXval[2] ) )
                ww=self.wave[qtq]/(1.+zabs)                 
                flux1=self.flux[qtq]
                spline = splrep(np.append(self.FXval[0],self.FXval[2]),np.append(self.FYval[0],self.FYval[2]),k=1)
                continuum = splev(ww,spline)
    
                # Check if it is an absorption or emission line
                if ((self.FYval[1] < self.FYval[0]) & (self.FYval[1] < self.FYval[2])):
                    ydata=1.- (flux1/continuum)
                    gmodel = Model(gaussian)
                    result = gmodel.fit(ydata, x=ww, amp=self.FYval[1]-self.FYval[1], cen=self.FXval[2], wid=0.5*(self.FXval[2]-self.FXval[0]))
                    Final_fit=(1.-result.best_fit)*continuum
                else:
                    ydata=(flux1/continuum)
                    gmodel = Model(gaussian)
                    result = gmodel.fit(ydata, x=ww, amp=self.FYval[1]-self.FYval[1], cen=self.FXval[2], wid=0.5*(self.FXval[2]-self.FXval[0]))
                    Final_fit=result.best_fit*continuum         
    
    
                print(result.fit_report())
                model_fit=self.ax.plot(ww, Final_fit, 'r-')
                plt.draw()
                print("Gaussian Fit")
    
                FXval=[]
                FYval=[]
         # Making sure any drawn line list remains drawn
        self.DrawLineList(self.label)
        plt.draw()
        self.fig.canvas.draw()

       

    def DrawLineList(self,label):
        del self.ax.lines[1:len(self.ax.lines)]
        self.ax.texts=[]
        self.label=label
        if label == 'None':
            line,=self.ax.plot([0],[0],'k--')
        else:
            from IGM import rb_setline as line        
            data=line.read_line_list(label)

            # Now make a smaller linelist within the current xaxes range.

            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()


            for i in range(0, len(data)):
                if ((data[i]['wrest']*(1.+self.zabs) >= np.double(xlim[0])) & (data[i]['wrest']*(1.+self.zabs) <= np.double(xlim[1]))):
                    xdata=[data[i]['wrest']*(1.+self.zabs),data[i]['wrest']*(1.+self.zabs)]
                    ydata=[0,2]
                    line,=self.ax.plot(xdata,ydata,'k--',)
                    ss=self.ax.transData.transform((0, .9))
                    tt=self.ax.text(xdata[0],0.85*ylim[1],data[i]['ion'],rotation=90) 
        plt.draw()              
    def set_redshift(self,text):
        self.zabs=np.double(text)
        self.DrawLineList(self.label)
        #xdata = self.wave/(1.+np.double(text))
        #self.spectrum.set_xdata(xdata)
        #self.ax.set_xlim(np.min(xdata), np.max(xdata))
        plt.draw()
        self.fig.canvas.draw()


    
def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(sqrt(2*pi)*wid)) * exp(-(x-cen)**2 /(2*wid**2))
    
    
    
def compute_EW(lam,flx,lam_lim,lam_ylim,err_flx,ax):
    qtq=np.where( ( lam >= lam_lim[0] ) & ( lam <= lam_lim[1] ) )
    ww=lam[qtq]                 
    flux1=flx[qtq]
    spline = splrep(lam_lim,lam_ylim,k=1)
    continuum = splev(ww,spline)
    ax.plot(ww,continuum,'k--')
    sig_flx1=err_flx[qtq]/continuum
    flux1=flux1/continuum

    EW=np.trapz(1.-flux1, x= ww)
    delw=np.double(ww[2]-ww[1])
    sig_w=delw*sig_flx1
    sig_wtot=np.sqrt(np.sum(sig_w**2.))
    return EW,sig_wtot