import matplotlib
matplotlib.use('TkAgg')
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
from pkg_resources import resource_filename
import PySimpleGUI as sg
from IGM import rb_setline as line       


class rb_plot_spec(object):

    def __init__(self,wave,flux,error,zabs=0.):
        self.wave=wave
        self.flux=flux
        self.smoothed_spectrum=flux
        self.error=error
        self.zabs=zabs
        self.label='None' # Initializing a label
        
        fig, ax = plt.subplots(1)
        self.fig=fig

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
        elif event.key=='Z':
            print('Testing to see if we can get a pop up window to work')
            self.set_redshift()
            self.ax.plot([event.xdata,event.xdata],[-10,10])

        elif event.key =='j':
            lambda_rest, LineList=self.identify_line_GUI()
            print(lambda_rest,LineList)
            self.zabs= (event.xdata -lambda_rest)/lambda_rest
            self.label=LineList
            self.DrawLineList(self.label)



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
            self.smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            self.specplot()
        #Unsmooth Spectrum
        elif event.key=='U':
            self.vel[0] -= 2
            if self.vel[0] <= 0:
                self.vel[0]=1;
            Filter_size=np.int(self.vel[0]) 
            self.smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            self.specplot()

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

            self.specplot()

            self.plot_keystroke(event)





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
            lineplot,=self.ax.plot([0],[0],'k--')
        else:                   
            data=line.read_line_list(label)

            # Now make a smaller linelist within the current xaxes range.

            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()


            for i in range(0, len(data)):
                if ((data[i]['wrest']*(1.+self.zabs) >= np.double(xlim[0])) & (data[i]['wrest']*(1.+self.zabs) <= np.double(xlim[1]))):
                    xdata=[data[i]['wrest']*(1.+self.zabs),data[i]['wrest']*(1.+self.zabs)]
                    ydata=[0,2]
                    lineplot,=self.ax.plot(xdata,ydata,'k--',)
                    ss=self.ax.transData.transform((0, .9))
                    tt=self.ax.text(xdata[0],0.85*ylim[1],data[i]['ion'],rotation=90) 
        plt.draw()

    def plot_keystroke(self,event):
        xdata=event.xdata
        ydata=event.ydata
        test,=self.ax.plot([xdata,0],[ydata,0],'ro',)


    def set_redshift(self):
        zabs,LineList=self.set_redshift_GUI()
        self.zabs=np.double(zabs)
        self.label=LineList
        self.DrawLineList(self.label)
        print('Target Redshfit set at : ' + np.str(self.zabs))
        #xdata = self.wave/(1.+np.double(text))
        #self.spectrum.set_xdata(xdata)
        #self.ax.set_xlim(np.min(xdata), np.max(xdata))
        plt.draw()
        self.fig.canvas.draw()

    def identify_line_GUI(self):
        if self.label == 'None':
            self.label='LLS Small'
        data=line.read_line_list(self.label)
        Transition_List=[]
        wavelist=[]
        for i in range(0,len(data)):
            Transition_List.append(data[i]['ion'])
            wavelist.append(data[i]['wrest'])

        layout = [
            [sg.Text('Please select the transition and LineList')],
            [sg.Listbox(values=Transition_List, size=(30, 6),key='_Transition_')],
            [sg.Text('LineList', size=(15, 1)),
                sg.Drop(values=(self.label,'LLS', 'LLS Small', 'DLA'),size=(15, 1), key='_Menu_')], 
            [sg.Button('Refresh'), sg.Button('Exit')]]


        window = sg.Window('Line Identification', layout,font=("Helvetica", 12))  


        while True:             # Event Loop  
            event, values = window.Read()  

            if event is None or event == 'Exit':  
                break
            self.label=values['_Menu_']  
            data=line.read_line_list(self.label)
            Transition_List=[]
            wavelist=[]
            for i in range(0,len(data)):
                Transition_List.append(data[i]['ion'])
                wavelist.append(data[i]['wrest'])
            window.Element('_Transition_').Update(Transition_List)  
        window.Close()



        Transition_List=np.array(Transition_List)
        wavelist=np.array(wavelist)

        Transition_rest=values['_Transition_']
        qq=np.where(Transition_List==Transition_rest)
        lambda_rest=wavelist[qq][0]
        LineList=values['_Menu_']
        return lambda_rest,LineList




    def set_redshift_GUI(self):
        layout = [
            [sg.Text('Please enter the desired redshift, LineList')],
            [sg.Text('Redshift', size=(15, 1)), sg.InputText(np.str(self.zabs))],
            [sg.Text('LineList', size=(15, 1)),
                sg.Drop(values=(self.label,'None','LLS', 'LLS Small', 'DLA'), auto_size_text=True)], 
            [sg.Submit(), sg.Cancel()]]

        window = sg.Window('Redshift Query', layout,font=("Helvetica", 12))   
        event, values = window.Read()
        window.Close()
        zabs=values[0]
        LineList=values[1]
        return zabs,LineList


    def specplot(self):
        ax=self.ax
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        ax.cla()
        ax.step(self.wave,self.smoothed_spectrum,'k-',lw=1,label='smooth')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        
        self.ax=ax








    
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