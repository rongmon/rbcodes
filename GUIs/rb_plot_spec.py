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
import ipdb
import pandas as pd

class rb_plot_spec(object):

    def __init__(self,wave,flux,error,zabs=0.):

        '''
               ---------------------------------------------------------------------------
        This is an interactive 1D spectrum viewer.
        The help scene activates by pressing h on the plot.


        The program only works properly if none of the toolbar buttons in the figure is activated. 
        It also needs pysimpleGUI code to be installed. 
        https://pysimplegui.readthedocs.io/en/latest/





        Useful Keystrokes:

            Keystrokes:
              
              r        :    Reset Spectrum and replot to default settings.
              h        :    Prints this help window.
              x or X   :    Set xmin, xmax
              b or t   :    Set ymin, ymax
              [ or ]   :    Pan left or right 
              s or S.  :    Smooth or Unsmooth spectra
              E        :    Two E keystrokes will compute rest frame equivalent width at a defined region
              F        :    Three keystrokes to fit a Gaussian profile. [Currently not drawing on the spectrum]

              #GUI ELEMENTS [WARNING UNSTABLE]
              Works with TkAGG backend and pysimplegui

              Z  :   pop up window to select absorber redshift and linelist
              j  :   pop up window to select a corresponding rest frame transition and linelist
              K  :   pop up window to select multiple absorber lines and plot them
              0  :   pop up window to select identified absorber list to show with 1d spectrum

              q     :    Quit Program.
         ---------------------------------------------------------------------------
        Written By:  Rongmon Bordoloi                                   August 2020.

        HEALTH WARNING: The GUI implementation is still in alpha version and is quite unstable.
        User must be careful to make sure that they exit individual GUIs first by pressing the correct button
        before closing the plot window.
        ''' 


        self.wave=wave
        self.flux=flux
        self.smoothed_spectrum=flux
        self.error=error
        self.zabs=zabs
        self.label='None' # Initializing a label
        
        fig, ax = plt.subplots(1,figsize=(20, 5))
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

        #create indentified line list
        # Very basic window.  Return values using auto numbered keys
        d={'zabs':[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],'List':['None','None','None','None','None','None','None','None','None','None'], 'color':['None','None','None','None','None','None','None','None','None','None']} 
        df=pd.DataFrame(data=d)    

        self.zabs_list=df

        plt.gcf().canvas.mpl_connect('key_press_event',self.ontype)
        plt.draw()
        plt.show()


    def ontype(self,event):
        zabs=np.double(0.)
        self.ax.draw_artist(self.ax) 
        # when the user hits 'r': clear the axes and plot the original spectrum
        if event.key=='r':
            self.ax.cla()
            self.ax.step(self.wave,self.flux,'k-',linewidth=1)
            self.ax.set_xlabel('Wavelength')
            self.ax.set_ylabel('Flux')
            xr=[np.min(self.wave),np.max(self.wave)]
            yr=[np.min(self.flux),np.max(self.flux)]
            self.ax.set_ylim([yr[0],yr[1]])
            self.ax.set_xlim([xr[0], xr[1]])
    
        # Set top y max
        elif event.key=='t':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_ylim([ylim[0],event.ydata])
            self.ax.set_xlim(xlim)
            plt.draw()
        elif event.key=='Z':
            print('Testing to see if we can get a pop up window to work')
            self.set_redshift()
            self.DrawLineList(self.label)

        elif event.key =='j':
            lambda_rest, LineList=self.identify_line_GUI()
            print(lambda_rest,LineList)
            self.zabs= (event.xdata -lambda_rest)/lambda_rest
            self.label=LineList
            self.DrawLineList(self.label)
            print('Target Redshfit set at : ' + np.str(self.zabs))

        # Manage and plot multiple absorbers
        #Clunky GUI to plot lines
        elif event.key =='K':
            self.manage_identified_absorbers()

        #Load a saved linelist
        elif event.key =='0':
            self.load_linelist_GUI()
            self.manage_identified_absorbers()



        # Set top y min
        elif event.key=='b':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_ylim([event.ydata,ylim[1]])
            self.ax.set_xlim(xlim)
            plt.draw()


        # Smooth spectrum
        elif event.key=='S':
            self.vel[0] += 2
            Filter_size=np.int(self.vel[0]) 
            self.smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            self.specplot()
            plt.draw()

        #Unsmooth Spectrum
        elif event.key=='U':
            self.vel[0] -= 2
            if self.vel[0] <= 0:
                self.vel[0]=1;
            Filter_size=np.int(self.vel[0]) 
            self.smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            self.specplot()
            plt.draw()

            # Set X max
        elif event.key=='X':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_xlim([xlim[0],event.xdata])
            self.ax.set_ylim(ylim)
            plt.draw()
        # Set x min
        elif event.key=='x':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_xlim([event.xdata,xlim[1]])
            self.ax.set_ylim(ylim)
            plt.draw()


        # Set pan spectrum
        elif event.key==']':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            delx=(xlim[1]-xlim[0])
            self.ax.set_xlim([xlim[1],xlim[1]+delx])
            self.ax.set_ylim(ylim)
            plt.draw()

        # Set pan spectrum
        elif event.key=='[':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            delx=(xlim[1]-xlim[0])
            self.ax.set_xlim([xlim[0]-delx,xlim[0]])
            self.ax.set_ylim(ylim)
            plt.draw()

        # Compute Equivalent Width between two points
        elif event.key=='E':
            #ekeycounts +=1
            self.lam_lim=np.append(self.lam_lim,event.xdata)
            self.lam_ylim=np.append(self.lam_ylim,event.ydata)

    
            # Keep running tab of all E clicks
            eclick=len(self.lam_lim);

            #self.specplot()

            self.ax.plot(event.xdata,event.ydata,'rs',ms=5,picker=5,label='EW_pt',markeredgecolor='k')
            #self.fig.canvas.draw()


            if eclick==2:
                # Check if the wave entries are monotonously increasing
                tab=self.lam_lim.argsort()


                EW,sig_EW,cont,wave_slice=self.compute_EW(self.wave/(1.+zabs),self.flux,self.lam_lim[tab],self.lam_ylim[tab],self.error)
                EW=np.array(EW)*1000.
                sig_EW=np.array(sig_EW)*1000.
                self.ax.plot(wave_slice,cont,'r--')
                #ipdb.set_trace()
                print('---------------------- Equivalent Width -------------------------------------')
                Wval='EW [mAA]: '+ '%.1f' % EW + ' +/- ' + '%.1f' % sig_EW
                self.ax.text(np.mean([self.lam_lim]),np.max(self.lam_ylim)+0.2,Wval, rotation=90,verticalalignment='bottom')
                print(Wval)
                print('---------------------------------------------------------------------------')
                self.lam_lim=[]
                self.lam_ylim=[]

            plt.draw()
            self.fig.canvas.draw()



    
                
        
        # Fit a Gaussian
        elif event.key=='F':
            self.FXval=np.append(self.FXval, event.xdata)
            self.FYval=np.append(self.FYval, event.ydata)
    
            fclick=len(self.FXval)    
            self.ax.plot(event.xdata,event.ydata,'rs',ms=5,picker=5,label='EW_pt',markeredgecolor='k')
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


            #If the user presses 'h': The help is printed on the screen
    elif event.key=='h':
        print(
        '''    
        ---------------------------------------------------------------------------
        This is an interactive 1D spectrum viewer.
        The help scene activates by pressing h on the plot.


        The program only works properly if none of the toolbar buttons in the figure is activated. 
        It also needs pysimpleGUI code to be installed. 
        https://pysimplegui.readthedocs.io/en/latest/





        Useful Keystrokes:

            Keystrokes:
              
              r        :    Reset Spectrum and replot to default settings.
              h        :    Prints this help window.
              x or X   :    Set xmin, xmax
              b or t   :    Set ymin, ymax
              [ or ]   :    Pan left or right 
              s or S.  :    Smooth or Unsmooth spectra
              E        :    Two E keystrokes will compute rest frame equivalent width at a defined region
              F        :    Three keystrokes to fit a Gaussian profile. [Currently not drawing on the spectrum]

              #GUI ELEMENTS [WARNING UNSTABLE]
              Works with TkAGG backend and pysimplegui

              Z  :   pop up window to select absorber redshift and linelist
              j  :   pop up window to select a corresponding rest frame transition and linelist
              K  :   pop up window to select multiple absorber lines and plot them
              0  :   pop up window to select identified absorber list to show with 1d spectrum

              q     :    Quit Program.
         ---------------------------------------------------------------------------
        Written By:  Rongmon Bordoloi                                   August 2020.

        HEALTH WARNING: The GUI implementation is still in alpha version and is quite unstable.
        User must be careful to make sure that they exit individual GUIs first by pressing the correct button
        before closing the plot window. 

        '''
        )



         # Making sure any drawn line list remains drawn
        self.DrawLineList(self.label)
        q=np.where(self.zabs_list['List'] != 'None')
        if len(q[0] >0): 
            self.draw_any_linelist()
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
                    ss=self.ax.transData.transform((0, .9))
                    ydata=[0,ylim[1]]
                    lineplot,=self.ax.plot(xdata,ydata,'k--',)                    
                    tt=self.ax.text(xdata[0],0.75*ylim[1],data[i]['ion']+' '+ np.str(self.zabs),rotation=90) 
        #print('Target Redshfit set at : ' + np.str(self.zabs))
        plt.draw()
    #This code will draw any linelist identified
    
    def draw_any_linelist(self):
        # Now make a smaller linelist within the current xaxes range.
        #plt.figure(self.fig.number)
        xlim=self.ax.get_xlim()
        ylim=self.ax.get_ylim()

        #Hardcoding that only 10 independent absorber systems can be plotted
        for i in range(0,10):
            if self.zabs_list['List'][i] != 'None':
                zabs=np.double(self.zabs_list['zabs'][i])
                linelist=self.zabs_list['List'][i]
                lineclr=self.zabs_list['color'][i]

                data=line.read_line_list(linelist)
                for i in range(0, len(data)):
                    if ((data[i]['wrest']*(1.+zabs) >= np.double(xlim[0])) & (data[i]['wrest']*(1.+zabs) <= np.double(xlim[1]))):
                        xdata=[data[i]['wrest']*(1.+zabs),data[i]['wrest']*(1.+zabs)]
                        ss=self.ax.transData.transform((0, .9))
                        ydata=[0,ylim[1]]
                        lineplot_list,=self.ax.plot(xdata,ydata,'--',color=lineclr)                    
                        tt_list=self.ax.text(xdata[0],0.75*ylim[1],data[i]['ion']+' '+ np.str('%.5f' %zabs),rotation=90) 
        #print('Target Redshfit set at : ' + np.str(self.zabs))
        plt.draw()

    def manage_identified_absorbers(self):
        sg.ChangeLookAndFeel('Dark')   
        col1=[ [sg.Text('1. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][0]),  size=(15, 1))],
               [sg.Text('2. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][1]),  size=(15, 1))],
               [sg.Text('3. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][2]),  size=(15, 1))],
               [sg.Text('4. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][3]),  size=(15, 1))],
               [sg.Text('5. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][4]),  size=(15, 1))],
               [sg.Text('6. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][5]),  size=(15, 1))],
               [sg.Text('7. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][6]),  size=(15, 1))],
               [sg.Text('8. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][7]),  size=(15, 1))],
               [sg.Text('9. zabs', size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][8]),  size=(15, 1))],
               [sg.Text('10. zabs',size=(5, 1)), sg.In(default_text=np.str(self.zabs_list['zabs'][9]),  size=(15, 1))]]

        col2= [[sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][0], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][1], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][2], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][3], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][4], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][5], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][6], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][7], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][8], size=(10,1))],
               [sg.Text('LineList', size=(5, 1)),sg.Spin(values=('None','LLS', 'LLS Small', 'DLA'), initial_value=self.zabs_list['List'][9], size=(10,1))]]


        col3=  [[sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][0] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][1] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][2] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][3] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][4] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][5] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][6] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][7] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][8] ,size=(5, 1))],
                [sg.Text('color', size=(5, 1)), sg.In(default_text=self.zabs_list['color'][9] ,size=(5, 1))]]
    



        layout = [[sg.Column(col1),sg.Column(col2),sg.Column(col3)], [sg.Submit(), sg.Exit(),sg.Button('Reset')]]


        window = sg.Window('Update Selected Absorbers', layout, font=("Helvetica", 12))

        while True:
            event, values = window.read()
            #update database
            for i in range(0,10):
                self.zabs_list.at[i, 'zabs'] = values[i]
                self.zabs_list.at[i, 'List'] = values[i+10]
                self.zabs_list.at[i, 'color'] = values[i+20]   
            print(self.zabs_list)
            self.specplot()
            self.draw_any_linelist()
            if event == sg.WIN_CLOSED or event == 'Exit': 
                break 
            elif event =='Reset':
                self.specplot()
        window.close()
 






    def plot_keystroke(self,event):
        xdata=event.xdata
        ydata=event.ydata
        test,=self.ax.plot(xdata,ydata,'r+',)
        print(xdata,ydata)
        plt.draw()
        self.fig.canvas.draw()



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
            [sg.Button('Reset'), sg.Button('Submit')]]


        window = sg.Window('Line Identification', layout,font=("Helvetica", 12))  


        while True:             # Event Loop  
            event, values = window.Read()  

            if event is None or event == 'Submit':  
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
        #ipdb.set_trace()
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
        ax=plt.gca()
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        ax.cla()
        ax.step(self.wave,self.smoothed_spectrum,'k-',lw=1,label='smooth')

        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        
        self.ax=ax


    def load_linelist_GUI(self):
        sg.ChangeLookAndFeel('Dark')

        event, values = sg.Window('Load Identifed Line list', [[sg.Text('Filename')], [sg.Input(), sg.FileBrowse()], [sg.OK(), sg.Cancel()] ]).read(close=True)
        tt=ascii.read(values[0])
        #count how many entries
        n_abs=len(tt['zabs'])
        # Now load the first 10 absorbers into the identifed linelist database

        if n_abs>10:
            for i in range(0,10):
                self.zabs_list.at[i, 'zabs'] = tt['zabs'][i]
                self.zabs_list.at[i, 'List'] = tt['List'][i]
                self.zabs_list.at[i, 'color'] = tt['color'][i]   

        else:
            for i in range(0,n_abs):
                self.zabs_list.at[i, 'zabs'] = tt['zabs'][i]
                self.zabs_list.at[i, 'List'] = tt['List'][i]
                self.zabs_list.at[i, 'color'] = tt['color'][i]   

 





    
    def compute_EW(self,lam,flx,lam_lim,lam_ylim,err_flx):
        qtq=np.where( ( lam >= lam_lim[0] ) & ( lam <= lam_lim[1] ) )
        ww=lam[qtq]                 
        flux1=flx[qtq]
        spline = splrep(lam_lim,lam_ylim,k=1)
        continuum = splev(ww,spline)
        #ax.plot(ww,continuum,'k--')
        sig_flx1=err_flx[qtq]/continuum
        flux1=flux1/continuum
    
        EW=np.trapz(1.-flux1, x= ww)
        delw=np.double(ww[2]-ww[1])
        sig_w=delw*sig_flx1
        sig_wtot=np.sqrt(np.sum(sig_w**2.))
        return EW,sig_wtot,continuum,ww








    
def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(sqrt(2*pi)*wid)) * exp(-(x-cen)**2 /(2*wid**2))
    
    