"""Interactive continuum fitter for jupyter notebook
"""
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
# Converting To velocities
from rbcodes.igm import rb_setline as line
from scipy.interpolate import splrep,splev
import sys
import os
from IPython.display import display


class interactive_cont(object):
    """This is an interactive continuum fitter for 1D spectrum.
    
        Attributes
        ----------
            wave :- numpy arary of wavelength
            flux :- numpy array of flux
            error:- numpy arrat of error [optional]
            xlim:- xrange optional
            **kwargs:- optional

        Returns
        -------
            Continuum fit: handfitted numpy array of continuum

        Example
        ------


            Useful Keystrokes:
    
                Mouse Clicks:
                
                    Left Click  : Select the median flux value within +/- 5 pixel from the x-coordinate.
                                  These points are used for the continuum fit.
                    Right Click : Delete the nearest continuum point.
    
                Keystrokes:
                  
                  b     :    Select a point for continuum fit at that exact (x,y) coordinate.
                  enter :    Perform a spline fit to data to create a continuum.
                  n     :    Show the normalized spectrum.
                  w     :    Only after pressing n: This will ourput the continuum. 
                  h     :    This Help screen.
                  r     :    Reset fit.
                  q     :    Quit Program.
             ---------------------------------------------------------------------------
            Written By:  Rongmon Bordoloi                                   July 13 2017.
    
    
            ----------------------------------------------------------------------------
            Note::The purpose of this code is to create a spline continuum fit from selected points. 
            The help scene activates by pressing h on the plot. 
            The program only works properly if none of the toolbar buttons in the figure is activated. 
            Basic code is taken from : http://www.ster.kuleuven.be/~pieterd/python/html/plotting/specnorm.html
            Heavily modified by Rongmon Bordoloi July 13/14 2017.
            Modified to add custom points and changed the look of the plots.
            Also added custom input options to read different formats. 
            Input file could be ascii, fits or pickle format
            Output will be in the same format as the input file. 
            Added help feature and graceful exit option. - Now pressing q will exit the program at any stage
            ---------------------------------------------------------------------------
        """
    def __init__(self,wave,flux,error=False,xlim=[-600.,600.],**kwargs):
        self.flux=flux
        self.error=error
        self.wave=wave
    
        self.fig, self.ax = plt.subplots()
        #This is where you feed in your velocity and flux to be fit
        spectrum, = self.ax.step(self.wave,self.flux,'b-',label='spectrum',linewidth=1)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')

        self.w = widgets.HTML()
        display(self.w)

        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        cid1 = self.fig.canvas.mpl_connect('key_press_event', self.onpress)
        cid2=self.fig.canvas.mpl_connect('pick_event',self.onpick)

        
        
    def onclick(self,event):
        # when none of the toolbar buttons is activated and the user clicks in the
        # plot somewhere, compute the median value of the spectrum in a 5angstrom
        # window around the x-coordinate of the clicked point. The y coordinate
        # of the clicked point is not important. Make sure the continuum points
        # `feel` it when it gets clicked, set the `feel-radius` (pickradius) to 5 points
        toolbar = plt.get_current_fig_manager().toolbar
        #self.w.value = 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
        self.w.value = f'button={event.button}, x={event.x}, y={event.y}, xdata={event.xdata}, ydata={event.ydata}'

        if event.button==1 and toolbar.mode=='':
            window = ((event.xdata-2.5)<=self.wave) & (self.wave<=(event.xdata+2.5))
            y = np.median(self.flux[window])
            self.ax.plot(event.xdata,y,'ro',ms=5,picker=5,label='cont_pnt',markeredgecolor='k')
            display(self.w)   
    def onpick(self,event):
        # when the user clicks right on a continuum point, remove it
        if event.mouseevent.button==3:
            if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
                event.artist.remove()
                display(self.w)
  
    def onpress(self,event):
        if event.key=='enter':
            cont_pnt_coord = []
            for artist in self.ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                    cont_pnt_coord.append(artist.get_data())
                elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                    artist.remove()
            cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
            sort_array = np.argsort(cont_pnt_coord[:,0])
            x,y = cont_pnt_coord[sort_array].T
            spline = splrep(x,y,k=3)
            continuum = splev(self.wave,spline)
            self.ax.plot(self.wave,continuum,'r-',lw=2,label='continuum')
            display(self.w)
    
        # when the user hits 'n' and a spline-continuum is fitted, normalise the
        # spectrum
        elif event.key=='n':
            continuum = None
            for artist in self.ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                    continuum = artist.get_data()[1]
                    break
            if continuum is not None:
                self.ax.cla()
                self.ax.step(self.wave,self.flux/continuum,'b-',label='normalised',linewidth=1)
                self.ax.step(self.wave,continuum,'r-',label='unnorm_cont',linewidth=0)            
                self.ax.plot([np.min(self.wave),np.max(self.wave)],[1,1],'k--')
                self.ax.set_xlim([np.min(self.wave),np.max(self.wave)])
                self.ax.xlabel('Wavelength')
                self.ax.ylabel('Relative Flux')
    
    
        # when the user hits 'r': clear the axes and plot the original spectrum
        elif event.key=='r':
            self.ax.cla()
            self.ax.step(self.wave,self.flux,'b-')
            display(self.w)
        # when the user hits 'b': selects a handpicked x,y value
        elif event.key=='b':
            self.ax.plot(event.xdata,event.ydata,'ro',ms=5,picker=5,label='cont_pnt',markeredgecolor='k')
            display(self.w)
        #If the user presses 'h': The help is printed on the screen
        elif event.key=='h':
            print(
            '''    
            ---------------------------------------------------------------------------
            This is an interactive continuum fitter for 1D spectrum.
            The purpose of this code is to create a spline continuum fit from selected points.
            The help scene activates by pressing h on the plot.
    
            The program only works properly if none of the toolbar buttons in the figure is activated. 
    
    
            Useful Keystrokes:
    
                Mouse Clicks:
                
                    Left Click  : Select the median flux value within +/- 5 pixel from the x-coordinate.
                                  These points are used for the continuum fit.
                    Right Click : Delete the nearest continuum point.
    
                Keystrokes:
                  
                  b     :    Select a point for continuum fit at that exact (x,y) coordinate.
                  enter :    Perform a spline fit to data to create a continuum.
                  n     :    Show the normalized spectrum.
                  w     :    Only after pressing n: This will ourput the continuum. 
                  h     :    This Help screen.
                  r     :    Reset fit.
                  q     :    Quit Program.
             ---------------------------------------------------------------------------
            Written By:  Rongmon Bordoloi                                   July 13 2017.
    
    
            ----------------------------------------------------------------------------
            
            Basic code is taken from : http://www.ster.kuleuven.be/~pieterd/python/html/plotting/specnorm.html
            Heavily modified by Rongmon Bordoloi July 13/14 2017.
            Modified to add custom points and changed the look of the plots.
            Also added custom input options to read different formats. 
            Input file could be ascii, fits or pickle format
            Output will be in the same format as the input file. 
            Added help feature and graceful exit option. - Now pressing q will exit the program at any stage
            ---------------------------------------------------------------------------
            '''
            )
    
    
        # At any time pressing q means graceful exit
        elif event.key=='q':
            self.w.value ='Good Bye!'
            plt.close()
                
    
        # when the user hits 'w': if the normalised spectrum exists, write it to a
        # file.
        elif event.key=='w':
            for artist in self.ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='unnorm_cont':#'normalised':
                    data = np.array(artist.get_data())
                    cont=(data.T[:,1])
                    self.cont=cont
                    self.w.value ='Final Continuum Chosen'
                    return self.cont
                    break
