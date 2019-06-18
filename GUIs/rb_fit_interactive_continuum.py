import matplotlib
matplotlib.use('QT5Agg')

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep,splev
import sys
import os

class rb_fit_interactive_continuum(object):

    def __init__(self,wave,flux,error):

        self.wave=wave
        self.flux=flux
        self.error=error
        
        spectrum, = plt.step(self.wave,self.flux,'b-',label='spectrum',linewidth=1)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')




        # Connect the different functions to the different events
        plt.gcf().canvas.mpl_connect('key_press_event',self.ontype)
        plt.gcf().canvas.mpl_connect('button_press_event',self.onclick)
        plt.gcf().canvas.mpl_connect('pick_event',self.onpick)
        plt.show() # show the window


    def onclick(self,event):
        # when none of the toolbar buttons is activated and the user clicks in the
        # plot somewhere, compute the median value of the spectrum in a 5angstrom
        # window around the x-coordinate of the clicked point. The y coordinate
        # of the clicked point is not important. Make sure the continuum points
        # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
        toolbar = plt.get_current_fig_manager().toolbar
        if event.button==1 and toolbar.mode=='':
            window = ((event.xdata-2.5)<=self.wave) & (self.wave<=(event.xdata+2.5))
            y = np.median(self.flux[window])
            plt.plot(event.xdata,y,'ro',ms=5,picker=5,label='cont_pnt',markeredgecolor='k')
        plt.draw()
    def onpick(self,event):
        # when the user clicks right on a continuum point, remove it
        if event.mouseevent.button==3:
            if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
                event.artist.remove()

    def ontype(self,event):
        #---------------------------------------------------------------------------
        # When the user hits enter:
        # 1. Cycle through the artists in the current axes. If it is a continuum
        #    point, remember its coordinates. If it is the fitted continuum from the
        #    previous step, remove it
        # 2. sort the continuum-point-array according to the x-values
        # 3. fit a spline and evaluate it in the wavelength points
        # 4. plot the continuum
        #
        # Original Code taken from : http://www.ster.kuleuven.be/~pieterd/python/html/plotting/specnorm.html
        # Modified by Rongmon Bordoloi July 13 2017.
        # Modified to add custom points and changed the look of the plots.
        # Also added custom input options to read different formats. 
        # Input file could be ascii, fits or pickle format
        # Output will be in the same format as the input file. 
        # Added help feature and graceful exit option. - Now pressing q will exit the program at any stage
        #---------------------------------------------------------------------------
        if event.key=='enter':
            cont_pnt_coord = []
            for artist in plt.gca().get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                    cont_pnt_coord.append(artist.get_data())
                elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                    artist.remove()
            cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
            sort_array = np.argsort(cont_pnt_coord[:,0])
            x,y = cont_pnt_coord[sort_array].T
            spline = splrep(x,y,k=3)
            continuum = splev(self.wave,spline)
            plt.plot(self.wave,continuum,'r-',lw=2,label='continuum')
    
        # when the user hits 'n' and a spline-continuum is fitted, normalise the
        # spectrum
        elif event.key=='n':
            continuum = None
            for artist in plt.gca().get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                    continuum = artist.get_data()[1]
                    break
            if continuum is not None:
                plt.cla()
                plt.step(self.wave,self.flux/continuum,'b-',label='normalised',linewidth=1)
                plt.step(self.wave,continuum,'r-',label='unnorm_cont',linewidth=0)            
                plt.plot([np.min(self.wave),np.max(self.wave)],[1,1],'k--')
                plt.xlim([np.min(self.wave),np.max(self.wave)])
                plt.xlabel('Wavelength')
                plt.ylabel('Relative Flux')
    
    
        # when the user hits 'r': clear the axes and plot the original spectrum
        elif event.key=='r':
            plt.cla()
            plt.step(self.wave,self.flux,'b-')
    
        # when the user hits 'b': selects a handpicked x,y value
        elif event.key=='b':
            plt.plot(event.xdata,event.ydata,'ro',ms=5,picker=5,label='cont_pnt',markeredgecolor='k')
            plt.draw()
    
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
            quit_index=0;
            for artist in plt.gca().get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                    quit_index=1
                if quit_index==1:
                    plt.close()
                    print('Interactive Contunuum Normalization Done.')
                    print('Hope you remembered to save the fit by pressing w!')
                    print('Good Bye!')
                    break
                else:
                    plt.close()
                    print('Quitting without normalizing. Moving along.....')
                    break
    
    
    
    
        # when the user hits 'w': if the normalised spectrum exists, write it to a
        # file.
        elif event.key=='w':
            for artist in plt.gca().get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='unnorm_cont':#'normalised':
                    data = np.array(artist.get_data())
                    cont=(data.T[:,1])
                    self.cont=cont
                    print('Final Continuum Chosen')
                    return self.cont
                    break
        plt.draw()
    
