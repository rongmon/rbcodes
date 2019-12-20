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
              w     :    Only after pressing n, write fitted continuum to file. 
                         [Output file saves wave, flux, error [if available], continuum]
                         [output file is in same format as input file with *_norm.* appended to the name.]
              h     :    This Help screen.
              r     :    Reset fit.
              q     :    Quit Program.
        ---------------------------------------------------------------------------
        Example:   Type in Terminal 

                        Case I:
                            > ipython rb_cont.py filename

                            Where filename could be file.fits, file.txt, file.dat, or file.p 
                            Fits and Pickle files should have dictionaries with keys = wave, flux and error [optional].
                            ascii files should have tab seperated columns with wave,flux and error [optional]
                            If it is none of these extentions then please specify as described below.


                        Case II:
                            > ipython rb_cont.py filename filetype

                            Where filename could be file.fits, file.txt, file.dat, or file.p 
                            Fits and Pickle files should have dictionaries with keys = wave, flux and error [optional].
                            ascii files should have tab seperated columns with wave,flux and error [optional]

                            filetype =  Type of file.
                                filetype could be written as : fits [for fits files.]
                                                             : ascii [for ascii files.]
                                                             : p for [pickle files.]

                        Case III: 

                        Add this path to your .cshrc file 
                            alias rb_cont   'ipython PATH_TO_THIS_FILE/rb_cont.py'

                            then this file can be run as 
                            > rb_cont filename
                                Or
                            > rb_cont filename filetype


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
        Added keyword xfits to read in xspecplot formatted files [RB 9.10.2017]
        Did minor syntax change to fix saving file issue [RB 06.17.2019]
        ---------------------------------------------------------------------------
'''
   
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep,splev
import pdb
import sys
import os

def onclick(event):
    # when none of the toolbar buttons is activated and the user clicks in the
    # plot somewhere, compute the median value of the spectrum in a 10angstrom
    # window around the x-coordinate of the clicked point. The y coordinate
    # of the clicked point is not important. Make sure the continuum points
    # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
    toolbar = plt.get_current_fig_manager().toolbar
    if event.button==1 and toolbar.mode=='':
        window = ((event.xdata-5)<=wave) & (wave<=(event.xdata+5))
        y = np.median(flux[window])
        plt.plot(event.xdata,y,'ro',ms=5,picker=5,label='cont_pnt',markeredgecolor='k')
    plt.draw()

def onpick(event):
    # when the user clicks right on a continuum point, remove it
    if event.mouseevent.button==3:
        if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
            event.artist.remove()

def ontype(event):
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
        continuum = splev(wave,spline)
        plt.plot(wave,continuum,'r-',lw=2,label='continuum')

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
            plt.step(wave,flux/continuum,'b-',label='normalised',linewidth=1)
            plt.step(wave,continuum,'r-',label='unnorm_cont',linewidth=0)            
            plt.plot([np.min(wave),np.max(wave)],[1,1],'k--')
            plt.xlim([np.min(wave),np.max(wave)])
            plt.title(filename)
            plt.xlabel('Wavelength')
            plt.ylabel('Relative Flux')


    # when the user hits 'r': clear the axes and plot the original spectrum
    elif event.key=='r':
        plt.cla()
        plt.step(wave,flux,'b-')

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
              w     :    Only after pressing n, write fitted continuum to file. 
                         [Output file saves wave, flux, error [if available], continuum]
                         [output file is in same format as input file with *_norm.* appended to the name.]
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
                outfilename=os.path.splitext(filename)[0]+'_norm.'+filetype
                if filetype=='ascii':
                    from astropy.table import Table, Column, MaskedColumn
                    if (len(tab)>=3):
                        table=Table([wave,flux,error,cont],names=['wave','flux','error','cont'])
                    else:
                        table=Table([wave,flux,cont],names=['wave','flux','cont'])
                    ascii.write(table,outfilename)
                if (filetype=='fits') | (filetype=='xfits'):
                    from astropy.table import Table, Column, MaskedColumn
                    if (len(tab)>=3):
                        table=Table([wave,flux,error,cont],names=['wave','flux','error','cont'])
                    else:
                        table=Table([wave,flux,cont],names=['wave','flux','cont'])
                    table.write(outfilename,format='fits')
                if filetype=='p':
                    table={}
                    table['wave']=wave
                    table['flux']=flux
                    table['cont']=cont
                    if (len(tab)>=3):
                        table['error']=error
                    pickle.dump(table, open( outfilename, "wb"), protocol=2 )
                print('Saved to file')
                break
    plt.draw()




if __name__ == "__main__":
    # Get the filename of the spectrum from the command line, and plot it
    filename = sys.argv[1]

    #Check if filetype is specified if not try to take what is given as extention to the file


    if (len(sys.argv) >2):
        filetype=sys.argv[2]
    else:
        #Take File Extention and try to match it
        tt=os.path.splitext(filename)[1]
        if (tt=='txt')| (tt=='dat'):
            filetype='ascii'
        else:
            filetype=tt[1:len(tt)]
    cwd=os.getcwd()
    print(cwd+'/'+filename)

    # Read in Files in differet formats
    if filetype=='ascii':
        from astropy.io import ascii
        dat=ascii.read(cwd+'/'+filename)
        tab=dat.keys()
        wave=np.array(dat[tab[0]])
        flux=np.array(dat[tab[1]])
        if (len(dat.keys())>=3):
            error=dat[tab[2]]
    elif filetype=='fits':
        from astropy.io import fits
        file=fits.open(cwd+'/'+filename)
        dat=file[1].data
        tab=dat.names
        wave=np.array(dat['wave'][0])
        flux=np.array(dat['flux'][0])
        if (len(tab)>=3):
            error=np.array(dat['error'][0])
    if filetype=='xfits':
        from astropy.io import fits
        hdu = fits.open(filename)
        wave = hdu['wavelength'].data
        flux = hdu['flux'].data
        error=hdu['error'].data
        hdu.close()
    elif filetype=='p':
        import pickle
        dat=pickle.load( open(cwd+'/'+filename, "rb" ))
        tab=dat.keys()
        wave=np.array(dat['wave'])
        flux=np.array(dat['flux'])
        if (len(tab)>=3):
            error=np.array(dat['error'])

    spectrum, = plt.step(wave,flux,'b-',label='spectrum',linewidth=1)
    plt.title(filename)
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')


    # Connect the different functions to the different events
    plt.gcf().canvas.mpl_connect('key_press_event',ontype)
    plt.gcf().canvas.mpl_connect('button_press_event',onclick)
    plt.gcf().canvas.mpl_connect('pick_event',onpick)
    plt.show() # show the window
