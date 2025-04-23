import matplotlib
matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep, splev
import sys
import os
import warnings

class rb_fit_interactive_continuum:
    """
    Interactive continuum fitter for 1D spectrum.

    This class provides an interactive matplotlib-based GUI for fitting
    continuum to spectral data. Users can select points by clicking on the
    plot, and fit a cubic spline through these points to create a continuum.

    Parameters
    ----------
    wave : array-like
        Wavelength values for the spectrum.
    flux : array-like
        Flux values for the spectrum.
    error : array-like
        Error values for the spectrum.

    Attributes
    ----------
    wave : ndarray
        Wavelength array.
    flux : ndarray
        Flux array.
    error : ndarray
        Error array.
    cont : ndarray or None
        The fitted continuum (set after fitting).
    norm_flag : int
        Flag indicating if normalization has been performed (0=no, 1=yes).

    Notes
    -----
    The interactive interface provides the following controls:

    Mouse Clicks:
        Left Click  : Select the median flux value within +/- 2.5 units from
                     the x-coordinate for continuum fitting.
        Right Click : Delete the nearest continuum point.

    Keystrokes:
        b     : Select a point for continuum fit at the exact (x,y) coordinate
                of the cursor.
        enter : Perform a spline fit to create a continuum.
        n     : Show the normalized spectrum.
        w     : After pressing 'n', this will save the continuum.
        h     : Display the help screen.
        r     : Reset fit.
        q     : Quit the interactive session.

    Examples
    --------
    >>> import numpy as np
    >>> wave = np.linspace(4000, 5000, 1000)
    >>> flux = np.ones_like(wave) + 0.1 * np.sin(wave/100)
    >>> error = 0.05 * np.ones_like(wave)
    >>> fitter = rb_fit_interactive_continuum(wave, flux, error)
    >>> # Interact with the plot window
    >>> # After pressing 'w', the continuum is accessible
    >>> continuum = fitter.cont
    """

    def __init__(self, wave, flux, error):
        """
        Initialize the interactive continuum fitter.

        Parameters
        ----------
        wave : array-like
            Wavelength values for the spectrum.
        flux : array-like
            Flux values for the spectrum.
        error : array-like
            Error values for the spectrum.

        Raises
        ------
        ValueError
            If input arrays have mismatched lengths or contain invalid values.
        """
        # Convert inputs to numpy arrays if they aren't already
        self.wave = np.asarray(wave)
        self.flux = np.asarray(flux)
        self.error = np.asarray(error)
        
        # Validate inputs
        if len(self.wave) != len(self.flux) or len(self.wave) != len(self.error):
            raise ValueError("Wavelength, flux, and error arrays must have the same length")
        
        if len(self.wave) == 0:
            raise ValueError("Input arrays cannot be empty")
            
        # Check for NaN or infinity
        if np.any(np.isnan(self.wave)) or np.any(np.isinf(self.wave)):
            warnings.warn("Wavelength array contains NaN or infinity values")
            
        if np.any(np.isnan(self.flux)) or np.any(np.isinf(self.flux)):
            warnings.warn("Flux array contains NaN or infinity values")
            
        if np.any(np.isnan(self.error)) or np.any(np.isinf(self.error)):
            warnings.warn("Error array contains NaN or infinity values")
        
        self.norm_flag = 0
        self.cont = None
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.step(self.wave, self.flux, 'b-', label='spectrum', linewidth=1)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title('Interactive Continuum Fitting')
        
        # Connect the different functions to the different events
        plt.gcf().canvas.mpl_connect('key_press_event', self.ontype)
        plt.gcf().canvas.mpl_connect('button_press_event', self.onclick)
        plt.gcf().canvas.mpl_connect('pick_event', self.onpick)
        plt.show()  # show the window

    def onclick(self, event):
        """
        Handle mouse click events for selecting continuum points.
        
        When a user left-clicks on the plot (and no toolbar buttons are active),
        this function computes the median flux value in a 5 unit window around
        the clicked x-coordinate and adds it as a continuum point.
        
        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            The mouse click event.
        """
        # Only proceed if the event has valid data coordinates
        if event.xdata is None or event.ydata is None:
            return
            
        toolbar = plt.get_current_fig_manager().toolbar
        if event.button == 1 and toolbar.mode == '':
            try:
                # Define a window around the clicked x-coordinate
                window = ((event.xdata - 2.5) <= self.wave) & (self.wave <= (event.xdata + 2.5))
                
                # Check if there are points in the window
                if not np.any(window):
                    warnings.warn("No data points found in the selected window")
                    return
                    
                # Compute the median flux in the window
                y = np.median(self.flux[window])
                
                # Plot the continuum point
                plt.plot(event.xdata, y, 'ro', ms=5, pickradius=5, 
                         label='cont_pnt', markeredgecolor='k', picker=True)
                plt.draw()
            except Exception as e:
                warnings.warn(f"Error selecting point: {str(e)}")

    def onpick(self, event):
        """
        Handle pick events for removing continuum points.
        
        When a user right-clicks on a continuum point, this function removes it.
        
        Parameters
        ----------
        event : matplotlib.backend_bases.PickEvent
            The pick event.
        """
        try:
            if event.mouseevent.button == 3:
                if hasattr(event.artist, 'get_label') and event.artist.get_label() == 'cont_pnt':
                    event.artist.remove()
                    plt.draw()
        except Exception as e:
            warnings.warn(f"Error removing point: {str(e)}")

    def ontype(self, event):
        """
        Handle keyboard events for controlling the fitting process.
        
        This function processes the following key commands:
        - 'enter': Fit a spline to the selected continuum points
        - 'n': Show the normalized spectrum
        - 'r': Reset the fit
        - 'b': Select a point at the exact cursor position
        - 'h': Display the help screen
        - 'q': Quit the interactive session
        - 'w': Save the continuum after normalization
        
        Parameters
        ----------
        event : matplotlib.backend_bases.KeyEvent
            The keyboard event.
        
        Returns
        -------
        ndarray or None
            The fitted continuum when 'w' is pressed after normalization,
            otherwise None.
        """
        if event.key == 'enter':
            try:
                # Collect continuum point coordinates
                cont_pnt_coord = []
                for artist in plt.gca().get_children():
                    if hasattr(artist, 'get_label') and artist.get_label() == 'cont_pnt':
                        cont_pnt_coord.append(artist.get_data())
                    elif hasattr(artist, 'get_label') and artist.get_label() == 'continuum':
                        artist.remove()
                
                # Check if there are enough points for a spline fit
                if len(cont_pnt_coord) < 4:
                    warnings.warn("At least 4 points are required for cubic spline fitting")
                    return
                
                # Process the coordinates
                cont_pnt_coord = np.array(cont_pnt_coord)[..., 0]
                sort_array = np.argsort(cont_pnt_coord[:, 0])
                x, y = cont_pnt_coord[sort_array].T
                
                # Fit a cubic spline
                spline = splrep(x, y, k=3)
                continuum = splev(self.wave, spline)
                plt.plot(self.wave, continuum, 'r-', lw=2, label='continuum')
                plt.draw()
            except Exception as e:
                warnings.warn(f"Error fitting continuum: {str(e)}")
                
        elif event.key == 'n':
            try:
                continuum = None
                for artist in plt.gca().get_children():
                    if hasattr(artist, 'get_label') and artist.get_label() == 'continuum':
                        continuum = artist.get_data()[1]
                        break
                
                if continuum is not None:
                    plt.cla()
                    plt.step(self.wave, self.flux / continuum, 'b-', label='normalised', linewidth=1)
                    plt.step(self.wave, continuum, 'r-', label='unnorm_cont', linewidth=0)
                    plt.plot([np.min(self.wave), np.max(self.wave)], [1, 1], 'k--')
                    plt.xlim([np.min(self.wave), np.max(self.wave)])
                    plt.xlabel('Wavelength')
                    plt.ylabel('Relative Flux')
                    plt.draw()
                else:
                    print("No continuum fitted. Press 'enter' to fit the continuum first.")
            except Exception as e:
                warnings.warn(f"Error normalizing spectrum: {str(e)}")
                
        elif event.key == 'r':
            try:
                plt.cla()
                plt.step(self.wave, self.flux, 'b-')
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.draw()
            except Exception as e:
                warnings.warn(f"Error resetting plot: {str(e)}")
                
        elif event.key == 'b':
            try:
                if event.xdata is not None and event.ydata is not None:
                    plt.plot(event.xdata, event.ydata, 'ro', ms=5, pickradius=5, 
                             label='cont_pnt', markeredgecolor='k', picker=True)
                    plt.draw()
            except Exception as e:
                warnings.warn(f"Error placing point: {str(e)}")
                
        elif event.key == 'h':
            print("""
            ---------------------------------------------------------------------------
            This is an interactive continuum fitter for 1D spectrum.
            The purpose of this code is to create a spline continuum fit from selected points.
            
            The program only works properly if none of the toolbar buttons in the figure is activated.
            
            Useful Keystrokes:
            
                Mouse Clicks:
                
                    Left Click  : Select the median flux value within +/- 2.5 units from
                                 the x-coordinate for continuum fitting.
                    Right Click : Delete the nearest continuum point.
                
                Keystrokes:
                  
                  b     : Select a point for continuum fit at the exact (x,y) coordinate
                         of the cursor.
                  enter : Perform a spline fit to create a continuum.
                  n     : Show the normalized spectrum.
                  w     : After pressing 'n', this will save the continuum.
                  h     : Display this help screen.
                  r     : Reset fit.
                  q     : Quit the interactive session.
            ---------------------------------------------------------------------------
            """)
                
        elif event.key == 'q':
            try:
                if self.norm_flag == 1:
                    plt.close()
                    print('Interactive Continuum Normalization Done.')
                    print('Hope you remembered to save the fit by pressing w!')
                    print('Good Bye!')
                else:
                    plt.close()
                    print('Quitting without normalizing. Moving along.....')
            except Exception as e:
                warnings.warn(f"Error quitting: {str(e)}")
                plt.close()
                
        elif event.key == 'w':
            try:
                for artist in plt.gca().get_children():
                    if hasattr(artist, 'get_label') and artist.get_label() == 'unnorm_cont':
                        data = np.array(artist.get_data())
                        cont = (data.T[:, 1])
                        self.cont = cont
                        self.norm_flag = 1
                        print('Final Continuum Chosen')
                        return self.cont
                
                if not self.norm_flag:
                    print("No normalized spectrum available. Press 'n' after fitting the continuum.")
            except Exception as e:
                warnings.warn(f"Error saving continuum: {str(e)}")
                
        plt.draw()
        
    def get_continuum(self):
        """
        Get the fitted continuum.
        
        Returns
        -------
        ndarray or None
            The fitted continuum if available, otherwise None.
        """
        return self.cont