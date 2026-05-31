# Project Documentation
[Back to Main Page](../main_readme.md)

## GUI Tools — Quick Reference

| Tool | Description | Docs |
|------|-------------|------|
| **launch_specgui** | Interactive `rb_spec` GUI — single spectrum and batch modes | [rb_spec docs](rb_spec/rb_spec.md) |
| **rb_multispec** | Multi-spectrum viewer with line ID, redshift overlay, and fitting | [multispec docs](multispec/multispec.md) |
| **rb_ifuview** | IFU datacube viewer — KCWI/MUSE/generic, spectra, moment maps, ds9 | [ifuview docs](ifuview/rb_ifuview.md) |
| **rb_zgui** | Redshift measurement GUI for JWST NIRCam/grism spectroscopy | [PDF tutorial](zgui/Tutorial_for_Emission_Line_Redshift_Estimator_GUI.pdf) |
| **rb_llsfitter** | GUI for Lyman Limit System column density fitting | [LLSFitter docs](LLSFitter/LLSFitter.md) |
| **interactive_continuum_fit** | Polynomial and spline continuum fitter with interactive masking | [docs](interactive_continuum_fit.md) |
| **rb_align** | Astrometry alignment — IFU cubes and 2D images to a reference frame | [rb_align docs](rb_align/rb_align.md) |

---

## rb_align

[Full Documentation](rb_align/rb_align.md)

`rb_align` aligns IFU datacubes (KCWI, MUSE, NIRSpec IFU, JWST/NIRISS) and 2D
images to a reference frame using interactive or automated source matching with
full WCS fitting via `astropy`. Works standalone, in batch pipelines, and
embedded in GUIs.

```python
import rbcodes.rb_align as rb_align
rb_align.help()                         # print all workflow examples
```

```python
from rbcodes.rb_align import wcs_align

c = wcs_align.from_file(reference='hst_ref.fits', targets=['cube.fits'],
                         input_type='ifu')
c.preprocess()
c.find_sources(strategy='interactive', stretch='zscale', save_catalog='sources.fits')
c.align()
c.qa(plot=True)
c.write_output()                        # → cube_wcsfix.fits
```

**Source strategies:** `interactive` | `gaia` | `dao` | `knots` | `cross_corr` | `batch` | `auto`

**Supported modes:** image↔image, IFU↔image, IFU↔IFU, batch (10+ exposures)

---

*Auto-generated documentation from docstrings*

## Modules

### interactive_continuum_fit

Interactive masking and continuum fitting tool for spectral analysis.

This module provides a PyQt5-based GUI for interactive selection of mask regions,
continuum fitting with both polynomial and spline methods.

### rb_interactive_mask

Interactive masking and continuum fitting tool for spectral analysis.

This module provides a PyQt5-based GUI for interactive selection of mask regions
and continuum fitting for spectroscopic data.

### launch_specgui

Unified launcher for rb_spec GUI applications.

This script replaces both the old launch_specgui.py and launch_batch.py files.
It can launch either the single-spectrum GUI or the batch processing GUI
depending on the command-line arguments provided.

Usage:
    python launch_specgui.py [filename] [-t filetype]    # Single spectrum mode (default)
    python launch_specgui.py -b [config_file]            # Batch processing mode
    python launch_specgui.py -v                          # Show version
    python launch_specgui.py --help                      # Show help

Examples:
    python launch_specgui.py spectrum.fits                    # Load single spectrum
    python launch_specgui.py spectrum.fits -t linetools       # Load with specific filetype
    python launch_specgui.py analysis.json                    # Load saved analysis
    python launch_specgui.py -b                               # Launch batch processor
    python launch_specgui.py -b batch_config.json             # Load batch configuration
    python launch_specgui.py -b batch_items.csv               # Import batch items from CSV

### interactive_cont_jn

Interactive continuum fitter for jupyter notebook

### PlotSpec_Integrated

Modules for PlotSpec_Integrated GUI

### rb_cont

Interactive continuum fitter for 1D spectra.

This script loads a spectrum from a file and launches an interactive
continuum fitting interface. The user can select points by clicking on
the plot, and fit a cubic spline through these points to create a continuum.
The normalized spectrum and fitted continuum can be saved to a file.

Example usage:
    python rb_cont.py filename [filetype]

Where:
    filename : Path to the file containing the spectrum
    filetype : (Optional) Type of file - 'fits', 'ascii', 'p' (pickle), or 'xfits'
               If not provided, it will be inferred from the file extension.

Controls:
    Mouse:
        Left Click  : Select the median flux value within +/- 2.5 units from 
                      the x-coordinate for continuum fitting.
        Right Click : Delete the nearest continuum point.
    
    Keyboard:
        b     : Select a point for continuum fit at the exact (x,y) coordinate 
                of the cursor.
        enter : Perform a spline fit to create a continuum.
        n     : Show the normalized spectrum.
        w     : After pressing 'n', this will save the continuum.
        h     : Display the help screen.
        r     : Reset fit.
        q     : Quit the interactive session.

### specgui.batch.panels.ew_range_editor

Separate module for EW range editing functionality.
This module provides an enhanced EW range editor with interactive features.

### specgui.batch.panels.batch_figure_generator

Batch figure generation module for rb_spec batch processing.

This module provides functionality to create publication-quality figures
from batch processing results, supporting both individual figures and
multi-page PDF layouts.

### specgui.batch.panels.spectrum_plotter

Spectrum plotting module for rb_spec GUI applications.

This module provides reusable plotting functions that generate matplotlib figures
from rb_spec objects. All plots are generated fresh from the rb_spec object state
to ensure consistency and avoid update artifacts.

### specgui.batch.panels.batch_selection_dialog

Batch Selection Dialog for choosing multiple systems for processing.

This dialog provides a familiar scrollable checklist interface with smart
selection buttons for choosing systems to process in batch mode.

### abstools.intervening_utils

Utilities for handling intervening absorption lines.
This module provides functions for reading and plotting
intervening absorption line data.

### abstools.Metal_Plot

Metal_Plot.py - Main module for the absorption line analysis toolbox.
Updated with improved signal-slot communication and error handling for greater stability.

### abstools.event_handler

Event handler module for the absorption line analysis toolbox.
This module handles all user interactions including mouse clicks,
key presses, and mouse motion events.
Improved with signal-slot communication and better error handling.

### abstools.config

Configuration module for the absorption line analysis toolbox.
This module centralizes constants, settings, and utility functions
used throughout the application.

### abstools.text_utils

Text utilities module for the absorption line analysis toolbox.
This module handles text formatting and display for measurements.

### abstools.equivalent_width

Equivalent Width module for the absorption line analysis toolbox.
This module handles the calculation of equivalent widths, column densities,
and other line measurements.
Enhanced with improved error handling and signal-slot communication.

### abstools.plotting

Plotting module for the absorption line analysis toolbox.
This module handles the plotting functionality for both continuum fitting
and normalized spectra.
Enhanced with improved error handling and signal-slot communication.

### abstools.abstools_launcher

AbsTools Launcher GUI

A PyQt5 GUI for setting up and launching the AbsTools absorption line analysis toolbox.
This launcher simplifies the process of loading spectra, setting redshift and lines to analyze,
and launching the Metal_Plot visualization tool.

Usage:
    python abstools_launcher.py

### abstools.simplified_abstools_launcher

Simplified AbsTools Launcher - Entry point for the Absorption Line Analysis Toolbox.

This script provides a convenient entry point for launching the AbsTools
wrapper GUI, with improved process isolation.

Usage:
    python simplified_abstools_launcher.py

### abstools.cleanup

Resource cleanup module for the absorption line analysis toolbox.
This module provides functions for safely cleaning up resources and preventing
segmentation faults when closing the application.

### abstools.Absorber

Absorber 

Inputs:
flux; wave; error; linelist; redshift; bin

1st. Asborber will bin the flux,wave and error to clean the data
2nd. will pull the actual lamd_rest from the atom.dat file for all lines
-->this (2nd) will return a dictionary of lamd_rest,ion_name,fval,gamma

3rd. Initialized entries for the Vstack plotting

### abstools.json_utils

JSON saving and loading functionality for the absorption line analysis toolbox.
This module provides functions to save and load analysis data in JSON format.

### abstools.ui_components

UI components module for the absorption line analysis toolbox.
This module handles the creation and management of UI elements and pages.
Now with JSON support added.
The main modification is to the SavePage class to add JSON save functionality.

## Classes

### MplCanvas (`interactive_continuum_fit`)

Matplotlib canvas for embedding in PyQt5.

### InteractiveContinuumFitWindow (`interactive_continuum_fit`)

Main window for interactive masking and continuum fitting.

### HelpDialog (`interactive_continuum_fit`)

Modal dialog to display HTML-formatted help content.

### ResultsDialog (`LLSFitter_GUI`)

Dialog to display detailed fit results

### MplCanvas (`rb_interactive_mask`)

Matplotlib canvas for embedding in PyQt5.

### InteractiveMaskWindow (`rb_interactive_mask`)

Main window for interactive masking and continuum fitting.

### interactive_cont (`interactive_cont_jn`)

This is an interactive continuum fitter for 1D spectrum.
    
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

### rb_fit_interactive_continuum (`rb_fit_interactive_continuum`)

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

### rb_spec (`rb_spec`)

A spectrum read into a class, spectrum will have following properties.

    Attributes 
    ----------
        wave: wavelength.
        flux: flux.
        error: error
        filename=filename and location
        filetype = False [default] : other options 
                 ascii, fits, HSLA, xfits, p [pickle], temp, and linetools [uses linetools.io routine for this]

        Optional: 
            All only valid for filetype=linetools option
            efil= errorfile [Default None]

    Returns
    -------
        This gives a rb_spec object with following attributes:

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
        self.N=  AOD column density
        self.N_e= AOD column density uncertainty
        self.logN=  log AOD column density
        self.logN_e= log AOD column density uncertainty
        self.Tau= Apparant optical depth as a function of velocity
        self.vel_centroid= EW weighted velocity centroid of the absorption line
        self.vel_disp=    1sigma velocity dispersion
        self.vel50_err = error on velocity centroid



        Written : Rongmon Bordoloi      April 2018
        Edit    : Rongmon Bordoloi      September 2018 Changed kwargs to be compatible to python 3   
        Edit    : Rongmon Bordoloi      Aug 2020: added linetools.io.readspec file
        Edit    : Rongmon Bordoloi      April 2021: added velocity centroid estimates
        Edit    : Rongmon Bordoloi      March 2022: Added more continuum fitting methods
        Edit    : Rongmon Bordoloi      April 2022: Added velocity centroid error
        Edit    : Rongmon Bordoloi      April 2022: Added different calling sequence to ingest numpy arrays directly. 
        Edit    : Rongmon Bordoloi      April 2022: Small updates to have all continuum fitting routines working
        Edit    : Rongmon Bordoloi      April 2022: Added plotting sliced spectrum option
        Edit    : Rongmon Bordoloi      September 2024: Added saving as json option, and loading the json dictionary as a rb_spec object.
        
        # WARNING: CALLING SEQUENCE HAS CHANGED SINCE APRIL 2022.
        CAREFULLY LOOK AT THE EXAMPLE BELOW

    Example
    -------
        import numpy as np
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt
        from rbcodes.GUIs.rb_spec import rb_spec as r 
        #List of absorber redshifts
        zabs=[0.511020,1.026311,1.564481]
        transition= 2796.3
        #Which absorber to analyze
        index=0
        filename='Quasar_Spectrum.fits'
        #Read in file
        s=r.from_file(filename,filetype='linetools')

        #-------------------------------------------------------------------------------
        # DETOUR --->
        #ALTERNATIVE 
        #IF YOU WANT TO DIRECTLY INJEST NUMPY ARRAYS DO THE FOLLOWING
        s=r.from_data(wave,flux,error)
        # HERE wave,flux,error are numpy arrays of wavelength,flux and error respectively.
        #-------------------------------------------------------------------------------
                
        #Shift spectra to rest frame
        s.shift_spec(zabs[index]);
        #Velocity window around transition
        xlim=[-1500,1500]
        
        #Slice Spectrum within that window
        s.slice_spec(transition,xlim[0],xlim[1],use_vel=True);
        
        #Fit continuum Mask the regions defined by velocity/ use_weight= True/False toggles if error is used for weighting the fit
        s.fit_continuum(mask=[-200,300,500,1100],domain=xlim,Legendre=3, use_weights=False)
        
        #-------------------------------------------------------------------------------
        # DETOUR 1--->
        #Alternative Fit continuum methods.
        #s.fit_continuum_ransac(window=149,mednorm=False)  
        
        #-------------------------------------------------------------------------------
        # DETOUR 2--->
        #Aternate continuum fitting method [interactive]
        s.fit_continuum(Interactive=True)
        #Aternate continuum fitting method [input prefit continuum]
        # Length of prefit continuum array = length of sliced spectrum
        s.fit_continuum(Legendre=False,prefit_cont=cont_arrary)
        #-------------------------------------------------------------------------------
        # DETOUR 3--->
        #Aternate continuum fitting method, using sigma clipping
        s.fit_continuum(domain=xlim,Legendre=3,sigma_clip=True)
        #-------------------------------------------------------------------------------
        
       
        
        #Compute EW
        #Compute equivalent width within a velocity window
        s.compute_EW(transition,vmin=-200.,vmax=360.);
        
        # Saving the analysis
        #--------------------
        # There are two options:
        # First method: save everything as a pickle file [default]
        s.save_slice('outfile.p')




        #---------------------
        # Second method: Saving information as a json file
        s.save_slice('outfile.json',file_format='json')
        


        #-----------------------------------------------------------
        #Loading the saved rb_spec object from the above two methods

        # Loading the rb_spec object back from the saved pickle file
        import pickle
        with open('outfile.p', 'rb') as f:
            # Load the pickled data
            sp_test = pickle.load(f)



        ----------------------------------------------------------
        # Loading the rb_spec object back from the saved json file
        from rbcodes.GUIs.rb_spec import load_rb_spec_object as r_load
        
        f='outfile.json'

        sp_test=r_load(f)






        #-------------------------------------------------------------------------------
        # Additional inspection routines        
        #plot the Full spectrum
        s.plot_spec()
        
        #Plot the sliced spectrum with the fitted continuum
        s.plot_slice()
    
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

### ToggleFrames (`zgui.gui_frame_io`)

Main class to search targeted frames

### GuessTransition (`zgui.guess_transition`)

Pop-up Widget to select ions

### FitsObj (`zgui.utils`)

Main object to save and traverse data between widgets

### Fits_2dAux (`zgui.utils`)

Auxiliary object to save and traverse non-important data

### MainWindow (`zgui.main`)

Simplified main window with proper positioning

### MessageBox (`multispecviewer.MessageBox`)

A widget that displays messages to the user with support for colored text.
    Designed to match the style of PlotSpec_Integrated.

### RedshiftInputWidget (`multispecviewer.RedshiftInputWidget`)

A standalone PyQt5 widget that provides UI elements for entering a redshift value
    and selecting a line list and color.
    
    Signals:
        submitted(float, str, str): Emitted when valid data is submitted
        linelist_changed(float, str, str): Emitted when linelist or color changes
        catalog_clicked(float, str, str): Emitted when catalog button is clicked

### IOManager (`multispecviewer.io_manager`)

Singleton class to handle all file I/O operations for MultispecViewer.
    Provides methods for loading and saving various file formats including
    JSON for combined data, and traditional text/CSV formats for backward compatibility.

### LineSelectionDialog (`multispecviewer.LineSelectionDialog`)

A non-modal dialog for selecting spectral lines from a list.
    Uses signals to communicate back with the parent window instead of blocking execution.

### SpectralPlot (`multispecviewer.multispec`)

A Matplotlib canvas for displaying spectral plots.
    Supports interactive key commands for adjusting axis limits and quitting the application.

### MainWindow (`multispecviewer.multispec`)

Main application window for displaying spectral plots.
    Includes a file selection button and a Matplotlib canvas.

### AbsorberManager (`multispecviewer.AbsorberManager`)

A widget for managing multiple absorber systems with their redshifts, line lists, and colors.
    Allows adding, removing, plotting, and hiding absorbers with a single checkbox.
    
    Note: To fully support this widget, the parent should implement:
    - plot_absorber_lines(row, z_abs, line_list, color) - Plot lines for this absorber
    - remove_absorber_lines(row) - Remove lines for this absorber
    - toggle_absorber_visibility(row, is_visible) - Toggle visibility of absorber lines
    - update_absorber_redshift(row, z_abs) - Update when redshift is changed manually

### RedshiftInputWidget (`multispecviewer.classic.RedshiftInputWidget_classic`)

A standalone PyQt5 widget that provides UI elements for entering a redshift value
    and selecting a line list (LLS or DLA).
    
    Signals:
        submitted(float, str): Emitted when valid data is submitted, containing the redshift value and selected line list
        linelist_changed(str): Emitted when linelist selection changes

### SpectralPlot (`multispecviewer.classic.multispec_classic`)

A Matplotlib canvas for displaying spectral plots.
    Supports interactive key commands for adjusting axis limits and quitting the application.

### MainWindow (`multispecviewer.classic.multispec_classic`)

Main application window for displaying spectral plots.
    Includes a file selection button and a Matplotlib canvas.

### LineSelectionDialog (`multispecviewer.classic.LineSelectionDialog_classic`)

A non-modal dialog for selecting spectral lines from a list.
    Uses signals to communicate back with the parent window instead of blocking execution.

### SpectrumController (`specgui.spectrum_controller`)

Controller class that manages the rb_spec instance and operations.

### RbSpecGUI (`specgui.main`)

Main application window for the rb_spec GUI.

### OutputPanel (`specgui.panels.output_panel`)

Panel for saving analysis results and generating output.

### AdvancedBatchDialog (`specgui.panels.advanced_batch_dialog`)

Advanced batch processing dialog with detailed visualization and controls.

### MatplotlibCanvas (`specgui.panels.measurement_panel`)

Canvas for matplotlib plots.

### MeasurementPanel (`specgui.panels.measurement_panel`)

Panel for computing equivalent width and column density.

### BatchPanel (`specgui.panels.batch_panel`)

Panel for batch processing multiple transitions or files.

### MatplotlibCanvas (`specgui.panels.continuum_panel`)

Simple matplotlib canvas for the spectrum preview.

### ContinuumPanel (`specgui.panels.continuum_panel`)

Panel that launches the existing interactive continuum fitter.

### MatplotlibCanvas (`specgui.panels.redshift_panel`)

Canvas for matplotlib plots.

### RedshiftPanel (`specgui.panels.redshift_panel`)

Panel for setting and applying the absorber redshift.

### InputPanel (`specgui.panels.input_panel`)

Panel for loading spectrum data from file or arrays.

### MatplotlibCanvas (`specgui.panels.transition_panel`)

Canvas for matplotlib plots.

### TransitionPanel (`specgui.panels.transition_panel`)

Panel for selecting a transition and slicing the spectrum.

### MasterBatchTable (`specgui.batch.master_batch_table`)

Master table using pandas DataFrame as single source of truth.

### BatchSpecGUI (`specgui.batch.batch_main`)

Main window for batch processing of absorption line spectra.

### BatchController (`specgui.batch.batch_controller`)

Controller for managing batch processing using master batch table.

### MatplotlibCanvas (`specgui.batch.panels.ew_range_editor`)

Canvas for matplotlib plots in the EW editor.

### BatchProcessor (`specgui.batch.panels.processing_panel`)

Worker object that processes batch items in a separate thread.

### ProcessingPanel (`specgui.batch.panels.processing_panel`)

Panel for configuring and running batch processing.

### ConfigurationPanel (`specgui.batch.panels.configuration_panel`)

Panel for configuring batch processing items using master table.

### ExportPanel (`specgui.batch.panels.export_panel`)

Panel for exporting batch processing results.

### BatchSelectionDialog (`specgui.batch.panels.batch_selection_dialog`)

Dialog for selecting multiple systems for batch processing.

### MatplotlibCanvas (`specgui.batch.panels.review_panel`)

Canvas for matplotlib plots.

### AxisLimitsDialog (`specgui.batch.panels.review_panel`)

Dialog for setting custom X/Y axis limits and analysis settings.

### ReviewPanel (`specgui.batch.panels.review_panel`)

Panel for reviewing batch processing results - with navigation and filtering.

### MainWindowSignals (`abstools.Metal_Plot`)

Signal definitions for the MainWindow class to enable robust signal-slot communication.

### AbsToolsCanvas (`abstools.Metal_Plot`)

Custom Matplotlib canvas for AbsTools with improved error handling.

### MainWindow (`abstools.Metal_Plot`)

Main application window for the absorption line analysis toolbox.
    Improved with signal-slot communication for better stability.

### Transitions (`abstools.Metal_Plot`)

Callable class to initialize and run the absorption line analysis application.
    Improved with better error handling and application lifecycle management.

### EventHandler (`abstools.event_handler`)

Class for handling all user interaction events in the application.

### EquivalentWidth (`abstools.equivalent_width`)

Class for calculating equivalent widths and related properties
    for absorption lines.

### Plotting (`abstools.plotting`)

Class for handling the plotting functionality of absorption line data.

### AbsToolsLauncher (`abstools.abstools_launcher`)

Main window for the AbsTools launcher GUI.

### ResourceCleanup (`abstools.cleanup`)

Class for safely cleaning up resources to prevent segmentation faults.

### NumPyJSONEncoder (`abstools.json_utils`)

Custom JSON encoder that handles NumPy arrays and other special types.

### HelpWindow (`abstools.ui_components`)

Window displaying help information for the application.

### SavePage (`abstools.ui_components`)

Dialog for saving analysis results in various formats.

### PageManager (`abstools.ui_components`)

Manager for handling page (tab) creation and initialization in the UI.

## Functions

### setup_ui() (`interactive_continuum_fit`)

Set up the user interface.

### show_help_dialog() (`interactive_continuum_fit`)

Open the help dialog.

### add_help_button_to_toolbar() (`interactive_continuum_fit`)

Add a help button to the matplotlib toolbar.

### toggle_method() (`interactive_continuum_fit`)

Toggle between polynomial and spline fitting methods.

### toggle_x_axis() (`interactive_continuum_fit`)

Toggle between velocity and wavelength x-axis.

### update_window_type() (`interactive_continuum_fit`)

Update the window type for median calculation.

### toggle_optimization() (`interactive_continuum_fit`)

Toggle optimization options.

### update_plots() (`interactive_continuum_fit`)

Update the plots with current data and masks.

### plot_masks() (`interactive_continuum_fit`)

Plot mask regions on the given axis.

### on_click() (`interactive_continuum_fit`)

Handle mouse click events.

### add_exact_spline_point() (`interactive_continuum_fit`)

Add a spline anchor point at the exact coordinates.

### add_median_spline_point() (`interactive_continuum_fit`)

Add a spline anchor point using the median flux in a window.

### remove_closest_spline_point() (`interactive_continuum_fit`)

Remove the spline point closest to the given coordinates.

### clear_spline_points() (`interactive_continuum_fit`)

Clear all spline points.

### handle_add_mask() (`interactive_continuum_fit`)

Handle adding mask regions.

### handle_remove_mask() (`interactive_continuum_fit`)

Handle removing mask regions with right clicks.

### merge_masks() (`interactive_continuum_fit`)

Merge overlapping mask regions.

### on_key_press() (`interactive_continuum_fit`)

Handle keyboard events.

### zoom_in() (`interactive_continuum_fit`)

Zoom in on the plot.

### zoom_out() (`interactive_continuum_fit`)

Zoom out on the plot.

### reset_zoom() (`interactive_continuum_fit`)

Reset the plot zoom to the domain limits.

### reset_masks() (`interactive_continuum_fit`)

Reset all masks.

### reset_all() (`interactive_continuum_fit`)

Reset everything: masks, spline points, and fitted continuum.

### undo_last_mask() (`interactive_continuum_fit`)

Remove the last mask pair.

### auto_mask() (`interactive_continuum_fit`)

Automatically generate masks using robust detection of spectral features,
        applicable in either wavelength or velocity space.
    
        Features:
        - Sigma clipping using MAD.
        - Error threshold masking.
        - Continuum fitting with fallback.
        - Minimum mask width and gap merging.
        - Applies in selected coordinate system: 'velocity' or 'wavelength'.

### manual_mask_entry() (`interactive_continuum_fit`)

Open a dialog for manual entry of mask values.

### fit_continuum() (`interactive_continuum_fit`)

Fit continuum using current parameters and masks.

### fit_polynomial_continuum() (`interactive_continuum_fit`)

Fit continuum using polynomial method with masked regions excluded.

### fit_spline_continuum() (`interactive_continuum_fit`)

Fit continuum using spline method with provided anchor points.

### check_unmasked_points() (`interactive_continuum_fit`)

Check if there are enough unmasked points for fitting.

### accept_results() (`interactive_continuum_fit`)

Accept results and close window.

### cancel() (`interactive_continuum_fit`)

Cancel and close without saving.

### closeEvent() (`interactive_continuum_fit`)

Handle window close event.

### launch_interactive_continuum_fit() (`interactive_continuum_fit`)

Launch the interactive continuum fitting GUI.
    
    Parameters
    ----------
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    velocity : array, optional
        Velocity array. If not provided, wavelength will be used as x-axis.
    existing_masks : list, optional
        Existing mask regions
    order : int, optional
        Initial polynomial order
    use_weights : bool, optional
        Whether to use weights
    domain : list, optional
        Velocity or wavelength domain limits
        
    Returns
    -------
    dict
        Results including masks, continuum, and fit parameters,
        or {'cancelled': True} if the user cancelled.

### launch_interactive_continuum_fit_dialog() (`interactive_continuum_fit`)

Launch the interactive continuum fitter as a modal dialog.
    
    This version uses QDialog.exec_() which works properly when called from
    an application that already has an event loop running.
    
    Parameters
    ----------
    Same parameters as launch_interactive_continuum_fit
    
    Returns
    -------
    dict
        Results including masks, continuum, and fit parameters,
        or {'cancelled': True} if the user cancelled.

### load_help_content() (`interactive_continuum_fit`)

Load the HTML help content from file.

### load_fallback_help() (`interactive_continuum_fit`)

Load a basic fallback help content if the file can't be found.

### set_default_regions() (`LLSFitter_GUI`)

Set the table to contain the default regions

### get_regions() (`LLSFitter_GUI`)

Get the current regions as a list of tuples

### validate_cell() (`LLSFitter_GUI`)

Validate that cell contents are numeric

### sort_regions() (`LLSFitter_GUI`)

Sort the continuum regions by minimum wavelength

### add_row() (`LLSFitter_GUI`)

Add a new empty row to the table and sort

### remove_selected_rows() (`LLSFitter_GUI`)

Remove the selected rows

### copy_to_clipboard() (`LLSFitter_GUI`)

Copy results to clipboard

### save_continuum_regions() (`LLSFitter_GUI`)

Save current continuum regions to a template file

### load_continuum_regions() (`LLSFitter_GUI`)

Load continuum regions from a template file

### suggest_continuum_regions() (`LLSFitter_GUI`)

Suggest continuum regions based on the spectrum

### create_menus() (`LLSFitter_GUI`)

Create application menus

### save_mcmc_posterior() (`LLSFitter_GUI`)

Save the MCMC posterior samples to a file

### save_results() (`LLSFitter_GUI`)

Save fit results to a JSON file

### load_results() (`LLSFitter_GUI`)

Load and display saved results from a JSON file

### _show_results_dialog() (`LLSFitter_GUI`)

Display loaded results in a dialog (view only)

### _restore_session() (`LLSFitter_GUI`)

Restore GUI fields from loaded results data

### export_current_plot() (`LLSFitter_GUI`)

Export the current plot to an image file

### view_detailed_results() (`LLSFitter_GUI`)

Show detailed results in a dialog

### show_about() (`LLSFitter_GUI`)

Show about dialog

### show_help() (`LLSFitter_GUI`)

Show help dialog

### submit_redshift() (`LLSFitter_GUI`)

Update the redshift when the submit button is clicked or Enter is pressed

### toggle_y_limits() (`LLSFitter_GUI`)

Enable or disable custom y-limit inputs based on checkbox state

### browse_file() (`LLSFitter_GUI`)

Open a file dialog to select a spectrum file

### load_spectrum() (`LLSFitter_GUI`)

Load the spectrum from the selected file

### get_current_continuum_regions() (`LLSFitter_GUI`)

Get the current continuum regions from the table

### update_lls_fitter_params() (`LLSFitter_GUI`)

Update the LLSFitter object with the current GUI parameters

### preview_continuum_regions() (`LLSFitter_GUI`)

Preview the continuum regions on the spectrum

### run_curve_fit() (`LLSFitter_GUI`)

Run the curve_fit process

### run_mcmc_fit() (`LLSFitter_GUI`)

Run the MCMC fitting process

### save_settings() (`LLSFitter_GUI`)

Save current settings

### load_settings() (`LLSFitter_GUI`)

Load saved settings

### closeEvent() (`LLSFitter_GUI`)

Handle close event to save settings

### setup_ui() (`rb_interactive_mask`)

Set up the user interface.

### toggle_mode() (`rb_interactive_mask`)

Toggle between add and remove mask modes.

### toggle_optimization() (`rb_interactive_mask`)

Toggle optimization options.

### update_plots() (`rb_interactive_mask`)

Update the plots with current data and masks.

### plot_masks() (`rb_interactive_mask`)

Plot mask regions on the given axis.

### on_click() (`rb_interactive_mask`)

Handle mouse click events.

### handle_add_mask() (`rb_interactive_mask`)

Handle adding mask regions.

### handle_remove_mask() (`rb_interactive_mask`)

Handle removing mask regions with right clicks.

### merge_masks() (`rb_interactive_mask`)

Merge overlapping mask regions.

### on_key_press() (`rb_interactive_mask`)

Handle keyboard events.

### reset_masks() (`rb_interactive_mask`)

Reset all masks.

### undo_last_mask() (`rb_interactive_mask`)

Remove the last mask pair.

### auto_mask() (`rb_interactive_mask`)

Automatically generate masks using improved detection of spectral features.
        
        This method uses a robust approach to detect both absorption and emission features,
        as well as regions with unreliable error estimates.

### manual_mask_entry() (`rb_interactive_mask`)

Open a dialog for manual entry of mask values.

### fit_continuum() (`rb_interactive_mask`)

Fit continuum using current parameters and masks.

### check_unmasked_points() (`rb_interactive_mask`)

Check if there are enough unmasked points for fitting.

### accept_results() (`rb_interactive_mask`)

Accept results and close window.

### cancel() (`rb_interactive_mask`)

Cancel and close without saving.

### closeEvent() (`rb_interactive_mask`)

Handle window close event.

### update_results_display() (`rb_interactive_mask`)

Update the text display with current fitting results.

### launch_interactive_mask() (`rb_interactive_mask`)

Launch the interactive masking GUI.
    
    Parameters
    ----------
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    velocity : array
        Velocity array
    existing_masks : list, optional
        Existing mask regions
    order : int, optional
        Initial polynomial order
    use_weights : bool, optional
        Whether to use weights
    domain : list, optional
        Velocity domain limits
        
    Returns
    -------
    dict
        Results including masks, continuum, and fit parameters,
        or {'cancelled': True} if the user cancelled.

### show_version() (`launch_specgui`)

Show version information.

### launch_single_gui() (`launch_specgui`)

Launch the single-spectrum GUI.

### launch_batch_gui() (`launch_specgui`)

Launch the batch processing GUI.

### main() (`launch_specgui`)

Main entry point for the unified launcher.

### launch() (`launch_specgui`)

Entry point for setup.py console scripts.

### onclick() (`PlotSpec_Integrated`)

Mouse click
            Left == 1; Right == 3

### specplot() (`PlotSpec_Integrated`)

Updates the plot with smoothed spectrum.
        This version simply updates the y-data of the existing line objects.

### clearstuff() (`PlotSpec_Integrated`)

clearstuff() ensures that all the plots from the previously drawn page are cleared before plotting the next page

### print_help() (`rb_cont`)

Print help information about the interactive continuum fitter.

### load_spectrum() (`rb_cont`)

Load a spectrum from a file.
    
    Parameters
    ----------
    filename : str
        Path to the file containing the spectrum.
    filetype : str, optional
        Type of file - 'fits', 'ascii', 'p' (pickle), or 'xfits'.
        If None, it will be inferred from the file extension.
    
    Returns
    -------
    tuple
        (wave, flux, error, metadata) where:
        - wave is the wavelength array
        - flux is the flux array
        - error is the error array
        - metadata is a dictionary with additional information
    
    Raises
    ------
    FileNotFoundError
        If the file doesn't exist.
    ValueError
        If the file type is not supported or the file doesn't contain valid data.
    ImportError
        If required modules are not available.

### save_spectrum() (`rb_cont`)

Save a spectrum with its fitted continuum to a file.
    
    Parameters
    ----------
    wave : array-like
        Wavelength array.
    flux : array-like
        Flux array.
    error : array-like
        Error array.
    continuum : array-like
        Fitted continuum array.
    filename : str
        Path to the original input file.
    filetype : str
        Type of file to save - 'fits', 'ascii', or 'p' (pickle).
    metadata : dict, optional
        Additional metadata to include in the file.
    
    Returns
    -------
    str
        Path to the saved file.
    
    Raises
    ------
    ValueError
        If the file type is not supported or if there's an error saving the file.
    ImportError
        If required modules are not available.
    
    Notes
    -----
    The output file will be saved with the same name as the input file but with
    '_norm' appended before the extension. For example, if the input file is
    'spectrum.fits', the output file will be 'spectrum_norm.fits'.

### main() (`rb_cont`)

Main function to run the interactive continuum fitting script.
    
    This function parses command-line arguments, loads the spectrum,
    launches the interactive fitter, and saves the results.

### __init__() (`rb_fit_interactive_continuum`)

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

### onclick() (`rb_fit_interactive_continuum`)

Handle mouse click events for selecting continuum points.
        
        When a user left-clicks on the plot (and no toolbar buttons are active),
        this function computes the median flux value in a 5 unit window around
        the clicked x-coordinate and adds it as a continuum point.
        
        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            The mouse click event.

### onpick() (`rb_fit_interactive_continuum`)

Handle pick events for removing continuum points.
        
        When a user right-clicks on a continuum point, this function removes it.
        
        Parameters
        ----------
        event : matplotlib.backend_bases.PickEvent
            The pick event.

### ontype() (`rb_fit_interactive_continuum`)

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

### get_continuum() (`rb_fit_interactive_continuum`)

Get the fitted continuum.
        
        Returns
        -------
        ndarray or None
            The fitted continuum if available, otherwise None.

### get_version() (`rb_spec`)

Return the current version of rb_spec.

### load_rb_spec_object() (`rb_spec`)

Load an rb_spec object from a JSON file and populate its attributes with precomputed values.

    This function reads a JSON file containing precomputed spectral data, initializes an 
    `rb_spec` object using key spectral arrays, and then dynamically sets all remaining 
    attributes from the JSON data. Lists in the JSON file are converted to `numpy` arrays 
    for consistency and efficient numerical operations.

    Parameters
    ----------
    filename : str
        Path to the JSON file containing the spectral data.
    verbose : bool, optional
        If True, prints a message when loading is complete. Default is True.

    Returns
    -------
    rb_spec
        An instance of the `rb_spec` class with all attributes loaded from the JSON file.
        Returns None if there is an error in loading the JSON file.

### __init__() (`rb_spec`)

Initializes a Spectrum object for absorption line analysis.
    
        Parameters:
        wave (array-like): Wavelength values.
        flux (array-like): Flux values.
        error (array-like): Error values.
        filename (str, optional): Filename if reading from a file. Default is False.
        
        Raises:
        ValueError: If input arrays have mismatched lengths or contain invalid values.

### version() (`rb_spec`)

Return the rb_spec version.

### from_file() (`rb_spec`)

Creates a Spectrum object from a file.
    
        Parameters:
        filename (str): Path to the file.
        filetype (str, optional): Type of file (e.g., 'ascii', 'fits', 'HSLA', etc.). If False, it is inferred.
        efil (str, optional): Error file if separate.
        kwargs: Additional arguments for specific file readers.
    
        Returns:
        Spectrum: An instance of the Spectrum class.
    
        Raises:
        FileNotFoundError: If the file does not exist.
        ImportError: If required modules cannot be imported.
        ValueError: If file type is not supported or file content is invalid.
        IOError: If there's an error reading the file.

### from_data() (`rb_spec`)

Creates a Spectrum object from given data arrays.
    
        Parameters:
        wave (array-like): Wavelength values.
        flux (array-like): Flux values.
        error (array-like): Error values.
    
        Returns:
        Spectrum: An instance of the Spectrum class.
    
        Raises:
        ValueError: If input arrays have mismatched lengths or contain invalid values.

### shift_spec() (`rb_spec`)

Shifts wavelength to absorber rest frame

### slice_spec() (`rb_spec`)

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

### fit_continuum_interactive() (`rb_spec`)

Launch an interactive GUI for continuum fitting with masking.
        
        This method opens a PyQt5 GUI that allows visual selection of mask regions
        and interactive continuum fitting.
        
        Parameters
        ----------
        mask : list, optional
            Initial mask regions to display. If not provided, uses existing masks.
        domain : list, optional
            Velocity domain limits [vmin, vmax] for fitting.
        order : int, optional
            Initial polynomial order for fitting.
        use_weights : bool, optional
            Whether to use flux errors as weights in fitting.
        
        Returns
        -------
        None
            Updates the object in-place with new masks and continuum fit.
        
        Notes
        -----
        The GUI allows:
        - Left-click pairs to add mask regions
        - Right-click to remove mask regions
        - Button to guess masks using sigma clipping
        - Adjustment of fitting parameters
        - Preview of the normalized spectrum
        
        See Also
        --------
        fit_continuum : Standard continuum fitting method

### fit_continuum() (`rb_spec`)

Fit continuum to the sliced spectrum using multiple methods.
        
        By default uses a Legendre polynomial fit. With Interactive=True,
        launches a GUI for interactive continuum fitting.
        
        Parameters
        ----------
        mask : list or False, optional
            Velocity ranges to mask during fitting, e.g., [vmin1, vmax1, vmin2, vmax2, ...].
            If False, no regions are masked.
        domain : list or False, optional
            Velocity domain limits [vmin, vmax] for fitting.
            If False, defaults to [-600, 600] km/s.
        Legendre : int or False, optional
            Order of Legendre polynomial to fit. If False, uses interactive fitting.
        Interactive : bool, optional
            If True, launches interactive continuum fitting GUI. Default is False.
        classic : bool, optional
            If True and Interactive=True, uses the classic GUI instead of the new one.
            Default is False (use new GUI).
        use_weights : bool, optional
            If True, uses flux errors as weights in fitting. Default is False.
        
        Other Parameters
        ----------------
        optimize_cont : bool, optional
            If True, uses BIC to determine the optimal polynomial order. Default is False.
        n_sigma : float, optional
            Sigma clipping threshold for outlier rejection. Default is 3.
        min_order : int, optional
            Minimum polynomial order to try when optimize_cont=True. Default is 1.
        max_order : int, optional
            Maximum polynomial order to try when optimize_cont=True. Default is 6.
        sigma_clip : bool, optional
            If True, uses sigma clipping during polynomial fitting. Default is False.
        prefit_cont : array, optional
            Predefined continuum array to use instead of fitting. Must match length of slice.
        
        Returns
        -------
        None
            Updates the object in-place with fitted continuum and normalized spectrum.

### fit_polynomial_ransac() (`rb_spec`)

Alternate continuum fitting, using ransac to fit polynomial
          degree: polynomial order

        residual_threshold : important RANSAC paramerter to identify inmask points to fit

        _n_bootstrap= number of bootstrap sampling to estimate fitting uncertainty [default = 100]

### fit_continuum_ransac() (`rb_spec`)

Alternate continuum fitting method. Does iterative ransac continumm fitting.

### plot_continuum_fit() (`rb_spec`)

Plot the fitted continuum and mark out masked regions if available.
        
        This function creates a two-panel plot showing:
        1. Original flux with fitted continuum and error
        2. Normalized spectrum with error
        
        Parameters
        ----------
        outfilename : str, optional
            If provided, the plot will be saved to this filename
        xlim : list, optional
            Velocity limits for x-axis [vmin, vmax]. If None, uses full range.
        mask_alpha : float, optional
            Transparency level for masked regions (default: 0.2)
            
        Returns
        -------
        matplotlib.figure.Figure
            The created figure object for further customization if needed
            
        Notes
        -----
        This method should be called after `fit_continuum` has been executed.
        
        Examples
        --------
        >>> spec.fit_continuum(mask=[-200, 300, 500, 1100], domain=[-1500, 1500], Legendre=3)
        >>> spec.plot_continuum_fit()
        >>> plt.show()  # Display the plot
        
        # Save the plot to a file
        >>> spec.plot_continuum_fit(outfilename='continuum_fit.png')
        
        # Customize the plot further
        >>> fig = spec.plot_continuum_fit()
        >>> fig.suptitle('My Custom Title')
        >>> plt.show()

### compute_EW() (`rb_spec`)

Computes rest frame equivalent width and column density for a desired atomic line.
        Around the species lam_cen and given vmin and vmax keyword values.

### plot_spec() (`rb_spec`)

Quick wrapper to call an interactive plotter for the full spectrum as given in input file.

### plot_slice() (`rb_spec`)

Quick wrapper to call an interactive plotter for the full spectrum as given in input file.

### save_slice() (`rb_spec`)

Saves the slice object for future processing.

        Parameters:
        -----------
        outfilename : str
            The file path to save the slice object.
        file_format : str, optional
            Format to save the object. Options: 'pickle' (default) or 'json'.

        Notes:
        ------
        - Pickle saves the entire object for later editing.
        - JSON saves only the output data, not the object, so it cannot be reloaded for editing.

### display_field_info() (`rb_spec`)

Display information about fields in the rb_spec object.
        
        Parameters
        ----------
        field : str, optional
            Specific field to display information about.
            If None, displays information about all fields.

### plot_doublet() (`rb_spec`)

Plot a given doublet defined by the lam1 and lam2 wavelength centers.

### ontype() (`zgui.spec_fit_gauss2d`)

Interactivae keyboard events
        Note:
            Always Update to help mannual for new keyboard events

### replot() (`zgui.spec_plot`)

Update existing line data instead of recreating lines

### _compute_distance() (`zgui.spec_plot`)

Compute the distance between the event xydata and selected point for Gaussian fitting

### clear_interactive_elements() (`zgui.spec_plot`)

Remove all interactive elements except base flux/error lines

### add_emission_marker() (`zgui.spec_plot`)

Add emission line marker and track it

### add_gaussian_fit() (`zgui.spec_plot`)

Add Gaussian fit curve with modern styling

### extract_1d() (`zgui.spec_plot`)

Extract 1d spectrum and error spectrum from input 2d spectrum

### plot_spec2d() (`zgui.spec_plot`)

Display 2D spec in top panel, 1D extraction in bottom panel
		self.flux2d, self.err2d - 2D spec info
		self.flux, self.error - 1D spec extraction

### replot2d() (`zgui.spec_plot`)

Update 1D spectrum with modern matplotlib practices

### replot2d_im() (`zgui.spec_plot`)

Re-plot smoothed/unsmoothed 2D spectrum

### ontype() (`zgui.spec_plot`)

Interactivae keyboard events
		Note:
			Always Update to help mannual for new keyboard events

### onclick() (`zgui.spec_plot`)

Mouse click
			Left == 1; Right == 3

### _on_sent_fitsobj() (`zgui.spec_plot`)

Receive FitsObj and store available frame sources.

### _on_frame_selected() (`zgui.spec_plot`)

Update the currently active frame for measurements.

		Called when user selects a different frame in the toolbar.

### _is_valid_number() (`zgui.menu_toolbars`)

Check if text can be converted to a float

### _colorbar_scale_toggled() (`zgui.menu_toolbars`)

Handle colorbar scale checkbox toggle

### _ensure_frame_columns() (`zgui.tableview_pandas`)

Ensure frame-specific columns exist for a given frame.
		Creates z_FRAMENAME and z_err_FRAMENAME columns if they don't exist.

### _on_context_menu() (`zgui.tableview_pandas`)

Handle right-click context menu on table.

### _set_primary_frame() (`zgui.tableview_pandas`)

Set a frame as the primary frame for a row.

		Updates:
		1. z and z_err columns with values from selected frame
		2. primary_frame column to track which frame is primary

### _build_wave() (`zgui.gui_io`)

Returns a NumPy array containing wavelength axis of a 2d specturm in Angstrom.
			Args:
				header (astropy.io.fits.Header): header that contains wavelength axis
				that is specified in 'CTYPE' keywords in NAXIS1 dimension.
			Returns:
				numpy.ndarray: Wavelength axis for this data.

### _detect_frame_sources() (`zgui.gui_io`)

Detect available frame sources from FITS file headers.

		Returns:
			list: Uppercase frame names found in FITS file (e.g., ['EMLINEA', 'EMLINEB'])

### setup_ui() (`zgui.main`)

Setup the user interface

### setup_connections() (`zgui.main`)

Setup signal connections between widgets

### launch_gui() (`zgui.main`)

Launch the zgui GUI application.
    
    Args:
        xspecio (bool): Whether to use XSpectrum1D for IO
        toggle_frames (bool): Enable frame toggling
        fitsfile (str): Optional FITS file to load on startup

### main() (`zgui.main`)

Entry point for the zgui command line interface

### prepare_absorber_object() (`multispecviewer.vStack`)

Prepare absorption line data for plotting
    
    Parameters:
    -----------
    z_abs : float
        Redshift of the absorber
    wave : array
        Wavelength array
    flux : array
        Flux array
    error : array
        Error array
    line_flg : str, optional
        Line list flag (default: 'LLS')
    vlim : list, optional
        Velocity limits (default: [-1000, 1000])
        
    Returns:
    --------
    dict
        Dictionary of ions data

### __init__() (`multispecviewer.vStack`)

Initialize vStack class for displaying velocity plots of absorption lines
        
        Parameters:
        -----------
        parent : object
            Parent object (usually SpectralPlot)
        wave : array
            Wavelength array
        flux : array
            Flux array
        error : array
            Error array
        line_flg : str
            Line list flag
        zabs : float, optional
            Redshift of the absorber (default: 0)
        vlim : list, optional
            Velocity limits (default: [-1000., 1000.])

### onkb() (`multispecviewer.vStack`)

Handle keyboard events
        
        Parameters:
        -----------
        event : matplotlib.backend_bases.KeyEvent
            Keyboard event

### vPlot() (`multispecviewer.vStack`)

Plot velocity panels
        
        Parameters:
        -----------
        ploti : int or array, optional
            Indices of panels to plot
        comment : str, optional
            Comment to add to plots
        yrange : list, optional
            Y-axis limits
        clearpage : bool, optional
            Whether to clear the page

### clearstuff() (`multispecviewer.vStack`)

Clear a specific panel
        
        Parameters:
        -----------
        i : int
            Panel index

### plotstuff() (`multispecviewer.vStack`)

Plot a specific velocity panel
        
        Parameters:
        -----------
        i : int
            Panel index
        comment : str, optional
            Comment to add to plot
        yrange : list, optional
            Y-axis limits

### plotText() (`multispecviewer.vStack`)

Return text description based on flag value
        
        Parameters:
        -----------
        flag : int, optional
            Flag value (default: 1)
            
        Returns:
        --------
        str
            Text description

### on_sent_message() (`multispecviewer.MessageBox`)

Send a message with optional color
        
        :param sent_message: Message text
        :param hexColor: Hex color code for the message (default white)

### append_message() (`multispecviewer.MessageBox`)

Append a message with optional color to the existing messages
        
        :param sent_message: Message text
        :param hexColor: Hex color code for the message (default white)

### clear() (`multispecviewer.MessageBox`)

Clear all messages from the text edit widget

### on_catalog_clicked() (`multispecviewer.RedshiftInputWidget`)

Handle catalog button click

### on_linelist_changed() (`multispecviewer.RedshiftInputWidget`)

Handle linelist selection change

### on_color_changed() (`multispecviewer.RedshiftInputWidget`)

Handle color selection change

### validate_and_submit() (`multispecviewer.RedshiftInputWidget`)

Validate inputs and emit the submitted signal if valid

### set_redshift() (`multispecviewer.RedshiftInputWidget`)

Sets the redshift value in the input field programmatically
        
        :param redshift: The redshift value to set

### display_examples() (`multispecviewer.rb_multispec`)

Display detailed usage examples for the rb_multispec tool.

### main() (`multispecviewer.rb_multispec`)

Wrapper script for rbcodes multispecviewer tools.
    Allows users to choose between the new version and the classic version.

### from_data() (`multispecviewer.rb_multispec`)

Launch MultispecViewer GUI with spectral data from arrays or spectrum objects.

### _prepare_spectra_list() (`multispecviewer.rb_multispec`)

Convert various input formats to a list of spectrum objects compatible with MultispecViewer.
    
    Parameters
    ----------
    spectrum_data : various
        Input spectrum data in various formats
        
    Returns
    -------
    list
        List of spectrum objects (rb_spectrum or XSpectrum1D)

### launch_empty() (`multispecviewer.rb_multispec`)

Launch empty MultispecViewer GUI.

### __new__() (`multispecviewer.io_manager`)

Ensure only one instance of IOManager exists (singleton pattern)

### _initialize() (`multispecviewer.io_manager`)

Initialize the IOManager with default settings

### set_message_box() (`multispecviewer.io_manager`)

Set the message box reference for displaying messages

### show_message() (`multispecviewer.io_manager`)

Display a message in the message box if available, otherwise print

### append_message() (`multispecviewer.io_manager`)

Append a message to the message box if available, otherwise print

### load_fits_files() (`multispecviewer.io_manager`)

Load FITS files with rb_spectrum fallback to XSpectrum1D.
        
        Args:
            file_paths (list, optional): List of file paths to load.
            
        Returns:
            list: List of spectrum objects (rb_spectrum or XSpectrum1D)

### save_line_list() (`multispecviewer.io_manager`)

Save a line list DataFrame to a file.
        
        Args:
            line_list (pd.DataFrame): DataFrame containing the line list
            file_path (str, optional): Path to save to. If None, use the last directory.
            format (str, optional): Format to save in ('txt', 'csv', 'json'). 
                                    If None, determine from file extension.
        
        Returns:
            bool: Success status
            str: Error message if any

### load_line_list() (`multispecviewer.io_manager`)

Load a line list from a file.
        
        Args:
            file_path (str, optional): Path to the file to load
            
        Returns:
            pd.DataFrame: The loaded line list
            str: Error message if any

### save_absorbers() (`multispecviewer.io_manager`)

Save absorbers DataFrame to a file.
        
        Args:
            absorbers_df (pd.DataFrame): DataFrame containing absorber data
            file_path (str, optional): Path to save to. If None, use the last directory.
            format (str, optional): Format to save in ('csv', 'json'). 
                                     If None, determine from file extension.
        
        Returns:
            bool: Success status
            str: Error message if any

### load_absorbers() (`multispecviewer.io_manager`)

Load absorber data from a file.
        
        Args:
            file_path (str, optional): Path to the file to load
            
        Returns:
            pd.DataFrame: The loaded absorber data
            str: Error message if any

### save_combined_data() (`multispecviewer.io_manager`)

Save line list, absorber data, and metadata to a single JSON file.
        
        Args:
            line_list (pd.DataFrame): Line list DataFrame
            absorbers_df (pd.DataFrame): Absorbers DataFrame
            spectrum_files (list, optional): List of spectrum filenames
            file_path (str, optional): Path to save to
            user_comment (str, optional): User-provided comment about the data
            metadata (dict, optional): Additional metadata to include
            
        Returns:
            bool: Success status
            str: Error message if any

### load_combined_data() (`multispecviewer.io_manager`)

Load combined data from a JSON file.
        
        Args:
            file_path (str, optional): Path to the file to load
            
        Returns:
            tuple: (line_list, absorbers_df, spectrum_files, metadata, error_message)
                   Any of these may be None if not found or if an error occurs

### convert_text_to_json() (`multispecviewer.io_manager`)

Convert text/CSV files to a combined JSON file.
        
        Args:
            line_list_path (str, optional): Path to line list text/CSV file
            absorbers_path (str, optional): Path to absorbers CSV file
            output_path (str, optional): Path for the output JSON file
            user_comment (str, optional): User-provided comment
            spectrum_files (list, optional): List of spectrum filenames
            metadata (dict, optional): Additional metadata
            
        Returns:
            bool: Success status
            str: Error message if any

### convert_json_to_text() (`multispecviewer.io_manager`)

Convert a combined JSON file to separate text/CSV files.
        
        Args:
            json_path (str, optional): Path to the JSON file
            output_dir (str, optional): Directory to save the output files
            
        Returns:
            tuple: (line_list_path, absorbers_path, metadata_path, error_message)
                   Any of these may be None if not applicable or if an error occurs

### get_user_comment_dialog() (`multispecviewer.io_manager`)

Create a dialog for entering metadata when saving files.

            Args:
                parent (QWidget, optional): Parent widget for the dialog
                previous_metadata (dict, optional): Dictionary with previously saved metadata to pre-fill

            Returns:
                tuple: (comment, metadata_dict, accepted)
                       comment: string with user's comment
                       metadata_dict: dictionary with additional metadata
                       accepted: boolean indicating if dialog was accepted

### integrated_save_data() (`multispecviewer.io_manager`)

Integrated method for saving data with format selection.
        Shows file dialog and optional metadata dialog.
        
        Args:
            parent (QWidget): Parent widget for dialogs
            line_list (pd.DataFrame, optional): Line list data
            absorbers_df (pd.DataFrame, optional): Absorber data
            spectrum_filenames (list, optional): List of spectrum filenames
            
        Returns:
            bool: Success status

### integrated_load_data() (`multispecviewer.io_manager`)

Integrated method for loading data with format detection.
        Shows file dialog and loads data accordingly.
        
        Args:
            parent (QWidget): Parent widget for dialogs
            
        Returns:
            tuple: (line_list, absorbers_df, error_message)

### initUI() (`multispecviewer.LineSelectionDialog`)

Initialize the user interface

### apply_theme() (`multispecviewer.LineSelectionDialog`)

Apply the dark theme styling consistent with other widgets

### on_selection_changed() (`multispecviewer.LineSelectionDialog`)

Handle item selection changes for immediate plotting

### on_select_clicked() (`multispecviewer.LineSelectionDialog`)

Handle select button click

### on_item_double_clicked() (`multispecviewer.LineSelectionDialog`)

Handle double-click on an item

### process_selection() (`multispecviewer.LineSelectionDialog`)

Process the selected item

### read_line_options() (`multispecviewer.utils`)

Read available line list options from the configuration file.
    
    Returns:
        list: List of line list option names. Returns default options if the
              configuration file cannot be read.

### show_help_dialog() (`multispecviewer.utils`)

Display a help dialog with keyboard shortcuts and usage information.
    
    Parameters:
    -----------
    parent : QWidget, optional
        Parent widget for the dialog

### reconcile_linelists() (`multispecviewer.utils`)

Reconcile multiple linelist files by merging duplicate entries based on line name and velocity separation.
    
    Parameters:
    -----------
    input_files : list
        List of file paths to linelist files (can be txt, csv, or json format)
    velocity_threshold : float, optional
        Maximum velocity difference in km/s to consider lines as duplicates (default: 20)
    output_file : str, optional
        Path to save the reconciled linelist (if None, doesn't save to file)
    create_absorber_df : bool, optional
        Whether to create an absorber DataFrame from unique redshifts (default: True)
        
    Returns:
    --------
    tuple
        (reconciled_linelist, absorber_df) where reconciled_linelist is a DataFrame with
        merged line entries and absorber_df is a DataFrame of unique absorber systems

### on_mouse_press() (`multispecviewer.multispec`)

Handles mouse press events. Right-click shows a menu of possible spectral lines.

### plot_spectra() (`multispecviewer.multispec`)

Plots the given list of spectra as subplots.
        
        :param spectra: List of spectra objects from linetools.XSpectrum1D

### on_key_press() (`multispecviewer.multispec`)

Handles key press events for interactive adjustments of the plots.

### reset_view() (`multispecviewer.multispec`)

Resets the view to original x and y limits and removes all lines.
        Called when 'r' key is pressed.

### replot() (`multispecviewer.multispec`)

Re-plot smoothed/unsmoothed spectrum

### plot_one_spec() (`multispecviewer.multispec`)

Plot a single spectrum panel in step mode with specific dark theme colors

### set_redshift_data() (`multispecviewer.multispec`)

Receives redshift and linelist data from the main window.

        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        :param color: str, the color to use for line plotting (default 'sky_blue')

### plot_redshift_lines() (`multispecviewer.multispec`)

Plots spectral lines based on the redshift and linelist.

### clear_redshift_lines() (`multispecviewer.multispec`)

Removes any previously plotted redshift lines.

### clear_quickid_lines() (`multispecviewer.multispec`)

Removes any previously plotted quick ID lines.

### check_lineid() (`multispecviewer.multispec`)

This method quickly draws some doublet/multiplet lines on the canvas for a quicklook
        
        :param wave0: The observed wavelength where the user clicked
        :param ionname: The ion to identify (e.g., 'CIV', 'MgII')
        :param yval: The y-position where the user clicked
        :param ax_index: The index of the active axis

### plot_absorber_lines() (`multispecviewer.multispec`)

Wrapper to call the canvas's plot_absorber_lines method

### remove_absorber_lines() (`multispecviewer.multispec`)

Wrapper to call the canvas's remove_absorber_lines method

### clear_all_absorber_lines() (`multispecviewer.multispec`)

Remove all absorber lines from the plot and reset all checkboxes in the absorber manager

### handle_redshift_submission() (`multispecviewer.multispec`)

Processes the redshift and linelist data submitted from the RedshiftInputWidget.
        Creates the appropriate lines on the plot.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        :param color: str, the selected color

### handle_linelist_changed() (`multispecviewer.multispec`)

Handle linelist or color selection change in the redshift widget
        
        :param redshift: float, the current redshift value
        :param linelist: str, the newly selected line list
        :param color: str, the selected color

### handle_catalog_clicked() (`multispecviewer.multispec`)

Handles the catalog button click. Adds the current redshift to the absorber manager.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list
        :param color: str, the selected color

### update_absorber_redshift() (`multispecviewer.multispec`)

Update the redshift value for an absorber and replot its lines.
        
        :param row: The row index in the absorber manager
        :param z_abs: The new redshift value

### select_fits_files() (`multispecviewer.multispec`)

Opens a file dialog for selecting FITS files and loads them for plotting.

### load_fits_files() (`multispecviewer.multispec`)

Loads and plots the selected FITS files using the IO Manager.

### handle_load_clicked() (`multispecviewer.multispec`)

Handle the Load button click event.
        Uses IOManager to load line list and absorber data.
        Provides option to append or overwrite existing data.

### handle_save_clicked() (`multispecviewer.multispec`)

Handle the Save button click event.
        Uses IOManager to save line list and absorber data.

### handle_show_clicked() (`multispecviewer.multispec`)

Handle the Show button click event.
        Displays/hides all identified lines from line_list on the plot.
        Toggle behavior: first click shows lines, second click hides them.

### display_line_list() (`multispecviewer.multispec`)

Display a sortable table of all identified lines with options to select and delete entries.

### handle_convert_clicked() (`multispecviewer.multispec`)

Handle conversion between file formats.

### populate_table_from_df() (`multispecviewer.AbsorberManager`)

Populate the table from an existing DataFrame

### initUI() (`multispecviewer.AbsorberManager`)

Initialize the user interface components

### _add_row_widgets() (`multispecviewer.AbsorberManager`)

Add widgets to a specific row in the table

### toggle_absorber_from_checkbox() (`multispecviewer.AbsorberManager`)

Handle checkbox state change

### remove_absorber_from_button() (`multispecviewer.AbsorberManager`)

Remove absorber using row index stored on the button

### _populate_row() (`multispecviewer.AbsorberManager`)

Populate a row with absorber data

### add_absorber() (`multispecviewer.AbsorberManager`)

Add a new absorber to the manager

### plot_absorber() (`multispecviewer.AbsorberManager`)

Plot the absorber at the specified row

### remove_absorber() (`multispecviewer.AbsorberManager`)

Remove the absorber at the specified row

### on_cell_changed() (`multispecviewer.AbsorberManager`)

Handle manual edits to cells

### get_absorber_count() (`multispecviewer.AbsorberManager`)

Return the number of defined absorbers

### get_absorber_data() (`multispecviewer.AbsorberManager`)

Get the data for a specific absorber

### get_all_absorber_data() (`multispecviewer.AbsorberManager`)

Get all absorber data as a DataFrame

### setup_dark_theme() (`multispecviewer.AbsorberManager`)

Apply a dark theme matching the main application's color scheme

### on_linelist_changed() (`multispecviewer.classic.RedshiftInputWidget_classic`)

Handle linelist selection change

### validate_and_submit() (`multispecviewer.classic.RedshiftInputWidget_classic`)

Validate inputs and emit the submitted signal if valid

### get_values() (`multispecviewer.classic.RedshiftInputWidget_classic`)

Get the current values from the inputs (can be used as an alternative to the signal)

### set_values() (`multispecviewer.classic.RedshiftInputWidget_classic`)

Set the values of the input fields

### set_redshift() (`multispecviewer.classic.RedshiftInputWidget_classic`)

Sets the redshift value in the input field programmatically
        
        :param redshift: The redshift value to set

### on_mouse_press() (`multispecviewer.classic.multispec_classic`)

Handles mouse press events. Right-click shows a menu of possible spectral lines.

### plot_spectra() (`multispecviewer.classic.multispec_classic`)

Plots the given list of spectra as subplots.
        
        :param spectra: List of spectra objects from linetools.XSpectrum1D

### on_key_press() (`multispecviewer.classic.multispec_classic`)

Handles key press events for interactive adjustments of the plots.

### reset_view() (`multispecviewer.classic.multispec_classic`)

Resets the view to original x and y limits and removes all lines.
        Called when 'r' key is pressed.

### replot() (`multispecviewer.classic.multispec_classic`)

Re-plot smoothed/unsmoothed spectrum

### plot_one_spec() (`multispecviewer.classic.multispec_classic`)

Plot a single spectrum panel

### set_redshift_data() (`multispecviewer.classic.multispec_classic`)

Receives redshift and linelist data from the main window.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list

### plot_redshift_lines() (`multispecviewer.classic.multispec_classic`)

Plots spectral lines based on the redshift and linelist.
        This method draws vertical lines at the expected wavelengths of 
        common spectral features adjusted by the redshift.

### clear_redshift_lines() (`multispecviewer.classic.multispec_classic`)

Removes any previously plotted redshift lines.

### clear_quickid_lines() (`multispecviewer.classic.multispec_classic`)

Removes any previously plotted quick ID lines.

### check_lineid() (`multispecviewer.classic.multispec_classic`)

This method quickly draws some doublet/multiplet lines on the canvas for a quicklook
        
        :param wave0: The observed wavelength where the user clicked
        :param ionname: The ion to identify (e.g., 'CIV', 'MgII')
        :param yval: The y-position where the user clicked
        :param ax_index: The index of the active axis

### handle_redshift_submission() (`multispecviewer.classic.multispec_classic`)

Processes the redshift and linelist data submitted by the RedshiftInputWidget
        and sends it to the SpectralPlot instance.
        
        :param redshift: float, the redshift value
        :param linelist: str, the selected line list

### select_fits_files() (`multispecviewer.classic.multispec_classic`)

Opens a file dialog for selecting FITS files and loads them for plotting.

### load_fits_files() (`multispecviewer.classic.multispec_classic`)

Loads and plots the selected FITS files using XSpectrum1D directly.
        Adds error spectrum as 5% of flux if none exists.

### initUI() (`multispecviewer.classic.LineSelectionDialog_classic`)

Initialize the user interface

### on_select_clicked() (`multispecviewer.classic.LineSelectionDialog_classic`)

Handle select button click

### on_item_double_clicked() (`multispecviewer.classic.LineSelectionDialog_classic`)

Handle double-click on an item

### process_selection() (`multispecviewer.classic.LineSelectionDialog_classic`)

Process the selected item

### reset_state() (`specgui.spectrum_controller`)

Reset all analysis-related attributes.

### load_from_file() (`specgui.spectrum_controller`)

Load spectrum from file using rb_spec.

### load_from_data() (`specgui.spectrum_controller`)

Load spectrum from numpy arrays.

### load_from_json() (`specgui.spectrum_controller`)

Load spectrum from saved JSON analysis.

### get_json_info() (`specgui.spectrum_controller`)

Get information from a loaded JSON file.

### has_spectrum() (`specgui.spectrum_controller`)

Check if a spectrum is loaded.

### get_spectrum_data() (`specgui.spectrum_controller`)

Get the current spectrum data for plotting.

### apply_redshift() (`specgui.spectrum_controller`)

Apply redshift to the spectrum using rb_spec.

### reset_after_redshift() (`specgui.spectrum_controller`)

Reset analysis after redshift change, preserving the loaded spectrum.

### reset_after_transition() (`specgui.spectrum_controller`)

Reset analysis after transition change, preserving spectrum and redshift.

### slice_spectrum() (`specgui.spectrum_controller`)

Slice the spectrum around a transition using rb_spec.

### get_sliced_data() (`specgui.spectrum_controller`)

Get the sliced spectrum data if available.

### get_transition_info() (`specgui.spectrum_controller`)

Get information about the currently selected transition.

### fit_continuum() (`specgui.spectrum_controller`)

Launch interactive continuum fitter as a dialog.

### get_continuum() (`specgui.spectrum_controller`)

Get the fitted continuum if available.

### compute_equivalent_width() (`specgui.spectrum_controller`)

Compute equivalent width using rb_spec's compute_EW method.

### has_continuum() (`specgui.spectrum_controller`)

Check if the spectrum has a fitted continuum.

### get_normalized_data() (`specgui.spectrum_controller`)

Get the normalized spectrum data.

### get_measurement_figure() (`specgui.spectrum_controller`)

Get the figure from rb_spec's compute_EW for display.

### save_analysis() (`specgui.spectrum_controller`)

Save the analysis to a file.

### export_continuum_plot() (`specgui.spectrum_controller`)

Export the continuum fit plot to a file.

### get_analysis_summary() (`specgui.spectrum_controller`)

Get a summary of the analysis for display.

### apply_flat_continuum() (`specgui.spectrum_controller`)

Apply a flat continuum (value = 1.0) to the spectrum.

### apply_current_continuum() (`specgui.spectrum_controller`)

Normalize flux/error by the *existing* continuum in self.spec.cont.

### export_measurement_plot() (`specgui.spectrum_controller`)

Export the equivalent width measurement plot directly from compute_EW.

### init_ui() (`specgui.main`)

Initialize the user interface.

### on_spectrum_loaded() (`specgui.main`)

Handle spectrum loaded signal.

### on_redshift_applied() (`specgui.main`)

Handle redshift applied signal.

### on_spectrum_sliced() (`specgui.main`)

Handle spectrum sliced signal.

### on_continuum_fitted() (`specgui.main`)

Handle continuum fitted signal.

### on_measurement_completed() (`specgui.main`)

Handle measurement completed signal.

### init_ui() (`specgui.panels.output_panel`)

Initialize the user interface.

### reset() (`specgui.panels.output_panel`)

Reset panel state when a new file is loaded.

### showEvent() (`specgui.panels.output_panel`)

Override the show event to update the filename when the panel becomes visible.

### update_default_filename() (`specgui.panels.output_panel`)

Update the default filename based on the current state.

### browse_save_path() (`specgui.panels.output_panel`)

Open file dialog to select save path.

### save_analysis() (`specgui.panels.output_panel`)

Save the analysis to a file.

### update_summary() (`specgui.panels.output_panel`)

Update the analysis summary display.

### export_plots() (`specgui.panels.output_panel`)

Export plots of the analysis using the same path as the JSON file.

### init_ui() (`specgui.panels.advanced_batch_dialog`)

Initialize the user interface.

### init_ui() (`specgui.panels.measurement_panel`)

Initialize the user interface.

### reset() (`specgui.panels.measurement_panel`)

Reset panel state when a new file is loaded.

### update_velocity_range() (`specgui.panels.measurement_panel`)

Update velocity range inputs if spectrum has changed and contains vmin/vmax values.

### toggle_interactive_selection() (`specgui.panels.measurement_panel`)

Enable or disable interactive velocity range selection.

### on_plot_click() (`specgui.panels.measurement_panel`)

Handle mouse clicks on the plot for setting vmin/vmax.

### on_plot_key_press() (`specgui.panels.measurement_panel`)

Handle key presses on the plot.

### compute_equivalent_width() (`specgui.panels.measurement_panel`)

Compute equivalent width using rb_spec's compute_EW method.

### update_results_table() (`specgui.panels.measurement_panel`)

Update the results table with the computation results.

### update_plot() (`specgui.panels.measurement_panel`)

Update the preview plot with normalized spectrum and measurement results.

### init_ui() (`specgui.panels.batch_panel`)

Initialize the user interface.

### setup_transitions_table() (`specgui.panels.batch_panel`)

Setup table for transitions batch mode.

### setup_files_table() (`specgui.panels.batch_panel`)

Setup table for files batch mode.

### toggle_batch_mode() (`specgui.panels.batch_panel`)

Switch between transitions and files batch modes.

### toggle_velocity_range() (`specgui.panels.batch_panel`)

Enable/disable custom velocity range inputs.

### load_from_file() (`specgui.panels.batch_panel`)

Load batch items from a file.

### load_from_json_file() (`specgui.panels.batch_panel`)

Load batch items from a JSON file containing rb_spec objects.

### load_from_excel_file() (`specgui.panels.batch_panel`)

Load batch items from an Excel file.

### load_from_csv_file() (`specgui.panels.batch_panel`)

Load batch items from a CSV/TSV file.

### _process_transitions_dataframe() (`specgui.panels.batch_panel`)

Process a dataframe containing transition information.

### _process_files_dataframe() (`specgui.panels.batch_panel`)

Process a dataframe containing file information.

### export_template() (`specgui.panels.batch_panel`)

Export a template file for batch input.

### add_batch_item() (`specgui.panels.batch_panel`)

Add an item to the batch processing list.

### add_transition_dialog() (`specgui.panels.batch_panel`)

Open dialog to add a transition to the batch.

### add_file_dialog() (`specgui.panels.batch_panel`)

Open dialog to add a file to the batch.

### remove_batch_item() (`specgui.panels.batch_panel`)

Remove the selected item from the batch.

### clear_batch_items() (`specgui.panels.batch_panel`)

Clear all batch items.

### run_batch() (`specgui.panels.batch_panel`)

Run the batch processing tasks.

### add_result_to_table() (`specgui.panels.batch_panel`)

Add a result to the results table.

### process_transition() (`specgui.panels.batch_panel`)

Process a transition batch item.

### process_file() (`specgui.panels.batch_panel`)

Process a file batch item.

### export_results() (`specgui.panels.batch_panel`)

Export the batch processing results.

### export_to_csv() (`specgui.panels.batch_panel`)

Export results to a CSV file.

### export_to_excel() (`specgui.panels.batch_panel`)

Export results to an Excel file.

### launch_advanced_mode() (`specgui.panels.batch_panel`)

Launch the advanced batch processing dialog.

### update_batch_table() (`specgui.panels.batch_panel`)

Update the batch table with current batch items.

### update_results_table() (`specgui.panels.batch_panel`)

Update the results table with current results.

### init_ui() (`specgui.panels.continuum_panel`)

Initialize the user interface.

### reset() (`specgui.panels.continuum_panel`)

Reset panel state when a new file is loaded.

### launch_fitter() (`specgui.panels.continuum_panel`)

Launch the interactive continuum fitter.

### update_plot() (`specgui.panels.continuum_panel`)

Update the preview plot with current spectrum and continuum.

### toggle_skip_continuum() (`specgui.panels.continuum_panel`)

Toggle between interactive fitter and flat continuum.

### apply_flat_continuum() (`specgui.panels.continuum_panel`)

Apply a flat continuum (equal to 1 everywhere).

### apply_current_continuum() (`specgui.panels.continuum_panel`)

Apply the currently stored continuum (no re-fitting).

### init_ui() (`specgui.panels.redshift_panel`)

Initialize the user interface.

### reset() (`specgui.panels.redshift_panel`)

Reset panel state when a new file is loaded.

### on_redshift_changed() (`specgui.panels.redshift_panel`)

Handle redshift spinbox value changes.

### reset_redshift() (`specgui.panels.redshift_panel`)

Reset the redshift and clear downstream analysis.

### on_common_redshift_selected() (`specgui.panels.redshift_panel`)

Handle selection from common redshifts dropdown.

### set_redshift() (`specgui.panels.redshift_panel`)

Set the redshift value (e.g., when loading from JSON).

### apply_redshift() (`specgui.panels.redshift_panel`)

Apply the redshift to the spectrum.

### update_plot() (`specgui.panels.redshift_panel`)

Update the spectrum plot.

### init_ui() (`specgui.panels.input_panel`)

Initialize the user interface.

### browse_file() (`specgui.panels.input_panel`)

Open file dialog to browse for spectrum file.

### load_from_file() (`specgui.panels.input_panel`)

Load spectrum from the selected file.

### init_ui() (`specgui.panels.transition_panel`)

Initialize the user interface.

### reset() (`specgui.panels.transition_panel`)

Reset panel state when a new file is loaded.

### reset_transition() (`specgui.panels.transition_panel`)

Reset the transition selection and clear downstream analysis.

### populate_transition_combo() (`specgui.panels.transition_panel`)

Populate the transition combo box with common transitions.

### on_transition_selected() (`specgui.panels.transition_panel`)

Handle selection from common transitions dropdown.

### on_transition_changed() (`specgui.panels.transition_panel`)

Handle direct transition wavelength changes.

### on_vmin_changed() (`specgui.panels.transition_panel`)

Handle velocity minimum changes.

### on_vmax_changed() (`specgui.panels.transition_panel`)

Handle velocity maximum changes.

### slice_spectrum() (`specgui.panels.transition_panel`)

Slice the spectrum around the selected transition.

### set_transition() (`specgui.panels.transition_panel`)

Set the transition wavelength and update UI (e.g., when loading from JSON).

### set_velocity_limits() (`specgui.panels.transition_panel`)

Set the velocity limits (e.g., when loading from JSON).

### update_plot() (`specgui.panels.transition_panel`)

Update the spectrum plot.

### _create_empty_dataframe() (`specgui.batch.master_batch_table`)

Create empty DataFrame with all required columns.

### add_item() (`specgui.batch.master_batch_table`)

Add a new batch item. Returns row index.

### remove_item() (`specgui.batch.master_batch_table`)

Remove a batch item by row index.

### get_item_count() (`specgui.batch.master_batch_table`)

Get the total number of items.

### get_all_items() (`specgui.batch.master_batch_table`)

Get all items as a list of row dictionaries.

### get_item() (`specgui.batch.master_batch_table`)

Get a single item by row index.

### update_template() (`specgui.batch.master_batch_table`)

Update template parameters for an item.

### update_analysis() (`specgui.batch.master_batch_table`)

Update analysis parameters for an item.

### update_results() (`specgui.batch.master_batch_table`)

Update results for an item.

### set_rb_spec_object() (`specgui.batch.master_batch_table`)

Store the rb_spec object for an item.

### get_rb_spec_object() (`specgui.batch.master_batch_table`)

Get the rb_spec object for an item.

### clear_all() (`specgui.batch.master_batch_table`)

Clear all items.

### validate_all() (`specgui.batch.master_batch_table`)

Validate all items. Returns (valid_indices, invalid_items_with_errors).

### get_items_by_status() (`specgui.batch.master_batch_table`)

Get row indices by processing status.

### to_dict() (`specgui.batch.master_batch_table`)

Convert entire table to dictionary for serialization.

### from_dict() (`specgui.batch.master_batch_table`)

Load table from dictionary.

### import_template_csv() (`specgui.batch.master_batch_table`)

Import template from CSV file.

### export_template_csv() (`specgui.batch.master_batch_table`)

Export just the template part as CSV.

### init_ui() (`specgui.batch.batch_main`)

Initialize the user interface.

### create_menus() (`specgui.batch.batch_main`)

Create the menu bar.

### create_status_bar() (`specgui.batch.batch_main`)

Create an enhanced status bar.

### connect_signals() (`specgui.batch.batch_main`)

Connect signals between panels and controller.

### update_tab_states() (`specgui.batch.batch_main`)

Update which tabs are enabled based on current state.

### on_configuration_changed() (`specgui.batch.batch_main`)

Handle configuration changes - updated for master table architecture.

### update_status() (`specgui.batch.batch_main`)

Update the status bar with a message.

### update_progress() (`specgui.batch.batch_main`)

Update the progress bar.

### on_batch_item_completed() (`specgui.batch.batch_main`)

Handle batch item completion.

### on_batch_item_failed() (`specgui.batch.batch_main`)

Handle batch item failure.

### on_batch_started() (`specgui.batch.batch_main`)

Handle batch processing started signal.

### on_batch_completed() (`specgui.batch.batch_main`)

Handle batch processing completed signal.

### on_export_completed() (`specgui.batch.batch_main`)

Handle export completed signal.

### on_tab_changed() (`specgui.batch.batch_main`)

Handle tab changes to provide guidance.

### update_window_title() (`specgui.batch.batch_main`)

Update the window title based on current state.

### new_configuration() (`specgui.batch.batch_main`)

Create a new batch configuration.

### open_configuration() (`specgui.batch.batch_main`)

Open an existing batch configuration.

### save_configuration() (`specgui.batch.batch_main`)

Save the current batch configuration.

### save_configuration_as() (`specgui.batch.batch_main`)

Save the batch configuration with a new name.

### create_csv_template() (`specgui.batch.batch_main`)

Create a CSV template for batch import.

### clear_all_items() (`specgui.batch.batch_main`)

Clear all batch items.

### validate_configuration() (`specgui.batch.batch_main`)

Validate the current batch configuration.

### run_batch_processing() (`specgui.batch.batch_main`)

Start batch processing.

### stop_batch_processing() (`specgui.batch.batch_main`)

Stop batch processing.

### show_about() (`specgui.batch.batch_main`)

Show about dialog.

### show_user_guide() (`specgui.batch.batch_main`)

Show user guide.

### closeEvent() (`specgui.batch.batch_main`)

Handle window close event.

### load_batch_configuration() (`specgui.batch.batch_controller`)

Load a batch configuration using hybrid format.

### _on_item_updated() (`specgui.batch.batch_controller`)

Handle item updates from master table.

### batch_items() (`specgui.batch.batch_controller`)

Get batch items in old format for backward compatibility.

### add_batch_item() (`specgui.batch.batch_controller`)

Add a batch item from dictionary (for backward compatibility).

### remove_batch_item() (`specgui.batch.batch_controller`)

Remove a batch item by row index.

### clear_all_items() (`specgui.batch.batch_controller`)

Clear all batch items.

### validate_all_items() (`specgui.batch.batch_controller`)

Validate all items in the master table.

### process_batch() (`specgui.batch.batch_controller`)

Process batch items using the master table.
        
        Parameters:
        -----------
        selected_row_indices : list or None
            If provided, only process these row indices. If None, process all items.
        
        Returns:
        --------
        success : bool
            True if batch processing completed successfully.
        results : list
            List of processing results in old format for backward compatibility.

### _process_single_item() (`specgui.batch.batch_controller`)

Process a single batch item using the master table.

### save_batch_configuration() (`specgui.batch.batch_controller`)

Save the current batch configuration using hybrid format.

### _recreate_rb_spec_objects() (`specgui.batch.batch_controller`)

Recreate rb_spec objects for items that have been processed.

### _recreate_single_rb_spec() (`specgui.batch.batch_controller`)

Recreate a single rb_spec object from saved parameters.

### export_template_csv() (`specgui.batch.batch_controller`)

Export just the template part as CSV.

### import_template_csv() (`specgui.batch.batch_controller`)

Import template from CSV.

### export_results_csv() (`specgui.batch.batch_controller`)

Export batch results to a CSV file with enhanced fields including velocity centroid and dispersion.

### export_error_log() (`specgui.batch.batch_controller`)

Export the error log to a file.

### get_item_count() (`specgui.batch.batch_controller`)

Get total number of items.

### get_valid_item_count() (`specgui.batch.batch_controller`)

Get number of valid items.

### get_completed_item_count() (`specgui.batch.batch_controller`)

Get number of completed items.

### get_error_item_count() (`specgui.batch.batch_controller`)

Get number of items with errors.

### edit_ew_range_dialog() (`specgui.batch.panels.ew_range_editor`)

Launch an enhanced EW range editor dialog.
    
    Parameters
    ----------
    item : object with template, analysis, results attributes
        The batch item to edit (now includes row_index)
    current_spec : rb_spec
        The current spectrum object
    controller : BatchController
        The batch controller for updates
    parent : QWidget
        Parent widget for the dialog
    
    Returns
    -------
    bool
        True if EW range was updated successfully, False otherwise

### _update_ew_range_in_master_table() (`specgui.batch.panels.ew_range_editor`)

Update EW range and recalculate measurements, ensuring master table consistency.

### _generate_spectrum_from_item() (`specgui.batch.panels.ew_range_editor`)

Generate rb_spec object from item parameters (local copy of the function).

### process_selected() (`specgui.batch.panels.processing_panel`)

Process only the selected batch items.

### process() (`specgui.batch.panels.processing_panel`)

Process batch items - kept for backward compatibility.

### cancel() (`specgui.batch.panels.processing_panel`)

Cancel batch processing.

### init_ui() (`specgui.batch.panels.processing_panel`)

Initialize the user interface.

### on_cont_method_changed() (`specgui.batch.panels.processing_panel`)

Handle changes to the continuum fitting method.

### update_process_options() (`specgui.batch.panels.processing_panel`)

Update the process options dropdown with current item counts.

### on_process_option_changed() (`specgui.batch.panels.processing_panel`)

Handle changes to the process option selection.

### set_selected_from_review() (`specgui.batch.panels.processing_panel`)

Set the items selected from the review panel.

### get_items_to_process() (`specgui.batch.panels.processing_panel`)

Get the list of row indices to process based on current selection.

### apply_settings_to_selected_items() (`specgui.batch.panels.processing_panel`)

Apply current UI settings to the selected items only.

### run_batch() (`specgui.batch.panels.processing_panel`)

Run batch processing on the selected items.

### _start_processing() (`specgui.batch.panels.processing_panel`)

Start the batch processing thread with selected items.

### _on_processing_finished() (`specgui.batch.panels.processing_panel`)

Handle processing finished signal.

### update_progress() (`specgui.batch.panels.processing_panel`)

Update the progress bar.

### on_item_completed() (`specgui.batch.panels.processing_panel`)

Handle item completed signal.

### on_item_failed() (`specgui.batch.panels.processing_panel`)

Handle item failed signal.

### cancel_processing() (`specgui.batch.panels.processing_panel`)

Cancel the current batch processing.

### init_ui() (`specgui.batch.panels.configuration_panel`)

Initialize the user interface.

### update_button_states() (`specgui.batch.panels.configuration_panel`)

Update button states based on table selection.

### edit_batch_item() (`specgui.batch.panels.configuration_panel`)

Edit the selected batch item.

### validate_ranges() (`specgui.batch.panels.configuration_panel`)

Validate that EW range is within slice range.

### import_multiple_json_files() (`specgui.batch.panels.configuration_panel`)

Import multiple rb_spec JSON files to create batch items.

### setup_table() (`specgui.batch.panels.configuration_panel`)

Setup table columns.

### add_batch_item() (`specgui.batch.panels.configuration_panel`)

Add an item to the batch processing list with smart defaults.

### _get_smart_defaults() (`specgui.batch.panels.configuration_panel`)

Get smart default values from existing table entries.

### refresh_table() (`specgui.batch.panels.configuration_panel`)

Refresh the table display from the master table.

### remove_batch_item() (`specgui.batch.panels.configuration_panel`)

Remove the selected item from the batch.

### clear_batch_items() (`specgui.batch.panels.configuration_panel`)

Clear all batch items.

### import_from_csv() (`specgui.batch.panels.configuration_panel`)

Import batch items from a CSV file.

### export_template_csv() (`specgui.batch.panels.configuration_panel`)

Export current configuration as CSV template.

### save_configuration() (`specgui.batch.panels.configuration_panel`)

Save the current batch configuration to a file.

### load_configuration() (`specgui.batch.panels.configuration_panel`)

Load a batch configuration from a file.

### create_csv_template() (`specgui.batch.panels.configuration_panel`)

Create an empty CSV template file.

### validate_configuration() (`specgui.batch.panels.configuration_panel`)

Validate the current batch configuration.

### get_item_count() (`specgui.batch.panels.configuration_panel`)

Get the current number of items.

### get_selected_item_ids() (`specgui.batch.panels.configuration_panel`)

Get IDs of currently selected items.

### create_individual_figures() (`specgui.batch.panels.batch_figure_generator`)

Create individual figure files for each batch item.
    
    Parameters
    ----------
    items : list
        List of batch items from master table
    output_directory : str
        Directory to save figures
    file_format : str
        Output format ('pdf' or 'png')
    master_table : MasterBatchTable
        Master table containing rb_spec objects
    
    Returns
    -------
    success_count : int
        Number of figures successfully created
    error_count : int
        Number of figures that failed to create

### create_multipage_pdf() (`specgui.batch.panels.batch_figure_generator`)

Create a multi-page PDF with 6 transitions per page.
    
    Parameters
    ----------
    items : list
        List of batch items from master table
    output_path : str
        Path for output PDF file
    master_table : MasterBatchTable
        Master table containing rb_spec objects
    
    Returns
    -------
    success : bool
        Whether the PDF was created successfully

### _get_spec_object() (`specgui.batch.panels.batch_figure_generator`)

Get rb_spec object - use existing one from master table.

### export_batch_figures() (`specgui.batch.panels.batch_figure_generator`)

Main function to export batch figures.
    
    Parameters
    ----------
    items : list
        List of batch items from master table
    output_directory : str
        Directory to save figures
    file_format : str
        Output format ('pdf' or 'png')
    figure_type : str
        Type of output ('individual' or 'multipage')
    master_table : MasterBatchTable
        Master table containing rb_spec objects
    
    Returns
    -------
    success : bool
        Whether the export was successful
    message : str
        Status message

### _get_clean_basename() (`specgui.batch.panels.batch_figure_generator`)

Get clean basename for filename generation.

### _create_single_figure() (`specgui.batch.panels.batch_figure_generator`)

Create a two-panel figure for a single transition.

### _create_multipage_figure() (`specgui.batch.panels.batch_figure_generator`)

Create a figure with up to 6 transitions (6 rows × 2 columns layout).

### _plot_flux_panel() (`specgui.batch.panels.batch_figure_generator`)

Plot flux panel with continuum and masked regions.

### _plot_normalized_panel() (`specgui.batch.panels.batch_figure_generator`)

Plot normalized flux panel with EW region and results.

### plot_spectrum_overview() (`specgui.batch.panels.spectrum_plotter`)

Create a comprehensive two-panel spectrum plot.
    
    Parameters
    ----------
    rb_spec : rb_spec object
        The spectrum object containing all data and analysis results
    item : object with template and analysis attributes
        Batch item containing metadata (filename, redshift, transition info, etc.)
    figure : matplotlib.figure.Figure, optional
        Existing figure to plot into. If None, creates new figure
    clear_figure : bool, optional
        Whether to clear the figure before plotting (default True)
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot

### _plot_flux_continuum_panel() (`specgui.batch.panels.spectrum_plotter`)

Plot the top panel: flux, error, continuum, and masked regions.

### _plot_normalized_panel() (`specgui.batch.panels.spectrum_plotter`)

Plot the bottom panel: normalized flux, EW region, and results.

### _add_results_text_box() (`specgui.batch.panels.spectrum_plotter`)

Add measurement results text box to the plot.

### plot_single_panel_normalized() (`specgui.batch.panels.spectrum_plotter`)

Create a single-panel normalized flux plot (for compact display).
    
    Parameters
    ----------
    rb_spec : rb_spec object
        The spectrum object containing all data and analysis results
    item : object with template and analysis attributes
        Batch item containing metadata
    figure : matplotlib.figure.Figure, optional
        Existing figure to plot into
    clear_figure : bool, optional
        Whether to clear the figure before plotting
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot

### plot_continuum_focus() (`specgui.batch.panels.spectrum_plotter`)

Create a continuum-focused plot showing flux + continuum fit details.
    
    This is useful for continuum editing workflows where you want to focus
    on the continuum fitting quality.

### init_ui() (`specgui.batch.panels.export_panel`)

Initialize the user interface.

### _update_defaults_on_change() (`specgui.batch.panels.export_panel`)

Update default directories when batch configuration changes - simplified.

### _is_csv_path_default() (`specgui.batch.panels.export_panel`)

Check if CSV path looks like our auto-generated name.

### _are_directories_default() (`specgui.batch.panels.export_panel`)

Check if directories haven't been manually changed by user.

### _update_csv_filename() (`specgui.batch.panels.export_panel`)

Update just the CSV filename, keeping the consistent directory.

### _update_default_directories() (`specgui.batch.panels.export_panel`)

Update the directory fields only.

### _set_default_directories() (`specgui.batch.panels.export_panel`)

Set default output directories and filenames - all use same directory.

### _get_default_directory() (`specgui.batch.panels.export_panel`)

Get the best default directory - simplified for consistency.

### _generate_smart_csv_filename() (`specgui.batch.panels.export_panel`)

Generate smart CSV filename based on batch content.

### _update_figure_type_options() (`specgui.batch.panels.export_panel`)

Update available figure type options based on selected format.

### browse_csv_path() (`specgui.batch.panels.export_panel`)

Browse for CSV file path.

### browse_json_dir() (`specgui.batch.panels.export_panel`)

Browse for JSON output directory.

### browse_figure_dir() (`specgui.batch.panels.export_panel`)

Browse for figure output directory.

### browse_error_path() (`specgui.batch.panels.export_panel`)

Browse for error log file path.

### export_csv() (`specgui.batch.panels.export_panel`)

Export results to a CSV file with enhanced fields.

### update_defaults_after_processing() (`specgui.batch.panels.export_panel`)

Update defaults after batch processing completes (use processing output directory).

### _export_enhanced_csv() (`specgui.batch.panels.export_panel`)

Export enhanced CSV with velocity centroid and dispersion.

### export_json() (`specgui.batch.panels.export_panel`)

Export individual JSON files for selected items.

### _export_individual_jsons() (`specgui.batch.panels.export_panel`)

Export individual JSON files for the given items.

### _resolve_filename_conflict() (`specgui.batch.panels.export_panel`)

Resolve filename conflicts by appending numbers.

### export_figures() (`specgui.batch.panels.export_panel`)

Export figures using the batch figure generator.

### export_error_log() (`specgui.batch.panels.export_panel`)

Export the error log to a file.

### init_ui() (`specgui.batch.panels.batch_selection_dialog`)

Initialize the user interface.

### populate_systems() (`specgui.batch.panels.batch_selection_dialog`)

Populate the systems table from master table.

### update_selection_counter() (`specgui.batch.panels.batch_selection_dialog`)

Update the selection counter label.

### select_failed_items() (`specgui.batch.panels.batch_selection_dialog`)

Select all items with error status.

### select_pending_items() (`specgui.batch.panels.batch_selection_dialog`)

Select all items ready for processing.

### select_complete_items() (`specgui.batch.panels.batch_selection_dialog`)

Select all completed items.

### select_all_items() (`specgui.batch.panels.batch_selection_dialog`)

Select all items.

### clear_selection() (`specgui.batch.panels.batch_selection_dialog`)

Clear all selections.

### invert_selection() (`specgui.batch.panels.batch_selection_dialog`)

Invert the current selection.

### preview_selection() (`specgui.batch.panels.batch_selection_dialog`)

Show a preview of the selected systems.

### get_selected_indices() (`specgui.batch.panels.batch_selection_dialog`)

Get the indices of selected systems.

### accept_selection() (`specgui.batch.panels.batch_selection_dialog`)

Accept the current selection.

### get_selection() (`specgui.batch.panels.batch_selection_dialog`)

Get the final selection after dialog closes.

### show_batch_selection_dialog() (`specgui.batch.panels.batch_selection_dialog`)

Convenience function to show the batch selection dialog.
    
    Returns:
        list: Selected system indices, or empty list if cancelled

### toggle_x_auto() (`specgui.batch.panels.review_panel`)

Toggle X-axis manual controls.

### toggle_y_upper_auto() (`specgui.batch.panels.review_panel`)

Toggle upper Y-axis manual controls.

### toggle_y_lower_auto() (`specgui.batch.panels.review_panel`)

Toggle lower Y-axis manual controls.

### toggle_snr_controls() (`specgui.batch.panels.review_panel`)

Toggle SNR bin size control.

### get_limits_and_settings() (`specgui.batch.panels.review_panel`)

Get the selected limits and analysis settings.

### init_ui() (`specgui.batch.panels.review_panel`)

Initialize UI with navigation controls and filtering.

### keyPressEvent() (`specgui.batch.panels.review_panel`)

Handle keyboard shortcuts.

### open_axis_limits_dialog() (`specgui.batch.panels.review_panel`)

Open dialog for setting custom axis limits and analysis settings.

### apply_custom_limits() (`specgui.batch.panels.review_panel`)

Apply custom axis limits to the current plot.

### toggle_interactive_ew() (`specgui.batch.panels.review_panel`)

Toggle interactive EW range selection mode.

### connect_ew_interactions() (`specgui.batch.panels.review_panel`)

Connect mouse and key events for interactive EW selection.

### disconnect_ew_interactions() (`specgui.batch.panels.review_panel`)

Disconnect interactive EW events.

### on_ew_plot_click() (`specgui.batch.panels.review_panel`)

Handle mouse clicks for EW range selection.

### update_ew_region_overlay() (`specgui.batch.panels.review_panel`)

Optimized update of just the EW region visualization.

### on_ew_plot_key() (`specgui.batch.panels.review_panel`)

Handle key presses in interactive EW mode.

### update_status_bar_info() (`specgui.batch.panels.review_panel`)

Send comprehensive info to main window status bar.

### on_filter_changed() (`specgui.batch.panels.review_panel`)

Handle filter radio button changes.

### get_filtered_items() (`specgui.batch.panels.review_panel`)

Get items based on current filter.

### get_overall_statistics() (`specgui.batch.panels.review_panel`)

Calculate overall statistics for all items.

### get_current_position_in_filter() (`specgui.batch.panels.review_panel`)

Get the current position within the filtered items.

### refresh_results_table() (`specgui.batch.panels.review_panel`)

Populate system dropdown from master table with current filter.

### on_system_selected() (`specgui.batch.panels.review_panel`)

Handle system selection from dropdown.

### update_navigation_state() (`specgui.batch.panels.review_panel`)

Update navigation button enabled states.

### navigate_to_first() (`specgui.batch.panels.review_panel`)

Navigate to first system in current filter.

### navigate_to_previous() (`specgui.batch.panels.review_panel`)

Navigate to previous system in current filter.

### navigate_to_next() (`specgui.batch.panels.review_panel`)

Navigate to next system in current filter.

### navigate_to_last() (`specgui.batch.panels.review_panel`)

Navigate to last system in current filter.

### open_batch_selection_dialog() (`specgui.batch.panels.review_panel`)

Open the batch selection dialog for choosing systems to process.

### display_item_spectrum() (`specgui.batch.panels.review_panel`)

Get existing spectrum from master table and copy it for editing.

### update_spectrum_plot() (`specgui.batch.panels.review_panel`)

Updated plotting method using the new plotting module.

### update_details() (`specgui.batch.panels.review_panel`)

Update details display with comprehensive system information.

### clear_preview() (`specgui.batch.panels.review_panel`)

Clear preview area.

### update_ew_range() (`specgui.batch.panels.review_panel`)

Update EW range using the same logic as the existing EW range editor.

### edit_continuum() (`specgui.batch.panels.review_panel`)

Launch continuum editor - using original working logic.

### apply_existing_continuum() (`specgui.batch.panels.review_panel`)

Apply the existing continuum from the spectrum file to the current item.

### _update_mark_failed_btn() (`specgui.batch.panels.review_panel`)

Update the Mark as Failed button label/style to reflect current item status.

### toggle_failed_status() (`specgui.batch.panels.review_panel`)

Toggle between manually_failed and needs_processing for the current item.

### grab_intervening_linelist() (`abstools.intervening_utils`)

Read intervening absorption line data from a file and filter based on wavelength range.
    
    Parameters:
    -----------
    filename : str
        Path to the file containing intervening line data
    z_gal : float
        Redshift of the host galaxy
    wrest_galaxy : float
        Rest wavelength of the line being analyzed
    wavelength : ndarray
        Array of wavelengths in the analysis window
        
    Returns:
    --------
    outlist : dict
        Dictionary containing filtered intervening line data

### plot_intervening_lines() (`abstools.intervening_utils`)

Plot intervening absorption lines on the given axes.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to plot on
    outlist : dict
        Dictionary containing intervening line data
    delv : float
        Velocity difference threshold
        
    Returns:
    --------
    None

### safe_draw() (`abstools.Metal_Plot`)

Safe drawing method with fallback

### connect_signals() (`abstools.Metal_Plot`)

Connect internal signals to their slots.

### show_error_message() (`abstools.Metal_Plot`)

Display an error message dialog.

### show_status_message() (`abstools.Metal_Plot`)

Display a message in the status bar.

### update_tab_display() (`abstools.Metal_Plot`)

Update display when tab changes.

### process_data_update() (`abstools.Metal_Plot`)

Process data update signals.

### setup_ui() (`abstools.Metal_Plot`)

Set up the user interface components.

### setup_event_connections() (`abstools.Metal_Plot`)

Set up event connections for mouse and keyboard input.

### getPage() (`abstools.Metal_Plot`)

Update page information when the tab changes.

### clean_exit() (`abstools.Metal_Plot`)

Properly clean up resources and exit the application.

### closeEvent() (`abstools.Metal_Plot`)

Handle application close event.

### NewTransition() (`abstools.Metal_Plot`)

Add a new transition to analyze.

### onsave() (`abstools.Metal_Plot`)

Open the save dialog.

### onload() (`abstools.Metal_Plot`)

Open the load dialog.

### opensub() (`abstools.Metal_Plot`)

Open the help window.

### onmotion() (`abstools.Metal_Plot`)

Handle mouse motion events.

### onpress() (`abstools.Metal_Plot`)

Handle key press events.

### onclick() (`abstools.Metal_Plot`)

Handle mouse click events.

### _setup_application() (`abstools.Metal_Plot`)

Set up the QApplication with proper error handling

### _run_application() (`abstools.Metal_Plot`)

Run the application event loop with error handling

### ensure_exit() (`abstools.Metal_Plot`)

Ensure the application fully exits without segmentation fault.

### on_motion() (`abstools.event_handler`)

Handle mouse motion events to track which axes is active and
        provide visual feedback.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all data and UI elements
        event : matplotlib.backend_bases.MouseEvent
            The mouse motion event
            
        Returns:
        --------
        None

### on_press() (`abstools.event_handler`)

Handle key press events for various functionality controls.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all data and UI elements
        event : matplotlib.backend_bases.KeyEvent
            The key press event
            
        Returns:
        --------
        None

### on_click() (`abstools.event_handler`)

Handle mouse click events for selecting velocity regions and limits.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all data and UI elements
        event : matplotlib.backend_bases.MouseEvent
            The mouse click event
            
        Returns:
        --------
        None

### shift2vel() (`abstools.config`)

Compute velocity shift given two redshifts.
    
    Parameters:
    -----------
    z1 : float
        Reference redshift (at rest)
    z2 : float
        Redshift for which relative velocity is computed
    rest_wavelength : float, optional
        Rest wavelength in Angstroms, default: 1215.67 (Lyman-alpha)
        
    Returns:
    --------
    vel : float
        Velocity in km/s

### plot_text() (`abstools.text_utils`)

Generate formatted text for displaying measurement results.
    
    Parameters:
    -----------
    parent : mainWindow instance
        The parent window containing all display settings
    line : dict
        Dictionary containing line measurement information
        
    Returns:
    --------
    None
        Updates the 'text' key in the line dictionary

### calculate() (`abstools.equivalent_width`)

Calculate equivalent width, column density, and other properties
        for an absorption line.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all the necessary data and UI elements
        page : int
            The page index in the tabbed interface
        ii : int
            Index of the ion to calculate properties for
        lims : list
            The velocity limits [vmin, vmax] for the integration
            
        Returns:
        --------
        bool
            True if calculation was successful, False otherwise

### calculate_all() (`abstools.equivalent_width`)

Calculate equivalent width for all transitions in all pages.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all the necessary data and UI elements
            
        Returns:
        --------
        tuple
            (successful_count, total_count) Number of successful calculations and total ions

### extract_measurements_table() (`abstools.extract_abstools_measurements`)

Extract measurements from an AbsTools JSON file into a pandas DataFrame
    
    Parameters:
    -----------
    json_file : str
        Path to the JSON file saved from AbsTools
        
    Returns:
    --------
    df : pandas.DataFrame
        DataFrame containing all measurements

### plot() (`abstools.plotting`)

Main plotting function for both continuum and normalized spectra.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window containing all the necessary data and UI elements
        ii : int
            Index of the ion to plot
        modify : bool, optional
            Whether to recalculate the continuum fit, default: False
        Print : bool, optional
            Whether to print the equivalent width text, default: False

### setup_ui() (`abstools.abstools_launcher`)

Set up the user interface components.

### browse_spectrum() (`abstools.abstools_launcher`)

Open a file dialog to select a spectrum file.

### browse_intervening() (`abstools.abstools_launcher`)

Open a file dialog to select an intervening absorbers file.

### browse_saved_intervening() (`abstools.abstools_launcher`)

Open a file dialog to select an intervening absorbers file for saved analysis.

### browse_analysis() (`abstools.abstools_launcher`)

Open a file dialog to select a saved analysis file.

### _sync_line_fields() (`abstools.abstools_launcher`)

Keep the hidden lines_edit and visible lines_display in sync.

### apply_preset() (`abstools.abstools_launcher`)

Apply the selected preset to the absorption lines field.

### update_status() (`abstools.abstools_launcher`)

Update the status text area with a new message.

### save_settings() (`abstools.abstools_launcher`)

Save current settings to QSettings.

### load_settings() (`abstools.abstools_launcher`)

Load saved settings from QSettings.

### clean_up_ui() (`abstools.abstools_launcher`)

Disconnect signals and prepare for clean exit.

### launch_abstools() (`abstools.abstools_launcher`)

Launch AbsTools with the configured settings.

### launch_metal_plot_separately() (`abstools.abstools_launcher`)

Launch Metal_Plot in a separate process to avoid event loop conflicts.
        
        Parameters:
        -----------
        data_file : str
            Path to the data file (either a temporary pickle file or a JSON file)
        is_json : bool
            Whether the file is a JSON file
        intervening_path : str, optional
            Path to intervening absorbers file

### _execute_script() (`abstools.abstools_launcher`)

Execute a Python script file in a separate process.

### run_abstools_new_spectrum() (`abstools.abstools_launcher`)

Run AbsTools with a new spectrum by creating an Absorber object and saving to a temporary file.
        Then launch Metal_Plot in a separate process.
        
        Parameters:
        -----------
        spectrum_path : str
            Path to the spectrum file
        format_name : str
            Format of the spectrum file
        redshift : float
            Redshift at which to perform analysis
        lines : list
            List of absorption line wavelengths to analyze
        window_limits : list
            List of [min, max] velocity window limits in km/s
        intervening_path : str, optional
            Path to intervening absorbers file

### run_abstools_load_analysis() (`abstools.abstools_launcher`)

Run AbsTools with a saved analysis file.
        Launch Metal_Plot in a separate process to avoid event loop conflicts.
        
        Parameters:
        -----------
        analysis_path : str
            Path to the saved analysis file
        intervening_path : str, optional
            Path to intervening absorbers file

### run_launcher() (`abstools.abstools_launcher`)

Run the AbsTools launcher application.

### run_as_separate_process() (`abstools.simplified_abstools_launcher`)

Launch AbsTools in a completely separate process.
    This function should be called from absorption_manager.

### run_gui() (`abstools.simplified_abstools_launcher`)

Run the AbsTools Launcher GUI.

### disconnect_matplotlib_events() (`abstools.cleanup`)

Safely disconnect all matplotlib event handlers.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window containing matplotlib connections

### clear_matplotlib_figures() (`abstools.cleanup`)

Safely clear all matplotlib figures.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window containing matplotlib figures

### disconnect_qt_signals() (`abstools.cleanup`)

Safely disconnect all Qt signals.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window containing Qt signals

### safe_exit() (`abstools.cleanup`)

Perform a safe exit of the application to prevent segmentation faults.
        
        Parameters:
        -----------
        window : MainWindow instance
            The main window to properly close

### load_abstools_json() (`abstools.batch_process_abstools`)

Load an AbsTools JSON file and return the ions dictionary

### batch_process_abstools_files() (`abstools.batch_process_abstools`)

Process multiple AbsTools JSON files matching the given pattern
    and compile their measurements
    
    Parameters:
    -----------
    json_pattern : str
        Glob pattern for JSON files, e.g., "data/*.json"
        
    Returns:
    --------
    all_measurements : pandas.DataFrame
        DataFrame containing measurements from all files

### plot_ew_comparison() (`abstools.batch_process_abstools`)

Create a comparison plot of EW measurements for selected ions
    across different spectra
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with measurements from batch_process_abstools_files
    ion_names : list, optional
        List of ion names to include in the plot, if None use all

### json_object_hook() (`abstools.json_utils`)

Custom object hook for JSON deserialization that handles special types.

### save_to_json() (`abstools.json_utils`)

Save the analysis data to a JSON file.
    
    Parameters:
    -----------
    data : dict
        The ions dictionary containing all analysis data
    filename : str
        Path to the output JSON file
        
    Returns:
    --------
    bool
        True if successful, False otherwise

### load_from_json() (`abstools.json_utils`)

Load analysis data from a JSON file.
    
    Parameters:
    -----------
    filename : str
        Path to the input JSON file
        
    Returns:
    --------
    dict or None
        The loaded ions dictionary, or None if loading failed

### add_json_support_to_save_page() (`abstools.json_utils`)

Code to add to SavePage class in ui_components.py to support JSON saving

### onjson() (`abstools.json_utils`)

Save analysis data as a JSON file

### load_from_json_example() (`abstools.json_utils`)

Example code to load analysis data from a JSON file in Metal_Plot.py

### setup_ui() (`abstools.ui_components`)

Set up the UI components for the save dialog

### onpdf() (`abstools.ui_components`)

Save plots as PDF files

### ontable() (`abstools.ui_components`)

Save measurement table to a file

### onpickle() (`abstools.ui_components`)

Save analysis data as a pickle file

### onjson() (`abstools.ui_components`)

Save analysis data as a JSON file

### add_page() (`abstools.ui_components`)

Add a new page (tab) to the parent window.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window to add a page to
            
        Returns:
        --------
        None

### initialize_page() (`abstools.ui_components`)

Initialize a new page with all necessary UI components.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window to initialize a page in
            
        Returns:
        --------
        None

### load_abstools_analysis() (`abstools.extract_abstools_data`)

Load an AbsTools JSON analysis file and return the data as a dictionary
    
    Parameters:
    -----------
    json_file : str
        Path to the JSON file saved from AbsTools
        
    Returns:
    --------
    ions : dict
        Dictionary containing all the analysis data

