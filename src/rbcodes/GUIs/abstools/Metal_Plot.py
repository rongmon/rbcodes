"""
Metal_Plot.py - Main module for the absorption line analysis toolbox.
Updated with improved signal-slot communication and error handling for greater stability.
"""
import matplotlib
matplotlib.use('Qt5Agg')  # Set this before importing any matplotlib modules
import sys
import os
import numpy as np
import pickle
import json
from matplotlib import rcParams
rcParams['lines.linewidth'] = .9
import webbrowser
from pkg_resources import resource_filename

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (
    QStyleFactory, QPushButton, QLineEdit, QMainWindow, 
    QInputDialog, QLabel, QMessageBox, QScrollBar, 
    QVBoxLayout, QHBoxLayout, QApplication, QFileDialog
)
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtCore import QTimer, pyqtSignal, QObject
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

#import rbcodes modules
from rbcodes.GUIs.abstools import Absorber
from rbcodes.GUIs.abstools.config import MAX_TABS, MAX_IONS_PER_TAB
from rbcodes.GUIs.abstools.plotting import Plotting
from rbcodes.GUIs.abstools.ui_components import HelpWindow, SavePage, PageManager
from rbcodes.GUIs.abstools.equivalent_width import EquivalentWidth
from rbcodes.GUIs.abstools.event_handler import EventHandler
from rbcodes.GUIs.abstools.cleanup import ResourceCleanup


# Try to import JSON utilities, but don't fail if they're not available
try:
    from rbcodes.GUIs.abstools.json_utils import load_from_json, save_to_json
    JSON_SUPPORT = True
except ImportError:
    JSON_SUPPORT = False
    print("JSON utilities not found. JSON loading support will be disabled.")


class MainWindowSignals(QObject):
    """
    Signal definitions for the MainWindow class to enable robust signal-slot communication.
    """
    data_updated = pyqtSignal(dict)            # General data update signal
    ion_selected = pyqtSignal(str, int)        # ion name, index
    continuum_fitted = pyqtSignal(object)      # continuum data
    ew_measured = pyqtSignal(str, float, float)  # ion name, EW value, EW error
    tab_changed = pyqtSignal(int)              # page index
    error_occurred = pyqtSignal(str)           # error message
    status_message = pyqtSignal(str)           # status message
    file_loaded = pyqtSignal(str)              # file path
    file_saved = pyqtSignal(str)               # file path


class AbsToolsCanvas(FigureCanvasQTAgg):
    """
    Custom Matplotlib canvas for AbsTools with improved error handling.
    """
    
    def __init__(self, parent=None, width=6, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)


        
        # Make canvas focusable
        self.fig.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.fig.canvas.setFocus()
        
        # Store parent for signal access
        self.parent = parent
        
    def safe_draw(self):
        """Safe drawing method with fallback"""
        try:
            self.draw_idle()
        except Exception as e:
            print(f"Error in draw_idle: {e}")
            try:
                self.draw()
            except Exception as e2:
                print(f"Critical drawing error: {e2}")
                if hasattr(self.parent, 'signals'):
                    self.parent.signals.error_occurred.emit(f"Failed to update display: {str(e2)}")


class MainWindow(QtWidgets.QTabWidget):
    """
    Main application window for the absorption line analysis toolbox.
    Improved with signal-slot communication for better stability.
    """
    
    def __init__(self, ions, parent=None, intervening=False):
        # Initialize signals
        self.signals = MainWindowSignals()
        
        # Check if this is a loaded file (complete ions dictionary with Target key)
        self.loaded_file = isinstance(ions, dict) and 'Target' in ions
        
        if self.loaded_file:
            # This is a pre-loaded complete ions dictionary
            self.ions = ions
            # Extract target properties from the loaded data
            self.z = ions['Target']['z']
            self.flux = ions['Target']['flux']
            self.wave = ions['Target']['wave']
            self.error = ions['Target']['error']
            # Get keys without the Target key
            self.keys = [key for key in ions.keys() if key != 'Target']
        else:
            # Extract full spectra properties from a new Absorber object
            self.z = ions['Target']['z']
            self.flux = ions['Target']['flux']
            self.wave = ions['Target']['wave']
            self.error = ions['Target']['error']
            
            # Initialize ions dictionary
            self.ions = ions
            self.keys = list(self.ions.keys())[:-1]  # Last item is the full target spectrum
        
        # Initialize parameters
        self.old_axes = None
        self.tab_names = ['Ions', 'Ions 2', 'Ions 3', 'Ions 4', 'Ions 5']
        self.wc = None  # Initialize overall parameter
        self.ylims = []
        self.vclim = None
        self.vclim_all = []
        self.name = None
        self.EWlim = [None, None]  # left, right
        self.event_button = None
        self.Manual_Mask = None
        self.pFlag = 1
        self.intervening = intervening
        self.save = False
        self.Lidx = None
        self.Ridx = None
        
        # Initialize parent class
        super(MainWindow, self).__init__(parent)
        
        # Initial page setup
        self.page = 0
        self.tabs = [QtWidgets.QWidget()]
        self.addTab(self.tabs[self.page], self.tab_names[self.page])
        # Create canvas first, then use its figure
        self.canvas = [AbsToolsCanvas(self, width=15, height=9, dpi=100)]
        self.figs = [self.canvas[0].fig]  # Use the figure from the canvas
                        
        self.nions = [len(self.keys)]
        if self.nions[0] > 6:
            self.nions[0] = 6
        
        # Set up focus for keyboard functionality
        self.setParent(parent)
        self.canvas[0].setFocusPolicy(QtCore.Qt.ClickFocus)  # Apply focus to canvas directly
        self.canvas[0].setFocus()
        
        # Set up UI layout
        self.setup_ui()
        
        # Initialize axes
        self.axesL = [list(range(MAX_IONS_PER_TAB))]
        self.axesR = [list(range(MAX_IONS_PER_TAB))]
        
        for ii in range(self.nions[0]):
            self.axesL[0][ii] = self.figs[0].add_subplot(MAX_IONS_PER_TAB, 2, 2 * ii + 1)
            self.axesR[0][ii] = self.figs[0].add_subplot(MAX_IONS_PER_TAB, 2, 2 * (ii + 1))
            self.figs[self.page].subplots_adjust(hspace=0.01)
            Plotting.plot(self, ii, modify=True)
        
        # Set up event connections
        self.setup_event_connections()
        
        # Connect signals
        self.connect_signals()
        
        # Set up additional pages if needed
        PageManager.add_page(self)
    
    def connect_signals(self):
        """Connect internal signals to their slots."""
        self.signals.error_occurred.connect(self.show_error_message)
        self.signals.status_message.connect(self.show_status_message)
        self.signals.tab_changed.connect(self.update_tab_display)
        self.signals.data_updated.connect(self.process_data_update)
        self.signals.file_loaded.connect(lambda path: self.statusBar().showMessage(f"Loaded file: {path}"))
        self.signals.file_saved.connect(lambda path: self.statusBar().showMessage(f"Saved file: {path}"))
    
    def show_error_message(self, message):
        """Display an error message dialog."""
        QMessageBox.warning(self, "Error", message)
    
    def show_status_message(self, message):
        """Display a message in the status bar."""
        if hasattr(self, 'status_bar'):
            self.status_bar.showMessage(message)
        else:
            print(message)  # Fallback
    
    def update_tab_display(self, page_index):
        """Update display when tab changes."""
        self.PageLabel.setText(f"Page: {page_index + 1}/{len(self.figs)}")
    
    def process_data_update(self, update_info):
        """Process data update signals."""
        action = update_info.get('action', '')
        
        if action == 'continuum_updated':
            ion_idx = update_info.get('ion_index')
            if ion_idx is not None and ion_idx < len(self.keys):
                Plotting.plot(self, ion_idx % 6, modify=False, Print=False)
        
        elif action == 'ew_measured':
            ion_name = update_info.get('ion_name')
            if ion_name in self.keys:
                ion_idx = self.keys.index(ion_name)
                Plotting.plot(self, ion_idx % 6, modify=False, Print=True)
    
    def setup_ui(self):
        """Set up the user interface components."""
        # Create control buttons
        self.status_bar = QtWidgets.QStatusBar()
        self.add_ion_button = QPushButton("Add Ion", self)
        self.add_ion_button.setGeometry(630, 30, 200, 30)
        self.add_ion_button.clicked.connect(lambda: self.NewTransition(self))
        
        self.openButton = QPushButton("Help", self)
        self.openButton.setGeometry(830, 30, 200, 30)
        self.openButton.clicked.connect(lambda: self.opensub(self))
        
        self.saveButton = QPushButton("Save", self)
        self.saveButton.setGeometry(430, 30, 200, 30)
        self.saveButton.clicked.connect(lambda: self.onsave(self))
        
        # Add load button
        self.loadButton = QPushButton("Load", self)
        self.loadButton.setGeometry(230, 30, 200, 30)
        self.loadButton.clicked.connect(lambda: self.onload(self))
        
        self.PageLabel = QtWidgets.QLabel("Page: " + str(self.currentIndex() + 1) + "/" + str(len(self.figs)), self)
        self.PageLabel.setStyleSheet("font: 16pt;color: white;background-color:QColor(53, 53, 53)")
        
            
        # Create layouts
        self.main_layout = QVBoxLayout()
        self.top_layout = QHBoxLayout()
        self.bot_layout = QHBoxLayout()
            
        self.spacerItem = QtWidgets.QSpacerItem(5, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        
        # Arrange widgets in layouts
        self.top_layout.addItem(self.spacerItem)
        self.top_layout.addWidget(self.loadButton)
        self.top_layout.addWidget(self.saveButton)
        self.top_layout.addWidget(self.add_ion_button)
        self.top_layout.addWidget(self.openButton)
        self.top_layout.addItem(self.spacerItem)
        
        self.bot_layout.addItem(self.spacerItem)
        self.bot_layout.addWidget(self.PageLabel)
        self.bot_layout.addItem(self.spacerItem)
        
        self.main_layout.addLayout(self.top_layout, stretch=1)
        self.main_layout.addWidget(self.canvas[0], stretch=8)
        self.main_layout.addLayout(self.bot_layout, stretch=1)
        # Add status bar to the bottom
        self.main_layout.addWidget(self.status_bar)
    
        self.tabs[self.page].setLayout(self.main_layout)
        
        self.tabs[self.page].setLayout(self.main_layout)
        
        # Set up page change handler
        self.currentChanged.connect(lambda: self.getPage())
    
    def setup_event_connections(self):
        """Set up event connections for mouse and keyboard input."""
        # Set up connectivity
        self.cid1 = self.figs[0].canvas.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.figs[0].canvas.mpl_connect("key_press_event", self.onpress)
        self.cid3 = self.figs[0].canvas.mpl_connect("motion_notify_event", self.onmotion)
    
    def getPage(self):
        """Update page information when the tab changes."""
        self.page = self.currentIndex()
        self.signals.tab_changed.emit(self.page)
    
    def clean_exit(self):
        """Properly clean up resources and exit the application."""
        try:
            print("Cleaning up resources before exit...")
            
            # Use the ResourceCleanup class for a thorough cleanup
            from cleanup import ResourceCleanup
            ResourceCleanup.safe_exit(self)
            
            # This line should never be reached due to os._exit in safe_exit
            import os
            os._exit(0)
            
        except Exception as e:
            print(f"Error during cleanup: {e}")
            # Force immediate exit even if cleanup fails
            import os
            os._exit(0)
        
    def closeEvent(self, event):
        """Handle application close event."""
        try:
            # Call our clean exit method
            self.clean_exit()
            
            # Accept the close event
            event.accept()
        except Exception as e:
            print(f"Error during closeEvent: {e}")
            self.signals.error_occurred.emit(f"Error during application exit: {str(e)}")
            event.accept()
    
    def NewTransition(self, parent):
        """Add a new transition to analyze."""
        try:
            # Get new transition wavelength from user
            new_line, ok = QInputDialog.getDouble(self, 'Add Line', 'Enter new transition:')
            
            if ok:
                # Create new absorber object for the line
                new_abs = Absorber.Absorber(self.z, self.wave, self.flux, self.error, [new_line])
                
                # Update the ions dictionary
                self.ions.update(new_abs.ions)
                self.ions['Target'] = self.ions.pop('Target')  # Moves Target to last index
                
                # Add new key to the keys list
                for key in list(new_abs.ions.keys())[:-1]:
                    self.keys.append(key)
                
                # Emit signal for data update
                self.signals.data_updated.emit({"action": "add_ion", "ion": key})
                
                # Determine whether to add to current page or create a new page
                if self.nions[self.page] < 6:
                    # Add to current page
                    ii = self.nions[self.page]
                    self.axesL[self.page][ii] = self.figs[self.page].add_subplot(MAX_IONS_PER_TAB, 2, 2 * ii + 1)
                    self.axesR[self.page][ii] = self.figs[self.page].add_subplot(MAX_IONS_PER_TAB, 2, 2 * (ii + 1))
                    self.figs[self.page].subplots_adjust(hspace=0.01)
                    self.nions[self.page] = self.nions[self.page] + 1
                    Plotting.plot(self, ii, modify=True)
                else:
                    # Need to add a new page
                    PageManager.add_page(self)
                
                self.signals.status_message.emit(f"Added new transition: {key}")
        
        except Exception as e:
            print(f"Error adding new transition: {e}")
            self.signals.error_occurred.emit(f"Failed to add transition: {str(e)}")
    
    def onsave(self, parent):
        """Open the save dialog."""
        try:
            self.savepg = SavePage(self)
            if self.savepg.closewin is None:
                return None
            else:
                self.savepg.show()
                self.signals.status_message.emit("Save dialog opened")
        except Exception as e:
            print(f"Error opening save dialog: {e}")
            self.signals.error_occurred.emit(f"Failed to open save dialog: {str(e)}")
    
    def onload(self, parent):
        """Open the load dialog."""
        try:
            # Check which format we support
            if JSON_SUPPORT:
                file_filter = "Analysis Files (*.p *.json);;Pickle Files (*.p);;JSON Files (*.json);;All Files (*)"
            else:
                file_filter = "Pickle Files (*.p);;All Files (*)"
                
            # Open file dialog
            file_path, _ = QFileDialog.getOpenFileName(
                self, "Load Analysis File", "", file_filter
            )
            
            if not file_path:
                return  # User cancelled
                
            if not os.path.exists(file_path):
                self.signals.error_occurred.emit(f"File not found: {file_path}")
                return
                
            # Load the file based on extension
            loaded_data = None
            
            if file_path.lower().endswith('.json'):
                if not JSON_SUPPORT:
                    self.signals.error_occurred.emit(
                        "JSON utilities not found. Please install the json_utils.py module."
                    )
                    return
                # Load from JSON
                loaded_data = load_from_json(file_path)
            
            elif file_path.lower().endswith('.p'):
                # Load from pickle
                try:
                    with open(file_path, 'rb') as f:
                        loaded_data = pickle.load(f)
                except Exception as e:
                    self.signals.error_occurred.emit(f"Failed to load pickle file: {str(e)}")
                    return
            else:
                self.signals.error_occurred.emit(
                    "Unrecognized file extension. Please use .p for pickle files or .json for JSON files."
                )
                return
                
            if loaded_data is None:
                self.signals.error_occurred.emit("Failed to load the analysis data.")
                return
                
            # Confirm before replacing current analysis
            reply = QMessageBox.question(
                self, "Confirm Load", 
                "Loading this file will replace your current analysis. Proceed?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            
            if reply == QMessageBox.Yes:
                # Close current window and open a new one with the loaded data
                self.closeEvent = lambda event: event.accept()  # Override closeEvent to prevent cleanup
                self.close()
                
                # Create new window with loaded data
                # Need to access the calling module to create a new Transitions instance
                from GUIs.abstools import Metal_Plot as M
                M.Transitions(loaded_data, intervening=self.intervening)
                
                self.signals.file_loaded.emit(file_path)
                
        except Exception as e:
            print(f"Error loading file: {e}")
            self.signals.error_occurred.emit(f"Failed to load file: {str(e)}")
    
    def opensub(self, parent):
        """Open the help window."""
        try:
            self.sub = HelpWindow()
            self.sub.show()
            self.signals.status_message.emit("Help window opened")
        except Exception as e:
            print(f"Error opening help window: {e}")
            self.signals.error_occurred.emit(f"Failed to open help window: {str(e)}")
    
    def onmotion(self, event):
        """Handle mouse motion events."""
        try:
            EventHandler.on_motion(self, event)
        except Exception as e:
            print(f"Error in motion handler: {e}")
            # We don't show errors for motion events to avoid overwhelming the user
    
    def onpress(self, event):
        """Handle key press events."""
        try:
            EventHandler.on_press(self, event)
        except Exception as e:
            print(f"Error in key press handler: {e}")
            self.signals.error_occurred.emit(f"Error processing keyboard command: {str(e)}")
    
    def onclick(self, event):
        """Handle mouse click events."""
        try:
            EventHandler.on_click(self, event)
        except Exception as e:
            print(f"Error in mouse click handler: {e}")
            self.signals.error_occurred.emit(f"Error processing mouse click: {str(e)}")
            self.vclim = None  # Reset in case of error
        

class Transitions:
    """
    Callable class to initialize and run the absorption line analysis application.
    Improved with better error handling and application lifecycle management.
    """
    def __init__(self, Abs, intervening=False):
        # Initialize QApplication if not already running
        app = self._setup_application()
        
        # Create main window
        main = MainWindow(Abs, intervening=intervening)
        main.resize(1400, 900)
        main.show()
        
        # Register cleanup handler
        app.aboutToQuit.connect(lambda: self.ensure_exit()) 
        
        # Connect signals for error handling
        if hasattr(main, 'signals'):
            main.signals.error_occurred.connect(
                lambda msg: QMessageBox.warning(main, "Error", msg)
            )
        
        # Start event loop only if we created a new app
        self._run_application(app)
    
    def _setup_application(self):
        """Set up the QApplication with proper error handling"""
        if not QtWidgets.QApplication.instance():
            app = QtWidgets.QApplication(sys.argv)
            app.setStyle("Fusion")
            
            # Set up dark theme
            palette = QPalette()
            palette.setColor(QPalette.Window, QColor(53, 53, 53))
            palette.setColor(QPalette.WindowText, QtCore.Qt.white)        
            palette.setColor(QPalette.Base, QColor(25, 25, 25))
            palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
            palette.setColor(QPalette.Button, QColor(53, 53, 53))
            palette.setColor(QPalette.ButtonText, QtCore.Qt.white)
            palette.setColor(QPalette.BrightText, QtCore.Qt.red)
            palette.setColor(QPalette.Link, QColor(42, 130, 218))
            palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
            palette.setColor(QPalette.Text, QtCore.Qt.white)
            
            app.setPalette(palette)
            return app
        else:
            return QtWidgets.QApplication.instance()
    
    def _run_application(self, app):
        """Run the application event loop with error handling"""
        # Ensure application quits properly when using 'q' key
        QtWidgets.QApplication.setQuitOnLastWindowClosed(True)
        
        # Check if this is a new app instance we created
        is_new_app = app == QtWidgets.QApplication.instance()
        
        if is_new_app:
            try:
                app.exec_()
                # Ensure we exit after app.exec_() returns
                sys.exit(0)
            except Exception as e:
                print(f"Error during application execution: {e}")
                sys.exit(1)
    
    def ensure_exit(self):
        """Ensure the application fully exits without segmentation fault."""
        try:
            print("Forcing immediate exit...")
            import os
            os._exit(0)  # Force immediate exit with no cleanup
        except:
            # This should never be reached, but just in case:
            import os
            os._exit(1)