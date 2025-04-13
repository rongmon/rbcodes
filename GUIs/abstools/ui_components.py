"""
UI components module for the absorption line analysis toolbox.
This module handles the creation and management of UI elements and pages.
Now with JSON support added.
The main modification is to the SavePage class to add JSON save functionality.
"""

import numpy as np
import pickle
import json
import webbrowser
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QPushButton, QLabel, QVBoxLayout, QHBoxLayout, QLineEdit, QInputDialog, QMessageBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from rbcodes.GUIs.abstools.config import HELP_TEXT, MAX_IONS_PER_TAB
from rbcodes.GUIs.abstools.plotting import Plotting

# Import the JSON utilities
try:
    from rbcodes.GUIs.abstools.json_utils import save_to_json, load_from_json
    JSON_SUPPORT = True
except ImportError:
    JSON_SUPPORT = False
    print("JSON utilities not found. JSON support will be disabled.")

class HelpWindow(QtWidgets.QWidget):
    """
    Window displaying help information for the application.
    """
    def __init__(self, parent=None):
        super(HelpWindow, self).__init__(parent)
        self.resize(400, 500)
        label = QtWidgets.QLabel(HELP_TEXT, self)

class SavePage(QtWidgets.QWidget):
    """
    Dialog for saving analysis results in various formats.
    """
    def __init__(self, parentvals, parent=None):
        super(SavePage, self).__init__(parent)
        self.resize(700, 450)  # Made taller to accommodate JSON controls
        self.closewin = False
        self.setup_ui(parentvals)
        
    def setup_ui(self, parentvals):
        """Set up the UI components for the save dialog"""
        from astropy.table import Table
        from astropy.io import ascii
        
        # Create a table of results for display and saving
        Table_e = Table()
        Table_e['Transitions'] = parentvals.keys
        EW = []
        EWsig = []
        N = []
        Nsig = []
        Vel = []
        EWlims_low = []
        EWlims_high = []
        
        # Check if all ions have been evaluated
        for ion in parentvals.keys:
            all_evaluated = True
            for items in parentvals.ions[ion]:
                if items in ['N', 'Nsig', 'EW', 'EWsig', 'med_vel'] and parentvals.ions[ion][items] is None:
                    all_evaluated = False
                    break
                    
            if not all_evaluated:
                # Show confirmation dialog for unevaluated ions
                reply = QtWidgets.QMessageBox.question(
                    self, 'Message', 
                    "Unevaluated ions, proceed to save?",
                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, 
                    QtWidgets.QMessageBox.No
                )
                
                if reply == QtWidgets.QMessageBox.Yes:
                    # Replace None values with NaN for saving
                    for ions in parentvals.keys:
                        for item in parentvals.ions[ions]:
                            if item in ['N', 'Nsig', 'EW', 'EWsig', 'med_vel'] and parentvals.ions[ions][item] is None:
                                parentvals.ions[ions][item] = float('nan')
                else:
                    self.closewin = None
                    return
                    
            # Add values to the table
            try:
                EW.append(round(parentvals.ions[ion]['EW'], 2) if parentvals.ions[ion]['EW'] is not None else float('nan'))
                EWsig.append(round(parentvals.ions[ion]['EWsig'], 2) if parentvals.ions[ion]['EWsig'] is not None else float('nan'))
                
                # Handle N values, checking for None and non-positive values
                if parentvals.ions[ion]['N'] is not None and parentvals.ions[ion]['N'] > 0:
                    N.append(round(np.log10(parentvals.ions[ion]['N']), 2))
                else:
                    N.append(float('nan'))
                    
                # Handle Nsig values
                if parentvals.ions[ion]['Nsig'] is not None and parentvals.ions[ion]['Nsig'] > 0:
                    Nsig.append(round(np.log10(parentvals.ions[ion]['Nsig']), 2))
                else:
                    Nsig.append(float('nan'))
                    
                Vel.append(round(parentvals.ions[ion]['med_vel'], 2) if parentvals.ions[ion]['med_vel'] is not None else float('nan'))
                
                # Handle velocity limits
                if parentvals.ions[ion]['EWlims'][0] is not None:
                    EWlims_low.append(round(parentvals.ions[ion]['EWlims'][0], 2))
                else:
                    EWlims_low.append(float('nan'))
                    
                if parentvals.ions[ion]['EWlims'][1] is not None:
                    EWlims_high.append(round(parentvals.ions[ion]['EWlims'][1], 2))
                else:
                    EWlims_high.append(float('nan'))
                    
            except Exception as e:
                print(f"Error processing ion {ion}: {e}")
                EW.append(float('nan'))
                EWsig.append(float('nan'))
                N.append(float('nan'))
                Nsig.append(float('nan'))
                Vel.append(float('nan'))
                EWlims_low.append(float('nan'))
                EWlims_high.append(float('nan'))
                
        # Add data to the table
        Table_e['EW'] = EW
        Table_e['EWsig'] = EWsig
        Table_e['Vmin'] = EWlims_low
        Table_e['Vmax'] = EWlims_high
        Table_e['N'] = N
        Table_e['Nsig'] = Nsig
        Table_e['Vel'] = Vel
        
        print(Table_e)
        
        # Create PDF save UI elements
        pdflabel = QLabel(r"Enter path and filename: (e.g. pathname\Ions.pdf)", self)
        pdflabel.setGeometry(100, 75, 400, 30)
        
        self.pdfline = QLineEdit(self)
        self.pdfline.setText(f'Spectrum_Analysis_z_{str(parentvals.z)}_Ions.pdf')
        self.pdfline.setGeometry(100, 100, 300, 30)
        
        self.pdfsave = QPushButton("Save PDF", self)
        self.pdfsave.setGeometry(410, 100, 200, 30)
        self.pdfsave.clicked.connect(lambda: self.onpdf(parentvals))
        
        # Create table save UI elements
        tablelabel = QLabel(r"Enter path and filename: (e.g. pathname\Table.dat)", self)
        tablelabel.setGeometry(100, 150, 400, 30)
        
        self.tableline = QLineEdit(self)
        self.tableline.setText(f"Spectrum_Analysis_z_{str(parentvals.z)}_Measurement_Table.dat")
        self.tableline.setGeometry(100, 175, 300, 30)
        
        self.tablesave = QPushButton(r"Save Table", self)
        self.tablesave.setGeometry(410, 175, 200, 30)
        self.tablesave.clicked.connect(lambda: self.ontable(parentvals, Table_e))
        
        # Create pickle save UI elements
        picklelabel = QLabel(r"Enter path and filename: (e.g. pathname\Table.p)", self)
        picklelabel.setGeometry(100, 225, 400, 30)
        
        self.pickleline = QLineEdit(self)
        self.pickleline.setText(f"Spectrum_Analysis_z_{str(parentvals.z)}.p")
        self.pickleline.setGeometry(100, 250, 300, 30)
        
        self.picklesave = QPushButton(r"Save Progress (Pickle)", self)
        self.picklesave.setGeometry(410, 250, 200, 30)
        self.picklesave.clicked.connect(lambda: self.onpickle(parentvals))
        
        # Create JSON save UI elements (New)
        if JSON_SUPPORT:
            jsonlabel = QLabel(r"Enter path and filename: (e.g. pathname\analysis.json)", self)
            jsonlabel.setGeometry(100, 300, 400, 30)
            
            self.jsonline = QLineEdit(self)
            self.jsonline.setText(f"Spectrum_Analysis_z_{str(parentvals.z)}.json")
            self.jsonline.setGeometry(100, 325, 300, 30)
            
            self.jsonsave = QPushButton(r"Save Progress (JSON)", self)
            self.jsonsave.setGeometry(410, 325, 200, 30)
            self.jsonsave.clicked.connect(lambda: self.onjson(parentvals))
            
            # Add note about JSON vs Pickle
            jsonnotelabel = QLabel(r"Note: JSON files are more portable but larger than Pickle files.", self)
            jsonnotelabel.setGeometry(100, 360, 500, 30)
            jsonnotelabel.setStyleSheet("font-style: italic; color: gray;")
    
    def onpdf(self, parentvals):
        """Save plots as PDF files"""
        try:
            # Need to eliminate highlighted axes
            if parentvals.old_axes:
                for pos in ['top', 'bottom', 'left', 'right']:
                    parentvals.old_axes.spines[pos].set_edgecolor('black')
                    parentvals.old_axes.spines[pos].set_linewidth(0.5)
                    
            figurefile = self.pdfline.text()
            i = 1
            for figures in parentvals.figs:
                figures.savefig(figurefile[:-4] + str(i) + '.pdf', bbox_inches='tight')
                i = i + 1
            self.pdfsave.setStyleSheet('background-color : green')
            parentvals.save = True
        except Exception as e:
            print(f"Error saving PDF: {e}")
            self.pdfsave.setStyleSheet('background-color : red')
            
    def ontable(self, parentvals, Table_e):
        """Save measurement table to a file"""
        try:
            from astropy.io import ascii
            
            savefile = self.tableline.text()
            ascii.write(Table_e, savefile, overwrite=True)
            self.tablesave.setStyleSheet("background-color : green")
            parentvals.save = True
        except Exception as e:
            print(f"Error saving table: {e}")
            self.tablesave.setStyleSheet("background-color : red")
            
    def onpickle(self, parentvals):
        """Save analysis data as a pickle file"""
        try:
            pfile = self.pickleline.text()
            with open(pfile, 'wb') as pklfile:
                pickle.dump(parentvals.ions, pklfile, protocol=pickle.HIGHEST_PROTOCOL)
            self.picklesave.setStyleSheet('background-color : green')
            parentvals.save = True
        except Exception as e:
            print(f"Error saving pickle file: {e}")
            self.picklesave.setStyleSheet('background-color : red')
            
    def onjson(self, parentvals):
        """Save analysis data as a JSON file"""
        if not JSON_SUPPORT:
            QMessageBox.warning(self, "JSON Support Not Available", 
                                "JSON utilities not found. Please make sure json_utils.py is in your path.")
            return
            
        try:
            json_file = self.jsonline.text()
            success = save_to_json(parentvals.ions, json_file)
            
            if success:
                self.jsonsave.setStyleSheet('background-color : green')
                parentvals.save = True
            else:
                self.jsonsave.setStyleSheet('background-color : red')
        except Exception as e:
            print(f"Error saving JSON file: {e}")
            self.jsonsave.setStyleSheet('background-color : red')


class PageManager:
    """
    Manager for handling page (tab) creation and initialization in the UI.
    """
    @staticmethod
    def add_page(parent):
        """
        Add a new page (tab) to the parent window.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window to add a page to
            
        Returns:
        --------
        None
        """
        # Check if we've already reached the maximum number of pages
        nall = len(parent.keys)
        
        if nall > 6 and len(parent.figs) < 2:
            parent.page = 1
            PageManager.initialize_page(parent)
            
            parent.cid4 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid5 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid6 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
            
        elif nall > 12 and len(parent.figs) < 3:
            parent.page = 2
            PageManager.initialize_page(parent)
            
            parent.cid7 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid8 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid9 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
            
        elif nall > 18 and len(parent.figs) < 4:
            parent.page = 3
            PageManager.initialize_page(parent)
            
            parent.cid10 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid11 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid12 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
            
        elif nall > 24 and len(parent.figs) < 5:
            parent.page = 4
            PageManager.initialize_page(parent)
            
            parent.cid13 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid14 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid15 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
            
    @staticmethod
    def initialize_page(parent):
        """
        Initialize a new page with all necessary UI components.
        
        Parameters:
        -----------
        parent : mainWindow instance
            The parent window to initialize a page in
            
        Returns:
        --------
        None
        """
        try:
            # Create tab, figure, and canvas
            parent.tabs.append(QtWidgets.QWidget())
            parent.addTab(parent.tabs[parent.page], parent.tab_names[parent.page])
            parent.figs.append(Figure())
            parent.canvas.append(FigureCanvasQTAgg(parent.figs[parent.page]))
            
            # Create UI controls
            add_ion_button = QPushButton("Add Ion", parent)
            add_ion_button.setGeometry(630, 30, 200, 30)
            add_ion_button.clicked.connect(lambda: parent.NewTransition(parent))
            
            open_button = QPushButton("Help", parent)
            open_button.setGeometry(830, 30, 200, 30)
            open_button.clicked.connect(lambda: parent.opensub(parent))
            
            save_button = QPushButton("Save", parent)
            save_button.setGeometry(430, 30, 200, 30)
            save_button.clicked.connect(lambda: parent.onsave(parent))
            
            page_label = QLabel(f"Page: {parent.page+1}/{len(parent.figs)}", parent)
            page_label.setStyleSheet("font: 16pt;color: white;background-color:QColor(53, 53, 53)")
            
            # Create layouts
            main_layout = QVBoxLayout()
            top_layout = QHBoxLayout()
            bot_layout = QHBoxLayout()
            
            spacer_item = QtWidgets.QSpacerItem(5, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
            
            # Add widgets to layouts
            top_layout.addItem(spacer_item)
            top_layout.addWidget(add_ion_button)
            top_layout.addWidget(save_button)
            top_layout.addWidget(open_button)
            top_layout.addItem(spacer_item)
            
            bot_layout.addItem(spacer_item)
            bot_layout.addWidget(page_label)
            bot_layout.addItem(spacer_item)
            
            main_layout.addLayout(top_layout, stretch=1)
            main_layout.addWidget(parent.canvas[parent.page], stretch=7)
            main_layout.addLayout(bot_layout, stretch=1)
            
            # Set the layout for the tab
            parent.tabs[parent.page].setLayout(main_layout)
            
            # Calculate number of ions on this page
            parent.nions.append(len(parent.keys) - 6 * parent.page)
            if parent.nions[parent.page] > MAX_IONS_PER_TAB:
                parent.nions[parent.page] = MAX_IONS_PER_TAB
                
            # Initialize left and right axes
            parent.axesL.append(list(range(MAX_IONS_PER_TAB)))
            parent.axesR.append(list(range(MAX_IONS_PER_TAB)))
            
            for ii in range(parent.nions[parent.page]):
                parent.axesL[parent.page][ii] = parent.figs[parent.page].add_subplot(MAX_IONS_PER_TAB, 2, 2 * ii + 1)
                parent.axesR[parent.page][ii] = parent.figs[parent.page].add_subplot(MAX_IONS_PER_TAB, 2, 2 * (ii + 1))
                parent.figs[parent.page].subplots_adjust(hspace=0.01)
                Plotting.plot(parent, ii, modify=True)
            
            # Set up focus for keyboard functionality
            parent.figs[parent.page].canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
            parent.figs[parent.page].canvas.setFocus()
            
        except Exception as e:
            print(f"Error initializing page: {e}")
            
        finally:
            # Force garbage collection to prevent memory leaks
            import gc
            gc.collect()
