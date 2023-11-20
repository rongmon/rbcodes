"""
Modules for PlotSpec_Integrated GUI
"""
import numpy as np
from scipy.interpolate import splrep,splev
from IGM import rb_setline as line   
from astropy.modeling import models, fitting    
import pdb
import sys
import os
from pathlib import Path
from utils import rb_utility as rt
import matplotlib as mpl
mpl.use('Qt5Agg')
mpl.rcParams['lines.linewidth'] = .9
clr=rt.rb_set_color()
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QComboBox,QHBoxLayout,QFileDialog,
    QLineEdit, QInputDialog,QListWidget, QVBoxLayout, QListWidgetItem,QLabel,QTableWidget,QGridLayout,QMessageBox,QBoxLayout,QDesktopWidget)

from PyQt5.QtGui import QPalette, QColor
from pkg_resources import resource_filename
import pandas as pd
import matplotlib.pyplot as plt 
from astropy.convolution import convolve, Box1DKernel
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from GUIs import guess_abs_line_vel_gui as g
from GUIs.abstools import Absorber as A 



HELP = '''
        MAIN GUI HELP:
        Left Widget is the absorber manager. Here you can add known absorbers, guess absorbers
        and plot or hide the transition lines. If an absorber is determined to be incorrect, removing 
        the abosrber will delete it from the manager and it will not exist in the output csv file once saving
        
        The main canvas is where all following keyboard events are tied to. If interacting with another 
        widget outside of the canvas, must reclick within the canvas to enable the keyboard events
        
        The active zabs manager (below canvas) will display what redshift is currently being used
        for identifying matching transitions. The Catalog button will automatically add this transition 
        to the absorber manager
        
        --------------Keyboard Events----------
        
        'r':   resets/clears the axes and replots the spectra
        'R':   Keeps the spectra active and will remove all lines/text from the canvas
        't':   Will restrict the ymax of the canvas to the users current mouse height
        'b':   Restricts ymin of the canvas to current mouse height
        'S':   Smoothes the spectra
        'U':   Unsmooth spectra
        'x':   Sets left x limit (xmin)
        'X':   Sets right x limit (xmax)
        'o':   Zoom out of x range
        ']':   Shifts canvas to the right
        '[':   Shifts canvas to the left
        'Y':   User can input their own y limits
        'E':   Two E keystrokes will compute rest frame equivalent width at a defined region
        'G':   Three keystrokes to fit a Gaussian profile. [Currently not drawing on the spectrum]

        'H':   Help Window
        'v':   Opens Separate Vstack GUI for the user to identify detected transitions
               Vstack commands will be discussed below : press shift+S TO EXITß
        
        'j':   Designed to be used by zooming into a small region of the spectra,
               finding an absortion region, put mouse in the middle of that region.
               'j' then opens a transition list window for the user to select which
               ion they believe the abosrber is located. Once a transition is clicked,
               the active zabs will display what the redshift should be based on the 
               mouse x location
               
        'Right Click' :  Same as j.

        -------Marking Doublets or Multiplets-------------

        'M':   MgII (2796,2803)
        'C':   CIV  (1548, 1550)
        'F':   FeII (2600,2586,2382)
        '6':   OVI  (1031, 1037)
        '4':   SiIV (1393, 1402)
        '8':   NeVIII (778, 770)
        '2':   Lyb/Lya
        '1':   Lya/Lyb
        
        --------Vstack GUI Keyboard Events--------
        Upon pressing 'v' in main canvas, a separate window will pop up
        
        '>':   Shifts page right (if more than one page)
        '<':   Shifts page left 
        'w':   Will change the transition flag between detection and Non-detection
        
        
        ------Notes-------
        Upon hitting 'v', the transition list will be saved. If redoing an analysis to correct
        or check transitions, it will identify user that the transition list has already been
        analyzed for the absorber. If continuing it will overwrite the previous results
        PLEASE NOTE, it will overwrite, not update the previous results.
        
        Save: Saving will generate its own file tree for files. Can rename a folder, but
            do not rename the files otherwise they will not load. Saving can be done without 
            evaluate linelists for each absorber or with partially evaluated line lists.
            Saving will access the current working directory, but can change the filepath to any 
            desired folder.
            
            The absorber manager will be saved as a .csv and the linelists saved as a .txt
            
        Load: To load, only give the parent folder, which has contents of the .csv and .txt files
              Loading will repopulate the previously confirmed absorber redshifts, colors, and linelist 
              used and will plot the detected absorbers as identified during the VStack GUI.
               
        

        '''

class mainWindow(QtWidgets.QMainWindow):#QtWidgets.QMainWindow
    
    def __init__(self,wave,flux,error,zabs=0,parent=None):
        
        self.c = 2.9979e5 #spl
        #initializing storing containers and variables to meet action requirements
        self.zabs_list = pd.DataFrame(data=None,columns = ['Zabs','list','color'])
        self.line_list = pd.DataFrame(data=None,columns = ['Name','Wave_obs','Zabs'])
        self.z_list = []
        #zabs_line_plot and text are where the mpl line objects and text objects are stored to be removed or hidden
        self.zabs_line_plot = []
        self.text = []
        self.linelist = []
        # identified_ lines and text are to store the evaluated detected absorber objects for plotting/hiding
        self.identified_lines = []
        self.identified_text = []
        #action variables to signify whether an event should be carried out
        self.hide = False
        self.row = None
        self.row_remove = False
        self.identified_line_active = False
        
        #spectrum properties
        self.wave=wave
        self.flux=flux
        self.smoothed_spectrum=flux
        self.error=error
        self.zabs=zabs
        self.lam_lim=[] #For compute_EW running tab
        self.lam_ylim=[]#For compute_EW running tab
        self.label='None' # Initializing a label
        
        #make longer color list
        clrlist=list(clr.keys())  

        self.combo_options =clrlist[1:]
        self.line_options = ['LLS','LLS Small','DLA','LBG','Gal','None']
        
        #---------------Initial page setup------------------# 
        super(mainWindow,self).__init__(parent)
        self.setWindowTitle('Absorber Idenification')
        
        #Main canvas widgets
        main_layout = QHBoxLayout()
        self.spectrum = Figure()
        self.ax = self.spectrum.add_subplot(111)
        self.ax.step(self.wave, self.error, '-',lw=0.5,color=clr['pale_red'],zorder=2)
        self.ax.step(self.wave, self.flux, '-',lw=0.5,color=clr['white'])
        self.init_xlims = [min(self.wave),max(self.wave)]
        self.canvas = FigureCanvasQTAgg(self.spectrum)
        toolbar = NavigationToolbar(self.canvas, self)
        
        #Zabs Manager Widgets
        self.abs_plot = manage_identified_absorbers(self)
        self.abs_plot.table.cellChanged.connect(self.cellchanged)
        
        self.manage_save = QPushButton("Save",self)
        self.manage_save.clicked.connect(lambda: self.saveCatalog_fn(self))
        
        self.manage_load = QPushButton("Load",self)
        self.manage_load.clicked.connect(lambda: self.LoadCatalog_fn(self))
        
        self.Identified_line_plot = QPushButton("Plot Identified Lines", self)
        self.Identified_line_plot.clicked.connect(lambda: Identified_plotter(self))
        
        #save layout (bottom of left panel)
        save_layout = QHBoxLayout()
        save_layout.addWidget(self.manage_save)
        save_layout.addWidget(self.manage_load)
        
        
        #left panel Main layout
        manage_layout = QVBoxLayout()
        manage_layout.addWidget(self.abs_plot.table)
        manage_layout.addLayout(save_layout)
        manage_layout.addWidget(self.Identified_line_plot)
        
        #canvas layout
        plot_layout = QtWidgets.QVBoxLayout()
        plot_layout.addWidget(toolbar,stretch=1)
        plot_layout.addWidget(self.canvas,stretch=5)
        self.plot_layout = plot_layout
        
        # active values widgets (bottom of main panel)
        self.active_zabs = QLineEdit(self)
        self.active_zabs.setText(str(self.zabs))
        self.active_zabs.textChanged[str].connect(lambda: self.zabs_changed())
        self.zabs_label = QLabel("Active Redshift",self)
        self.line_label = QLabel("Active LineList",self)
        self.combo_lines = QComboBox()
        for items in self.line_options:
            self.combo_lines.addItem(items)
        self.main_linelist = self.combo_lines.currentText()

        # Shows the active redshift and transition
        self.combo_color_main = QComboBox()
        for items in self.combo_options:
            self.combo_color_main.addItem(items)
        self.color = self.combo_color_main.currentText()
        color_label = QLabel('Color',self)
        
        #Message Window for user
        self.message_window = QLabel("Message Window")
        self.message_window.setStyleSheet('background-color : black')
        
        #active values layout (bottom of canvas layout)
        active_elem_layout = QtWidgets.QFormLayout()
        active_elem_layout.addRow(self.zabs_label,self.active_zabs)
        active_elem_layout.addRow(self.line_label,self.combo_lines)
        active_elem_layout.addRow(color_label,self.combo_color_main)
        
        #Catalog to connect Guessed line to the zabs manager:
        catalog = QPushButton("Catalog",self)
        catalog.clicked.connect(lambda: self.update_manager())#Catalog(self))
        plot = QPushButton("Plot",self)
        plot.clicked.connect(lambda: Redshift_Guess(self))
        refresh = QPushButton("Refresh",self)
        refresh.clicked.connect(lambda: self.Refreshed(self))
        
        #Spacer is to reduce the size of Active redshift and transition widgets
        spacer = QHBoxLayout()
        plot_cat = QVBoxLayout()
        spacerItem = QtWidgets.QSpacerItem(100, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        spacer.addWidget(self.message_window,stretch=1)
        spacer.addLayout(active_elem_layout,stretch=1)
        plot_cat.addWidget(plot)
        plot_cat.addWidget(catalog)
        plot_cat.addWidget(refresh)
        spacer.addLayout(plot_cat)
        spacer.addItem(spacerItem)
        spacer.addItem(spacerItem)
        spacer.addItem(spacerItem)
        plot_layout.addLayout(spacer,stretch=1)

        
        main_layout.addLayout(manage_layout,28)
        main_layout.addLayout(plot_layout,80)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QtWidgets.QWidget()
        widget.setLayout(main_layout)
        self.setCentralWidget(widget)
        self.show()
        #--------------------------end of layouts/widgets initialization------------#
        
        #configure canvas axes and viewing window
        self.ax.set_xlabel('Wavelength')
        self.ax.set_ylabel('Flux')
        xr=[min(self.wave),max(self.wave)]
        yr=[0.,np.median(flux)*2.5]
        self.ax.set_ylim(yr)
        self.ax.set_xlim(xr)
        self.init_ylims = self.ax.get_ylim()
        #---------------------------------------------------

        self.vel=np.array([1.])
        self.lam_lim=[]
        self.lam_ylim=[]
        self.FXval=[]
        self.FYval=[]


        
        #connect keyboard events to the canvas
        self.setParent(parent)
        self.spectrum.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.spectrum.canvas.setFocus()
        self.cid = self.spectrum.canvas.mpl_connect('key_press_event',self.ontype)
        self.cid_m = self.spectrum.canvas.mpl_connect('button_press_event', self.onclick)

        
        try:
            self.abs_plot.table.cellChanged.connect(self.cellchanged)
        except:
            pass

    def onclick(self, event):
        '''Mouse click
            Left == 1; Right == 3
        '''
        if event.button == 3:
            #Manual mode
            self.xdata = event.xdata
            self.manT = Manual_Transition(self)
            self.manT.show()

        
    #keyboard events
    def ontype(self,event):
        zabs=np.double(0.)
        # when the user hits 'r': clear the axes (except flux and error) and return to the initial x&y lims
        if event.key=='r':
            #del self.ax.lines[2:]
            self.lam_lim=[]
            self.lam_ylim=[]
            self.FXval=[]
            self.FYval=[]

            #self.ax.texts = []
            #while self.ax.texts:
            #    self.ax.texts.pop()
            #while self.ax.collections:
            #    self.ax.collections.pop()

            self.ax.set_ylim(self.init_ylims)
            self.ax.set_xlim(self.init_xlims)
            self.spectrum.canvas.draw()
            
        #another refresh to keep the current flux values but remove the plotted lines
        elif event.key == 'R':
            del self.ax.lines[2:]
            self.lam_lim=[]
            self.lam_ylim=[]
            self.FXval=[]
            self.FYval=[]
            try:
                for ii in self.text[-1]:
                    ii.remove()
            except: 
                pass
            #self.ax.texts = []
            while self.ax.texts:
                self.ax.texts.pop()
            while self.ax.collections:
                self.ax.collections.pop()


            # Give initial axes limits
#             self.ax.set_ylim(self.init_ylims)
#             self.ax.set_xlim(self.init_xlims)
            
            if self.identified_line_active == True:
                self.identified_line_active = False
                self.Identified_line_plot.setStyleSheet('background-color : QColor(53, 53, 53)')
                self.message_window.setText(" ")
            self.spectrum.canvas.draw()
#         # Set top y max
        elif event.key=='t':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_ylim([ylim[0],event.ydata])
            self.ax.set_xlim(xlim)
            self.spectrum.canvas.draw() 
#         # Set top y min
        elif event.key=='b':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_ylim([event.ydata,ylim[1]])
            self.ax.set_xlim(xlim)
            self.spectrum.canvas.draw() 
#         # Smooth spectrum
        elif event.key=='S':
            self.vel[0] += 2
            Filter_size=np.int(self.vel[0]) 
            self.smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))
            self.smoothed_error =convolve(self.error, Box1DKernel(Filter_size))
            self.specplot()
            self.spectrum.canvas.draw()  
#         #Unsmooth Spectrum
        elif event.key=='U':
            self.vel[0] -= 2
            if self.vel[0] <= 0:
                self.vel[0]=1;
            Filter_size=np.int(self.vel[0]) 
            self.smoothed_spectrum =convolve(self.flux, Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            self.smoothed_error =convolve(self.error, Box1DKernel(Filter_size))
            self.specplot()
    
        # Set X max
        elif event.key=='X':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_xlim([xlim[0],event.xdata])
            self.ax.set_ylim(ylim)
            self.spectrum.canvas.draw() 
        # Set x min
        elif event.key=='x':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            self.ax.set_xlim([event.xdata,xlim[1]])
            self.ax.set_ylim(ylim)
            self.spectrum.canvas.draw() 

        # Set pan spectrum
        elif event.key==']':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            delx=(xlim[1]-xlim[0])
            self.ax.set_xlim([xlim[1],xlim[1]+delx])
            self.ax.set_ylim(ylim)
            self.spectrum.canvas.draw() 
        # Set pan spectrum
        elif event.key=='[':
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            delx=(xlim[1]-xlim[0])
            self.ax.set_xlim([xlim[0]-delx,xlim[0]])
            self.ax.set_ylim(ylim)
            self.spectrum.canvas.draw() 
        #zoom out of xrange
        elif (event.key == 'o'):
            xlim=self.ax.get_xlim()
            ylim=self.ax.get_ylim()
            xcen = (xlim[0]+xlim[1])/2.0
            delx   = xlim[1] - xcen
            self.ax.set_xlim([xcen - 1.5*delx,xcen + 1.5*delx])
            self.spectrum.canvas.draw() 

        #guess line CIV
        elif (event.key == 'C'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'CIV',yval)
        #guess line MgII
        elif (event.key == 'M'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'MgII',yval)
        #guess line FeII
        elif (event.key == 'F'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'FeII',yval)

        #guess line OVI
        elif (event.key == '6'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'OVI',yval)

        #guess line SiIV
        elif (event.key == '4'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'SiIV',yval)

        #guess line NeVIII
        elif (event.key == '8'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'NeVIII',yval)

        #guess line Lyb
        elif (event.key == '2'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'Lyb',yval)

        #guess line Lya
        elif (event.key == '1'):
            wave0 = event.xdata
            yval=event.ydata        
            self.check_lineid(wave0,'Lya',yval)


        elif event.key == 'Y':
            Windowname='Manual y-Limits'
            instruction='Input range (e.g. 0.,2.)'
            ylim, ok = QInputDialog.getText(self,Windowname,instruction)
            if ok:
                ylimit = ylim.split(',')
                ylimit = np.array(ylimit).astype('float32')
                self.ax.set_ylim(ylimit)
                self.spectrum.canvas.draw()
        elif ((event.key == 'h') or (event.key =='H')):
            self.help = HelpWindow()
            self.help.show()
            
        elif event.key =='v':
            # first check to see if we are in the main window or the vstack window
            if not isinstance(self, vStack):
                # check if absorber linelist has already been catologed.
                if self.zabs in self.line_list.Zabs.tolist():
                    
                    #if so, ask user if they would like to re-eval the results
                    buttonReply = QMessageBox.question(self,"Reevaluate" ,"Current Zabs LineList already evaluated: Reevaluate and overwrite?",
                                                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                    if buttonReply == QMessageBox.Yes:
                        self.ion_selection = vStack(self,self.wave,self.flux,self.error,self.label,zabs=self.zabs)
                        
                #otherwise, proceed without manual consent
                else:
                    self.ion_selection = vStack(self,self.wave,self.flux,self.error,self.label,zabs=self.zabs)
                


        elif event.key =='j':
            self.xdata = event.xdata
            self.manT = Manual_Transition(self)
            self.manT.show()
            

        # Compute Equivalent Width between two points
        elif event.key=='E':
            #ekeycounts +=1
            self.lam_lim=np.append(self.lam_lim,event.xdata)
            self.lam_ylim=np.append(self.lam_ylim,event.ydata)

    
            # Keep running tab of all E clicks
            eclick=len(self.lam_lim);


            self.ax.plot(event.xdata,event.ydata,'rs',ms=5,label='EW_pt',markeredgecolor='k')
            #self.fig.canvas.draw()

            if eclick==2:
                # Check if the wave entries are monotonously increasing
                tab=self.lam_lim.argsort()
                
                #calculate and display column density if EW measurement intersects a line
                if len(self.ax.lines) > 2:
                    fval = []
                    name = []
                    wrest = []
                    try:
                        for item in self.data:
                            if ((item['wrest']*(1+self.zabs)>self.lam_lim[tab][0]) and (item['wrest']*(1+self.zabs)<self.lam_lim[tab][1])):
                                fval.append(item['fval'])
                                name.append(item['ion'])
                                wrest.append(item['wrest'])
                    except:
                        pass
                    #if only one line is intersected, use that transitions values
                    if (len(fval)>0 and len(fval)<2):
                        vel_lims = (np.asarray(self.lam_lim[tab])-wrest[0]*(1+self.zabs))*self.c/(wrest[0]*(1+self.zabs))
                        EW,sig_EW,cont,wave_slice,N,Nerr = self.compute_EW(self.wave/(1.+self.zabs),self.flux,self.lam_lim[tab]/(1+self.zabs),self.lam_ylim[tab],self.error,fval[0],wrest[0])
                        EW=np.array(EW)*1000.
                        sig_EW=np.array(sig_EW)*1000.
                        self.ax.plot(wave_slice,cont,'r--')
                        Wval = 'EW [mAA]: '+ '%.1f' % EW + ' +/- ' + '%.1f' % sig_EW + '\n\n'
                        logNerr = str(np.round(np.log10((N+Nerr)/(N-Nerr)/2),2))
                        formatName = '$N_{'+name[0]+'}$'
                        logN = Wval + 'log ' + formatName + ' [$cm^{-2}$]:'  + str(np.round(np.log10(N),2)) + ' +/- ' +logNerr
                        self.ax.text(np.mean([self.lam_lim])+.05,np.max(self.lam_ylim)+0.2,logN, rotation=90,verticalalignment='bottom')
                        self.message_window.setText(Wval + 'N '+name[0] + '[1/cm^2]'+': '+str(np.round(np.log10(N),2)) + ' +/- ' +logNerr)
                        self.lam_lim=[]
                        self.lam_ylim=[]
                        
                    #otherwise, ask user which transition they would like to use
                    elif len(fval)>1:
                        self.tab = tab
                        self.fvals = fval
                        self.fval_wrest = wrest
                        self.fval_name = name
                        self.popup = Popup_col_density(self)
                        self.popup.show()
                        
                    # if no lines are intersected only use EW
                    else:
                        EW,sig_EW,cont,wave_slice=self.compute_EW(self.wave/(1.+self.zabs),self.flux,self.lam_lim[tab]/(1+self.zabs),self.lam_ylim[tab],self.error)
                        EW=np.array(EW)*1000.
                        sig_EW=np.array(sig_EW)*1000.
                        self.ax.plot(wave_slice,cont,'r--')

                        Wval='EW [mAA]: '+ '%.1f' % EW + ' +/- ' + '%.1f' % sig_EW
                        self.ax.text(np.mean([self.lam_lim]),np.max(self.lam_ylim)+0.2,Wval, rotation=90,verticalalignment='bottom')
                        self.message_window.setText(Wval)
                        self.lam_lim=[]
                        self.lam_ylim=[]
                        
                #if no transitions plotted only evaluate EW 
                else:
                    EW,sig_EW,cont,wave_slice=self.compute_EW(self.wave/(1.+self.zabs),self.flux,self.lam_lim[tab]/(1+self.zabs),self.lam_ylim[tab],self.error)
                    EW=np.array(EW)*1000.
                    sig_EW=np.array(sig_EW)*1000.
                    self.ax.plot(wave_slice,cont,'r--')

                    Wval='EW [mAA]: '+ '%.1f' % EW + ' +/- ' + '%.1f' % sig_EW
                    self.ax.text(np.mean([self.lam_lim]),np.max(self.lam_ylim)+0.2,Wval, rotation=90,verticalalignment='bottom')
                    self.message_window.setText(Wval)
                    self.lam_lim=[]
                    self.lam_ylim=[]

                
            self.spectrum.canvas.draw() 

         # Fit a Gaussian
        elif event.key=='G':

            self.FXval=np.append(self.FXval, event.xdata)
            self.FYval=np.append(self.FYval, event.ydata)
    
            fclick=len(self.FXval)    
            self.ax.plot(event.xdata,event.ydata,'rs',ms=5,label='EW_pt',markeredgecolor='k')
            self.spectrum.canvas.draw() 

            

            #Start Fitting
            if fclick==3:
                # Fit the data using a Gaussian
                g_init = models.Gaussian1D(amplitude=np.double(self.FYval[1]), mean=np.double(self.FXval[1]), stddev=0.5*(np.double(self.FXval[2])-np.double(self.FXval[0])))
                fit_g = fitting.LevMarLSQFitter()

                # First fit a quick continuum
                qtq=np.where( ( (self.wave >= self.FXval[0] ) & ( (self.wave <= self.FXval[2] ) )))
                ww=self.wave[qtq]                 
                flux1=self.flux[qtq]
                spline = splrep(np.append(self.FXval[0],self.FXval[2]),np.append(self.FYval[0],self.FYval[2]),k=1)
                continuum = splev(ww,spline)
    
                # Check if it is an absorption or emission line
                if ((self.FYval[1] < self.FYval[0]) & (self.FYval[1] < self.FYval[2])):
                    ydata=1.- (flux1/continuum)
                    g = fit_g(g_init, ww, ydata)
                    Final_fit=(1.-g(ww))*continuum
                else:
                    ydata=(flux1/continuum)-1
                    g = fit_g(g_init, ww, ydata)
                    Final_fit=(1.+g(ww))*continuum         
    
                model_fit=self.ax.plot(ww, Final_fit, 'r-')        
                values0=' Amp: '+'%.3f' %  g.parameters[0] 
                values1=' Center: '+'%.3f' %  g.parameters[1] 
                values2=' Sig: '+'%.3f' %  g.parameters[2] 


                self.message1=self.ax.text(np.mean(ww),np.max(self.FYval)+0.2,values0, rotation=90,verticalalignment='bottom')
                self.message2=self.ax.text(np.mean(ww)+np.std(ww),np.max(self.FYval)+0.2,values1, rotation=90,verticalalignment='bottom')
                self.message3=self.ax.text(np.mean(ww)-np.std(ww),np.max(self.FYval)+0.2,values2, rotation=90,verticalalignment='bottom')



                self.message_window.setText(values0+'\n'+values1+'\n'+values2+'\n\n')

                self.FXval=[]
                self.FYval=[]

            self.spectrum.canvas.draw() 

        elif event.key=='q' or event.key=='Q':
            # sys.exit()
            # self.manT.close()
            # app.exit(0)
            # QtWidgets.QApplication.instance().exit(0)
            # self.close()
            # self.destroy() # ha, kills the window but doesn't exit the application :(
            # self.closeEvent(self) # nope
            
            # app = QtWidgets.QApplication.instance()
            # sys.exit(app.exec_())

            # the only way I have found to quit the application without segfaulting
            os._exit(0) # though not advised as os._exit() DOESN'T call cleanup functions, flush of stdio buffers, etc.
            # this PyQT crashing on exit has been around since 2004 (at least): https://bugzilla.redhat.com/show_bug.cgi?id=133375


    
    def compute_EW(self,lam,flx,lam_lim,lam_ylim,err_flx,fval = -99,wrest=0):
        qtq=np.where( ( lam >= lam_lim[0] ) & ( lam <= lam_lim[1] ) )
        ww=lam[qtq]  

        flux1=flx[qtq]
        spline = splrep(lam_lim,lam_ylim,k=1)
        continuum = splev(ww,spline)
        #ax.plot(ww,continuum,'k--')
        sig_flx1=err_flx[qtq]/continuum
        flux1=flux1/continuum
        
        #Removing any NAN values 
        sq=np.isnan(flux1);
        tmp_flx=sig_flx1[sq]
        flux1[sq]=tmp_flx

        #clip the spectrum. If the flux is less than 0+N*sigma, then we're saturated. Clip the flux array(to avoid inifinite optical depth) and set the saturated flag
        sat_lim=-.01
        q=np.where(flux1<=sat_lim);
        tmp_flx=sig_flx1[q]
        flux1[q]=tmp_flx
        q=np.where(flux1<=0.);
        tmp_flx=sig_flx1[q]+0.01
        flux1[q]=tmp_flx;
        
        
        EW=np.trapz(1.-flux1, x= ww)
        delw=np.double(ww[2]-ww[1])
        sig_w=delw*sig_flx1
        sig_wtot=np.sqrt(np.sum(sig_w**2.))
        
        
        if fval != -99:
            lambda_r=(lam/(1.+self.zabs))[qtq]
            vel = ((lam-wrest*(1.0 + self.zabs))*self.c/(wrest*(1.0 + self.zabs)))[qtq]
            f0=fval
        #compute apparent optical depth
            Tau_a =np.log(1./flux1);
            #compute the median optical depth weighted velcity.
            Tau50=np.cumsum(Tau_a)/np.max(Tau_a)
            vel50=np.interp(0.5,Tau50,vel)
            del_vel_j=np.diff(vel);
            del_vel_j=np.append([del_vel_j[0]],del_vel_j)
            # Column density per pixel as a function of velocity
            nv = Tau_a/((2.654e-15)*f0*lambda_r);# in units cm^-2 / (km s^-1), SS91 
            n = nv* del_vel_j# column density per bin obtained by multiplying differential Nv by bin width 
            tauerr = err_flx[qtq]/flux1;
            nerr = (tauerr/((2.654e-15)*f0*lambda_r))*del_vel_j; 
            col = np.sum(n);
            colerr = np.sum((nerr)**2.)**0.5;
            ww = ww*(1+self.zabs)
            return EW,sig_wtot,continuum,ww,col,colerr
        ww = ww*(1+self.zabs)
        return EW,sig_wtot,continuum,ww



            
    #Update manager brings items from the active zabs layout (below canvas) to an entry in the zabs_manager
    def update_manager(self):
        self.color = self.combo_color_main.currentText()
        try: linelist = self.manT.combo_ll.currentText()
        except: linelist =self.combo_lines.currentText() 
        new_row = pd.Series(data = {'Zabs': self.zabs, 'list': linelist, 'color': self.color})
        self.zabs_list=self.zabs_list.append(new_row,ignore_index=True)
        
        self.zabs_line_plot.append(self.temp_plots)
        self.text.append(self.tt_temp)
        
        #if you're filling a final slot within the zabs manager, create another slot
        if self.abs_plot.table.rowCount() == self.zabs_list.shape[0]:
            iterations = self.abs_plot.table.rowCount()
            for ii in range(iterations):
                if self.abs_plot.table.item(ii,1) == None:
                    #set z
                    self.abs_plot.table.setItem(ii,1,QtWidgets.QTableWidgetItem(str(np.round(self.zabs,4))))
                    #set linelist
                    self.abs_plot.table.cellWidget(ii,0).setCurrentIndex(self.combo_lines.currentIndex())
                    #set color
                    self.abs_plot.table.cellWidget(ii,2).setCurrentIndex(self.combo_color_main.currentIndex())
            Catalog(self)
            try: self.manT.close()
            except: pass
#             self.active_ion.setText(' ')
            self.spectrum.canvas.draw()

        #else populate the first empty slot without creating a new
        else:
            iterations = self.abs_plot.table.rowCount()
            for ii in range(iterations):
                if self.abs_plot.table.item(ii,1) == None:
                    self.abs_plot.table.setItem(ii,1,QtWidgets.QTableWidgetItem(str(np.round(self.zabs,4))))
                    #set linelist
                    self.abs_plot.table.cellWidget(ii,0).setCurrentIndex(self.combo_lines.currentIndex())
                    #set color
                    self.abs_plot.table.cellWidget(ii,2).setCurrentIndex(self.combo_color_main.currentIndex())
                    self.abs_plot.table.resizeColumnsToContents()
                    try: self.manT.close()
                    except: pass
#                     self.active_ion.setText(' ')

                    self.spectrum.canvas.draw()
                    break
                    
    def Refreshed(self,parent):
        self.ax.set_ylim(self.init_ylims)
        self.ax.set_xlim(self.init_xlims)
        self.lam_lim=[]
        self.lam_ylim=[]
        self.FXval=[]
        self.FYval=[]

        
        if self.identified_line_active == True:
            self.identified_line_active = False
            self.Identified_line_plot.setStyleSheet('background-color : QColor(53, 53, 53)')
            self.message_window.setText(" ")
        if len(self.ax.lines)>2:
            del self.ax.lines[2:]
            #self.ax.texts = []
            while self.ax.texts:
                self.ax.texts.pop()
            while self.ax.collections:
                self.ax.collections.pop()

            try:
                for ii in self.text[-1]:
                    ii.remove()
            except:
                pass

            #may need condition to change hide color to gray
            self.spectrum.canvas.draw()
        
    
    def cellchanged(self):
        col = self.abs_plot.table.currentColumn()
        row = self.abs_plot.table.currentRow()
        try: text = self.abs_plot.table.currentItem().text()
        except: pass
        try: self.abs_plot.table.cellWidget(row,col).setText(text)
        except: pass
        try: self.active_zabs.setText(text); self.zabs = np.double(text)
        except: pass
    def zabs_changed(self):
        try: self.zabs = np.double(self.active_zabs.text())
        except: self.active_zabs.setText("Please input numerical redshift")
            
    def saveCatalog_fn(self,parent):
        self.saving = SaveCatalog(parent)
        if self.saving.continueSaving == True:
            self.saving.show()
        else:
            self.saving.close()
            
    def LoadCatalog_fn(self,parent):
        self.loading = LoadCatalog(parent)
        self.loading.show()

    #This code will quickly draw some doublet lines on the canvas for a quicklook
    def check_lineid(self, wave0, ionname,yval):
        if (ionname == 'CIV'):
            wave1 = wave0 * 1550.77845 / 1548.2049           
            z = wave0 / 1548.2049 - 1
            print(f"CIV: z = {z}")
            self.message_window.setText(f"CIV: z = {z}")

        elif (ionname == 'MgII'):
            wave1 = wave0 * 2803.5314853 / 2796.3542699
            z = wave0 / 2796.354 - 1
            print(f"MgII: z = {z}")
            self.message_window.setText(f"MgII: z = {z}")

        elif (ionname == 'FeII'):
            wave1 = wave0 * 2586.6495659 / 2600.1724835 
            wave2=wave0 *2382.7641781/ 2600.1724835 
            z = wave0 / 2600.1724835 - 1
            print(f"FeII: z = {z}")
            self.message_window.setText(f"FeII: z = {z}")

        elif (ionname == 'OVI'):
            wave1 = wave0 * 1037.6167 / 1031.9261
            z = wave0 / 1031.9261 - 1
            print(f"OVI: z = {z}")
            self.message_window.setText(f"OVI: z = {z}")

        elif (ionname == 'NeVIII'):
            wave1 = wave0 * 780.324 / 770.409
            z = wave0 / 770.409 - 1
            print(f"NeVIII: z = {z}")
            self.message_window.setText(f"NeVIII: z = {z}")

        elif (ionname == 'SiIV'):
            wave1 = wave0 *  1402.77291 / 1393.76018
            z = wave0 / 1393.76018 - 1
            print(f"SiIV: z = {z}")
            self.message_window.setText(f"SiIV: z = {z}")

        elif (ionname == 'SiIV'):
            wave1 = wave0 *  1402.77291 / 1393.76018
            z = wave0 / 1393.76018 - 1
            print(f"SiIV: z = {z}")
            self.message_window.setText(f"SiIV: z = {z}")

        elif (ionname == 'Lyb'):
            wave1 = wave0 *  1215.6701 / 1025.7223
            z = wave0 / 1025.7223 - 1
            print(f"HI: z = {z}")
            self.message_window.setText(f"HI: z = {z}")

        elif (ionname == 'Lya'):
            wave1 = wave0 *   1025.7223/ 1215.6701
            z = wave0 / 1215.6701 - 1
            print(f"HI: z = {z}")
            self.message_window.setText(f"HI: z = {z}")


  
               
        self.ax.plot([wave0,wave0,wave1,wave1],[yval,yval+.5,yval+.5,yval],color='r')
        
        if ionname == 'FeII':
            self.ax.text(0.5*(wave0+wave2),yval+0.7,ionname +' z: ' +str(np.round(z,4)) ,rotation=90,verticalalignment='bottom') 
            self.ax.plot([wave1,wave1,wave2,wave2],[yval,yval+.5,yval+.5,yval],color='r')
        else:
            self.ax.text(0.5*(wave0+wave1),yval+0.7,ionname +' z: ' +str(np.round(z,4)) ,rotation=90,verticalalignment='bottom') 


        self.canvas.draw()


    #spec plot is connected to the smoothing/unsmoothing events
    #flux is always stored as ax.lines[0], error as ax.lines[1]. So we can add two new S/US plots
    #replace lines[0] and [1], then delete the last two lines added to the ax.lines list to keep the reference as [1] and [0]
    def specplot(self):
        ax=self.spectrum.gca()
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        replace_flux = ax.step(self.wave,self.smoothed_spectrum,'-',lw=0.5,label='smooth',color=clr['white'])
        self.ax.lines[1] = replace_flux[0]
        del self.ax.lines[-1]
        replace_error = ax.step(self.wave,self.smoothed_error,'-',lw=0.5,label='smooth',color=clr['pale_red'],zorder=2)
        self.ax.lines[0] = replace_error[0]
        del self.ax.lines[-1]
        
        self.spectrum.canvas.draw() 
        

    # Read and draw linelist    
    def DrawLineList(self,label,color='white',remove = False,hide =False,Man_Transition = False):
        self.active_zabs.setText(str(self.zabs))
        self.throw_away = False
        
        #this is for 'Z' functionality deletes all lines and plots a new line
        if ((remove == True) & (self.row_remove == False) & (self.hide == False)):
            del self.ax.lines[2:len(self.ax.lines)]
            #self.ax.texts = []
            while self.ax.texts:
                self.ax.texts.pop()
            while self.ax.collections:
                self.ax.collections.pop()
            self.throw_away = True
            
        #HIDE PROCEDURE
        if ((remove == True) & (self.hide == True)):
            for ii in self.text[self.row][:]:
                ii.remove()
            for ii in self.zabs_line_plot[self.row][:]:
                ii.remove()
            
            self.zabs_line_plot[self.row] = []
            self.text[self.row] = []
            self.spectrum.canvas.draw()
            return
        
        # REMOVE PROCEDURE
        elif ((remove == True) & (self.row_remove == True)):
            #remove lines from canvas
            for ii in self.text[self.row][:]:
                try: ii.remove()
                except: pass

            for ii in self.zabs_line_plot[self.row][:]:
                try: ii.remove()
                except: pass
                
            #delete line identifiers
            if len(self.text)>0:
                del self.text[self.row]
            try: 
                del self.z_list[self.row]
            except:
                if len(self.z_list) < 1:
                    pass
                else:
                    print('copy error and cntrl f!')
            del self.zabs_line_plot[self.row]
            self.spectrum.canvas.draw()
            return
        
        #PLOT PROCEDURE (otherwise continue plotting)
        else:
            self.label=label
            linecolor=clr[color]
            data=line.read_line_list(label)
            self.data = data
            xlim=self.init_xlims
            ylim=self.ax.get_ylim()
            
            tt_temp = []
            temp_plots = []
            for i in range(0, len(data)):
                if ((data[i]['wrest']*(1.+self.zabs) >= np.double(xlim[0])) & (data[i]['wrest']*(1.+self.zabs) <= np.double(xlim[1]))):
                    xdata=[data[i]['wrest']*(1.+self.zabs),data[i]['wrest']*(1.+self.zabs)]
                    ss=self.ax.transData.transform((0, .9))
                    ydata=[0,ylim[1]]
                    lineplot,=self.ax.plot(xdata,ydata,'--',color=linecolor)                    
                    tt=self.ax.text(xdata[0],0.75*ylim[1],data[i]['ion']+'  z='+ np.str(self.zabs),rotation=90)
                    
                    #append text and plot artist objects
                    tt_temp.append(tt)
                    temp_plots.append(lineplot)


            self.tt_temp = tt_temp
            self.temp_lines = len(tt_temp)#uneccessary do use additional id
            self.temp_plots = temp_plots
            #if first zabs or new zabs

            if ((self.zabs not in self.z_list or len(self.zabs_line_plot)<1) and (self.throw_away == False) and (Man_Transition == False)):
                self.zabs_line_plot.append(temp_plots)
                self.text.append(tt_temp)
                self.z_list.append(self.zabs)
            
            #if replotting a zabs, then reindex its lines
            elif ((self.zabs in self.z_list) and (self.throw_away == False) and (Man_Transition == False)): #maybe need the zabs update option to make sure nothing happens in the 'Z' event
                self.zabs_line_plot[self.row] = temp_plots
                self.text[self.row] = tt_temp
            #need to store new text objects     
            if self.hide == True:
                self.text[self.row] = tt_temp
                self.zabs_line_plot[self.row] = temp_plots
        self.spectrum.canvas.draw()


    def get_linelist(self):
        items = ( "LLS", "LLS Small", "DLA","LBG","Gala","None")

        item, ok = QInputDialog.getItem(self, "Select Linelist", 
            "Line Lists", items, 0, False)            
        
        if ok and item:
            self.label=item


class Redshift_Guess:
    def __init__(self,parent):
        parent.zabs = np.double(parent.active_zabs.text())
        color = parent.combo_color_main.currentText()
        label = parent.combo_lines.currentText()
        parent.DrawLineList(label,color=color,remove = True,hide =False)

        
class manage_identified_absorbers(QWidget):
    def __init__(self,parent):
        super().__init__()
        self.resize(300,250)
        self.layout = QHBoxLayout()
        self.table = QTableWidget()
        self.setWindowTitle('Zabs Manager')
        
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(('Line Lists','z','Color','Plot','Remove','Hide'))
        self.table.setRowCount(2)
        #make longer color list
        clrlist=list(clr.keys()) 
        self.combo_options = clrlist[1:]#['yellow','orange','red','green','white']
        self.line_options = ['LLS','LLS Small','DLA','LBG','Gal','None']
        
        for i in range(2):
            combo = QComboBox()
            for items in self.combo_options:
                combo.addItem(items)
            self.table.setCellWidget(i,2,combo)
        for i in range(2):
            combo = QComboBox()
            self.plotbut = QPushButton("Plot",self)
            self.plotbut.clicked.connect(lambda: self.plot(parent))
            
            self.removebut = QPushButton("Remove",self)
            #overwrites objectName
            self.removebut.setObjectName(str(i))
            self.removebut.clicked.connect(lambda: self.remove(parent))
            
            self.hidebut = QPushButton("Hide",self)
            self.hidebut.clicked.connect(lambda: self.hide(parent))
            
            for items in self.line_options:
                combo.addItem(items)
            self.table.setCellWidget(i,0,combo)
            self.table.setCellWidget(i,3,self.plotbut)
            self.table.setCellWidget(i,4,self.removebut)
            self.table.setCellWidget(i,5,self.hidebut)

        self.table.cellChanged.connect(self.cellchanged)
        self.layout.addWidget(self.table)
        self.setLayout(self.layout)
        self.identified_line=None
        self.table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.table.resizeColumnsToContents()

    def plot(self,parent):
        row,column = self.get_index(self)
        parent.row = row
        linelist = parent.abs_plot.table.cellWidget(row,0).currentText()
        z = np.double(parent.abs_plot.table.item(row,1).text())
        parent.zabs = z
        color = parent.abs_plot.table.cellWidget(row,2).currentText()

        new_row = pd.Series(data = {'Zabs': z, 'list': linelist, 'color': color})
        if parent.zabs_list.shape[0]< 1:
            parent.zabs_list=parent.zabs_list.append(new_row,ignore_index=True)
    
        elif parent.zabs not in parent.zabs_list.Zabs.to_numpy():
            parent.zabs_list=parent.zabs_list.append(new_row,ignore_index=True)
#         elif parent.zabs in parent.zabs_list.Zabs.to_numpy():

        parent.hide = True #update text objects for removal/hiding
        parent.DrawLineList(linelist,color,remove = False)
        parent.hide = False
        
        #ensure entered redshift is visible
        self.table.resizeColumnsToContents()
        
        #if the user has hidden the plot, turn hide background back to gray
        self.table.cellWidget(row,5).setStyleSheet('background-color : QColor(53, 53, 53)')
        #need to write a function for making the table entries the below code is redundant
        if self.table.rowCount() == parent.zabs_list.shape[0]:
            new_row = self.table.rowCount()
            self.table.setRowCount(self.table.rowCount()+1)
            #self.table.insertRow(self.table.rowCount())
            
            combo_color = QComboBox()
            for items in self.combo_options:
                combo_color.addItem(items)
            self.table.setCellWidget(new_row,2,combo_color)
            
            self.plotbut = QPushButton("Plot",self)
            self.plotbut.clicked.connect(lambda: self.plot(parent))
            
            self.removebut = QPushButton("Remove",self)

#             self.removebut.setObjectName(str(i))
            self.removebut.clicked.connect(lambda: self.remove(parent))
            
            self.hidebut = QPushButton("Hide",self)
            self.hidebut.clicked.connect(lambda: self.hide(parent))
            combo = QComboBox()
            for items in self.line_options:
                combo.addItem(items)
            self.table.setCellWidget(new_row,0,combo)
            self.table.setCellWidget(new_row,3,self.plotbut)
            self.table.setCellWidget(new_row,4,self.removebut)
            self.table.setCellWidget(new_row,5,self.hidebut)
            self.table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
            
            
    
            
    def hide(self,parent):
        row,column = self.get_index(self)#parent.abs_plot.table.currentRow()
        parent.row = row
        linelist = parent.abs_plot.table.cellWidget(row,0).currentText()
        color = parent.abs_plot.table.cellWidget(row,2).currentText()
        parent.hide = True
        parent.DrawLineList(linelist,color,remove = True)
        parent.hide = False
        parent.abs_plot.table.cellWidget(row,5).setStyleSheet('background-color : green')
        
    def remove(self,parent):
        #get row of deleted zabs values
        row,column = self.get_index(self)
        parent.row = row
        
        #need to remove the objects stored in canvas for plot/hide functionality with condition row_remove==True
        parent.row_remove = True
        color = parent.abs_plot.table.cellWidget(row,2).currentText()
        parent.DrawLineList('label',color,remove = True)
        parent.row_remove = False
        
        #remove from the line_list manager if user has stored
        zabs = float(parent.abs_plot.table.item(row,1).text())
        if zabs in parent.line_list.Zabs.tolist():
            index = parent.line_list[parent.line_list['Zabs'] == zabs].index
            parent.line_list = parent.line_list.drop(index,inplace=False)

        #remove from the table widget
        self.table.removeRow(row)
        
        #remove from the zabs_manager
        parent.zabs_list = parent.zabs_list.drop(parent.zabs_list.index[row])
        

    def cellchanged(self):
        col = self.table.currentColumn()
        row = self.table.currentRow()
        try: text = self.table.currentItem().text()
        except: pass
    def get_index(self,parent):
        button = QtWidgets.QApplication.focusWidget()
        index = self.table.indexAt(button.pos())
        if index.isValid():
            row = index.row()
            column = index.column()
        return row,column
    
        



class Manual_Transition(QWidget):
    def __init__(self,parent):
        super().__init__()
        self.resize(200,900)
        self.layout = QVBoxLayout()
        self.line_options = ['LLS','LLS Small','DLA','LBG','Gal','None']
        self.combo_ll = QComboBox()
        for items in self.line_options:
            self.combo_ll.addItem(items)
        self.layout.addWidget(self.combo_ll)
        self.layout.setAlignment(QtCore.Qt.AlignTop)
        self.setLayout(self.layout)
        self.combo_ll.currentIndexChanged.connect(lambda: self.line_change(parent))
    
        
        #Obtain Absorber
        self.Transitions = QListWidget()
        data = line.read_line_list('LLS')
        self.wavelist = []
        for ii in range(len(data)):
            self.Transitions.addItem(data[ii]['ion'])
            self.wavelist.append(data[ii]['wrest'])
        self.Transitions.itemClicked.connect(lambda: self.transition_change(parent))
            
        self.layout.addWidget(self.Transitions)
        
        #Need to obtain linelist
    def line_change(self,parent):
        parent.label = self.combo_ll.currentText()
        data = line.read_line_list(self.combo_ll.currentText())
        self.Transitions.clear()
        
        for ii in range(len(data)):
            self.Transitions.addItem(data[ii]['ion'])
            
    def transition_change(self,parent):
        
        del parent.ax.lines[2:]
        #parent.ax.texts = []
        while parent.ax.texts:
            parent.ax.texts.pop()
        while parent.ax.collections:
            parent.ax.collections.pop()

        
        parent.label = self.combo_ll.currentText()
#         parent.active_ion.setText(self.Transitions.currentItem().text())
        lambda_rest = self.wavelist[self.Transitions.currentRow()]
        parent.lambda_rest = lambda_rest
        
        parent.zabs= np.round((parent.xdata -lambda_rest)/lambda_rest,4)
        parent.active_zabs.setText(str(parent.zabs))
        parent.DrawLineList(parent.label,color=parent.combo_color_main.currentText(),remove=False,Man_Transition=True)


        
class Catalog:
    def __init__(self,parent):
        try:self.table = parent.abs_plot.table
        except: self.table=parent.table
        new_row = self.table.rowCount()
        self.table.setRowCount(self.table.rowCount()+1)


        combo_color = QComboBox()
        for items in parent.combo_options:
            combo_color.addItem(items)
        self.table.setCellWidget(new_row,2,combo_color)

        self.plotbut = QPushButton("Plot",self.table)
        self.plotbut.clicked.connect(lambda: manage_identified_absorbers.plot(parent.abs_plot,parent))

        self.removebut = QPushButton("Remove",self.table)

        self.removebut.clicked.connect(lambda: manage_identified_absorbers.remove(parent.abs_plot,parent))

        self.hidebut = QPushButton("Hide",self.table)
        self.hidebut.clicked.connect(lambda: manage_identified_absorbers.hide(parent.abs_plot,parent))
        combo = QComboBox()
        for items in parent.line_options:
            combo.addItem(items)
        self.table.setCellWidget(new_row,0,combo)
        self.table.setCellWidget(new_row,3,self.plotbut)
        self.table.setCellWidget(new_row,4,self.removebut)
        self.table.setCellWidget(new_row,5,self.hidebut)
        self.table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.table.resizeColumnsToContents()
        parent.spectrum.canvas.focusWidget()
        try:parent.manT.close()
        except: pass

        
class SaveCatalog(QWidget):
    def __init__(self,parent):
        super().__init__()
        self.resize(700,200)
        
        #set a timer after setting background green. Then close save window
        def onsave(self,parent,method):
            if method == 0:
                file = self.line.text()
                parent_dir = Path(file).parent
                if parent_dir.is_dir() == False:
                    parent_dir.mkdir(exist_ok=True)
                parent.zabs_list.to_csv(file)
                self.savebut.setStyleSheet('background-color : green')
            elif method == 1:
                file = self.line1.text()
                parent_dir = Path(file).parent
                if parent_dir.is_dir() == False:
                    parent_dir.mkdir(exist_ok=True)
                parent.line_list.to_csv(file, sep=' ')
                self.savebut1.setStyleSheet('background-color : green')
            else:
                directory = Path(self.line2.text())
                if directory.is_dir() == False:
                    directory.mkdir(exist_ok=True)

                if ((self.LineLists == 'All') or (self.LineLists == 'Partial')):
                    parent.line_list.to_csv(str((directory/'Identified_LineList.txt').resolve()), sep=' ')


                file = directory / 'Absorber_Catalog.csv'
                parent.zabs_list.to_csv(str(file.resolve()))
                self.savebut2.setStyleSheet('background-color : green')
                
        self.continueSaving = True
#         layout = QHBoxLayout()
        
        #if no linelists have been identified
        if parent.line_list.shape[0] == 0:
            buttonReply = QMessageBox.question(self,"Missing Linelist" ,"No Linelists have been obtained: Proceed?",QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if buttonReply == QMessageBox.Yes:
                self.LineLists = 'None' 
            else: self.continueSaving = False
        
        #if more absorbers in zabs manager than unique zabs in linelist manager
        if (parent.zabs_list.shape[0] != len(np.unique(parent.line_list.Zabs.tolist()))) and (self.continueSaving==True) and (parent.line_list.shape[0] != 0):
            buttonReply = QMessageBox.question(self,"Missing Linelist" ,"Missing linelists for the cataloged absorbers, continue?",QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            
            if buttonReply == QMessageBox.Yes:
                self.LineLists = 'Partial'
            else: self.continueSaving = False
        else:
            self.LineLists = 'All'
            
        if self.continueSaving == True:
            directory = Path.cwd()
            directory = directory /'SpecPlot_Projects'
            
            #layouts
            main = QVBoxLayout()
            layout = QHBoxLayout()
            layout1 = QHBoxLayout()
            layout2 = QHBoxLayout()
            layouts = [layout,layout1,layout2]
            
            #labels
            label = QLabel('Absorber Catolog Save (formats accepted: .csv)')
            label1 = QLabel('Linelist Save (formats accepted: .txt, .log)')
            label2 = QLabel('Directory Save (enter parent folder, filename defaults: Absorber_Catolog.csv, Identified_Linelist.txt)')
            labels = [label,label1,label2]
            
            #line edits
            self.line = QLineEdit(self)
            self.line1 = QLineEdit(self)
            self.line2 = QLineEdit(self)
            self.line.setText(str((directory / 'Absorber_Catalog.csv').resolve()))
            self.line1.setText(str((directory / 'Identified_Linelist.txt').resolve()))
            self.line2.setText(str(directory.resolve()))
            lines = [self.line,self.line1,self.line2]
            
            #save buttons
            self.savebut = QPushButton("Save",self)
            self.savebut.clicked.connect(lambda: onsave(self,parent,0))
            self.savebut1 = QPushButton("Save",self)
            self.savebut1.clicked.connect(lambda: onsave(self,parent,1))
            self.savebut2 = QPushButton("Save",self)
            self.savebut2.clicked.connect(lambda: onsave(self,parent,2))
            savebuts = [self.savebut,self.savebut1,self.savebut2]
            
            #create layout
            for ii in range(3):
                main.addWidget(labels[ii])
                layouts[ii].addWidget(lines[ii])
                layouts[ii].addWidget(savebuts[ii])
                main.addLayout(layouts[ii])
                

            
            self.setLayout(main)
            
        #Tell user which absorbers dont have a linelist by marking the comboboxes in red background and close save window
        elif self.continueSaving == False: #and ((self.LineLists == 'Partial') or (self.LineLists == 'None')):
            line_zabs = np.unique(parent.line_list.Zabs.tolist())
            manager_zabs = parent.zabs_list.Zabs.tolist()
            
            # for zabs in zabs manager that arent in line_list
            if len(line_zabs>0):
                for z in manager_zabs:
                    if z in line_zabs:
                        pass
                    else:
                        index = parent.zabs_list.Zabs[parent.zabs_list.Zabs == z].index[0]
                        parent.abs_plot.table.cellWidget(index,0).setStyleSheet('background-color : red')
                        parent.abs_plot.table.cellWidget(index,2).setStyleSheet('background-color : red')
    
                    #setcellWidgetbackground color
                self.close()
            else:
                index = parent.zabs_list.Zabs[parent.zabs_list.Zabs == parent.zabs_list.Zabs.tolist()[0]].index[0]
                parent.abs_plot.table.cellWidget(index,0).setStyleSheet('background-color : red')
                parent.abs_plot.table.cellWidget(index,2).setStyleSheet('background-color : red')
                self.close()
        #else close save window
        else:
            self.close()

class LoadCatalog(QWidget):
    def __init__(self,parent):
        super().__init__()
        main_layout = QVBoxLayout(self)
        self.setLayout(main_layout)
        #lay1= absorbercat.csv load; lay2= linelist.txt load; lay3=directory containing both load
        layout = QHBoxLayout(self)
        layout1 = QHBoxLayout(self)
        layout2 = QHBoxLayout(self)
        
        self.label = QLabel("Enter Zabs Manager dir+file (Absorber_Catalog.csv is default)")
        self.entry = QLineEdit(self)
        self.browse = QPushButton("Browse")
        self.browse.clicked.connect(lambda: self.browsefiles(parent,0))
        
        self.label1 = QLabel("Enter Linelist dir+file (Identified_Linelist.txt is default)")
        self.entry1 = QLineEdit(self)
        self.browse1 = QPushButton("Browse")
        self.browse1.clicked.connect(lambda: self.browsefiles(parent,1))
        
        self.label2 = QLabel("Enter Directory containing both")
        self.entry2 = QLineEdit(self)
        self.browse2 = QPushButton("Browse")
        self.browse2.clicked.connect(lambda: self.browsefiles(parent,2))
        
        #Add widgets to layout
        labs = [self.label,self.label1,self.label2]
        entrys = [self.entry,self.entry1,self.entry2]
        browses = [self.browse,self.browse1,self.browse2]
        lays = [layout,layout1,layout2]
        for ii in range(len(labs)):
            main_layout.addWidget(labs[ii])
            lays[ii].addWidget(entrys[ii])
#             lays[ii].addWidget(buttons[ii])
            lays[ii].addWidget(browses[ii])
            main_layout.addLayout(lays[ii])
        
        
        
        #setAutoDefault registers 'enter' to "load" function
        #this should correspond to which of the boxes are filled
        #self.button.setAutoDefault(True)
#         self.entry.returnPressed.connect(self.button.click)
                                
    def browsefiles(self,parent,method):
        if method == 0: # get csv
            fname,_ = QFileDialog.getOpenFileName(self, 'Open file', os.getcwd(),"CSV files (*csv)")
            parent.zabs_list = pd.read_csv(fname)
            parent.zabs_list = parent.zabs_list[parent.zabs_list.keys()[1:]]
            #Populate zabs list and hide plots
            self.populate_zabs_manager(parent)
            
            
        elif method == 1: #get txt
            fname,_ = QFileDialog.getOpenFileName(self, 'Open file', os.getcwd(),"Text files (*txt *.log)")
            parent.line_list = pd.read_csv(fname,sep=' ')
            parent.line_list = parent.line_list[parent.line_list.keys()[1:]]
            self.plot_identified_lines(parent)
            
        else:
            dirs = QFileDialog.getExistingDirectory(self,"Open a folder",os.getenv("HOME"),QFileDialog.ShowDirsOnly)
            for files in os.listdir(dirs):
            #read in zabs manager
                if files.endswith('.csv'):
                    fq_filename = os.path.join(dirs,files)
                    parent.zabs_list = pd.read_csv(fq_filename)
                    parent.zabs_list = parent.zabs_list[parent.zabs_list.keys()[1:]]
                    self.populate_zabs_manager(parent)
                #read in linelist manager
                elif files.endswith('.txt'):
                    fq_filename = os.path.join(dirs,files)
                    parent.line_list = pd.read_csv(fq_filename,sep=' ')
                    parent.line_list = parent.line_list[parent.line_list.keys()[1:]]
                    self.plot_identified_lines(parent)

            self.close()
    def populate_zabs_manager(self,parent):
        zabs = parent.zabs_list['Zabs'].values.tolist()
        self.green_index = []
        for ii in range(parent.zabs_list.shape[0]):       
            if ii >= parent.abs_plot.table.rowCount():
                #Catalog class create the new zabs_manager rows
                Catalog(parent)
                
            #Fill the newly created rows with saved values
            parent.row = ii
            parent.zabs = zabs[ii]
            parent.abs_plot.table.setItem(ii,1,QtWidgets.QTableWidgetItem(str(zabs[ii])))
            parent.abs_plot.table.cellWidget(ii,2).setCurrentText(parent.zabs_list['color'].iloc[ii])
            parent.abs_plot.table.cellWidget(ii,0).setCurrentText(parent.zabs_list['list'].iloc[ii])
            parent.DrawLineList(parent.zabs_list['list'].iloc[ii],color = parent.zabs_list['color'].iloc[ii])
            if ii > 0:
                parent.hide = True
                parent.DrawLineList(parent.zabs_list['list'].iloc[ii],parent.zabs_list['color'].iloc[ii],remove = True)
                parent.hide = False
                parent.abs_plot.table.cellWidget(ii,5).setStyleSheet('background-color : green')
                self.green_index.append(ii)
        Catalog(parent)
    
    def plot_identified_lines(self,parent):
        #if zabs_manager already loaded
        if parent.zabs_list.shape[0]>0:
            # clear plotted lines/text from zabs_manager such that the new lines for plot/remove/hide are same color as linelist
            del parent.ax.lines[2:]; parent.zabs_line_plot = []
            #parent.ax.texts = [];
            while parent.ax.texts:
                parent.ax.texts.pop()
            while parent.ax.collections:
                parent.ax.collections.pop()

            parent.texts = []

            for z in parent.zabs_list.Zabs.tolist():
                index = parent.line_list[parent.line_list['Zabs'] == z].index
                color = parent.zabs_list.color[parent.zabs_list.Zabs == z].values[0]
                ylim=parent.ax.get_ylim()
                #for all lines at that redshift
                for i in index:
                    xdata = [parent.line_list.loc[i].Wave_obs,parent.line_list.loc[i].Wave_obs]
                    #ylow = np.interp(xdata[0],parent.wave,parent.flux)+.75
                    lineplot,=parent.ax.plot(xdata,[1.5*ylim[0],0.75*ylim[1]],'--',color=clr[color])
                    tt = parent.ax.text(xdata[0],0.75*ylim[1],parent.line_list.loc[i].Name+'  z='+ np.str(parent.line_list.loc[i].Zabs),rotation=90)
                    parent.identified_lines.append(lineplot)
                    parent.identified_text.append(tt)
                parent.text = [[]]*len(parent.zabs_list.Zabs.tolist())
                parent.zabs_line_plot = [[]]*len(parent.zabs_list.Zabs.tolist())
                
            #if linelist has absorbers that are not cataloged
            if len(np.unique(parent.line_list.Zabs.tolist())) > len(parent.zabs_list.Zabs.tolist()):
                zabs_list = parent.zabs_list.Zabs.tolist()
                line_list_zabs = np.unique(parent.line_list.Zabs.tolist())
                for z in line_list_zabs:
                    if z not in zabs_list:
                        index = parent.line_list[parent.line_list['Zabs'] == z].index
                        ylim=parent.ax.get_ylim()
                        for i in index:
                            xdata = [parent.line_list.loc[i].Wave_obs,parent.line_list.loc[i].Wave_obs]
                            #ylow = np.interp(xdata[0],parent.wave,parent.flux)+.75
                            lineplot,=parent.ax.plot(xdata,[1.5*ylim[0],0.75*ylim[1]],'--',color='y')
                            tt = parent.ax.text(xdata[0],0.75*ylim[1],parent.line_list.loc[i].Name+'  z='+ np.str(parent.line_list.loc[i].Zabs),rotation=90)
                            parent.identified_lines.append(lineplot)
                            parent.identified_text.append(tt)

        else:

            tt_all = []
            temp_plot = []
            ylim=parent.ax.get_ylim()
            for i in range(parent.line_list.shape[0]):
                xdata = [parent.line_list.loc[i].Wave_obs,parent.line_list.loc[i].Wave_obs]
                #ylow = np.interp(xdata[0],parent.wave,parent.flux)+.75
                lineplot,=parent.ax.plot(xdata,[1.5*ylim[0],0.75*ylim[1]],'--',color='y')
                tt = parent.ax.text(xdata[0],0.75*ylim[1],parent.line_list.loc[i].Name+'  z='+ np.str(parent.line_list.loc[i].Zabs),rotation=90)
                parent.identified_lines.append(lineplot)
                parent.identified_text.append(tt)
                
        try:
            for ii in self.green_index:
                parent.abs_plot.table.cellWidget(ii,5).setStyleSheet('background-color : QColor(53, 53, 53)')
        except:
            pass
        
        parent.Identified_line_plot.setStyleSheet('background-color : green')
        parent.message_window.setText(" ")
        parent.identified_line_active = True
        parent.spectrum.canvas.draw()
   
                
class HelpWindow(QtWidgets.QWidget):
    def __init__(self,parent=None):
        super(HelpWindow, self).__init__(parent)
        self.resize(500,850)
        label = QtWidgets.QLabel(HELP,self)
        
        
class Identified_plotter:
    def __init__(self,parent):
        if parent.line_list.shape[0]==0:
            parent.Identified_line_plot.setStyleSheet('background-color : red')
            parent.message_window.setText("No transitions have been identified")
        else:
            #if active remove plots and set indicator back to gray
            if parent.identified_line_active == True:
                for ii in parent.identified_lines:
                    ii.remove()
                for ii in parent.identified_text:
                    ii.remove()
                parent.identified_text = []; parent.identified_lines = []
                parent.Identified_line_plot.setStyleSheet('background-color : QColor(53, 53, 53)')
                parent.identified_line_active = False
            #else 
            else:
                parent.identified_line_active = True
                if parent.zabs_list.shape[0]>0:
                    for z in parent.zabs_list.Zabs.tolist():
                        index = parent.line_list[parent.line_list['Zabs'] == z].index
                        color = parent.zabs_list.color[parent.zabs_list.Zabs == z].values[0]
                        ylim=parent.ax.get_ylim()
                        for i in index:
                            xdata = [parent.line_list.loc[i].Wave_obs,parent.line_list.loc[i].Wave_obs]
                            #ylow = np.interp(xdata[0],parent.wave,parent.flux)+.75
                            lineplot,=parent.ax.plot(xdata,[1.5*ylim[0],0.75*ylim[1]],'--',color=clr[color])
                            tt = parent.ax.text(xdata[0],0.75*ylim[1],parent.line_list.loc[i].Name+'  z='+ np.str(parent.line_list.loc[i].Zabs),rotation=90)
                            parent.identified_lines.append(lineplot)
                            parent.identified_text.append(tt)

                
                    #if linelist has absorbers that are not cataloged
                    if len(np.unique(parent.line_list.Zabs.tolist())) > len(parent.zabs_list.Zabs.tolist()):
                        zabs_list = parent.zabs_list.Zabs.tolist()
                        line_list_zabs = np.unique(parent.line_list.Zabs.tolist())
                        for z in line_list_zabs:
                            if z not in zabs_list:
                                index = parent.line_list[parent.line_list['Zabs'] == z].index
                                ylim=parent.ax.get_ylim()
                                for i in index:
                                    xdata = [parent.line_list.loc[i].Wave_obs,parent.line_list.loc[i].Wave_obs]
                                    #ylow = np.interp(xdata[0],parent.wave,parent.flux)+.75
                                    lineplot,=parent.ax.plot(xdata,[1.5*ylim[0],0.75*ylim[1]],'--',color='y')
                                    tt = parent.ax.text(xdata[0],0.75*ylim[1],parent.line_list.loc[i].Name+'  z='+ np.str(parent.line_list.loc[i].Zabs),rotation=90)
                                    parent.identified_lines.append(lineplot)
                                    parent.identified_text.append(tt)
                else:
                    ylim=parent.ax.get_ylim()
                    for i in range(parent.line_list.shape[0]):
                        xdata = [parent.line_list.loc[i].Wave_obs,parent.line_list.loc[i].Wave_obs]
                        #ylow = np.interp(xdata[0],parent.wave,parent.flux)+.75
                        lineplot,=parent.ax.plot(xdata,[1.5*ylim[0],0.75*ylim[1]],'--',color='y')
                        tt = parent.ax.text(xdata[0],0.75*ylim[1],parent.line_list.loc[i].Name+'  z='+ np.str(parent.line_list.loc[i].Zabs),rotation=90)
                        parent.identified_lines.append(lineplot)
                        parent.identified_text.append(tt)
                                    
                parent.Identified_line_plot.setStyleSheet('background-color : green')        
            parent.spectrum.canvas.draw()
                
                
                
#---------------for Vstack---------------------#

def prepare_absorber_object(z_abs,wave,flux,error,line_flg='LLS',vlim=[-1000,1000]):
    
    # Read the full linelist
    data=line.read_line_list(line_flg)
    wavelist=[]
    for i in range(0,len(data)):
        wavelist.append(data[i]['wrest'])
        
        
    wavelist=np.array(wavelist)

    
    # select the lines within the wavelength range only
    q= np.where((wavelist > np.min(wave)/(1.+z_abs)) &  (wavelist < (np.max(wave)/(1.+z_abs))));
    
    # Total transitions visible within the wavelength window
    nTot= len(q[0]);
    
    wavelist_selected=wavelist[q]
    
    absys=A.Absorber(z_abs,wave,flux,error,list(wavelist_selected),window_lim=vlim,nofrills=True)   
    
    return absys.ions
    

class vStack:
    def __init__(self,parent,wave,flux,error,line_flg,zabs=0,vlim=[-1000.,1000.]):
        self.parent = parent
        self.parent_canvas = parent.canvas
        self.zabs=zabs
        self.vlim=vlim
        self.ions=prepare_absorber_object(zabs,wave,flux,error,line_flg=line_flg)
        #-----full spectra properties---------#
        self.z = self.ions['Target']['z']; self.flux = self.ions['Target']['flux']
        self.wave = self.ions['Target']['wave']; self.error = self.ions['Target']['error']
        
               
        self.keys = list(self.ions.keys())[:-1] # last item is the full target spectrum
        
        self.nions = np.int(len(self.keys))
        #Flag to know if it is a detection or not
        #Set everything by default to non-detection
        for i in (self.keys):
            self.ions[i]['flag']=0
        
        
        #-----Sorting out how many pages are needed---------#

        
        self.page=1
        self.plotppage=12
        self.nrow=int(self.plotppage/3)
        self.ncol=int((self.plotppage)/self.nrow)
        self.npages=int((self.nions//self.plotppage))
        # calculate number of pages needed
        if self.npages % (self.plotppage) > 0:
            self.npages += 1
        
        
        fig= Figure()#figure(figsize=(12,8))
        self.fig=fig
        self.axes=list(range(self.plotppage))
        for i in range(self.plotppage):
            self.axes[i]=self.fig.add_subplot(self.nrow,self.ncol,i+1)
        self.axes=np.array(self.axes)
        self.vPlot()
        self.fig.subplots_adjust(hspace=0)

        self.canvas = FigureCanvasQTAgg(self.fig)
        
        self.fig.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.fig.canvas.setFocus()
        self.cid = self.fig.canvas.mpl_connect('key_press_event',self.onkb)
        
        parent.plot_layout.removeWidget(parent.canvas)
        plt.close(self.parent.spectrum)
        parent.plot_layout.insertWidget(1,self.canvas)#FigureCanvasQTAgg(self.fig))


    def onkb(self,event):
        #set up custom y-limit
        if event.key=='Y':
            if event.inaxes in self.axes:
                i=np.where(event.inaxes==self.axes)[0][0]+self.plotppage*(self.page-1)
                Windowname='Manual y-Limits'
                instruction_text='Input range (e.g. 0.,2.)'
                #temp=input_txt_dlg(Windowname,instruction_text)
                #print(temp)
                #yrangetext=temp.filename
                #yrange = yrangetext.split(',')
                #yrange = np.array(yrange).astype('float32')
                ylim, ok = QInputDialog.getText(self.parent,Windowname,instruction_text)
                if ok:
                    ylimit = ylim.split(',')
                    ylimit = np.array(ylimit).astype('float32')
                    self.vPlot(ploti=i,yrange=[ylimit[0],ylimit[1]])
                
        #page right
        elif event.key=='>':
            self.page+=1
            if self.page>self.npages: 
                self.page=1
            self.vPlot()
        #page left
        elif event.key=='<':
            self.page-=1
            if self.page<1: 
                self.page=self.npages
            self.vPlot()
            
        #Toggle between detection-non-detection or blended    
        elif event.key =='w': #Detected,non-detected, blended-detection 
            if event.inaxes in self.axes:
                i=np.where(event.inaxes==self.axes)[0][0]+self.plotppage*(self.page-1)
                #set up a dumb toggling cycle
                temp_flag=self.ions[self.keys[i]]['flag']+1

                if temp_flag==0:
                    temp_flag =1
                elif temp_flag==1:
                    temp_flag==2
                else:
                    temp_flag=0
                self.ions[self.keys[i]]['flag']= temp_flag#
                self.vPlot(ploti=i,comment=False)
                
        #save linelist
        elif event.key=='S':
            #reinsert the primary spectrum canvas to the mainlayout for keyboard functionality
            self.parent.plot_layout.removeWidget(self.canvas)
            self.parent.plot_layout.insertWidget(1,self.parent.canvas)
            
            #Need to check if it is reevaluating a linelist, if so delete all lines with same 'zabs' value
            if self.parent.line_list.shape[0] > 0:
                zabs_list = self.parent.line_list.Zabs.tolist()
                if self.parent.zabs in zabs_list:
                    index = self.parent.line_list[self.parent.line_list['Zabs'] == self.parent.zabs].index
                    self.parent.line_list = self.parent.line_list.drop(index,inplace=False)

            keys = list(self.ions.keys())
            # based on line evaluation, add lines to the overall zabs manager linelist ('flag=1' is detected lines to add)
            for key in keys[:-1]:
                if self.ions[key]['flag'] == 1:
                    wave_obs = self.ions[key]['lam_0_z']
                    name = self.ions[key]['name']
                    zabs = self.ions['Target']['z']
                    new_row = pd.Series(data = {'Name': name, 'Wave_obs': wave_obs, 'Zabs': zabs})
                    self.parent.line_list=self.parent.line_list.append(new_row,ignore_index=True)
                    
            #lets keep zabs_list sorted by ascending zabs
            self.parent.line_list = self.parent.line_list.sort_values(by='Zabs')
            
            #close to canvas so spectrum canvas is again visible
            self.canvas.close()
            
            #after quitting save function, lineList and color will be red for unsaved lines, check and replace if saving a missing linelist
            
            #first need to find row
            row = self.parent.zabs_list.Zabs[self.parent.zabs_list.Zabs == self.parent.zabs].index[0]
            
            # red is #ff0000, if not red it will have no palette so this will need to be a try/pass
            try:
                if self.parent.abs_plot.table.cellWidget(row,0).palette().color(QtGui.QPalette.Background).name() == '#ff0000':
                    self.parent.abs_plot.table.cellWidget(row,0).setStyleSheet('background-color : QColor(53, 53, 53)')
                    self.parent.abs_plot.table.cellWidget(row,2).setStyleSheet('background-color : QColor(53, 53, 53)')
            except:
                pass
    def vPlot(self,ploti=None,comment=False,yrange=None):#spec,i=0):
        # global axesR
        if ploti is None:
            ploti=np.arange(self.plotppage*(self.page-1),min(self.plotppage*self.page,self.nions))
        else:
            ploti=[ploti]

        # clear the axes if needed
        for i in range(self.plotppage):
            if i not in ploti:
                self.clearstuff(i)
        # plot the things
        for i in ploti:
            self.plotstuff(i,comment=comment,yrange=yrange)

        self.fig.canvas.draw()

    def clearstuff(self, i):
        """ clearstuff() ensures that all the plots from the previously drawn page are cleared before plotting the next page"""
        ax=self.axes[i % self.plotppage]
        ax.clear()
        ax.set_axis_off()


        
    def plotstuff(self,i,comment=False,yrange=False):
        
        ax=self.axes[i % self.plotppage]
        #---------------Define variables for readability--------------#
        vel = self.ions[self.keys[i]]['vel']
        wave = self.ions[self.keys[i]]['wave']
        error = self.ions[self.keys[i]]['error']
        flux = self.ions[self.keys[i]]['flux']
        name = self.ions[self.keys[i]]['name']
        window_lim = self.ions[self.keys[i]]['window_lim']
        flag=self.ions[self.keys[i]]['flag']
        f0 = self.ions[self.keys[i]]['f']

        ax.clear()
        
        # Set a title giving the page number
        if i % self.plotppage ==0:
            ax.set_title('Page '+ str(self.page)+' of ' + str(self.npages),color=clr['teal'])


        ax.step(vel,flux/np.nanmean(flux),where='mid',color=clr['teal'])
        ax.step(vel,error/np.nanmean(flux),where='mid',color=clr['orange2'])
        
        ax.axhline(1,color=clr['light_gray'],linestyle='dotted')
        ax.axvline(0,color=clr['light_gray'],linestyle='dotted')
        ax.text(x=0.05, y=0.815, s=name, fontsize=10, transform=ax.transAxes,color=clr['red'])
        ax.text(x=0.75, y=0.815, s='f0: '+str(f0), fontsize=10, transform=ax.transAxes,color=clr['red'])
        
        if comment != False:
            ax.text(x=0.85, y=0.815, s=comment, fontsize=12, transform=ax.transAxes,color=clr['teal'])
        if yrange != False:
            ax.set_ylim(yrange)

        
        if flag is not None: #Display some measurement
            textout=self.plotText(flag=flag)
            if flag==1:
                textcolor=clr['yellow']
            else:
                textcolor=clr['light_gray']
            ax.text(x=0.05, y=0.01, s=textout, fontsize=12, transform=ax.transAxes,color=textcolor)
    
    def plotText(self,flag=1):
        if flag==1:
            text='Detection'       
        elif flag==0:
            text='Non-Detection'      
        elif flag==2:
            text ='Blended-detection'
        return text

#Initial inputs and callable class to run proram        
class input_txt_dlg(QWidget):
    def __init__(self,Windowname,instruction,default_text='test'):
        app = QApplication(sys.argv)
        main = popup_windows(Windowname,instruction,default_text=default_text)
        self.filename=main.filename
        main.show()
        sys.exit(app.exec_())        
        
        
class Popup_col_density(QWidget):
    #creates popup displaying the multiple transitions encompassed by EW
    def __init__(self,parent):
        super().__init__()
        self.list = QListWidget(self)
        self.setWindowTitle('Select Transition')
        self.resize(250,100)
        for ii in range(len(parent.fvals)):
            self.list.addItem(parent.fval_name[ii])
        
        self.list.itemDoubleClicked.connect(lambda: self.ion_selected(parent))
    
    #
    def ion_selected(self,parent):
        item_idx = self.list.currentRow()
        vel_lims = (np.asarray(parent.lam_lim[parent.tab])-parent.fval_wrest[item_idx]*(1+parent.zabs))*parent.c/(parent.fval_wrest[item_idx]*(1+parent.zabs))
        EW,sig_EW,cont,wave_slice,N,Nerr = parent.compute_EW(parent.wave/(1.+parent.zabs),parent.flux,parent.lam_lim[parent.tab]/(1+parent.zabs),parent.lam_ylim[parent.tab],parent.error,parent.fvals[item_idx],parent.fval_wrest[item_idx])
        EW=np.array(EW)*1000.
        sig_EW=np.array(sig_EW)*1000.
        parent.ax.plot(wave_slice,cont,'r--')
        Wval = 'EW [mAA]: '+ '%.1f' % EW + ' +/- ' + '%.1f' % sig_EW + '\n'
        logNerr = str(np.round(np.log10((N+Nerr)/(N-Nerr)/2),2))
        formatName = '$N_{'+parent.fval_name[item_idx]+'}$'
        logN = Wval + 'log ' + ' [$cm^{-2}$]: ' +str(np.round(np.log10(N),2)) + ' +/- ' +logNerr
        parent.ax.text(np.mean([parent.lam_lim])+.05,np.max(parent.lam_ylim)+0.2,logN, rotation=90,verticalalignment='bottom')
        parent.message_window.setText(logN)
        parent.lam_lim=[]
        parent.lam_ylim=[]
        parent.spectrum.canvas.draw()
        self.close()
        
     
        
#--------------------initialization of Application-----------------
class rb_plotspec():
    
    def __init__(self,wave,flux,error,zabs=0.):


        if not QtWidgets.QApplication.instance():
            app = QtWidgets.QApplication(sys.argv)
            app.setStyle("Fusion")

            # Now use a palette to switch to dark colors:
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

        else:
            app = QtWidgets.QApplication.instance() 



        #app = QtWidgets.QApplication(sys.argv)
         # Force the style to be the same on all OSs:
        main = mainWindow(wave,flux,error)
        main.resize(1700,900)
        
        #Center app on initialization
        qr = main.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        main.move(qr.topLeft())
        
        main.show()
        QtWidgets.QApplication.setQuitOnLastWindowClosed(True)
        app.exec_()
        app.quit()

