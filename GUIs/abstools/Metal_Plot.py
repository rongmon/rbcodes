import numpy as np
from IGM import compute_EW
from PyQt5 import QtCore, QtGui, QtWidgets
from GUIs.abstools import Absorber
from PyQt5.QtWidgets import QStyleFactory, QPushButton,QLineEdit,QMainWindow,QInputDialog,QLabel,QMessageBox,QScrollBar,QVBoxLayout,QHBoxLayout
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import pickle
import numpy.polynomial.legendre as L
import sys
from astropy.io import ascii
from utils import rb_utility as rt
from matplotlib import rcParams
rcParams['lines.linewidth'] = .9
clr=rt.rb_set_color()
import webbrowser
from pkg_resources import resource_filename

from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QPalette, QColor

from astropy.table import Table
from astropy.io import ascii

HELP =  '''
        ---------------------------------------------------------------------------
        This is an interactive 1D absorption line measurement toolbox.
        This allows for interactive continuum fitting and equivalent width measurement 
        of CGM/IGM/ISM absorption lines.

       Screen Layout:
            LHS/RHS = Left Hand Side/Right Hand Side
            LMB/RMB = Left Mouse Button/Right Mouse Button
            
            LHS shows raw spectrum with overlaid legendre poly for continuum fitting
            ---grayed regions indicate masked regions
            
            RHS shows normalized spectrum with velocity limits
        ------------------------------Mouse Clicks------------------------------------

            
        Useful Mouse Clicks:
            
            LHS LMB     : Add wavelengths within region set by two clicks to continuum fit.
            LHS RMB     : Remove wavelengths from continuum fit.
            RHS LMB     : Set lower velocity limit
            RHS RMB     : Set upper velocity limit
            
        Useful Keystrokes:            

            v           : place mouse on desired subplot
                                   LHS: manually enter regions to mask continuum 
                                   RHS: manually enter EW intergration limits
            
            V (RHS only): Use active subplot velocity limits for all RHS plots
            
            Up arrow    : Increase Polynomial Order [default 4]
            Down arrow  : Decrease Polynomial Order [default 4]

            m           : Measure EW/N for active subplot
            M           : Measure EW/N for ALL subplots
            
            1/2/0 (RHS only): flag absorber as
                              (0) positive detection
                              (1) upper limit 
                              (2) lower limit

            t (RHS only): Cycle text printed on absorbers. Displays logN, or EW
            ------------------------------------------------------------------------------
            Each tab displays up to 6 transitions. There are maximum 5 tabs allowed.
            This limits the total number of transitions that can be simultaneously analyized to 30.



 

            Written By: Sean Clark, Rongmon Bordoloi [2021]

                    '''



def shift2vel(z1,z2):
    # z1 at rest
    # z2 for which relative velocity is computed
    spl=2.9979e5;  #speed of light
    vel=((1215.67*(1.+z2))-(1215.67*(1.0 + z1)))*spl/(1215.67*(1.0 + z1))#(self.wrest-str['wave']*(1.0 + 0.))*spl/(str['wave']*(1.0 + 0.))
    return vel

def grab_intervening_linelist(filename,z_gal,wrest_galaxy,wavelength):
    import pandas as pd 
    #s=ascii.read(filename)
    s=pd.read_csv(filename,sep=' ')
    ion=s['Name'].values#s['col1']
    wobs=s['Wave_obs'].values #s['col3']
    zobs=s['Zabs'].values#s['col4']
    wrest=wobs/(1.+zobs)#s['col2']

    spl=2.9979e5;  #speed of light

    #compute relative velocity difference with the host galaxy
    delv=shift2vel(z_gal,zobs)
    # Now select only lines within the given wavelength range
    window_max=np.max(wavelength)#*(1.+z_gal)
    window_min=np.min(wavelength)#*(1.+z_gal)
    q=np.where( (wobs >= window_min) & (wobs <= window_max))#& ( np.abs(delv) >= 200.) )
    outlist={}
    if np.sum(q)>0:
        #If there are lines create a new list
        outlist['ion']=ion[q]
        outlist['wrest']=wrest[q]
        outlist['wobs']=wobs[q]
        outlist['zobs']=zobs[q]
        wobs_small=wobs[q]
        outlist['number']=len(wobs_small)
        vel=np.zeros((len(wobs_small),))
        outlist['delv']=delv[q]
        #Now computing velocity for each
        for i in range(0,len(wobs_small)):
            vel[i] = (outlist['wobs'][i]-wrest_galaxy*(1.0 + z_gal))*spl/(wrest_galaxy*(1.0 + z_gal))
        outlist['vel']=vel
    else:
        outlist['number']=0

    return outlist



def plot_intervening_lines(ax,outlist,delv):
    #Plot all intervening lines. 
    # Other lines associated with current absorber = blue
    # All other intervening lines = Red
    if outlist['number']>0:
        #print('There are Intervening Lines!')
        vellist=outlist['vel']
        relative_vel_z=outlist['delv']
        for index in range(0,outlist['number']):
            if (np.abs(vellist[index]) <delv):
                #print(vellist[index])

                if np.abs(relative_vel_z[index]) >200.:
                    color =clr['pale_red']
                else:
                    color='b'
                ax.text(vellist[index],1.05, np.str(outlist['ion'][index])+' '+ np.str('%.0f' % outlist['wrest'][index]),
                    fontsize=8,rotation=90, rotation_mode='anchor',color=color)
                ax.text(vellist[index]+50.,1.05, 'z = '+np.str('%.3f' % outlist['zobs'][index]),
                    fontsize=8,rotation=90, rotation_mode='anchor',color=color)



class mainWindow(QtWidgets.QTabWidget):
    
    def __init__(self,ions, parent=None,intervening=False):
        #-----full spectra properties---------#
        self.z = ions['Target']['z']; self.flux = ions['Target']['flux']
        self.wave = ions['Target']['wave']; self.error = ions['Target']['error']
        #------------------------------------#
        self.old_axes = None
        self.tab_names = ['Ions', 'Ions 2','Ions 3', 'Ions 4', 'Ions 5']
        self.ions = ions
        self.keys = list(self.ions.keys())[:-1] # last item is the full target spectrum
        self.wc = None #initialize overall parameter
        self.ylims = []
        self.vclim = None
        self.vclim_all = []
        self.name = None
        self.EWlim = [None,None] #left,right
        self.event_button = None
        self.Manual_Mask = None
        self.pFlag = 1
        self.intervening=intervening
        self.save = False
        super(mainWindow,self).__init__(parent)
        
#---------------Initial page setup------------------# 
        self.page = 0
        self.tabs = [QtWidgets.QWidget()]
        self.addTab(self.tabs[self.page],self.tab_names[self.page])
        self.figs = [Figure()]
        self.canvas = [FigureCanvasQTAgg(self.figs[self.page])]

        self.nions = [len(self.keys)]
        if self.nions[0] > 6: self.nions[0]=6
        self.page=0
        #Need to set click focus for keyboard functionality
        self.setParent(parent)
        self.figs[0].canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.figs[0].canvas.setFocus()
        
        
        
        #LAYOUT
        #Top layer is a horizontal layout with help|save|add_ion; middle is canvas; bottom is Page number
        self.add_ion_button = QPushButton("Add Ion",self)
        self.add_ion_button.setGeometry(630,30,200,30)
        self.add_ion_button.clicked.connect(lambda: self.NewTransition(self))
        
        self.openButton = QPushButton("Help",  self)
        self.openButton.setGeometry(830,30,200,30)
        self.openButton.clicked.connect(lambda: self.opensub(self))
        
        self.saveButton = QPushButton("Save",  self)
        self.saveButton.setGeometry(430,30,200,30)
        self.saveButton.clicked.connect(lambda: self.onsave(self))
        
        self.PageLabel = QtWidgets.QLabel("Page: " + str(self.currentIndex()+1)+"/"+str(len(self.figs)),self)
        self.PageLabel.setStyleSheet("font: 16pt;color: white;background-color:QColor(53, 53, 53)")
        
        self.main_layout = QVBoxLayout()
        self.top_layout = QHBoxLayout()
        self.bot_layout = QHBoxLayout()
        
        self.spacerItem = QtWidgets.QSpacerItem(5, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        
        self.top_layout.addItem(self.spacerItem)
        self.top_layout.addWidget(self.add_ion_button)
        self.top_layout.addWidget(self.saveButton)
        self.top_layout.addWidget(self.openButton)
        self.top_layout.addItem(self.spacerItem)
        
        self.bot_layout.addItem(self.spacerItem)
        self.bot_layout.addWidget(self.PageLabel)
        self.bot_layout.addItem(self.spacerItem)
        
        self.main_layout.addLayout(self.top_layout,stretch=1)
        self.main_layout.addWidget(self.canvas[0],stretch=8)
        self.main_layout.addLayout(self.bot_layout,stretch=0.5)
        self.tabs[self.page].setLayout(self.main_layout)
        
        
        #initializing left and right axes
        self.axesL = [list(range(6))]; self.axesR = [list(range(6))]
        for ii in range(self.nions[0]):
            self.axesL[0][ii] = self.figs[0].add_subplot(6,2,2*ii+1)
            self.axesR[0][ii] = self.figs[0].add_subplot(6,2,2*(ii+1))
            self.figs[self.page].subplots_adjust(hspace=0.01)
            Plotting(self,ii,modify=True)
            
         # Set up connectivity
        self.cid1 = self.figs[0].canvas.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.figs[0].canvas.mpl_connect("key_press_event", self.onpress)
        self.cid3 = self.figs[0].canvas.mpl_connect("motion_notify_event",self.onmotion)
        
#----------------Setup for Additional pages-----------#  
        AddPage(self)

                       
#---------------------Save Button/function----------------# 
        

#--------------------------------------------------------------#

#         self.PageLabel = QtWidgets.QLabel("Page: " + str(self.currentIndex()+1)+"/"+str(len(self.figs)),self)
#         self.PageLabel.setStyleSheet("font: 16pt;color: black;background-color:white")
#         self.PageLabel.setGeometry(630,850,200,30)
        
        
        
        def getPage(self):
            self.page = self.currentIndex()
            self.PageLabel.setText("Page: " + str(self.page+1)+"/"+str(len(self.figs)))
        self.currentChanged.connect(lambda: getPage(self))
            
#-------------------Add Ion Button------------------------------# 
    def NewTransition(self,parent):
    #will run back through absorber class to identify line, obtain slice window, 
        new_line,ok3 = QInputDialog.getDouble(self,'Add Line','Enter new transition:')#,decimals=4)
        if ok3:
            #add new line
            new_abs = Absorber.Absorber(self.z,self.wave,self.flux,self.error,[new_line])
            #update initial dictionary
            self.ions.update(new_abs.ions); self.ions['Target'] = self.ions.pop('Target')# moves Target to last index
            for key in list(new_abs.ions.keys())[:-1]: self.keys.append(key)

#                 self.page = len(self.nions) - 1 #finds current page max
            if self.nions[self.page] < 6: #proceed with filling page
                ii = self.nions[self.page]
                self.axesL[self.page][ii] = self.figs[self.page].add_subplot(6,2,2*(ii)+1)
                self.axesR[self.page][ii] = self.figs[self.page].add_subplot(6,2,2*(ii+1))
                self.figs[self.page].subplots_adjust(hspace=0.01)
                self.nions[self.page] = self.nions[self.page]+1
                Plotting(self,ii,modify=True)

            else: # need to add a new page
                AddPage(self)
                    
#         self.add_ion_button = QPushButton("Add Ion",self)
#         self.add_ion_button.setGeometry(630,30,200,30)
#         self.add_ion_button.clicked.connect(lambda: NewTransition(self))
        
#--------------------------------------------------------------------# 


#-------HELP------#
    def onsave(self,parent):
            self.savepg = SavePage(self)
            if self.savepg.closewin == None:
                return None
            else:
                self.savepg.show()
                
    def opensub(self,parent):
        self.sub = HelpWindow()
        self.sub.show()

        #html
        webbrowser.open_new_tab("C:\\Users\\Sean's PC\\Python Modules\\Help.html")
            
#         self.openButton = QPushButton("Help",  self)
#         self.openButton.setGeometry(830,30,200,30)
#         self.openButton.clicked.connect(lambda: opensub(self))

        
        
    def onmotion(self,event):
        #self.PageLabel.setText("Page: " + str(self.currentIndex()+1)+"/"+str(len(self.figs)))
        
        if event.xdata != None and event.ydata != None:
            #pdb.set_trace()
            
#             for qq in range(len(self.axesL)):
#                 if (event.inaxes in self.axesL[qq]) | (event.inaxes in self.axesR[qq]) :
#                     self.page = qq

            self.page = np.where((np.asarray(self.axesL)==event.inaxes)|(np.asarray(self.axesR)==event.inaxes))[0][0]
            if len(np.where(np.asarray(self.axesL[self.page]) == event.inaxes)[0]) == 0:
                self.Lidx = None
                self.Ridx = np.where(np.asarray(self.axesR[self.page])==event.inaxes)[0][0]
                
                if self.old_axes  and (self.old_axes != self.axesR[self.page][self.Ridx]):
                    for pos in ['top','bottom','left','right']:
                        self.old_axes.spines[pos].set_edgecolor('black')
                        self.old_axes.spines[pos].set_linewidth(0.5)
                    self.figs[self.page].canvas.draw()
                if self.old_axes != self.axesR[self.page][self.Ridx]:
                    for pos in ['top','bottom','left','right']:
                        self.axesR[self.page][self.Ridx].spines[pos].set_edgecolor('#01DF01')
                        self.axesR[self.page][self.Ridx].spines[pos].set_linewidth(2)
                    self.figs[self.page].canvas.draw()
                    self.old_axes = self.axesR[self.page][self.Ridx]
                
            else:
                self.Ridx = None
                self.Lidx = np.where(np.asarray(self.axesL[self.page]) == event.inaxes)[0][0]
                if (self.old_axes != None) and (self.old_axes != self.axesL[self.page][self.Lidx]):
                    for pos in ['top','bottom','left','right']:
                        self.old_axes.spines[pos].set_edgecolor('black')
                        self.old_axes.spines[pos].set_linewidth(0.5)
                        
                if self.old_axes != self.axesL[self.page][self.Lidx]:
                        for pos in ['top','bottom','left','right']:
                            self.axesL[self.page][self.Lidx].spines[pos].set_edgecolor('#01DF01')
                            self.axesL[self.page][self.Lidx].spines[pos].set_linewidth(2)
                        self.figs[self.page].canvas.draw()
                        self.old_axes = self.axesL[self.page][self.Lidx]



#----------------------key button events-----------------------------#            
    
    '''on press is used to reduce the order for the polyfit, toggle measurement displays, measure properties, and assist selecting EW vel bounds'''
    def onpress(self,event):
        if event.key == 'v':
            if self.old_axes in self.axesL[self.page]:
                mask_reg,ok = QInputDialog.getText(self,'Manual Mask Limits','Input Region to Mask (e.g. 200,250)')

                if ok:
                    key_idx = (self.page*6)+self.Lidx
                    vel = self.ions[self.keys[key_idx]]['vel']
                    wc = self.ions[self.keys[key_idx]]['wc']

                    mask = mask_reg.split(',')
                    mask = np.array(mask).astype('float32')

                    wc=((vel<mask[0]) | (vel>mask[1])) & wc
                    self.ions[self.keys[key_idx]]['wc'] = wc
                    Plotting(self,self.Lidx,modify=True)
                    
            elif self.old_axes in self.axesR[self.page]:
                integ_lims,ok = QInputDialog.getText(self,'Manual EW Limits','Input integration region (eg -100,100)')
                if ok:
                    key_idx = (self.page*6)+self.Ridx
                    integ_lims = integ_lims.split(',')
                    integ_lims = np.array(integ_lims).astype('float32')

                    self.ions[self.keys[key_idx]]['EWlims'][0] = integ_lims[0]
                    self.ions[self.keys[key_idx]]['EWlims'][1] = integ_lims[1]
                    Plotting(self,self.Ridx,modify=False,Print=False)

        
        #right arrow key (directional) increases poly order 
        if event.key == 'up':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx
                
                self.ions[self.keys[key_idx]]['order'] = self.ions[self.keys[key_idx]]['order']+1
                Plotting(self,self.Lidx,modify=True)
            else:
                print('click on a left transition window first')
                
        #reduce polynomial        
        if event.key == 'down':
            if self.Lidx is not None:
                key_idx = self.page*6+self.Lidx
#                 if self.page == 0:
                self.ions[self.keys[key_idx]]['order'] = self.ions[self.keys[key_idx]]['order']-1
                Plotting(self,self.Lidx,modify=True)
            else:
                print('click on a transition window first')
                
        #evaluate last limit subplot
        if event.key == 'm':
            key_idx = self.page*6+self.Ridx
            EW(self,self.page,self.Ridx,self.ions[self.keys[key_idx]]['EWlims'])

                
        #evaluate all ions with limits on page
        if event.key == 'M':
            for jj in range(len(self.figs)):
                self.page = jj
                for ii in range(self.nions[self.page]):
                    EW(self,self.page,ii,self.ions[self.keys[ii+self.page*6]]['EWlims'])

            
        
        #use same range for all ion Velocity limits and directly measured by following with clicks bounds on a single subplot
        if event.key == 'V':
            EWlims = self.ions[self.keys[self.Ridx+6*self.page]]['EWlims']
            for jj in range(len(self.figs)):
                self.page = jj
                for ii in range(self.nions[self.page]):

                    self.ions[self.keys[ii+self.page*6]]['EWlims'] = EWlims
                    Plotting(self,ii,modify=False,Print=False)

        if event.key in ['0','1','2']: #detection, upperlimit, lower limit
            key_idx = self.page*6 +self.Ridx
            self.ions[self.keys[key_idx]]['flag'] = int(event.key)
            Plotting(self,self.Ridx,modify=False,Print=True)
            
        if event.key == 't':#toggle the EW/N display
            self.pFlag = (self.pFlag+1)%3
            for jj in range(len(self.figs)):
                self.page = jj
                for ii in range(self.nions[self.page]):
                    Plotting(self,ii,modify=False,Print=True)
            
    

        
#------------------------------click button events----------------------------#        
            
    def onclick(self,event):              
         
        if event.button in [1,3]:
            '''Left hand side is for fitting'''
                        
                # if not first canvas, +6 per pae is needed to appropraitely index the ions
            if self.Lidx is not None:

                key_idx = (self.page*6)+self.Lidx
                vel = self.ions[self.keys[key_idx]]['vel']
                name = self.ions[self.keys[key_idx]]['name']

                if self.vclim is None:
                    self.vclim = [event.xdata]
                    self.axesL[self.page][self.Lidx].plot(event.xdata,event.ydata,'ro',ms=5)
                    if self.page == 0: self.figs[self.page].canvas.draw()
                    else: self.figs[self.page].canvas.draw()


                else:
                    vclim = np.sort(np.append(self.vclim,event.xdata))
                    self.vclim = None
                    wc = self.ions[self.keys[key_idx]]['wc']

                    if event.button == 1:#if left click add to masks
                        wc=((vel<vclim[0]) | (vel>vclim[1])) & wc
                    else: #remove mask
                        wc=((vel>vclim[0]) & (vel<vclim[1])) | wc


                    #update wc for plotting
                    self.ions[self.keys[key_idx]]['wc'] = wc

                    #replot
                    Plotting(self,self.Lidx,modify=True)
                        
                        
                    
                '''Right hand side for picking velocity limits for EW measurements'''
            check = []; 
                        
                
            if self.Ridx is not None:
                key_idx = (self.page*6)+self.Ridx

                    #if left click then define leftward vel limit
                if event.button == 1:
                    self.EWlim[0] = event.xdata#used for plotting all with same range 'V' command
                    self.ions[self.keys[key_idx]]['EWlims'][0] = event.xdata
                    #plot selected limits
                    Plotting(self,self.Ridx,modify=False,Print=False)

                    #if right click define rightward vel limit
                elif event.button == 3:
                    self.EWlim[1] = event.xdata
                    self.ions[self.keys[key_idx]]['EWlims'][1] = event.xdata
                    Plotting(self,self.Ridx,modify=False,Print=False)

                           

class Plotting:
    def __init__(self,parent,ii,modify=False,Print=False):
        
        # following if statement is to determine which canvas it is on for proper indexing (page1/page2)
        key_idx = ii+6*parent.page
        
        #---------------Define variables for readability--------------#
        vel = parent.ions[parent.keys[key_idx]]['vel']
        wave = parent.ions[parent.keys[key_idx]]['wave']
        error = parent.ions[parent.keys[key_idx]]['error']
        flux = parent.ions[parent.keys[key_idx]]['flux']
        weight = parent.ions[parent.keys[key_idx]]['weight']
        name = parent.ions[parent.keys[key_idx]]['name']
        #L = parent.ions[parent.keys[key_idx]].L
        wc =np.array(parent.ions[parent.keys[key_idx]]['wc']&(error != 0)) #error != 0 is a bad pixel mask 
        cont = parent.ions[parent.keys[key_idx]]['cont']
        window_lim = parent.ions[parent.keys[key_idx]]['window_lim']
        order = parent.ions[parent.keys[key_idx]]['order']
        EWlims = parent.ions[parent.keys[key_idx]]['EWlims']
        lam_0=parent.ions[parent.keys[key_idx]]['lam_0']
        fvals=parent.ions[parent.keys[key_idx]]['f']


        #--------------------------------------------------------------#
        if Print == False:
            #if modify == True, adjustments on LHS(continuum) are needed otherwise, only replot the RHS
            if modify == True: 
                #re-eval continuum
                parent.ions[parent.keys[key_idx]]['pco']=L.Legendre.fit(wave[wc],flux[wc],order,w=weight[wc])
                parent.ions[parent.keys[key_idx]]['cont'] = parent.ions[parent.keys[key_idx]]['pco'](wave)
                cont = parent.ions[parent.keys[key_idx]]['cont']
                
                #gray_idx is to avoid line plotted through spectra from discontinuity in flux/vel/err from wc mask.
                #uses np.diff to find where wc (true/false array) changes and plots the line segments instead of the full line with missing values
                if wc[0] == True:
                    gray_idx = np.where(np.diff(wc,prepend=np.nan))[0][1:]
                else:
                    gray_idx = np.where(np.diff(wc,prepend=np.nan))[0]
                    
                parent.axesL[parent.page][ii].clear()
                parent.axesL[parent.page][ii].step(vel,flux,color='k',where='mid');parent.axesL[parent.page][ii].step(vel,error,color='r',where='mid')
                parent.axesL[parent.page][ii].step(vel,cont,color='b',where='mid')
                for zz in range(int(len(gray_idx)/2)):
                    #gray_idx = gray_idx-1
                    parent.axesL[parent.page][ii].step(vel[gray_idx[zz*2]:gray_idx[2*zz+1]],flux[gray_idx[2*zz]:gray_idx[2*zz+1]],vel[gray_idx[2*zz]:gray_idx[2*zz+1]],error[gray_idx[2*zz]:gray_idx[2*zz+1]],where='mid',color='lightgray',linewidth=1.3,alpha=1)

            #clear axes to redraw modifications
            parent.axesR[parent.page][ii].clear()


            #replot right (flux; error)
            parent.axesR[parent.page][ii].step(vel,flux/cont,color='k',where='mid');parent.axesR[parent.page][ii].step(vel,error/cont,color='r',where='mid')
            parent.axesR[parent.page][ii].axhline(y=1,xmin=window_lim[0],xmax=window_lim[1],ls='--',c='b')
            

            #clear y ticks and label plots
            parent.axesL[parent.page][ii].set_yticks([]); #parent.axesR[i][ii].set_yticks([])
            parent.axesL[parent.page][ii].set_ylabel(name); parent.axesR[parent.page][ii].set_ylabel(name)


            #set axes bounds
            parent.axesL[parent.page][ii].set_xlim(window_lim); parent.axesR[parent.page][ii].set_xlim(window_lim)
            parent.axesR[parent.page][ii].set_ylim([0,2.2])

            #set x ticks only on bottom row
            if ii != parent.nions[parent.page] - 1:
                parent.axesL[parent.page][ii].set_xticks([]);parent.axesR[parent.page][ii].set_xticks([])
            else:
                if ii > 0: #ensures if a new subplot is added the xticks will be removed
                    parent.axesL[parent.page][ii-1].set_xticks([]);parent.axesR[parent.page][ii-1].set_xticks([])
                parent.axesL[parent.page][ii].set_xlabel('Velocity (km/s)')
                parent.axesR[parent.page][ii].set_xlabel('Velocity (km/s)')
            #plot column titles
            if ii == 0:
                parent.axesL[parent.page][0].set_title('Continuum Fitter')
                parent.axesR[parent.page][0].set_title('Normalized Spectra')

            #plot vel = 0 line
            parent.axesL[parent.page][ii].axvline(0,ymin=0,ymax=parent.axesL[parent.page][ii].get_ylim()[1],ls='--',c='k')
            parent.axesR[parent.page][ii].axvline(0,ymin=0,ymax=parent.axesR[parent.page][ii].get_ylim()[1],ls='--',c='k')

            #Plot EW velocity limit
            if EWlims[0] is not None: parent.axesR[parent.page][ii].axvline(EWlims[0],ymin=0,ymax=2.5,ls='--',c='b')
            if EWlims[1] is not None: parent.axesR[parent.page][ii].axvline(EWlims[1],ymin=0,ymax=2.5,ls='--',c='r')

            #### NOW PLOT INTERVENING LINES IF given
            if parent.intervening != False:
                outlist=grab_intervening_linelist(parent.intervening,np.double(parent.z),lam_0,wave)
                plot_intervening_lines( parent.axesR[parent.page][ii],outlist,np.max(vel))

            #plot fvals
            #xloc = parent.axesR[parent.page][ii].get_xlim()[1]
            #yloc = parent.axesR[parent.page][ii].get_ylim()[1]
            parent.axesR[parent.page][ii].text(0.85,0.85,'f: '+np.str('%.4f' % fvals),transform=parent.axesR[parent.page][ii].transAxes,color=clr['teal'])
            
            #redraw MUST BE LAST ITEM IN LIST
            parent.figs[parent.page].canvas.draw()

        #plot EW measurements on frame
        if Print == True:
            if parent.ions[parent.keys[key_idx]]['EW_text'] is not None:
                parent.ions[parent.keys[key_idx]]['EW_text'].remove()

            plotText(parent,parent.ions[parent.keys[key_idx]])
            text = parent.ions[parent.keys[key_idx]]['text']
            xlim = parent.ions[parent.keys[key_idx]]['EWlims'][0]
            EWtoggle = parent.axesR[parent.page][ii].text(.05,0.85,text, transform=parent.axesR[parent.page][ii].transAxes)
            parent.ions[parent.keys[key_idx]]['EW_text'] = EWtoggle
            parent.figs[parent.page].canvas.draw()
            
                
        
                
            
#Calculates N,EW,V,          
class EW:
    def __init__(self,parent,page,ii,lims):

        #determine which page is being accessed
        key_idx = ii+6*parent.page
            
        #define variables for readability
        vel = parent.ions[parent.keys[key_idx]]['vel']
        wave = parent.ions[parent.keys[key_idx]]['wave']
        error = parent.ions[parent.keys[key_idx]]['error']
        flux = parent.ions[parent.keys[key_idx]]['flux']
        name = parent.ions[parent.keys[key_idx]]['name']
        zabs = parent.ions[parent.keys[key_idx]]['z']
        f0 = parent.ions[parent.keys[key_idx]]['f']
        lam_0 = parent.ions[parent.keys[key_idx]]['lam_0']
        cont = parent.ions[parent.keys[key_idx]]['cont']

        #compute EW/N/med_vel
        output = compute_EW.compute_EW(wave,flux/cont,lam_0,lims,error/cont,plot=False,zabs=zabs,f0=f0)
        #save variables in ion's respective dictionary
        parent.ions[parent.keys[key_idx]]['N'] = output['col']
        parent.ions[parent.keys[key_idx]]['Nsig'] = output['colerr']
        parent.ions[parent.keys[key_idx]]['EW'] = output['ew_tot']*1000
        parent.ions[parent.keys[key_idx]]['EWsig'] = output['err_ew_tot']*1000
        parent.ions[parent.keys[key_idx]]['med_vel'] = output['med_vel']
        parent.ions[parent.keys[key_idx]]['EWlims'] = lims
        
        #plot EW on page
        Plotting(parent,ii,modify=False,Print=True)

        
# plots measurements, with toggle will display EW/N 
class plotText:
    def __init__(self,parent,line):
        EW_det_text= np.str('%.0f' % line['EW']) + ' $\pm$ ' + np.str('%.0f' % line['EWsig']) + ' m$\AA$'
        EW_limit_text="<{:.0f} m$\AA$".format(2.*line['EWsig']) #+  ' m$\AA$'
        logN_det_text= np.str('%.2f' % np.log10(line['N'])) +' $\pm$ ' + np.str('%.3f' % (np.log10(line['N']+line['Nsig']) - np.log10(line['N']))) + ' /cm$^2$'

        #line.flag is the line specific upper/lower/detections
        #pflag is the toggle button for which to show
        # Measurement: 
        if line['flag']==0:
            if parent.pFlag==0: text=""
            elif parent.pFlag==1: text=EW_det_text
            elif parent.pFlag==2: text=logN_det_text

        # Upper Limit
        elif line['flag']==1:
            if parent.pFlag==0: text=""
            elif parent.pFlag==1: text=EW_limit_text
            elif parent.pFlag==2: text="<{:.2f}".format(np.log10(line['Nsig']))

        elif line['flag']==2:
            if parent.pFlag==0: text=""
            elif parent.pFlag==1: text=EW_det_text
            elif parent.pFlag==2: text=">{:.2f}".format(np.log10(line['N']))
        line['text'] = text
        
class AddPage:
    
    def __init__(self,parent,):
        #what happens initiall if you run it and ionlist == 18?
        #and need a condition to not run if page has already been initialized
        nall = len(parent.keys)

        if ((nall>6)&(len(parent.figs)<2)):
            parent.page = 1 
            
            initialize(parent)
            
            parent.cid4 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid5 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid6 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
        elif ((nall>12)&(len(parent.figs)<3)):
            parent.page = 2
            initialize(parent)
            
            parent.cid7 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid8 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid9 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
        elif ((nall>18)&(len(parent.figs)<4)):
            parent.page = 3
            
            initialize(parent)
            
            parent.cid10 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid11 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid12 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
        elif ((nall>24)&(len(parent.figs)<5)):
            parent.page = 4
            
            self.initialize(parent)
            
            parent.cid13 = parent.figs[parent.page].canvas.mpl_connect("button_press_event", parent.onclick)
            parent.cid14 = parent.figs[parent.page].canvas.mpl_connect("key_press_event", parent.onpress)
            parent.cid15 = parent.figs[parent.page].canvas.mpl_connect("motion_notify_event", parent.onmotion)
            
class initialize:
    def __init__(self,parent):
        parent.tabs.append(QtWidgets.QWidget())
        parent.addTab(parent.tabs[parent.page], parent.tab_names[parent.page])
        parent.figs.append(Figure())
        parent.canvas.append(FigureCanvasQTAgg(parent.figs[parent.page]))
        layout = QtWidgets.QVBoxLayout()
#         layout.addWidget(parent.canvas[parent.page])
#         parent.tabs[parent.page].setLayout(layout)
        
        self.add_ion_button = QPushButton("Add Ion",parent)
        self.add_ion_button.setGeometry(630,30,200,30)
        self.add_ion_button.clicked.connect(lambda: parent.NewTransition(parent))
        
        self.openButton = QPushButton("Help",  parent)
        self.openButton.setGeometry(830,30,200,30)
        self.openButton.clicked.connect(lambda: parent.opensub(parent))
        
        self.saveButton = QPushButton("Save",  parent)
        self.saveButton.setGeometry(430,30,200,30)
        self.saveButton.clicked.connect(lambda: parent.onsave(parent))
        
        self.PageLabel = QtWidgets.QLabel("Page: " + str(parent.page+1)+"/"+str(len(parent.figs)),parent)
        self.PageLabel.setStyleSheet("font: 16pt;color: white;background-color:QColor(53, 53, 53)")
        
        self.main_layout = QVBoxLayout()
        self.top_layout = QHBoxLayout()
        self.bot_layout = QHBoxLayout()
        
        self.spacerItem = QtWidgets.QSpacerItem(5, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.top_layout.addItem(self.spacerItem)
        self.top_layout.addWidget(self.add_ion_button)
        self.top_layout.addWidget(self.saveButton)
        self.top_layout.addWidget(self.openButton)
        self.top_layout.addItem(self.spacerItem)
        
        self.bot_layout.addItem(self.spacerItem)
        self.bot_layout.addWidget(self.PageLabel)
        self.bot_layout.addItem(self.spacerItem)
                                   
        self.main_layout.addLayout(self.top_layout,stretch=1)
        self.main_layout.addWidget(parent.canvas[parent.page],stretch=7)
        self.main_layout.addLayout(self.bot_layout,stretch=0.5)
        parent.tabs[1].setLayout(self.main_layout)
        
        
        #len(parent.keys/6 -len(parent.page)) #can you make more general? you ave to
        parent.nions.append(len(parent.keys)-6*parent.page)
        
        #initializing left and right axes
        parent.axesL.append(list(range(6))); parent.axesR.append(list(range(6)))                   
        for ii in range(parent.nions[parent.page]):
            parent.axesL[parent.page][ii] = parent.figs[parent.page].add_subplot(6,2,2*ii+1)
            parent.axesR[parent.page][ii] = parent.figs[parent.page].add_subplot(6,2,2*(ii+1))
            parent.figs[parent.page].subplots_adjust(hspace=0.01)
            Plotting(parent,ii,modify=True)


        # Set up connectivity and Need to set click focus for keyboard functionality
        parent.figs[parent.page].canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        parent.figs[parent.page].canvas.setFocus()            



class HelpWindow(QtWidgets.QWidget):

    def __init__(self,parent=None):
        super(HelpWindow, self).__init__(parent)
        self.resize(400,500)
        label = QtWidgets.QLabel(HELP,self)
        
class SavePage(QtWidgets.QWidget):
    def __init__(self,parentvals,parent=None):
        super(SavePage, self).__init__(parent)
        self.resize(700,400)
        
        def onpdf(self,parentvals):
            # need to eliminate highlighted axes
            for pos in ['top','bottom','left','right']:
                parentvals.old_axes.spines[pos].set_edgecolor('black')
                parentvals.old_axes.spines[pos].set_linewidth(0.5)
            figurefile = self.pdfline.text()
            i = 1
            for figures in parentvals.figs:
                figures.savefig(figurefile[:-4]+str(i)+'.pdf',bbox_inches='tight')
                i = i+1
            self.pdfsave.setStyleSheet('background-color : green')
            parentvals.save = True
        def ontable(self,parentvals,Table_e):
            savefile = self.tableline.text()
            ascii.write(Table_e,savefile,overwrite=True)
            #how to display to user it has been saved?
            self.tablesave.setStyleSheet("background-color : green")
            parentvals.save = True
            
        def onpickle(self,parentvals):
            pfile = self.pickleline.text()
            with open(pfile,'wb') as pklfile:
                pickle.dump(parentvals.ions,pklfile,protocol=pickle.HIGHEST_PROTOCOL)
            self.picklesave.setStyleSheet('background-color : green')
            parentvals.save = True
        # want to print out table to the user regardless of whether it is saved
        Table_e = Table()
        Table_e['Transitions'] = parentvals.keys
        EW=[];EWsig=[];N=[];Nsig=[];Vel = []
        EWlims_low = []; EWlims_high = []
        for ion in parentvals.keys:
            for items in parentvals.ions[ion]:
                if np.any(parentvals.ions[ion][items] == None):
                    reply = QMessageBox.question(self, 'Message', "Unevaluated ions, proceed to save?",
                                            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                    if reply == QMessageBox.Yes:
                        for ions in parentvals.keys:
                            for item in parentvals.ions[ions]:
                                if np.any(parentvals.ions[ions][item] == None):
                                    parentvals.ions[ions][item] = np.NaN
                    else:
                        self.closewin = None
                        return self.closewin
                else:
                    self.closewin = False
            EW.append(np.round(parentvals.ions[ion]['EW'],2))
            EWsig.append(np.round(parentvals.ions[ion]['EWsig'],2))
            N.append(np.round(np.log10(parentvals.ions[ion]['N']),2))
            Nsig.append(np.round(np.log10(parentvals.ions[ion]['Nsig']),2))
            Vel.append(np.round(parentvals.ions[ion]['med_vel'],2))
            EWlims_low.append(np.round(parentvals.ions[ion]['EWlims'][0],2))
            EWlims_high.append(np.round(parentvals.ions[ion]['EWlims'][1],2))
        Table_e['EW'] = EW; Table_e['EWsig'] = EWsig; Table_e['Vmin']=EWlims_low
        Table_e['Vmax'] = EWlims_high
        Table_e['N'] = N; Table_e['Nsig'] = Nsig; Table_e['Vel']=Vel
        print(Table_e)
        
        #pdf save
        pdflabel = QLabel("Enter path and filename: (e.g. pathname\Ions.pdf)",self)
        pdflabel.setGeometry(100,100,400,30)
        
        self.pdfline = QLineEdit(self)
        self.pdfline.setText('Spectrum_Analysis_z_'+str(parentvals.z)+'_Ions.pdf')
        self.pdfline.setGeometry(100,125,300,30)
        
        self.pdfsave = QPushButton("save PDF",self)
        self.pdfsave.setGeometry(410,125,200,30)
        self.pdfsave.clicked.connect(lambda: onpdf(self,parentvals))
        #table save
        tablelabel = QLabel("Enter path and filename: (e.g. pathname\Table.dat)",self)
        tablelabel.setGeometry(100,175,400,30)
        
        self.tableline = QLineEdit(self)
        self.tableline.setText("Spectrum_Analysis_z_"+str(parentvals.z)+"_Measurement_Table.dat")
        self.tableline.setGeometry(100,200,300,30)
        
        self.tablesave = QPushButton("Save Table",self)
        self.tablesave.setGeometry(410,200,200,30)
        self.tablesave.clicked.connect(lambda: ontable(self,parentvals,Table_e))
        
        #pickle save
        picklelabel = QLabel("Enter path and filename: (e.g. pathname\Table.p)",self)
        picklelabel.setGeometry(100,250,400,30)
        
        self.pickleline = QLineEdit(self)
        self.pickleline.setText("Spectrum_Analysis_z_"+str(parentvals.z)+".p")
        self.pickleline.setGeometry(100,275,300,30)
        
        self.picklesave = QPushButton("Save Progress",self)
        self.picklesave.setGeometry(410,275,200,30)
        self.picklesave.clicked.connect(lambda: onpickle(self,parentvals))

        
#Initial inputs and callable class to run proram        
class Transitions:
    def __init__(self,Abs,intervening=False):
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
        #app.setStyle("Fusion")

        # Now use a palette to switch to dark colors:
        #palette = QPalette()
        #palette.setColor(QPalette.Window, QColor(53, 53, 53))
        #palette.setColor(QPalette.WindowText, QtCore.Qt.white)        
        #palette.setColor(QPalette.Base, QColor(25, 25, 25))
        #palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        #palette.setColor(QPalette.Button, QColor(53, 53, 53))
        #palette.setColor(QPalette.ButtonText, QtCore.Qt.white)
        #palette.setColor(QPalette.BrightText, QtCore.Qt.red)
        #palette.setColor(QPalette.Link, QColor(42, 130, 218))
        #palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        #palette.setColor(QPalette.Text, QtCore.Qt.white)

        #app.setPalette(palette)
        main = mainWindow(Abs,intervening=intervening)
        main.resize(1400,900)

        main.show()
        #app.exec_()

        QtWidgets.QApplication.setQuitOnLastWindowClosed(True)
        app.exec_()
        app.quit()
        #sys.exit()
