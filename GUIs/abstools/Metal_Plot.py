import numpy as np
from IGM import compute_EW
from PyQt5 import QtCore, QtGui, QtWidgets
from GUIs.abstools import Absorber
from PyQt5.QtWidgets import QStyleFactory, QPushButton,QLineEdit,QMainWindow,QInputDialog
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
clr=rt.rb_set_color()

HELP = '''
            The LHS shows the spectrum with a Legengre polynomial continuum fit overlaid.
                Points excluded from continuum fit are gray.
            The RHS shows the normalized spectrum and velocity limits for measurements.

            Mouse Clicks
            LHS LMB: Add wavelengths within region set by two clicks to continuum fit.
            LHS RMB: Remove wavelengths from continuum fit.
            RHS LMB: Set lower velocity limit
            RHS RMB: Set upper velocity limit

            Keyboard Commands
            Q: quit
            left arrow: Increase Polynomial Order [default 4]
            right arrow: Decrease Polynomial Order [default 4]
                --uses the last modified/ selected LHS spectra
                ---if no modifications, click on the spectra of which the order of fit should be changed

            V: Initializes for all transitions to use the same velocity limits for integration
                --Must be followed by clicking RHS on any spectra to define the velocity limits
            m: Measure EW/N for current subplot axes
            M: If all subplot axes on page have defined limits, measures EW/N for all ions on page

            1/2/0 : While on RHS, flag absorber as (1) upper limit , (2) lower limit, or (0) neither.
            t : Cycle text printed on absorbers. Displays logN, or EW
                    Display logN and EW as detection or limit based on flag above


            If Specified it can plot any intervening absorber list.
                    '''



def shift2vel(z1,z2):
    # z1 at rest
    # z2 for which relative velocity is computed
    spl=2.9979e5;  #speed of light
    vel=((1215.67*(1.+z2))-(1215.67*(1.0 + z1)))*spl/(1215.67*(1.0 + z1))#(self.wrest-str['wave']*(1.0 + 0.))*spl/(str['wave']*(1.0 + 0.))
    return vel

def grab_intervening_linelist(filename,z_gal,wrest_galaxy,wavelength):
    s=ascii.read(filename)
    ion=s['col1']
    wrest=s['col2']
    wobs=s['col3']
    zobs=s['col4']
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
        print('There are Intervening Lines!')
        vellist=outlist['vel']
        relative_vel_z=outlist['delv']
        for index in range(0,outlist['number']):
            if (np.abs(vellist[index]) <delv):
                print(vellist[index])

                if np.abs(relative_vel_z[index]) >200.:
                    color =clr['pale_red']
                else:
                    color='b'
                ax.text(vellist[index],1.05, np.str(outlist['ion'][index])+' '+ np.str(outlist['wrest'][index]),
                    fontsize=4,rotation=90, rotation_mode='anchor',color=color)
                ax.text(vellist[index]+25.,1.05, 'z = '+np.str('%.3f' % outlist['zobs'][index]),
                    fontsize=4,rotation=90, rotation_mode='anchor',color=color)



class mainWindow(QtWidgets.QTabWidget):
    
    def __init__(self,ions, parent=None,intervening=False):
    
        #-----full spectra properties---------#
        self.z = ions['Target']['z']; self.flux = ions['Target']['flux']
        self.wave = ions['Target']['wave']; self.error = ions['Target']['error']
        #------------------------------------#
        
        self.ions = ions
        self.keys = list(self.ions.keys())[:-1] # last item is the full target spectrum
        self.wc = None #initialize overall parameter
        self.ylims = []
        self.vclim = None
        self.vclim_all = []
        self.name = None
        self.EWlim = [None,None] #left,right
        self.event_button = None
        self.pFlag = 1
        self.intervening=intervening
        super(mainWindow,self).__init__(parent)
        
#---------------Initial page setup------------------#        
        self.tab1 = QtWidgets.QWidget()
        self.addTab(self.tab1,'First Ions')
        self.fig = Figure()
        self.canvas = FigureCanvasQTAgg(self.fig)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        self.tab1.setLayout(layout)
        self.nions = len(self.keys)
        if self.nions > 6: self.nions=6
        self.page=0
        #Need to set click focus for keyboard functionality
        self.setParent(parent)
        self.fig.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
        self.fig.canvas.setFocus()
        
        #initializing left and right axes
        self.axesL = [list(range(6)),list(range(6))]; self.axesR = [list(range(6)),list(range(6))]
        for ii in range(self.nions):
            #line = Abs.ions[self.keys[ii]]
            self.axesL[0][ii] = self.fig.add_subplot(6,2,2*ii+1)
            self.axesR[0][ii] = self.fig.add_subplot(6,2,2*(ii+1))
            self.fig.subplots_adjust(hspace=0.01)
            Plotting(self,ii)
            
         # Set up connectivity
        self.cid1 = self.fig.canvas.mpl_connect("button_press_event", self.onclick)
        self.cid2 = self.fig.canvas.mpl_connect("key_press_event", self.onpress)
        
#----------------Setup for second page-----------#        
        if len(self.keys) > 6:
            self.page=1
            self.fig2 = Figure()
            self.tab2 = QtWidgets.QWidget()
            self.addTab(self.tab2,'Additional Ions')
            self.canvas2 = FigureCanvasQTAgg(self.fig2)
            layout2 = QtWidgets.QVBoxLayout()
            layout2.addWidget(self.canvas2)
            self.tab2.setLayout(layout2)
            self.nions =6
            self.nions2 = len(self.keys)-6
            self.fig2.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
            self.fig2.canvas.setFocus()
            for ii in range(self.nions2):
                self.axesL[1][ii] = self.fig2.add_subplot(6,2,2*ii+1)
                self.axesR[1][ii] = self.fig2.add_subplot(6,2,2*(ii+1))
                self.fig2.subplots_adjust(hspace=0.01)
                Plotting(self,ii)
                
            self.fig2.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
            self.fig2.canvas.setFocus()
            self.cid3 = self.fig2.canvas.mpl_connect("button_press_event", self.onclick)
            self.cid4 = self.fig2.canvas.mpl_connect("key_press_event", self.onpress)
                       
#---------------------Save Button/function----------------# 
        def save(self):
            from astropy.table import Table
            from astropy.io import ascii
            Table_e = Table()
            Table_e['Transitions'] = self.keys
            EW=[];EWsig=[];N=[];Nsig=[];Vel = []
            EWlims_low = []; EWlims_high = []
            for ion in self.keys:
                EW.append(np.round(self.ions[ion]['EW']*1000,2))
                EWsig.append(np.round(self.ions[ion]['EWsig']*1000,2))
                N.append(np.round(np.log10(self.ions[ion]['N']),2))
                Nsig.append(np.round(np.log10(self.ions[ion]['Nsig']),2))
                Vel.append(np.round(self.ions[ion]['med_vel'],2))
                EWlims_low.append(np.round(self.ions[ion]['EWlims'][0],2))
                EWlims_high.append(np.round(self.ions[ion]['EWlims'][1],2))
            Table_e['EW'] = EW; Table_e['EWsig'] = EWsig; Table_e['Lower Limit']=EWlims_low
            Table_e['Upper Limit'] = EWlims_high
            Table_e['N'] = N; Table_e['Nsig'] = Nsig; Table_e['Vel']=Vel
         
            pfile,ok = QInputDialog.getText(self,'Input Pickle Path','Input path to save pickle file: ')
            Tablefile,ok2 = QInputDialog.getText(self,'Save Values','Input path to save tabulated measurements (.dat): ')
            if ok:
                #saved = fileDat(self,self.ions)
                with open(pfile,'wb') as pklfile:
                    pickle.dump(self.ions,pklfile,protocol=pickle.HIGHEST_PROTOCOL)
            else:
                print('Pickle File Not Saved. Data Will Be Lost')
            if ok2:
                ascii.write(Table_e, Tablefile, overwrite=True)
                print(Table_e)
            else:
                print(Table_e)

        button = QPushButton("SAVE",self)
        button.setGeometry(500,30,200,30)
        button.clicked.connect(lambda: save(self))

#--------------------------------------------------------------#
            
#-------------------Add Ion Button------------------------------# 
        def NewTransition(self):
        #will run back through absorber class to identify line, obtain slice window, 
            new_line,ok3 = QInputDialog.getDouble(self,'Add Line','Enter new transition:')#,decimals=4)
            if ok3:
                #add new line
                new_abs = Absorber.Absorber(self.z,self.wave,self.flux,self.error,[new_line])
                #update initial dictionary
                self.ions.update(new_abs.ions); self.ions['Target'] = self.ions.pop('Target')# moves Target to last index

                if self.nions < 6:
                    self.page=0; ii = self.nions; 
                    for key in list(new_abs.ions.keys())[:-1]: self.keys.append(key)
                    self.axesL[0][ii] = self.fig.add_subplot(6,2,2*(ii)+1)
                    self.axesR[0][ii] = self.fig.add_subplot(6,2,2*(ii+1))
                    self.fig.subplots_adjust(hspace=0.01)
                    self.nions = self.nions+1
                    Plotting(self,ii)
                elif self.nions2<6:
                    self.page=1; ii = self.nions2; 
                    for key in list(new_abs.ions.keys())[:-1]: self.keys.append(key)
                    self.axesL[1][ii] = self.fig2.add_subplot(6,2,2*ii+1)
                    self.axesR[1][ii] = self.fig2.add_subplot(6,2,2*(ii+1))
                    self.fig2.subplots_adjust(hspace=0.01)
                    self.nions2 = self.nions2+1
                    Plotting(self,ii)
                    

        
        add_ion_button = QPushButton("Add Ion",self)
        add_ion_button.setGeometry(700,30,200,30)
        add_ion_button.clicked.connect(lambda: NewTransition(self))

#--------------------------------------------------------------------#  


#----------------------key button events-----------------------------#            
        
    '''on press is used to reduce the order for the polyfit, toggle measurement displays, measure properties, and assist selecting EW vel bounds'''
    
    def onpress(self,event):
        #right arrow key (directional) increases poly order 
        if event.key == 'right':
            if self.Lidx is not None: 
                if self.page == 0:
                    self.ions[self.keys[self.Lidx]]['order'] = self.ions[self.keys[self.Lidx]]['order']+1
                    Plotting(self,self.Lidx,modify=True)
                elif self.page == 1:
                    self.ions[self.keys[self.Lidx+6]]['order'] = self.ions[self.keys[self.Lidx+6]]['order']+1
                    Plotting(self,self.Lidx,modify=True)
            else:
                print('click on a left transition window first')
                
        #reduce polynomial        
        if event.key == 'left':
            if self.Lidx is not None: #need a "try" statement?
                if self.page == 0:
                    self.ions[self.keys[self.Lidx]]['order'] = self.ions[self.keys[self.Lidx]]['order']-1
                    Plotting(self,self.Lidx,modify=True)
                elif self.page == 1:
                    self.ions[self.keys[self.Lidx+6]]['order'] = self.ions[self.keys[self.Lidx+6]]['order']-1
                    Plotting(self,self.Lidx,modify=True)
            else:
                print('click on a transition window first')
                
        #evaluate last limit subplot
        if event.key == 'm':
            print(self.page)
            if self.page == 0:
                EW(self,self.page,self.Ridx,self.ions[self.keys[self.Ridx]]['EWlims'])
            elif self.page == 1:
                EW(self,self.page,self.Ridx,self.ions[self.keys[self.Ridx+6]]['EWlims'])
                
        #evaluate all ions with limits on page
        if event.key == 'M':
            print(self.page)
            if self.page == 0:
                for ii in range(self.nions):
                    EW(self,self.page,ii,self.ions[self.keys[ii]]['EWlims'])
            elif self.page == 1:
                for ii in range(self.nions2):
                    EW(self,self.page,ii,self.ions[self.keys[ii+6]]['EWlims'])
            
        
        #use same range for all ion Velocity limits and directly measured by following with clicks bounds on a single subplot
        if event.key == 'V':
            self.event_button = event.key
        if event.key in ['0','1','2']: #detection, upperlimit, lower limit
            if self.page == 0: key_idx = self.Ridx
            if self.page == 1: key_idx = self.Ridx+6
            self.ions[self.keys[key_idx]]['flag'] = int(event.key)
            Plotting(self,self.Ridx,modify=False,Print=True)
            
        if event.key == 't':#toggle the EW/N display
            self.pFlag = (self.pFlag+1)%3
            if self.page==0: loop = self.nions
            if self.page==1: loop = self.nions2
            for ii in range(loop):
                Plotting(self,ii,modify=False,Print=True)
            
            #replot with toggled index
        if event.key == '?':
            print(HELP)

        
#------------------------------click button events----------------------------#        
    def onclick(self, event):
        #double click to select an axes
        if event.dblclick:
            if event.inaxes in self.axesL[0]:
                self.page = 0
                for ii in range(len(self.axesL[self.page])):
                    if (self.axesL[self.page][ii]==event.inaxes): self.Lidx = ii
            elif event.inaxes in self.axesL[1]:
                self.page = 1
                for ii in range(len(self.axesL[self.page])):
                    if (self.axesL[self.page][ii]==event.inaxes): self.Lidx = ii
                        
            elif event.inaxes in self.axesR[0]:
                self.page = 0
                for ii in range(len(self.axesR[self.page])):
                    if (self.axesR[self.page][ii]==event.inaxes): self.Ridx = ii
            elif event.inaxes in self.axesR[1]:
                self.page = 1
                for ii in range(len(self.axesR[self.page])):
                    if (self.axesR[self.page][ii]==event.inaxes): self.Ridx = ii
            
                    
         
        if event.button in [1,3]:
        
            '''Left hand side is for fitting'''
            if event.inaxes is not None and ((event.inaxes in self.axesL[0]) or (event.inaxes in self.axesL[1])):
                self.Lidx = None
                if event.inaxes in self.axesL[0]: 
                    self.page = 0
                elif event.inaxes in self.axesL[1]:
                    self.page = 1
                    
                #Find the specific subplot
                for ii in range(len(self.axesL[self.page])):
                    if (self.axesL[self.page][ii]==event.inaxes): self.Lidx = ii
                        
                # if not first canvas, +6 is needed to appropraitely index the ions
                if self.Lidx is not None:
                    if self.page == 0: key_idx = self.Lidx
                    else: key_idx = self.Lidx+6
                        
                    vel = self.ions[self.keys[key_idx]]['vel']
                    name = self.ions[self.keys[key_idx]]['name']

                    if self.vclim is None:
                        self.vclim = [event.xdata]
                        self.axesL[self.page][self.Lidx].plot(event.xdata,event.ydata,'ro',ms=5)
                        if self.page == 0: self.fig.canvas.draw()
                        else: self.fig2.canvas.draw()
                        
                        
                    else:
                        vclim = np.sort(np.append(self.vclim,event.xdata))
                        #self.vclim_all.append(vclim)
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
                        
                        
                else:
                    print('Not in Axes: Reclick')
                    
                '''Right hand side for picking velocity limits for EW measurements'''
            if event.inaxes is not None and ((event.inaxes in self.axesR[0]) or (event.inaxes in self.axesR[1])):
                self.Ridx = None
                
                #identify page
                if event.inaxes in self.axesR[0]:
                    self.page = 0
                elif event.inaxes in self.axesR[1]:
                    self.page = 1

                # get subplot index
                for ii in range(len(self.axesR[self.page])):
                    if (self.axesR[self.page][ii]==event.inaxes): self.Ridx = ii
                        
                if self.Ridx is not None:
                    if self.page == 0: key_idx = self.Ridx
                    else: key_idx = self.Ridx+6
                        #if left click then define leftward vel limit
                    if event.button == 1:
                        self.EWlim[0] = event.xdata#used for plotting all with same range 'V' command
                        self.ions[self.keys[key_idx]]['EWlims'][0] = event.xdata
                        #plot selected limits
                        if self.page == 0: Plotting(self,self.Ridx,modify=False,Print=False)
                        elif self.page == 1: Plotting(self,self.Ridx,modify=False,Print=False)

                        #if right click define rightward vel limit
                    elif event.button == 3:
                        self.EWlim[1] = event.xdata
                        self.ions[self.keys[key_idx]]['EWlims'][1] = event.xdata
                        if self.page == 0: Plotting(self,self.Ridx,modify=False,Print=False)
                        elif self.page == 1: Plotting(self,self.Ridx,modify=False,Print=False)

                        
                        '''If V has been selected, the clicked limits will apply to All subplot frames
                            and the EW will be measured for all ions in page'''
                    if (self.event_button == 'V')&(None not in self.EWlim):
                        if self.page == 0: loop = self.nions
                        if self.page == 1: loop = self.nions2
                        for ii in range(loop):
                            self.ions[self.keys[ii]]['EWlims'][0] = self.EWlim[0]
                            self.ions[self.keys[ii]]['EWlims'][1] = self.EWlim[1]
                            
                            Plotting(self,ii,modify=False,Print=False)

                        self.event_button = None
                        self.EWlim = [None,None]
                            
                           

class Plotting:
    def __init__(self,parent,ii,modify=False,Print=False):
        
        # following if statement is to determine which canvas it is on for proper indexing (page1/page2)
        if parent.page == 0:
            key_idx = ii
            i = 0
        else:
            key_idx = ii+6
            i = 1
        
        #---------------Define variables for readability--------------#
        vel = parent.ions[parent.keys[key_idx]]['vel']
        wave = parent.ions[parent.keys[key_idx]]['wave']
        error = parent.ions[parent.keys[key_idx]]['error']
        flux = parent.ions[parent.keys[key_idx]]['flux']
        weight = parent.ions[parent.keys[key_idx]]['weight']
        name = parent.ions[parent.keys[key_idx]]['name']
        #L = parent.ions[parent.keys[key_idx]].L
        wc =np.array(parent.ions[parent.keys[key_idx]]['wc']) 
        cont = parent.ions[parent.keys[key_idx]]['cont']
        window_lim = parent.ions[parent.keys[key_idx]]['window_lim']
        order = parent.ions[parent.keys[key_idx]]['order']
        EWlims = parent.ions[parent.keys[key_idx]]['EWlims']
        lam_0=parent.ions[parent.keys[key_idx]]['lam_0']

        #--------------------------------------------------------------#
        if Print == False:
            if modify == True: #modify is whether or not the continuum is being adjusted
                #re-eval continuum
                parent.ions[parent.keys[key_idx]]['pco']=L.Legendre.fit(wave[wc],flux[wc],order,w=weight[wc])
                parent.ions[parent.keys[key_idx]]['cont'] = parent.ions[parent.keys[key_idx]]['pco'](wave)
                cont = parent.ions[parent.keys[key_idx]]['cont']
                gray_idx = np.where(np.diff(wc,prepend=np.nan))[0][1:]

            #clear axes to redraw modifications
            parent.axesL[i][ii].clear()
            parent.axesR[i][ii].clear()

            #replot left (flux; error)
            parent.axesL[i][ii].step(vel,flux,color='k',where='mid');parent.axesL[i][ii].step(vel,error,color='r',where='mid')
            parent.axesL[i][ii].step(vel,cont,color='b',where='mid')
            #replot right (flux; error)
            parent.axesR[i][ii].step(vel,flux/cont,color='k',where='mid');parent.axesR[i][ii].step(vel,error/cont,color='b',where='mid')
            parent.axesR[i][ii].axhline(y=1,xmin=window_lim[0],xmax=window_lim[1],ls='--',c='b')
            
            #plotting grayed spectra
            if modify == True:
                for zz in range(int(len(gray_idx)/2)):
                    parent.axesL[i][ii].step(vel[gray_idx[zz*2]:gray_idx[2*zz+1]],flux[gray_idx[2*zz]:gray_idx[2*zz+1]],vel[gray_idx[2*zz]:gray_idx[2*zz+1]],error[gray_idx[2*zz]:gray_idx[2*zz+1]],color='gray',alpha=1)


            #clear y ticks and label plots
            parent.axesL[i][ii].set_yticks([]); #parent.axesR[i][ii].set_yticks([])
            parent.axesL[i][ii].set_ylabel(name); parent.axesR[i][ii].set_ylabel(name)


            #set axes bounds
            parent.axesL[i][ii].set_xlim(window_lim); parent.axesR[i][ii].set_xlim(window_lim)
            parent.axesR[i][ii].set_ylim([0,2.2])

            #set x ticks only on bottom row
            if ii != len(parent.axesL[i]) - 1:
                parent.axesL[i][ii].set_xticks([]);parent.axesR[i][ii].set_xticks([])
            #plot column titles
            if ii == 0:
                parent.axesL[i][0].set_title('Continuum Fitter')
                parent.axesR[i][0].set_title('Normalized Spectra')

            #plot vel = 0 line
            parent.axesL[i][ii].axvline(0,ymin=0,ymax=parent.axesL[i][ii].get_ylim()[1],ls='--',c='k')
            parent.axesR[i][ii].axvline(0,ymin=0,ymax=parent.axesR[i][ii].get_ylim()[1],ls='--',c='k')

            #Plot EW velocity limit
            if EWlims[0] is not None: parent.axesR[i][ii].axvline(EWlims[0],ymin=0,ymax=2.5,ls='--',c='b')
            if EWlims[1] is not None: parent.axesR[i][ii].axvline(EWlims[1],ymin=0,ymax=2.5,ls='--',c='r')

            #### NOW PLOT INTERVENING LINES IF given
            if parent.intervening != False:
                outlist=grab_intervening_linelist(parent.intervening,np.double(parent.z),lam_0,wave)
                plot_intervening_lines( parent.axesR[i][ii],outlist,np.max(vel))


            #redraw MUST BE LAST ITEM IN LIST
            if parent.page == 0: parent.fig.canvas.draw()
            else: parent.fig2.canvas.draw()

        #plot EW measurements on frame
        if Print == True:
            if parent.ions[parent.keys[key_idx]]['EW_text'] is not None:
                parent.ions[parent.keys[key_idx]]['EW_text'].remove()

            plotText(parent,parent.ions[parent.keys[key_idx]])
            text = parent.ions[parent.keys[key_idx]]['text']
            if parent.page == 0:
                xlim = parent.ions[parent.keys[key_idx]]['EWlims'][0]
                EWtoggle = parent.axesR[i][ii].text(-850,1.8,text)
                parent.ions[parent.keys[key_idx]]['EW_text'] = EWtoggle
                #parent.axesR[i][ii].text(xlim,ylim,'EW:'+str(np.round(parent.ions[parent.keys[key_idx]].EW*1000,2)),color='blue')
                parent.fig.canvas.draw()
            else:
                xlim = parent.ions[parent.keys[key_idx]]['EWlims'][0]
                EWtoggle = parent.axesR[i][ii].text(-850,1.8,text)
                parent.ions[parent.keys[key_idx]]['EW_text'] = EWtoggle
                
                #parent.axesR[i][ii].text(xlim,ylim,'EW:'+str(np.round(parent.ions[parent.keys[key_idx]].EW*1000,2)),color='blue')
                parent.fig2.canvas.draw()
                
        
                
            
                    
class EW:
    
    def __init__(self,parent,page,ii,lims):

        #determine which page is being accessed
        if page == 0: key_idx = ii
        else: key_idx = ii+6
            
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
        if page == 0: Plotting(parent,ii,modify=False,Print=True)
        else: Plotting(parent,ii,modify=False,Print=True)
            

        

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
        
    
class Transitions:
    def __init__(self,Abs,intervening=False):
        app = QtWidgets.QApplication(sys.argv)
        main = mainWindow(Abs,intervening=intervening)
        main.resize(1400,900)

        main.show()
        sys.exit(app.exec_())
