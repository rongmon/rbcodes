import numpy as np
from rbcodes.igm import rb_setline as line       
from rbcodes.utils import rb_utility as rt
from rbcodes.gui.abstools import Absorber as A

import matplotlib as mpl
#mpl.use('TkAgg')
mpl.rcParams['lines.linewidth'] = .9
clr=rt.rb_set_color()
import matplotlib.pyplot as plt 
plt.style.use('dark_background')

import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, 
QLineEdit, QInputDialog)


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
    

class vStack(object):
    def __init__(self,wave,flux,error,line_flg,zabs=0,vlim=[-1000.,1000.]):  

        self.zabs=zabs
        self.vlim=vlim
        self.ions=prepare_absorber_object(zabs,wave,flux,error,line_flg='LLS')
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
        self.npages=int((self.nions/self.plotppage))
        
        
        fig=plt.figure(figsize=(12,8))
        self.fig=fig
        axes=list(range(self.plotppage))
        for i in range(self.plotppage):
            axes[i]=fig.add_subplot(self.nrow,self.ncol,i+1)
        self.axes=np.array(axes)
        self.vPlot()
        #fig.canvas.mpl_connect('button_press_event',self.onclick)
        fig.canvas.mpl_connect('key_press_event',self.onkb)
        #fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0)

        print('Press ? for help')
        plt.show()

    def onkb(self,event):
        # global axesL,axesR,spec,pfile
        if event.key=='Q':
            plt.close()
            sys.exit()
            # print 'Ending setv'
        #set up custom y-limit
        elif event.key=='Y':
            if event.inaxes in self.axes:
                i=np.where(event.inaxes==self.axes)[0][0]+self.plotppage*(self.page-1)
                Windowname='Manual y-Limits'
                instruction_text='Input range (e.g. 0.,2.)'
                temp=input_txt_dlg(Windowname,instruction_text)
                yrangetext=temp.filename

                yrange = yrangetext.split(',')
                yrange = np.array(yrange).astype('float32')
                self.vPlot(ploti=i,yrange=[yrange[0],yrange[1]])



        elif event.key=='>':
            self.page+=1
            if self.page>self.npages: 
                self.page=1
            # print self.page
            self.vPlot()
        elif event.key=='<':
            self.page-=1
            if self.page<1: 
                self.page=self.npages
                print(self.page)
            # print self.page
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
                #print(self.ions[self.keys[i]]['flag'])
                self.vPlot(ploti=i,comment=False)
        elif event.key=='S':
            self.Print_Selected_LineList()

        elif event.key=='?':
            print('''

            This GUI allows user to select if an absorption line is present or not
            Keyboard Commands
            Q: quit
            >: next page
            <: previous page
            S: Save the List of identified absorbers. [Still need to be implemented]
            w : Toggle between (1) Detection , (2) Non-Detection.
            Y : Set custom y-limit
            ''')

    def vPlot(self,ploti=None,comment=False,yrange=None):#spec,i=0):
        # global axesR
        # for line in self.transitions[self.page : self.page+self.plotppage]:
        if ploti is None:
            ploti=np.arange(self.plotppage*(self.page-1),min(self.plotppage*self.page,self.nions))
        else:
            ploti=[ploti]
            # ploti=range(self.plotppage)+self.plotppage*(self.page-1)
        # pdb.set_trace()
        for i in ploti:
            self.plotstuff(i,comment=comment,yrange=yrange)
        plt.draw()

        
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
            #print(flag)
            if flag==1:
                textcolor=clr['yellow']
            else:
                textcolor=clr['light_gray']
            ax.text(x=0.05, y=0.01, s=textout, fontsize=12, transform=ax.transAxes,color=textcolor)
        plt.draw()

    
    
    def plotText(self,flag=1):
        if flag==1:
            text='Detection'       
        elif flag==0:
            text='Non-Detection'      
        elif flag==2:
            text ='Blended-detection'
        return text
       
    # Executes on button press in Print_list.
    def Print_Selected_LineList(self):

        pfile='LineList_Identified.log'
        print('-----------------------------')

        print('                             ')
        print('                             ')

        #uname=input('Saving to '+pfile+'. Return to confirm, or enter a new name or c to cancel: ' )
        
        #if uname=='c':
        #    return
        #elif uname!='':
        #    pfile=uname
        Windowname='Save LineList'
        instruction='Enter filename:'

        temp=input_txt_dlg(Windowname,instruction,default_text=pfile)
        pfile=temp.filename
        print('-----------------------------')

        print('Saving to: '+pfile)
        #import pandas as pd
        from astropy.io import ascii
        from astropy.table import Table


        # Create a pandas dataframe for objects marked as detections or blends
        #loop through ions and extract all the identified lines
        lst1 = []
        lst2=[]
        lst3=[]
        for i in (self.keys):
            #name = self.ions[self.keys[i]]['name']
            if self.ions[i]['flag']==1:
                #print(self.ions[i]['name']+' '+np.str('%.6f' % self.ions[i]['lam_0_z'])+' '+np.str(self.ions[i]['z']))
                zip1=self.ions[i]['name']
                zip2=np.str('%.6f' % self.ions[i]['lam_0_z'])#self.ions[i]['lam_0_z']#
                zip3=self.ions[i]['z']#np.str(self.ions[i]['z'])
                lst1.append(zip1)
                lst2.append(zip2)
                lst3.append(zip3)


        df = Table([lst1, lst2, lst3], names=('Name', 'wave_obs','zabs'))
        ascii.write(df, pfile, delimiter=' ',overwrite=True)
        

class popup_windows(QWidget):
    def __init__(self,Windowname,instruction,default_text='test'):
        super().__init__()
        filename=self.initUI(Windowname,instruction,default_text=default_text)

        self.filename=filename

    def initUI(self,Windowname,instruction,default_text='test'):

        text, ok = QInputDialog.getText(self, Windowname, instruction)
        #self.tableline.setText("Spectrum_Analysis_z_"+str(parentvals.z)+"_Measurement_Table.dat")

        if ok:
            self.filename=text
        return text



#Initial inputs and callable class to run proram        
class input_txt_dlg:
    def __init__(self,Windowname,instruction,default_text='test'):
        app = QApplication(sys.argv)
        main = popup_windows(Windowname,instruction,default_text=default_text)
        self.filename=main.filename
        #main.show()
        #sys.exit(app.exec_())

        
