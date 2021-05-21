import numpy as np
from IGM import rb_setline as line       
from utils import rb_utility as rt
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = .9
clr=rt.rb_set_color()
from GUIs.abstools import Absorber as A

import matplotlib.pyplot as plt 
plt.style.use('dark_background')
import sys

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
        	self.ions[i]['flag']=2
        
        
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
            # print 'Ending setv'
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
        elif event.key in ['1','2','3']: #Detected,non-detected, blended-detection            
            if event.inaxes in self.axes:
                i=np.where(event.inaxes==self.axes)[0][0]+self.plotppage*(self.page-1)
                self.ions[self.keys[i]]['flag']=int(event.key)
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
            1/2/3 : Flag absorber as (1) Detection , (2) Non-Detection, or (3) Blended Detectin.
            ''')

    def vPlot(self,ploti=None,comment=False):#spec,i=0):
        # global axesR
        # for line in self.transitions[self.page : self.page+self.plotppage]:
        if ploti is None:
            ploti=np.arange(self.plotppage*(self.page-1),min(self.plotppage*self.page,self.nions))
        else:
            ploti=[ploti]
            # ploti=range(self.plotppage)+self.plotppage*(self.page-1)
        # pdb.set_trace()
        for i in ploti:
            self.plotstuff(i,comment=comment)
        plt.draw()

        
    def plotstuff(self,i,comment=False):
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
        
        ax.axhline(1)
        ax.axvline(0,color='k',linestyle='dotted')
        ax.text(x=0.05, y=0.815, s=name, fontsize=10, transform=ax.transAxes,color=clr['red'])
        ax.text(x=0.75, y=0.815, s='f0: '+str(f0), fontsize=10, transform=ax.transAxes,color=clr['red'])
        
        if comment != False:
            ax.text(x=0.85, y=0.815, s=comment, fontsize=12, transform=ax.transAxes,color=clr['teal'])
            

        
        if flag is not None: #Display some measurement
            text=self.plotText(flag=flag)
            if flag==1:
            	textcolor=clr['yellow']
            else:
            	textcolor=clr['light_gray']
            ax.text(x=0.05, y=0.01, s=text, fontsize=12, transform=ax.transAxes,color=textcolor)
        plt.draw()

    
    
    def plotText(self,flag=1):
        if flag==1:
            text='Deteaction'       # Upper Limit
        elif flag==2:
            text='Non-Detection'       # Upper Limit
        # Lower Limit [only for logN]
        elif flag==3:
            text ='Blended-detection'
        return text
       
    # Executes on button press in Print_list.
    def Print_Selected_LineList(self):

        pfile='LineList_Identified.log'
        #uname=input('Saving to '+pfile+'. Return to confirm, or enter a new name or c to cancel: ' )
        
        #if uname=='c':
        #    return
        #elif uname!='':
        #    pfile=uname
        #import pandas as pd

        # Create a pandas dataframe for objects marked as detections or blends
        #loop through ions and extract all the identified lines

        for i in (self.keys):
        	#name = self.ions[self.keys[i]]['name']
        	if self.ions[i]['flag']==1:
        		print(self.ions[i]['name']+' '+np.str('%.6f' % self.ions[i]['lam_0_z'])+' '+np.str(self.ions[i]['z']))
        
    #logfile=get(handles.logfile,'String');
    #logname=[handles.PathName  logfile];
    #fileID = fopen(logname,'a+');
    #%fprintf(fileID,'%6s %12s\n','x','exp(x)');
    #%fprintf(fileID,'%6.2f %12.8f\n',A);
    
    
    #display('Start Printing Line Lists')
   #z_abs=str2double(get(handles.z_guess,'String'));
   #     for ind=1:length(handles.l_str.species)
   #         ss=get(handles.hcb(ind));
   #         if ss.Value==1
   #             tab=[ss.String  '  '  num2str(handles.l_str.lambda(ind)*(1.+z_abs))  ' '  num2str(z_abs)];
   #             formatspec=['%'  num2str(length(tab))  's\n'];
   #            fprintf(formatspec,tab);
   #             fprintf(fileID,formatspec,tab);
   #             
   #         end
   #     end
   #fclose(fileID);  

