import matplotlib
matplotlib.use('Qt5Agg')
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from astropy.convolution import convolve, Box1DKernel
import pdb
import sys
import os


class Spec_Inspect(object):

    def __init__(self,two_d_spec):
        self.two_d_spec=two_d_spec
        self.active_two_d_spec=two_d_spec
        fig, ax = plt.subplots(2, 1, sharex=True)
        self.fig=fig
        self.ax=ax
        self.pix=np.array([1.])
        #Sum in dispersion direction to do initial selection
        temp=np.sum(two_d_spec,axis=1)
        temp_cumsum=np.cumsum(temp)/np.sum(temp)
        xlist=np.arange(0,len(temp_cumsum),1)
        self.active_1d_spec=extract_1_d(self.active_two_d_spec)

        self.extration_y=[int(np.interp(0.05, temp_cumsum,xlist)), int(np.interp(0.95, temp_cumsum,xlist))]
        self.temp_extraction_y=[]
        print(self.extration_y)
        self.active_1d_spec=extract_1_d(self.active_two_d_spec[self.extration_y[0]:self.extration_y[1],:])
        self.master_plotter()

        

        # Connect the different functions to the different events
        self.fig.canvas.mpl_connect('key_press_event', self.ontype)
        #plt.gcf().canvas.mpl_connect('button_press_event',self.onclick)
        #plt.gcf().canvas.mpl_connect('pick_event',self.onpick)
        plt.show() # show the window


    def master_plotter(self,one_d_only=False):

        self.ax[1].cla()
        
        if one_d_only==False:
            self.ax[0].cla()
            im = self.ax[0].imshow(self.two_d_spec,origin = 'lower', vmin = -10, vmax = 65)
            self.fig.colorbar(im, ax=self.ax[0], label='Interactive colorbar',location='top')

        sp=self.ax[1].plot(self.active_1d_spec)
        xlim=self.ax[0].get_xlim()
        self.extr_y_min=self.ax[0].hlines(self.extration_y[0],xlim[0],xlim[1],colors='r', linestyles='dashed',label='ext_pt_min')
        self.extr_y_max=self.ax[0].hlines(self.extration_y[1],xlim[0],xlim[1],colors='r', linestyles='dashed',label='ext_pt_min')

        self.ax[0].set_aspect('auto')




    def ontype(self,event):    

        if event.key=='c':
            #Figure out the min max of extraction box
            self.ax[0].plot(event.xdata,event.ydata,'r+')
            plt.draw()

            self.temp_extraction_y=np.append(self.temp_extraction_y,event.ydata)



            print(event.xdata,event.ydata,event.key,event.x,event.y)
            print(len(self.temp_extraction_y))

            if len(self.temp_extraction_y)==2:
                #First remove previous extraction window lines HOW?

                ext_min_y=int(np.round(min(self.temp_extraction_y)))
                ext_max_y=int(np.round(max(self.temp_extraction_y)))
                xlim=self.ax[0].get_xlim()
                self.ax[0].hlines(ext_min_y,xlim[0],xlim[1],colors='r', linestyles='dashed',label='ext_pt_min')
                self.ax[0].hlines(ext_max_y,xlim[0],xlim[1],colors='r', linestyles='dashed',label='ext_pt_max')
                #print(self.extration_y,event.x,event.y)
                self.active_two_d_spec=self.active_two_d_spec[ext_min_y:ext_max_y,:]
                self.active_1d_spec=extract_1_d(self.active_two_d_spec)
                self.master_plotter(one_d_only=True)
    
                self.temp_extraction_y=[]
            plt.draw()

        #Reset Everything
        elif event.key=='r':
            self.active_two_d_spec=self.two_d_spec
            self.active_1d_spec=extract_1_d(self.active_two_d_spec)
            self.master_plotter()
            plt.draw()

        # Set top y max
        elif event.key=='t':
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()
            self.ax[1].set_ylim([ylim[0],event.ydata])
            self.ax[1].set_xlim(xlim)
            plt.draw()
        # Set top y min
        elif event.key=='b':
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()
            self.ax[1].set_ylim([event.ydata,ylim[1]])
            self.ax[1].set_xlim(xlim)
            plt.draw()


        # Smooth spectrum
        elif event.key=='S':
            self.pix[0] += 2
            Filter_size=np.int(self.pix[0]) 
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()
            self.active_1d_spec =convolve(extract_1_d(self.active_two_d_spec), Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            self.master_plotter(one_d_only=True)
            self.ax[1].set_ylim(ylim)
            self.ax[1].set_xlim(xlim)

            plt.draw()
        #Unsmooth Spectrum
        elif event.key=='U':
            self.pix[0] -= 2
            if self.pix[0] <= 0:
                self.pix[0]=1;
            Filter_size=np.int(self.pix[0]) 
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()

            self.active_1d_spec =convolve(extract_1_d(self.active_two_d_spec), Box1DKernel(Filter_size))#medfilt(flux,np.int(Filter_size))
            self.master_plotter(one_d_only=True)
            self.ax[1].set_ylim(ylim)
            self.ax[1].set_xlim(xlim)

            plt.draw()

            # Set X max
        elif event.key=='X':
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()
            self.ax[1].set_xlim([xlim[0],event.xdata])
            self.ax[1].set_ylim(ylim)
            plt.draw()
        # Set x min
        elif event.key=='x':
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()
            self.ax[1].set_xlim([event.xdata,xlim[1]])
            self.ax[1].set_ylim(ylim)
            plt.draw()

        # Set pan spectrum
        elif event.key==']':
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()
            delx=(xlim[1]-xlim[0])
            self.ax[1].set_xlim([xlim[1],xlim[1]+delx])
            self.ax[1].set_ylim(ylim)
            plt.draw()

        # Set pan spectrum
        elif event.key=='[':
            xlim=self.ax[1].get_xlim()
            ylim=self.ax[1].get_ylim()
            delx=(xlim[1]-xlim[0])
            self.ax[1].set_xlim([xlim[0]-delx,xlim[0]])
            self.ax[1].set_ylim(ylim)
            plt.draw()

    

def extract_1_d(input_2d_spec):
    return np.sum(input_2d_spec,axis=0)


  
