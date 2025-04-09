""" GUI for multiple line EW measurement."""

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from rbcodes.igm import rb_setline as rs 
import pdb


class vStack(object):
    def __init__(self,wave,flux,error,wrestlist,zabs=0,vlim=[-600.,600.]):  
        self.zabs=zabs
        self.vlim=vlim
        self.wave=wave/(1.+zabs)
        self.page=1
        self.plotppage=6
        self.wrestlist=wrestlist
        #pdb.set_trace()
        self.npages=int(np.ceil(len(self.wrestlist)/self.plotppage))
        fig=plt.figure(figsize=(12,8))
        self.fig=fig
        axesR=list(range(self.plotppage))
        axesL=list(range(self.plotppage))
        for i in range(self.plotppage):
            axesL[i]=fig.add_subplot(6,2,2*i+1)
            axesR[i]=fig.add_subplot(6,2,2*(i+1))
        self.axesL=np.array(axesL)
        self.axesR=np.array(axesR)
        axesL=np.array(axesL)
        axesR=np.array(axesR)
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
            self.page=(self.page)%(self.npages)+1
            # print self.page
            self.vPlot()
        elif event.key=='<':
            self.page-=1
            if self.page==0: self.page=self.npages
            # print self.page
            self.vPlot()
        elif event.key=='s':
            pass

    def vPlot(self,ploti=None,left=True,right=True):
        if ploti is None:
            if self.page==self.npages: #turn extra axes off:
                for i in np.arange(len(self.wrestlist),self.plotppage*self.page)-(self.npages-1)*self.plotppage:
                # for i in np.arange(self.plotppage*self.page-self.nlines,self.plotppage): #+((self.plotppage*self.page) % self.nlines):
                    self.axesL[i].set_visible(False)
                    self.axesR[i].set_visible(False)
                    ax.text(1,1,'test:'+ np.str(12))
            ploti=np.arange(self.plotppage*(self.page-1),min(self.plotppage*self.page,len(self.wrestlist)))
        else:
            ploti=[ploti]
            # ploti=range(self.plotppage)+self.plotppage*(self.page-1)
        # pdb.set_trace()
        #for i in ploti:
            #if right: self.plotR(i)
            #if left: self.plotL(i)
        plt.draw()


