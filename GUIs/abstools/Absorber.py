


'''
Absorber 

Inputs:
flux; wave; error; linelist; redshift; bin

1st. Asborber will bin the flux,wave and error to clean the data
2nd. will pull the actual lamd_rest from the atom.dat file for all lines
-->this (2nd) will return a dictionary of lamd_rest,ion_name,fval,gamma

3rd. Initialized entries for the Vstack plotting
'''


from IGM import rb_setline as rb_setline
import numpy as np
import numpy.polynomial.legendre as L
c =  299792.458

class Absorber:
    
    #defining variables to be used in the transition plotting GUI
    def Transition(self,ion_dict,line_dat,wave,flux,error,z,mask,window_lim,nofrills=False):
            # Edit RB May21, 2020: added nofrills keyword to toggle between the continuum property. 
            #        Added to avoid issues of calling Absorber class very near to the edge of the detector.
            #        Does not apply when calling abstools.

            # VARIABLE INITIALIZATION
            ion_dict['f'] = line_dat['fval']
            ion_dict['lam_0'] = line_dat['wave']
            ion_dict['name'] = line_dat['name']
            ion_dict['gamma'] = line_dat['gamma']
            ion_dict['z'] = z
            ion_dict['window_lim'] = window_lim 

            '''Shifting to vel-space centered on lam_o'''
            ion_dict['lam_0_z'] = ion_dict['lam_0']*(1+z)
            ion_dict['vel'] = (wave-ion_dict['lam_0_z'])/ion_dict['lam_0_z']*c

            '''create window for flux,wave,error based on max and min velocity'''
            window = (ion_dict['vel']>window_lim[0]) & (ion_dict['vel']<window_lim[1])
            ion_dict['flux'] = flux[window]; ion_dict['wave']=wave[window]
            ion_dict['error'] = error[window]; ion_dict['vel'] = ion_dict['vel'][window]

            '''Initial Polyfit assuming a masked region of -200:200 and polyorder=4
            cont= continuum, pco= polynomial coefficients for cont fitting; weight= parameter to fit polys
            order = order of polynomial

            lets also give each ion object the Legendre function for ease of use during plotting'''

            if nofrills==False:
                ion_dict['wc'] = ((ion_dict['vel']<mask[0])|(ion_dict['vel']>mask[1]))
                ion_dict['weight'] = 1/(ion_dict['error']**2)
                ion_dict['order'] = 4 #order of poly fit
                ion_dict['pco']=L.Legendre.fit(ion_dict['wave'][ion_dict['wc']],ion_dict['flux'][ion_dict['wc']],ion_dict['order'],w=ion_dict['weight'][ion_dict['wc']])
                ion_dict['cont'] = ion_dict['pco'](ion_dict['wave'])
                #resett the wc so the plotter can mask as many absorbers as needed for fixing
                ion_dict['wc'] = [True]*len(ion_dict['wc'])


            '''Property initializations:'''
            ion_dict['N']=None; ion_dict['Nsig']=None; ion_dict['EW']=None; ion_dict['EWsig']=None
            ion_dict['med_vel'] = None; ion_dict['EWlims'] = [mask[0],mask[1]]; ion_dict['flag'] = 0
            #for text
            ion_dict['EW_text'] = None
    
    
    
    def __init__(self,z,wave,flux,error,lines=None,mask_init=[-200,200],window_lim=[-1000,1000],load_file = False,nofrills=False):
        mask = mask_init
        self.z =z
        self.ions = {}

        if lines:
            for line in lines:
                line_dat = rb_setline.rb_setline(line,method='closest')
                
                #if using Transition class uncomment below line. Also comment transition def, while uncommenting transition class, comment out lines 80-82
                #self.ions[line_dat['name']] =Transition(line_dat,wave,flux,error,self.z,mask_init,window_lim)
        
                self.ions[line_dat['name']] = {}
                ion_dict = self.ions[line_dat['name']]
                self.Transition(ion_dict,line_dat,wave,flux,error,z,mask,window_lim,nofrills=nofrills)
                                
            #last dictionary item for full spectra data                    
            self.ions['Target'] = {}
            self.ions['Target']['flux'] = flux
            self.ions['Target']['wave'] = wave
            self.ions['Target']['error'] = error
            self.ions['Target']['z'] = z
        else:
            print('Input Linelist and rerun')
            
            
        
        
            
# '''Initializing the properties of each specific transition'''       
# class Transition(object):
#     def __init__(self,line_dat,wave,flux,error,z,mask,window_lim):
#         #line_dat has properties lam_o, fval,name,gamma
#         self.gamma = line_dat['gamma']
#         self.f = line_dat['fval']
#         self.lam_0 = line_dat['wave']
#         self.name = line_dat['name']
#         self.z = z
#         self.window_lim = window_lim 
        
#         '''Shifting to vel-space centered on lam_o'''
#         self.lam_0_z = self.lam_0*(1+z)
#         self.vel = (wave-self.lam_0_z)/self.lam_0_z*c
        
#         '''create window for flux,wave,error based on max and min velocity'''
#         vmin,vmax = -1000,1000
#         window = (self.vel>vmin) & (self.vel<vmax)
#         self.flux = flux[window]; self.wave=wave[window]
#         self.error = error[window]; self.vel = self.vel[window]
        
#         '''Initial Polyfit assuming a masked region of -200:200 and polyorder=4
#         cont= continuum, pco= polynomial coefficients for cont fitting; weight= parameter to fit polys
#         order = order of polynomial
        
#         lets also give each ion object the Legendre function for ease of use during plotting'''
#         #self.L = L
#         #wc = ((self.vel<mask[0])|(self.vel>mask[1]))
#         self.wc = ((self.vel<mask[0])|(self.vel>mask[1]))
#         self.weight = 1/(self.error**2)
#         self.order = 4 #order of poly fit
#         self.pco=L.Legendre.fit(self.wave[self.wc],self.flux[self.wc],self.order,w=self.weight[self.wc])
#         self.cont = self.pco(self.wave)
#         #resett the wc so the plotter can mask as many absorbers as needed for fixing
#         self.wc = [True]*len(self.wc)
        

#         '''Property initializations:'''
#         self.N=None; self.Nsig=None; self.EW=None; self.EWsig=None
#         self.med_vel = None; self.EWlims = [None,None]; self.flag = 0
#         #for text
#         self.EW_text = None
