""" Read in atomic line information for a given or approximate rest frame  wavelength."""

from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
from astropy.io import ascii
from pkg_resources import resource_filename


'''
 Function to read in atomic line information for a given rest frame  wavelength.
                           Or 
 For the line matching the closest wavelength. 

 Input :
           lambda_rest :-  Rest Frame wavelength (in \AA) of the line to match
            method     :-    'closest' ->  If set will match the closest line.
                              'Exact'  ->  If set will match the exact wavelength.
 
 Output:    dic :- Dictionary with fval,lambda and species name.

 Example:   str=rb_setline(2796.3,'closest')


Written By: Rongmon Bordoloi                Jan 2018, Python 2.7
Edit:       Rongmon Bordoloi                            Sep 2018, Depreciated kwargs to be compatible with python 3
'''
def rb_setline(lambda_rest,method,linelist='atom'):
    
    #if kwargs.has_key('linelist'):
    #   linelist=kwargs['linelist']
    #else:
    #   linelist='LLS'
    
    line_str=read_line_list(linelist)
    wavelist=np.zeros((len(line_str),))
    name = np.empty(len(line_str), dtype='object')
    fval=np.zeros((len(line_str),))
    if linelist=='atom':
        gamma=np.zeros((len(line_str),))


    for i in range(0,len(wavelist)):
        wavelist[i]=np.double(line_str[i]['wrest'])
        fval[i]=np.float(line_str[i]['fval'])
        name[i]=np.str(line_str[i]['ion'])
        if linelist=='atom':
            gamma[i]=np.str(line_str[i]['gamma'])

    if method=='Exact':
        q= np.where( (np.abs(lambda_rest-wavelist) < 1e-3))
        if linelist=='atom':
            outstr={'wave':wavelist[q],'fval':fval[q],'name':name[q],'gamma':gamma[q]}
        else:
            outstr={'wave':wavelist[q],'fval':fval[q],'name':name[q]}


    elif method=='closest':
        idx=(np.abs(lambda_rest-wavelist)).argmin()
        if linelist=='atom':
            outstr={'wave':wavelist[idx],'fval':fval[idx],'name':name[idx],'gamma':gamma[idx]}  

        else:

            outstr={'wave':wavelist[idx],'fval':fval[idx],'name':name[idx]} 
    else:
        raise NameError('Specify the matching method, closest or Exact')

    return outstr



def read_line_list(label):
    

    if label=='atom':
        filename=resource_filename('IGM','lines/atom_full.dat')
    elif label == 'LLS':
        filename=resource_filename('IGM','lines/lls.lst')
    elif label == 'LLS Small':
        filename=resource_filename('IGM','lines/lls_sub.lst')
    elif label == 'DLA':
        filename=resource_filename('IGM','lines/dla.lst')
    else:
        print('Give Correct LineList')

    data = []

    if label=='atom':

        s=ascii.read(filename)

        for line in range(0,len(s['col1'])):
            source = {}
            source['wrest'] = float(s['col2'][line])
            source['ion'] = s['col1'][line]+' '+np.str(np.int(s['col2'][line]))
            source['fval']=float(s['col3'][line])
            source['gamma']=float(s['col4'][line])

            data.append(source)


    else:       
        f=open(filename,'r')
        header1 = f.readline()
        for line in f:
            line = line.strip()
            columns = line.split()
            source = {}
            source['wrest'] = float(columns[0])
            source['ion'] = columns[1]+' '+columns[2]
            source['fval']=float(columns[3])
            data.append(source)


    return data
