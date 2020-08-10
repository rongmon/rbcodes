"""  Several utility functions. (1) show for loop progress (2) load color list"""
#!/usr/bin/python
from numpy import double
import sys
from time import time


def rb_perccount(jj,maxjj):
    """
    -------------------------------------------------------------------
    This code reports the percentage of a for loop done to the screen.
    **** NOTE *** This code currently only works on command line.
	Need to figure out how to make it work in notebook environment.

	Calling Sequence: rb_perccount(I,Imax)
    Input:
            I = Current iteration between 0 and Imax
	     Imax = Maximum number of iterations
	
    Output: 
            Percentage of the for loop done.
            Shows the amount of time elapsed during the job.
	
    WARNING: Do not print anything to the screen between calls of this function!
	
       Written by Rongmon Bordoloi Nov 15 2016
	--------------------------------------------------------------------
    """
    global _start_time
    jj=double(jj)
    maxjj=double(maxjj)
    CURSOR_UP_ONE = '\x1b[1A'
    ERASE_LINE = '\x1b[2K'
    if (jj ==0.0):
        pc_done=0.0
        _start_time=time()
        elapsed=double(_start_time-(_start_time))
    else:
        end_time=time()
        elapsed=double(end_time-(_start_time))
        pc_done  =  ((jj+1.)/(maxjj+1.)) * 100.

    elapsed_str = format_interval(elapsed)
    remainTime = format_interval(elapsed/(jj+1.)*(maxjj-jj))
    print ('{:s}'.format('Percentage Done: ')+'{:5.2f}'.format(pc_done)+'{:s}'.format('%    Time Elapsed: ') +elapsed_str+'{:s}'.format('   Time Remaining: ') +remainTime)
    sys.stdout.write(CURSOR_UP_ONE)
    #if (jj != maxjj):
    #    sys.stdout.write(ERASE_LINE)

def format_interval(t):
    """
        Formats a number of seconds as a clock time, [H:]MM:SS
        
        Parameters
        ----------
        t  : int
        Number of seconds.
        Returns
        -------
        out  : str
        [H:]MM:SS
        """
    mins, s = divmod(int(t), 60)
    h, m = divmod(mins, 60)
    if h:
        return '{0:d}:{1:02d}:{2:02d}'.format(h, m, s)
    else:
        return '{0:02d}:{1:02d}'.format(m, s)

   
def rb_set_color():
   clr={}
   clr['black']=[  0.,   0.,   0.]
   clr['white'] =[ 255./ 255., 255./ 255., 255./ 255.]
   clr['red'] =[ 255./ 255.,   0./ 255.,   0./ 255.]
   clr['green'] =[   0., 255./ 255.,   0.]
   clr['blue'] =[   0.,   0., 255./ 255.]
   clr['dark_orange']=[1., 0.5, 0.2]
   clr['orange'] =[ 230./ 255., 159./ 255.,   0./ 255.]
   clr['sky_blue'] =[  86./ 255., 180./ 255., 233./ 255.]
   clr['bluish_green'] =[   0./ 255., 158./ 255., 115./ 255.]
   clr['yellow'] =[ 240./ 255., 228./ 255.,  66./ 255.]
   clr['blue2'] =[   0./ 255., 114./ 255., 178./ 255.]
   clr['vermillion'] =[ 213./ 255.,  94./ 255.,   0.]
   clr['reddish_purple'] =[ 204./ 255., 121./ 255., 167./ 255.] 
   clr['cream'] =[ 248./ 255., 248./ 255., 248./ 255.]
   clr['cyan'] =[   0./ 255., 255./ 255., 255./ 255.]
   clr['light_lime_green'] =[ 153./ 255., 255./ 255., 153./ 255.]
   clr['pale_lime_green'] =[ 200./ 255., 255./ 255., 200./ 255.]
   clr['purple_wordle'] =[ 204./ 255.,   0., 255./ 255.]
   clr['light_purple'] =[ 231./ 255., 113./ 255., 255./ 255.]
   clr['orange2'] =[ 217./ 255., 110./ 255.,   0./ 255.]
   clr['light_orange'] =[ 233./ 255., 172./ 255.,  55./ 255.]
   clr['teal'] =[  62./ 255., 133./ 255., 181./ 255.]
   clr['pale_red'] =[ 255./ 255., 118./ 255., 110./ 255.]
   clr['pale_cyan'] =[ 133./ 255., 249./ 255., 255./ 255.]
   clr['dark_red'] =[ 100./ 255.,   0.,   0.]
   clr['dark_green'] =[   0., 100./ 255.,   0.] 
   clr['dark_blue'] =[   0.,   0., 100./ 255.]
   clr['gray']=[.5,.5,.5]
   clr['light_gray']=[.8,.8,.8]

   return clr
