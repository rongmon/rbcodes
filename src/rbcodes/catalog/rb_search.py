import numpy as np
def cone_search(ra_center,dec_center,ra_list,dec_list,angular_scale):
    """
    Function to do a cone search around any (ra,dec) pointing with respect
    to a list of ra,dec entries.

    Parameters
    ----------
    
        ra_center     = RA of center (Degrees)
        dec_center    = DEC of center (Degrees) 
        ra_list       = RA list of input catalogue (Degrees)
        dec_list      = DEC list of input catalogue (Degrees) 
        angular_scale = angular search radius in arcsec
               
    
    Returns
    -------
        out:- struture containing RA,DEC and logical operator identifying the objects
    
    Written by R.B.  Mar 11 2013
    -------------------------------------------------------
    """
    deg2rad=np.pi/180
    rad2deg=180./np.pi

    ra_center=deg2rad*np.array(ra_center)
    dec_center=deg2rad*np.array(dec_center)
    ra_list=deg2rad*np.array(ra_list)
    dec_list=deg2rad*np.array(dec_list)

    theta =np.arccos(np.sin(dec_center)*np.sin(dec_list) + np.cos(dec_center)*np.cos(dec_list)*np.cos(ra_center-ra_list))
 
    sq= np.where(theta*rad2deg*3600<=angular_scale)
    
    if sq[0].size>0:
        outstr={'ra':ra_list[sq],'dec':dec_list[sq],'index':sq,'angle':theta[sq]*rad2deg*3600}
    else:
        outstr={'ra':-99,'dec':-99,'index':sq,'angle':-99}
    return outstr
