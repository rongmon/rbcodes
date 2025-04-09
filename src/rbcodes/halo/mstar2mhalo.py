import numpy as np
def mstar2mhalo(mstar,z):
    """
     Function to compute the halo mass and virial radius of a galaxy.
     We use the Moster 2010 relation to compute these values.
     Moster et al 2010. The Astrophysical Journal, 710:903?923, 2010 February 20
    
    
    Input:-     mstar:- stellar mass of the galaxy under consideration
                                                   (e.g.  10^(10.4) )
                   z :-  Redshift of the galaxy
    
    Returns:- 
        m200:-   halo mass according to the paper
       r200:-   halo virial radius according to NFW profile with c=.2
       Written :- R.B., 2018
    ------------------------------------------------------------------------
    """
    if z <= 0.1:
        m_M_0=0.02820;
        M1=10.**(11.884);
        beta=1.057;
        gamma=0.556;
    elif (z>0.1) & (z <= 0.5):
        m_M_0=0.0254;
        M1=10.**(11.95);
        beta=1.37;
        gamma=0.55;
    elif (z>0.5) & (z <= 0.7):
        m_M_0=0.0215;
        M1=10.**(11.93);
        beta=1.18;
        gamma=0.48;
    elif (z>0.7) & (z <= 0.9):
        m_M_0=0.0142;
        M1=10.**(11.98);
        beta=0.91;
        gamma=0.43;
    elif (z>0.9) & (z <= 1.1):
        m_M_0=0.0175;
        M1=10.**(12.05);
        beta=1.66;
        gamma=0.52;
    elif (z>1.1) & (z <= 1.5):
        m_M_0=0.0110;
        M1=10.**(12.15);
        beta=1.29;
        gamma=0.41;
    elif (z>1.5) & (z <= 1.8):
        m_M_0=0.0116;
        M1=10.**(12.28);
        beta=1.53;
        gamma=0.41;
    elif (z>1.8) & (z <= 2.5):
        m_M_0=0.0130;
        M1=10.**(12.22);
        beta=0.9;
        gamma=0.30;
    elif (z>2.5) & (z <= 3.5):
        m_M_0=0.0130;
        M1=10.**(12.21);
        beta=0.90;
        gamma=0.46;
    else:
        print('specify correct redshift range')
            

    M=10.**(np.linspace(10.,14.5,100)) #;   % Halo Mass

    mstarfit = M * 2.* m_M_0*( ((M/M1)**(-beta)) + (M/M1)**(gamma))**(-1.);
    
    index=np.argsort(mstarfit)
    
    M_sort=M[index]
    mstar_sort=mstarfit[index]
    
    m200=np.interp(mstar,mstar_sort,M_sort);
    
    
    
    # NFW virial radius 
    c=.2;    
    delta_vir=200.;
    omegam=0.238;
    rho0=1.4876862e+11;
    msun=1.98892e33;
    
    deltac = delta_vir / 3.*c**3./ (np.log10(1.+c) - c / (1.+c));
    rho_crit = rho0 * omegam * (1.+z)**3 ;
    r200 = ( m200/delta_vir/rho_crit/(4.*np.pi/3.) )** 0.33333 * 1000. 
    
        
    return (m200,r200)
    