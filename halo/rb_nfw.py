import numpy as np
import scipy.integrate as integrate
def rb_nfw(m200,c,z):
    """
    Function to compute a NFW profile.
    Velocity Dispersion equation taken from Hoeft M.; Mucket J. P. & Gottlober, S, 2004, ApJ 602,1
    http://adsabs.harvard.edu/cgi-bin/bib_query?2004ApJ...602..162H


    Input :-  
       m200 :-   Halo mass
       c    :-   NFW concentration paramter
       z    :- redshift


    Returns :-

       A bunch of stuff
       





    
    """

    #Setting up cosmology
    rho0=1.4876862e+11;
    omegam=0.238000;
    msun=1.98892e+33;
    delta_vir=200.;
    G=6.6730003e-08;
    kmpsToCmps = 1.0*10.**(5.);
    Rvir=200.;
    kpc2cm=3.086*10.**(21);
    
    deltac = (delta_vir/3.)*( (c**3.)/( np.log(1.+c) - (c / (1.+c))));
    rho_crit =rho0*omegam*(1.+z)**3.;
    r200 =(m200/delta_vir / rho_crit / (4.*np.pi/3.) )**0.33333 * 1000. ;
    v200 = ((6.67e-8 * m200 * msun / (r200* 3.086*10.**(21.)) )**0.5)/1e5 ;
    
    r =np.linspace(1.,3.*r200,500);  # kpc
    rs = r200 / c; 
    ss=(((r/rs)*(1.+(r/rs))**2.)*1000.**3);
    rho = (rho_crit * deltac)/(ss); 
    M_r = 4.*np.pi* integrate.cumtrapz((r**2)*rho, r,initial=0.)
    
    x = r/r200 ;
    tab=1./x*(np.log(1.+c*x)-c*x/(1.+c*x))/(np.log(1.+c)-c/(1.+c));
    vcirc = v200*(tab)**0.5 ;
    maxvcirc =  np.max(vcirc) ;
    q=np.where((vcirc == np.max(vcirc)));
    maxvcircr = r[q];
    
    
    # Now compute V_Esc as per nfw.pro Binney & Tremaine  equation 2.31
    Phi_new = r * 0.0;
    vesc = r * 0.0 ;
    for ir in range(2,len(r)-4):
        term1 = (np.trapz(rho[0:ir]*(r[0:ir]**2.),x=r[0:ir])/(r[ir]))* msun; 
        term2 = np.trapz(rho[ir:len(r)]*r[ir:len(r)],x=r[ir:len(r)])*msun;        
        Phi_new[ir] = -4. *np.pi*6.67e-8*(term1 + term2)/3.086e21 ;
        vesc[ir] = ((2. * np.abs(Phi_new[ir]))**0.5) / 1e5 ; # See Binney & Tremaine (2-22) 
    

    # Chage Units to do velocity dispersion calculations
    rcm=r*kpc2cm;

    #M_r in gram
    M_r_gram=M_r*msun;

    Phi=G*integrate.cumtrapz((M_r_gram/rcm**(2)),rcm,initial=0);
    
    Phi=Phi*(1./((1e5)**2.));#%km^2/s^2
    Phi_out=np.max(Phi);

    k=0.41;
    a=0.29;

    sig = np.sqrt(a *(( Phi/Phi_out)**(k))*(Phi_out -Phi));
    
    nfw={}
    qqqt=np.where((vesc==0.))
    vesc[qqqt]=1e-99

    nfw["m200"]=m200;
    nfw["c"]=c;
    nfw["r200"]=r200;
    nfw["v200"]=v200;
    nfw["maxvcirc"]=maxvcirc;
    nfw["maxvcircr"]=maxvcircr;
    nfw["r"]=r;
    nfw["rho"]=rho;
    nfw["vcirc"]=vcirc;
    nfw["M_r"]=M_r;
    nfw["sig_v"]=sig;
    nfw["vesc"]=vesc;
    
    return nfw