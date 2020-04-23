import os
import numpy as np
import scipy.integrate
from scipy.interpolate import interp1d
from scipy.special import jv
from collections import OrderedDict

def ReadParameterFile(df):
    """
    Return cosmological and simulation parameters
    """
    #
    # Read each line of the parameter file
    #
    fr = open(df, "r")
    lines = fr.readlines()
    fr.close()
    #
    # Parse for relevent parameters
    #
    cparams = { }
    sparams = { }
    for i in range(len(lines)):
        cline = str.split(lines[i])
        if "#" not in lines[i] and len(cline) > 1:
            if cline[0] == "OMEGA_CDM":
                cparams["omega_c"] = float(cline[-1])
            
            if cline[0] == "DEUT":
                cparams["omega_b"] = float(cline[-1])
            
            if cline[0] == "OMEGA_NU":
                cparams["omega_nu"] = float(cline[-1])
            
            if cline[0] == "HUBBLE":
                cparams["h"] = float(cline[-1])
            
            if cline[0] == "SS8":
                cparams["sigma8"] = float(cline[-1])
           
            if cline[0] == "CDM_MASS":
                cparams["mass"] = float(cline[-1])
            
            if cline[0] == "SCATTER":
                cparams["scatter"] = float(cline[-1])

            if cline[0] == "NS":
                cparams["n_s"] = float(cline[-1])
            
            if cline[0] == "Z_FIN":
                sparams["z_final"] = float(cline[-1])
            
            if cline[0] == "TRANS":
                cparams["transfer"] = cline[-1]
            
            if cline[0] == "Z_IN":
                cparams["z_initial"] = float(cline[-1])
                sparams["z_initial"] = float(cline[-1])
            
            if cline[0] == "NSTEPS":
                sparams["n_steps"] = int(cline[-1])
            
            if cline[0] == "ALPHA":
                sparams["alpha"] = float(cline[-1])
            
            if cline[0] == "RL":
                sparams["box_size"] = float(cline[-1])
            
            if cline[0] == "NP_1":
                sparams["np1"] = int(cline[-1])
            
            if cline[0] == "NP_2":
                sparams["np2"] = int(cline[-1])
            
            if cline[0] == "SPEC_1":
                sparams["spec_1"] = int(cline[-1])
            
            if cline[0] == "SPEC_2":
                sparams["spec_2"] = int(cline[-1])
            
            if cline[0] == "GLASS_START_1":
                sparams["glass_1"] = cline[-1]
            
            if cline[0] == "GLASS_START_2":
                sparams["glass_2"] = cline[-1]
            
            if cline[0] == "NG":
                sparams["ngrid"] = int(cline[-1])
            
            if cline[0] == "FULL_ALIVE_DUMP":
                sparams["snapshots"] = ParseSnapshots(cline)
            
            if cline[0] == "T_CMB":
                cparams["t_cmb"] = float(cline[-1])
            
            if cline[0] == "N_EFF_MASSLESS":
                cparams["n_eff_massless"] = float(cline[-1])
            
            if cline[0] == "N_EFF_MASSIVE":
                cparams["n_eff_massive"] = float(cline[-1])
    
            if cline[0] == "Z_PR":
                td = []
                for z in cline:
                    if z=="Z_PR":
                        continue
                    td.append(float(z))
                sparams["z_print"] = td

            if cline[0] == "FREQUENCY":
                sparams["frequency"] = float(cline[-1])
            
            if cline[0] == "GAUSSIAN_WIDTH":
                sparams["width"] = float(cline[-1])
            
            if cline[0] == "INITIAL_SHIFT":
                sparams["shift"] = float(cline[-1])
            
            if cline[0] == "DURATION":
                sparams["duration"] = float(cline[-1])

    #
    # Some adjustments
    #
    # Change omega_b from omega_b*h^2 to omega_b
    cparams["omega_b"] = cparams["omega_b"] / cparams["h"]**2
    # Combined baryon and cdm
    cparams["omega_cb"] = cparams["omega_c"] + cparams["omega_b"]
    # Total omega_m
    cparams["omega_m"] = cparams["omega_c"] + cparams["omega_b"] + cparams["omega_nu"]
    # Neutrino and other contributions 
    cparams["omega_n"] = 0.
    cparams["omega_k"] = 0.
    cparams["omega_r"] = 2.471e-5/cparams["h"]**2 * (cparams["t_cmb"]/2.725)**4
    cparams["f_nu_massless"] = (7./8.) * (4./11.)**(4./3.) * cparams["n_eff_massless"]
    cparams["omega_n_massless"] = cparams["f_nu_massless"] * cparams["omega_r"]
    cparams["f_nu_massive"]     = (7./8.) * (4./11.)**(4./3.) * cparams["n_eff_massive"]
    if cparams["n_eff_massive"] > 0.:
        cparams["m_nu"] = 93.14*cparams["h"]**2 * cparams["omega_n"] / cparams["n_eff_massive"]
    else:
        cparams["m_nu"] = 0.
    cparams["omega_l"] = 1. - (cparams["omega_cb"] + cparams["omega_r"] + cparams["omega_n_massless"] + cparams["omega_n"] + cparams["omega_k"])
    
    return sparams, cparams

def CreateRedshiftDict(directory):
    """
    Parse the directory and create a dictionary of redshift:file
    """
    d={}
    f=[]
    files = os.listdir(directory)
    for name in files:
        if ".pk." in name:
            f.append(name)
        elif ".tf." in name:
            f.append(name)
    for name in f:
        zsplit = name.split("z")
        split = zsplit[0].split(".")
        if ".ini" in name:
            redshift=-1
        else:
            redshift=float(zsplit[-1])
        d[redshift]=name
    dd = OrderedDict(sorted(d.items(), key=lambda x: x[0], reverse=True))

    return dd

def CreateStepDict(directory):
    """
    Parse the directory and create a dictionary of redshift:file
    """
    d={}
    f=[]
    files = os.listdir(directory)
    for name in files:
        if ".1D." in name:
            f.append(name)
    for name in f:
        split = name.split(".")
        if ".ini" in name:
            step=0
        else:
            step = int(split[-2])+1
        d[step]=name
    dd = OrderedDict(sorted(d.items(), key=lambda x: x[0], reverse=False))

    return dd

def ConvertStepToRedshift(snaps, zi, zf, ns, zpr):
    """
    Converts simulation snapshot numbers to their corresponding redshift.
    """
    # Convert initial/final redshift to scale factor
    ai = 1./(1.+zi)
    af = 1./(1.+zf)

    apr = [1.0/(1.0 + z) for z in zpr]

    a = [ai]
    da_default = (af - ai)/ns
    da_previous = 0.0
    step = 0.0
    just_adjusted = False
    ind = 0
    da_previous = da_default
    for i in range(ns):
 
        step = da_default
        
        if(just_adjusted):
            step = 2.0*da_default - da_previous
            just_adjusted=False
        if (apr[ind] - a[-1]) < step:
            step = apr[ind]-a[-1]
        
        a.append(a[-1]+step)
        da_previous = step
        if(a[-1]==apr[ind]):
            ind+=1
            just_adjusted=True

    final = [1.0/ac - 1 for ac in a]

    zs = [final[shot] for shot in snaps]
    
    return zs

def log_interpolate(xx,yy,kind='linear'):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interp1d(logx, logy, kind=kind)
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

def GrowthFactorClassic(z, params):
    """
    Computes the cosmological growth factor at individual redshift z.
    """
    z_primoridal = 100000.
    
    x1 = 1./(1. + z_primoridal)
    x2 = 1./(1. + z)
    
    x  = np.array([x1, x2])
    y1 = np.array([x1, 0.])
    
    y2 = scipy.integrate.odeint(d_gf_cdmrnu, y1, x, args=(params,))[1]

    dplus = y2[0]
    ddot  = y2[1]
    
    x1 = 1./(1. + z_primoridal)
    x2 = 1.
    x  = np.array([x1, x2])
    y1 = np.array([x1, 0.])
    y2 = scipy.integrate.odeint(d_gf_cdmrnu, y1, x, args=(params,))[1]

    gf = dplus/y2[0]
    gd = ddot/y2[0]

    return gf

def GrowthFactorQuantum(z, params, k):
    """
    Computes the cosmological growth factor at individual redshift z.
    """
    z_primoridal = 100000.
    
    x1 = 1./(1. + z_primoridal)
    x2 = 1./(1. + z)
    
    x  = np.array([x1, x2])
    y1 = np.array([x1, 0.])
    
    y2 = scipy.integrate.odeint(delta_deriv, y1, x, args=(params,k))[1]

    dplus = y2[0]
    ddot  = y2[1]
    
    x1 = 1./(1. + z_primoridal)
    x2 = 1.
    x  = np.array([x1, x2])
    y1 = np.array([x1, 0.])
    y2 = scipy.integrate.odeint(delta_deriv, y1, x, args=(params,k))[1]

    gf = dplus/y2[0]
    gd = ddot/y2[0]

    return gf

#def Compute_Jeans(a, params):
#    NU=1.91715234e-26
#    alpha=NU/params["h"]
#    k1=(6*a*params["omega_m"])**(.25)*(params["mass"]/alpha)**(1/2)
#    return k1

def d_gf_cdmrnu(y, a, params):
    """
    Growth factor ODE to be integrated (see Adrian's code for more details).
    """

    dydx = 0.*y
    H = H_cdmrnu(params, a)
    dydx[0] = y[1]/a/H
    dydx[1] = -2.*y[1]/a + 1.5*params["omega_cb"]*y[0]/(H*a**4)
    
    return dydx

def delta_deriv(y,a, params, k):
    NU=1.91715234e-26
    alpha=NU*params["h"]
    kappa1=(k**4)*(1/6.0*params["omega_cb"])*(alpha/params["mass"])**2 #units of Mpc^4

    H = H_cdmrnu(params, a)
    dydx=0.*y

    dydx[0] = y[1]/a/H
    dydx[1] = -2.*y[1]/a + 1.5*params["omega_cb"]*y[0]/(H*a**4)*(1-kappa1/a)
    #dydx[1] = -(1/adot**2)*(add + 2.*adot**2/a)*y[1] +1.5*params["omega_cb"]*y[0]/(adot*adot*a**3)*(1-kappa1/a)

    return dydx

def ModeGrowthQ(k, cparams, zfin=0):
    z_primordial = 100000. #very very old
    
    a_primordial = 1./(1. + z_primordial)
    af = 1./(1.0 + zfin)
    
    sigma8=cparams["sigma8"]
    ns=cparams["n_s"]
    
    tm=KlypinHoltzmann(k, cparams)

    karray=np.logspace(-2, 2, 200)
    tmarray=KlypinHoltzmann(karray, cparams)

    s2=ComputeSigma2(karray, tmarray, ns)
    snorm=sigma8**2/s2


    #s2=ComputeSigma2(k, tm, ns)
    #snorm=sigma8**2/s2

    pk = snorm*k**ns*tm**2 #pk normalized to today
    pk*=a_primordial**2 #power spectrum scaled to primordial

    delta = np.sqrt(pk)#complex(np.random.normal(), np.random.normal());

    x=np.array([a_primordial, af])
    y1=np.array([delta, 0])

    del_k = scipy.integrate.solve_ivp(lambda a,y: delta_deriv(a,y,cparams,k), x, y1, method="BDF") 

    a=del_k.t
    delta=del_k.y[0]

    return a, delta

def ModeGrowthQIC(k, amp, cparams, zin, zfin=0):
    
    ain = 1.0/(1+zin)
    afin=1.0/(1+zfin)

    x=np.array([ain, afin])
    y1=np.array([amp, 0])

    del_k = scipy.integrate.solve_ivp(lambda a,y: delta_deriv(a,y,cparams,k), x, y1, method="BDF") 

    a=del_k.t
    delta=del_k.y[0]

    return a, delta

def ModeGrowthClassic(k, cparams, zfin=0):
    z_primordial = 100000. #very very old
    
    a_primordial = 1./(1. + z_primordial)
    af = 1./(1.0 + zfin)
    
    sigma8=cparams["sigma8"]
    ns=cparams["n_s"]
    
    tm=KlypinHoltzmann(k, cparams)

    karray=np.logspace(-2, 2, 200)
    tmarray=KlypinHoltzmann(karray, cparams)

    s2=ComputeSigma2(karray, tmarray, ns)
    snorm=sigma8**2/s2
    
    pk = snorm*k**ns*tm**2 #pk normalized to today
    pk*=a_primordial**2 #power spectrum scaled to primordial

    delta = np.sqrt(pk)#complex(np.random.normal(), np.random.normal());

    x=np.array([a_primordial, af])
    y1=np.array([delta, 0])

    del_k = scipy.integrate.solve_ivp(lambda a,y: d_gf_cdmrnu(y,a,cparams), x, y1, method="BDF") 
    
    a=del_k.t
    delta=del_k.y[0]

    return a,delta

def H_cdmrnu(params, a):
    """
    Hubble factor as a function of a.
    """

    om_r  = (1. + params["f_nu_massless"]) * params["omega_r"]/a**4
    om_cb = params["omega_cb"]/a**3
    om_nm = Omega_nu_massive_a(params, a)
    om_lm = params["omega_l"]
    om_k  = params["omega_k"]/a**2

    return np.sqrt(om_r + om_cb + om_nm + om_lm + om_k)

def addot(params, a):
    Hratio = H_cdmrnu(params,a)
    addot = (1./2)*(Hratio**2)*a*(3*(1./Hratio)**2 - 1.0)
    return addot

def Omega_nu_massive_a(params, a):
    """
    Chosen to correspond to either its corresponding radiation term or matter term -- whichever is larger.
    """

    mat = params["omega_n"]/a**3
    rad = params["f_nu_massive"]*params["omega_r"]/a**4

    return (mat>=rad)*mat + (rad>mat)*rad

def AxionPowerSpectrum(df, cparams, snorm):

    ns=cparams["n_s"]
    k, tm = np.genfromtxt(df, comments="#", dtype="float64", unpack=True, usecols=(0,6))

    pk=snorm*k**ns*tm**2

    return k,pk

def LinearTheoryPowerSpectrum(df, cparams, redshift):
    """
    Compute the total matter linear theory power spectrum
    """
    sigma8=cparams["sigma8"]
    ns=cparams["n_s"]
    #
    # Read total matter transfer function 
    #

    k, tm = np.genfromtxt(df, comments="#", dtype="float64", unpack=True, usecols=(0,6))

    if(cparams["transfer"]=="KH"):
        k=np.logspace(-2, 2, 200)
        tm=KlypinHoltzmann(k, cparams)

    #
    # Compute scale dependent (or not) growth function
    #

    #
    # Compute sigma(8)^2 with the current normalization and then renormalize
    # 

    s2 = ComputeSigma2(k, tm, ns)
    
    snorm = sigma8**2 / s2
    #
    # Linear theory power spectrum at z = 0 with the proper normalization 
    #
    print("snorm {}".format(snorm))
    pk = snorm * k**ns * tm**2

    #
    # Scale this by the linear growth factor
    #

    Dz = GrowthFactorClassic(redshift, cparams)
    pk *= Dz**2

    return k, pk

def ComputeSigma2(k, tk, ns):
    """
    Computes sigma(8)^2: sigma(8)^2 = 1/(2pi^2) int_0^infty k^2 P(k) W(k*8)^2 dk
    """

    assert k.shape == tk.shape

    # Compute the integrand
    y = (1./2./np.pi**2) * k**2 * k**ns * tk**2 * Wfilter(8*k)**2

    # Integrate
    s = np.trapz(y, k)

    return s

def ComputeSnorm(s8, ns, df):
    k, tm = np.genfromtxt(df, comments="#", dtype="float64", unpack=True, usecols=(0,6))

    s2 = ComputeSigma2(k, tm, ns)

    snorm = s8**2/s2

    return snorm

def Wfilter(x):
    """
    Fourier transform of a top-hat filter in real space.
    """

    w = (3./x**3) * (np.sin(x) - x*np.cos(x))

    return w

def KlypinHoltzmann(k, params):
    Omega_m=params["omega_m"]
    Omega_b=params["omega_b"]
    h=params["h"]
    akh1=((46.9*Omega_m*h*h)**0.670)*(1.0+(Omega_m*h*h*32.1)**(-0.532))
    akh2=((12.0*Omega_m*h*h)**0.424)*(1.0+(Omega_m*h*h*45.0)**(-0.582))
    alpha=((akh1)**(-1.0*Omega_b/Omega_m))*(akh2)**((-1.0*Omega_b/Omega_m)**3.0)
    kh_tmp = Omega_m*h*np.sqrt(alpha)*(1.0-Omega_b/Omega_m)**.6
    
    tt=(2.728/2.7)**2
    qkh=k*tt/kh_tmp

    tf=np.log(1.0+2.34*qkh)/(2.34*qkh)*(1.0+13.0*qkh+(10.5*qkh)**2.0 + (10.4*qkh)**3.0 + (6.51*qkh)**4.0)**(-0.25)

    return tf

def BesselGrowth(z, params, k):
    a = 1./(1. + z)
    NU=1.91715234e-26
    alpha=NU*params["h"]
    mass = params["mass"]
    
    arg = (1.0/np.sqrt(params["omega_m"]*a))*(alpha/mass)*k**2
    #arg = (alpha/mass)*k**2/np.sqrt(a)

    tmp =  a**(-0.25)*jv(-5.0/2, arg)

    arg = (1.0/np.sqrt(params["omega_m"]*1.0))*(alpha/mass)*k**2
    
    gf = tmp/(jv(-5.0/2.0, arg))
    
    return gf


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
