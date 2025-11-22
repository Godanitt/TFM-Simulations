

import numpy as np

""" 
Modelo de pablo para el ajuste del cociente
"""


# ------------------------------------------------------------------------------
# Datos del paper: 
# Observation of strong wavelength-shifting in the argon-tetraﬂuoromethane system 
# DOI 10.3389/fdest.2023.1282854

K_Ar3rd_to_CF4_plus_star    = 49    # ns-1
K_Ar3rd_to_Ar               = 4.1   # ns-1
n                           = 1     # ns # (??)

# Santorilli and Diego et al doi 10.1140/epjc/s10052-021-09375-3 (en principio igual que el de Pablo Amoedo)
tau_3rd                      = 5.02               

# Realmente  es el cociente K_Ar_dblstar_to_CF3_star/K_Ar_dblstar_to_Ar_star = 36.5, pero a efectos prácticos es lo mismo asumir esto
K_Ar_dblstar_to_CF3_star    = 36.5  # ns-1
K_Ar_dblstar_to_Ar_star     = 1     # ns-1

# ------------------------------------------------------------------------------

def Pgamma_CF3_refined(f_cf4, Pgamma_CF3_star_dir, P_Ar_dblstar,n,alpha,beta):
    f_cf4       = np.asarray(f_cf4, dtype=float)
    denom       = f_cf4*K_Ar_dblstar_to_CF3_star + ((1 - f_cf4)) * K_Ar_dblstar_to_Ar_star
    frac        = f_cf4*K_Ar_dblstar_to_CF3_star / denom
    
    return  Pgamma_CF3_star_dir * alpha +  P_Ar_dblstar * frac * beta


def Pgamma_CF4_refined(f_cf4, Pgamma_CF4_plus_star_dir, P_Ar_3rd, n,kNoQ,kdis):
    
    f_cf4       = np.asarray(f_cf4, dtype=float)
    
    
    denom       = (1 / tau_3rd) + f_cf4 * n * (K_Ar3rd_to_CF4_plus_star + K_Ar3rd_to_Ar)
    numer       = f_cf4 * n * K_Ar3rd_to_CF4_plus_star
    frac        = np.where(denom == 0, 0, numer / denom)  # evitar divisiones por cero
    
    relajacion  = 1/(n*f_cf4*kNoQ+1) # NoQ -> No Autoquenching (si autoquenchea emite)
    disociacion = n*f_cf4/(n*f_cf4+kdis)
    
    return  Pgamma_CF4_plus_star_dir * relajacion * disociacion +  P_Ar_3rd * frac


def Pgamma_Ar3rd_refined(f_cf4, P_Ar_3rd, n):
    
    f_cf4       = np.asarray(f_cf4, dtype=float)
    denom       = (1 / tau_3rd) + f_cf4 * n * (K_Ar3rd_to_CF4_plus_star + K_Ar3rd_to_Ar)
    numer       = 1 / tau_3rd
    frac        = np.where(denom == 0, 0, numer / denom)
    
    return  P_Ar_3rd * frac * 0.4866





def Pgamma_vis_Cociente(f_cf4, Pgamma_CF3_star_dir, P_Ar_dblstar,n,values0,alpha,beta):
    f_cf40,Pgamma_CF3_star_dir0,P_Ar_dblstar0,n0=values0 
    a=Pgamma_CF3_refined(f_cf4, Pgamma_CF3_star_dir,P_Ar_dblstar,n,alpha,beta)
    b=Pgamma_CF3_refined(f_cf40, Pgamma_CF3_star_dir0,P_Ar_dblstar0,n0,alpha,beta)    
    return  a/b


def Pgamma_UV_Cociente(f_cf4, Pgamma_CF4_plus_star_dir, P_Ar_3rd, n,values0,kcool,kdis):
    f_cf40,Pgamma_CF4_plus_star_dir0,P_Ar_3rd0,n0=values0
    #n0=n
    a=Pgamma_CF4_refined(f_cf4, Pgamma_CF4_plus_star_dir, P_Ar_3rd, n,kcool,kdis)
    b=Pgamma_CF4_refined(f_cf40, Pgamma_CF4_plus_star_dir0, P_Ar_3rd0, n0,kcool,kdis)
    c=Pgamma_Ar3rd_refined(f_cf4, P_Ar_3rd, b)
    d=Pgamma_Ar3rd_refined(f_cf40, P_Ar_3rd0, n0)
    return  (a+c)/(b+d)

