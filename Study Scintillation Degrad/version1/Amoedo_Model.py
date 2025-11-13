import numpy as np

""" 
Modelo de pablo puro.
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

def Pgamma_CF3(f_cf4, Pgamma_CF3_star_dir, P_Ar_dblstar):
    f_cf4 = np.asarray(f_cf4, dtype=float)
    denom = K_Ar_dblstar_to_CF3_star*f_cf4 + ((1 - f_cf4)) * K_Ar_dblstar_to_Ar_star
    frac = K_Ar_dblstar_to_CF3_star*f_cf4 / denom
    return Pgamma_CF3_star_dir + P_Ar_dblstar * frac


def Pgamma_CF4(f_cf4, Pgamma_CF4_plus_star_dir, P_Ar_3rd, n):
    f_cf4 = np.asarray(f_cf4, dtype=float)
    denom = (1 / tau_3rd) + f_cf4 * n * (K_Ar3rd_to_CF4_plus_star + K_Ar3rd_to_Ar)
    numer = f_cf4 * n * K_Ar3rd_to_CF4_plus_star
    # evitar divisiones por cero
    frac = np.where(denom == 0, 0, numer / denom)
    return Pgamma_CF4_plus_star_dir + P_Ar_3rd * frac


def Pgamma_Ar3rd(f_cf4, P_Ar_3rd, n):
    f_cf4 = np.asarray(f_cf4, dtype=float)
    denom = (1 / tau_3rd) + f_cf4 * n * (K_Ar3rd_to_CF4_plus_star + K_Ar3rd_to_Ar)
    numer = 1 / tau_3rd
    frac = np.where(denom == 0, 0, numer / denom)
    return P_Ar_3rd * frac

