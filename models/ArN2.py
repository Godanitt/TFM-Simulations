
import numpy as np



# % de CF4 en Ar
cf4_pct = np.array([0, 1.0, 2.0, 5.0, 10, 20, 30, 50, 75, 100])/100

# Potencial de ionización (según la columna Ar/CF4)
ion_pot = np.array([26.4, 26.7, 26.9, 27.4, 28.1, 29.4, 30.2, 31.7, 33.0, 34.3])

def ion_potential(f):
    f_cf4 = np.asarray(f, dtype=float)
    W=np.interp(f_cf4,cf4_pct,ion_pot)
    return W


def theory_yield_N2_uv(x, degrad_data, fN2, n, activate_components = False):
    fN2 = np.asarray(fN2, dtype=float)

    concentration = degrad_data["concentration"]
    P_N2 = degrad_data["N2_star"].to_numpy()
    P_Ar_Star = degrad_data["Ar_Star"].to_numpy()
    P_Ar_dbleStar = degrad_data["Ar_dbleStar"].to_numpy()

    N           = x[0]
    p_Star      = x[1]
    K1          = x[2]
    K2          = x[3]
    K3          = x[4]
    tau_emision = x[5]
    p_dbleStar  = x[6]
    K4          = x[7]
    K5          = x[8]

    if len(fN2)>len(P_N2):
        P_N2 = np.interp(fN2,concentration,P_N2)
        P_Ar_Star =  np.interp(fN2,concentration,P_Ar_Star)
        P_Ar_dbleStar =  np.interp(fN2,concentration,P_Ar_dbleStar)

    denom = n * fN2 * K1 + (n**2) * (1.0 - fN2) * K2 
    frac1  = np.where(denom == 0, 0.0, K1 * n * fN2 / denom)
     

    denom = 1 + fN2 * n * tau_emision + (1.0 - fN2) * n * K3
    frac2  = np.where(denom == 0, 0.0, 1 / denom)

    denom = (fN2) * n * K4 + 1/30 + (1-fN2) * n * K5
    frac3  = np.where(denom == 0, 0.0, (1/30)/ denom)
    
    W = 1 # (1/ion_potential(fN2))

    if activate_components:
        return (W * N * (P_N2 + (P_Ar_Star + P_Ar_dbleStar * p_dbleStar * frac3) * p_Star * frac1) * frac2,
                W * N * P_N2 * frac2)
    else:
        return W * N * (P_N2 + (P_Ar_Star + P_Ar_dbleStar * p_dbleStar * frac3) * p_Star * frac1) * frac2

############################
## VERSION ANTIGUA

def _theory_yield_N2_uv(x, degrad_data, fN2, n, activate_components = False):
    fN2 = np.asarray(fN2, dtype=float)

    concentration = degrad_data["concentration"]
    P_N2 = degrad_data["N2_star"].to_numpy()
    P_Ar_Star = degrad_data["Ar_Star"].to_numpy()
    P_Ar_dbleStar = degrad_data["Ar_dbleStar"].to_numpy()

    N           = x[0]
    p_Star      = x[1]
    K1          = x[2]
    K2          = x[3]
    K3          = x[4]
    tau_emision = x[5]
    p_dbleStar  = x[6]
    K4          = x[7]
    K5          = x[8]

    if len(fN2)>len(P_N2):
        P_N2 = np.interp(fN2,concentration,P_N2)
        P_Ar_Star =  np.interp(fN2,concentration,P_Ar_Star)
        P_Ar_dbleStar =  np.interp(fN2,concentration,P_Ar_dbleStar)

    denom = n * fN2 * K1 + (n**2) * (1.0 - fN2) * K2 
    frac1  = np.where(denom == 0, 0.0, K1 * n * fN2 / denom)
     

    denom = 1 + fN2 * n * tau_emision + (1.0 - fN2) * n * K3
    frac2  = np.where(denom == 0, 0.0, 1 / denom)

    denom = (fN2) * n * K4 + 1/30 + (1-fN2) * n * K5
    frac3  = np.where(denom == 0, 0.0, (1/30)/ denom)
    
    W = 1 # (1/ion_potential(fN2))

    if activate_components:
        return (W * N * (P_N2 + (P_Ar_Star + P_Ar_dbleStar * frac3) * frac1) * frac2,
                W * N * P_N2 * frac2)
    else:
        return W * N * (P_N2 + (P_Ar_Star + P_Ar_dbleStar * p_dbleStar * frac3) * p_Star * frac1) * frac2
