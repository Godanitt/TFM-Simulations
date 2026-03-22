
import numpy as np
from scipy.interpolate import interp1d


from scipy.interpolate import PchipInterpolator


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
    W = 1 

    concentration = degrad_data["concentration"]
    Pob_N2 = degrad_data["N2_star"].to_numpy()
    Pob_Ar_meta = degrad_data["Ar_meta"].to_numpy()
    Pob_Ar_res = degrad_data["Ar_res"].to_numpy()
    Pob_Ar_dbleStar = degrad_data["Ar_dbleStar"].to_numpy()

    cols = ["N2_star", "Ar_meta", "Ar_res", "Ar_dbleStar"]
    Y = degrad_data[cols].to_numpy()   # shape: (n_puntos, 4)

    if len(fN2) > len(concentration):
        # Por si acaso, ordena x e y
        idx = np.argsort(concentration)
        xn = concentration[idx]
        y = Y[idx]

        interp = PchipInterpolator(xn, y, axis=0)
        Y_new = interp(fN2)
    else:
        Y_new = Y

    Pob_N2, Pob_Ar_meta, Pob_Ar_res, Pob_Ar_dbleStar = Y_new.T

    Nnorm               = x[0]

    tau_N2              = x[1]
    K_N2_Q_N2           = x[2]
    K_N2_Q_Ar           = x[3]

    P_N2                = x[4]

    K_ArMeta_Q_N2c      = x[5]
    K_ArMeta_Q_N2b      = x[6]
    K_ArMeta_Q_2Ar      = x[7]

    P_Ar2               = x[8]
    K_Ar2_Q_N2          = x[9]     
    tau_Ar2             = x[10]     

    K_ArRes_Q_N2c       = x[5]
    K_ArRes_Q_N2b       = x[6]  

    P_Ar22               = x[11]
    tau_meta_Ar2         = x[12]  
    K_ArRes_Q_2Ar       = x[13]     

    P_Ar_dbleStar_1     = x[14]
    P_Ar_dbleStar_2     = x[15]         

    tau_Ar_dbleStar     = x[16]
    K_Ar_dbleStar_Q_Ar  = x[17]       

    frac_1 = (1/tau_N2) / (1/tau_N2 + n * fN2 * K_N2_Q_N2 + n * (1 - fN2) * K_N2_Q_Ar)

    factor_N2 = frac_1

    frac_2 = (K_ArMeta_Q_N2b * fN2 * n) / ((K_ArMeta_Q_N2b + K_ArMeta_Q_N2c) * fN2 * n + (K_ArMeta_Q_2Ar * (1-fN2) * n**2))
    frac_3 = (K_ArMeta_Q_2Ar * (1-fN2) * n**2) / ((K_ArMeta_Q_N2b + K_ArMeta_Q_N2c) * fN2 * n + (K_ArMeta_Q_2Ar * (1-fN2) * n**2))
    frac_4 = (K_Ar2_Q_N2 * fN2 * n) / (K_Ar2_Q_N2 * fN2 * n + 1/tau_Ar2)

    factor_Ar_meta = frac_2 + P_Ar2 * frac_3 * frac_4

    frac_5 = (K_ArRes_Q_N2c * fN2 * n) / ((K_ArRes_Q_N2b+K_ArRes_Q_N2c) * fN2 * n + K_ArRes_Q_2Ar * (1-fN2) * n**2 )
    frac_6 = (K_ArRes_Q_2Ar * (1-fN2) * n**2) / ((K_ArMeta_Q_N2b + K_ArMeta_Q_N2c) * fN2 * n + (K_ArRes_Q_2Ar * (1-fN2) * n**2))
    frac_7 = (K_Ar2_Q_N2 * fN2 * n) / (K_Ar2_Q_N2 * fN2 * n + 1/tau_meta_Ar2)

    factor_Ar_res = frac_5 + P_Ar22 * frac_6 * frac_7

    frac_8 = (n * fN2) / (n * fN2 + tau_Ar_dbleStar + n * (1-fN2) * K_Ar_dbleStar_Q_Ar) 
    frac_9 = (n * (1-fN2)) / (n * fN2 * (1/K_Ar_dbleStar_Q_Ar) + tau_Ar_dbleStar  + n * (1-fN2) ) 

    factor_Ar_dbleStar_meta =  P_Ar_dbleStar_1 * frac_8 +  P_Ar_dbleStar_2 * frac_9 * (0.9 * factor_Ar_meta + 0.1 * factor_Ar_res)



    if activate_components:
        return (W * Nnorm * factor_N2 * (Pob_N2 * P_N2 + Pob_Ar_meta * factor_Ar_meta + Pob_Ar_res * factor_Ar_res + Pob_Ar_dbleStar * factor_Ar_dbleStar_meta),
                W * Nnorm * factor_N2 * (Pob_N2 * P_N2),
                W * Nnorm * factor_N2 * ( Pob_Ar_meta * factor_Ar_meta),
                W * Nnorm * factor_N2 * (Pob_Ar_res * factor_Ar_res),
                W * Nnorm * factor_N2 * (Pob_Ar_dbleStar * factor_Ar_dbleStar_meta)
        )
    else:
        return  W * Nnorm * factor_N2 * (Pob_N2 * P_N2 + Pob_Ar_meta * factor_Ar_meta + Pob_Ar_res * factor_Ar_res + Pob_Ar_dbleStar * factor_Ar_dbleStar_meta) 



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
