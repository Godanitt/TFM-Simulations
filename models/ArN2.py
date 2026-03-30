
import numpy as np
from scipy.interpolate import interp1d


from scipy.interpolate import PchipInterpolator

def W_ArN2(xN2, WAr=26.4, WN2=34.8):
    return 1.0 / ((1.0-xN2)/WAr + xN2/WN2)


def theory_yield_N2_uv(x, degrad_data, fN2, n, activate_components = False):
    fN2 = np.asarray(fN2, dtype=float)
    W = W_ArN2(fN2)

    concentration = degrad_data["concentration"]
    Pob_N2 = degrad_data["N2_star"].to_numpy()
    Pob_Ar_meta = degrad_data["Ar_meta"].to_numpy()
    Pob_Ar_res = degrad_data["Ar_res"].to_numpy()
    Pob_Ar_dbleStar = degrad_data["Ar_dbleStar"].to_numpy()

    cols = ["N2_star", "Ar_meta", "Ar_res", "Ar_dbleStar"]
    Y = degrad_data[cols].to_numpy()   # shape: (n_puntos, 4)
    
    Y = np.asarray(Y, dtype=float)
    fN2 = np.asarray(fN2, dtype=float)

    if len(fN2) > len(concentration):
        conc = np.asarray(concentration, dtype=float)
        yvals = np.asarray(Y, dtype=float)

        idx = np.argsort(conc)
        conc_sorted = conc[idx]
        y_sorted = yvals[idx]

        interp = PchipInterpolator(conc_sorted, y_sorted)
        Y_interp = interp(fN2)
    else:
        Y_interp = Y

    Pob_N2, Pob_Ar_meta, Pob_Ar_res, Pob_Ar_dbleStar = Y_interp.T

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
        return ((1/W)* Nnorm * factor_N2 * (Pob_N2 * P_N2 + Pob_Ar_meta * factor_Ar_meta + Pob_Ar_res * factor_Ar_res + Pob_Ar_dbleStar * factor_Ar_dbleStar_meta),
                (1/W)* Nnorm * factor_N2 * (Pob_N2 * P_N2),
                (1/W)* Nnorm * factor_N2 * ( Pob_Ar_meta * factor_Ar_meta),
                (1/W)* Nnorm * factor_N2 * (Pob_Ar_res * factor_Ar_res),
                (1/W)* Nnorm * factor_N2 * (Pob_Ar_dbleStar * factor_Ar_dbleStar_meta)
        )
    else:
        return  (1/W)* Nnorm * factor_N2 * (Pob_N2 * P_N2 + Pob_Ar_meta * factor_Ar_meta + Pob_Ar_res * factor_Ar_res + Pob_Ar_dbleStar * factor_Ar_dbleStar_meta) 



############################
## VERSION ANTIGUA

def _theory_yield_N2_uv(x, degrad_data, fN2, n, activate_components = False):
    fN2 = np.asarray(fN2, dtype=float)

    W = W_ArN2(fN2)# (1/ion_potential(fN2))


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
    
    if activate_components:
        return ((1/W) * N * (P_N2 + (P_Ar_Star + P_Ar_dbleStar * frac3) * frac1) * frac2,
                (1/W) * N * P_N2 * frac2)
    else:
        return (1/W)* N * (P_N2 + (P_Ar_Star + P_Ar_dbleStar * p_dbleStar * frac3) * p_Star * frac1) * frac2
