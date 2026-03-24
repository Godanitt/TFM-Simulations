
import numpy as np
import numpy as np
from scipy.interpolate import interp1d


from scipy.interpolate import PchipInterpolator



# Santorilli and Diego et al doi 10.1140/epjc/s10052-021-09375-3 (en principio igual que el de Pablo Amoedo)
tau_3rd                      = 5.02               
tercer_continuo              = 0.4866

# % de CF4 en Ar
cf4_pct = np.array([0, 1.0, 2.0, 5.0, 10, 20, 30, 50, 75, 100])/100

# Potencial de ionización (según la columna Ar/CF4)
ion_pot = np.array([26.4, 26.7, 26.9, 27.4, 28.1, 29.4, 30.2, 31.7, 33.0, 34.3])

def ion_potential(f):
    f_cf4 = np.asarray(f, dtype=float)
    W=np.interp(f_cf4,cf4_pct,ion_pot)
    return W

def theory_yield_vis(x, degrad_data, fCF4, n,activate_components = False):
    f_cf4 = np.asarray(fCF4, dtype=float)

    concentration = degrad_data["concentration"]

    cols = ["CF3", "Ar_dbleStar", "CF4", "Ar_3rd"]
    Y = degrad_data[cols].to_numpy()   # shape: (n_puntos, 4)

    if len(f_cf4) > len(concentration):
        # Por si acaso, ordena x e y
        idx = np.argsort(concentration)
        xn = concentration[idx]
        y = Y[idx]

        interp = PchipInterpolator(xn, y, axis=0)
        Y_new = interp(f_cf4)
    else:
        Y_new = Y

    P_CF3, P_Ar_dbleStar, P_CF4, P_Ar_3rd = Y_new.T


    N           = x[0]
    p_CF3       = x[1]
    p_DbleStar  = x[2]
    K           = x[3]
    
    K2           = x[9]

    denom = n * f_cf4 * K2 + n * (1.0 - f_cf4) * K + 1/30
    frac  = np.where(denom == 0, 0.0, K2 * n * f_cf4 / denom)
    
    

    # OJO: aquí faltaba un "*" en tu ejemplo: p_CF3(P_CF3 + ...) → p_CF3 * (...)
    return (1/ion_potential(f_cf4))*N*(p_CF3 * P_CF3 + frac * p_DbleStar * P_Ar_dbleStar)


def theory_yield_uv(x, degrad_data, fCF4, n, activate_components = False):
    f_cf4 = np.asarray(fCF4, dtype=float)


    concentration = degrad_data["concentration"]

    cols = ["CF3", "Ar_dbleStar", "CF4", "Ar_3rd"]
    Y = degrad_data[cols].to_numpy()   # shape: (n_puntos, 4)

    if len(f_cf4) > len(concentration):
        # Por si acaso, ordena x e y
        idx = np.argsort(concentration)
        xn = concentration[idx]
        y = Y[idx]

        interp = PchipInterpolator(xn, y, axis=0)
        Y_new = interp(f_cf4)
    else:
        Y_new = Y

    P_CF3, P_Ar_dbleStar, P_CF4, P_Ar_3rd = Y_new.T

    N      = x[0]
    K1     = x[4]
    K2     = x[5]
    p_CF3  = x[6]
    K3     = x[7]
    K4     = x[8]

    # frac1 = nf / (nf + K1)
    numer = f_cf4 * n
    denom = f_cf4 * n + K1
    frac1 = np.where(denom == 0, 0.0, numer / denom)

    # frac2 = 1 / (1 + K2 n f_cf4)
    numer = 1.0
    denom = 1.0 + K2 * n * f_cf4
    frac2 = np.where(denom == 0, 0.0, numer / denom)

    # frac3
    denom = (1.0 / tau_3rd) + f_cf4 * n * (K3)
    numer = f_cf4 * n * K3 
    frac3 = np.where(denom == 0, 0.0, numer / denom)

    # frac4
    denom = (1.0 / tau_3rd) + f_cf4 * n * (K3)
    numer = 1.0 / tau_3rd
    frac4 = np.where(denom == 0, 0.0, numer / denom)

    return (1/ion_potential(f_cf4))* N * ((frac1 * frac2) * (p_CF3 * P_CF4 + frac3 * P_Ar_3rd * K4)
        + tercer_continuo * frac4 * P_Ar_3rd )
