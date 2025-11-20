


import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.optimize as opt
import os 
import dill
import sys

PARENT_DIR = os.path.abspath(os.path.join(".."))    # proyecto/
sys.path.append(PARENT_DIR)
from Amoedo_Model_DivisionFit import Pgamma_UV_Cociente, Pgamma_vis_Cociente,Pgamma_CF3_refined,Pgamma_CF4_refined,Pgamma_Ar3rd_refined

""" 
Ajuste del cociente con datos de Degrad con ajuste lineal
"""

#############################################################################################################
######################## Lectura de las poblaciones de Degrad ################################################

# === Carpeta donde están los pickles ===
DATA_DIR = os.path.join("..", "pickle_data")
yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

DATA_DIR = os.path.join("..", "pickle_data")
pCF4         = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
pCF3         = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF3.pkl"))
pArDbleStar  = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_dbleStar.pkl"))
pAr3rd       = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_3rd.pkl"))

######################## Lectura de las poblaciones de Degrad ################################################

fCF4            =  pCF4["fCF4"].to_numpy()

name_CF4            =  pCF4.columns.to_numpy()[1::2]
name_CF3            =  pCF3.columns.to_numpy()[1::2]
name_Ar_dbleStar    =  pArDbleStar.columns.to_numpy()[1::2]
name_Ar_3rd         =  pAr3rd.columns.to_numpy()[1::2]

name_yield_uv       =  yield_uv.columns.to_numpy()
name_yield_vis      =  yield_vis.columns.to_numpy()

###################### Inicialización de variables de interés ################################################

fCF4_real = yield_uv["fCF4 real"].to_numpy()
err_fCF4_real = yield_uv["Err fCF4 real"]


###################### Inicializamos los dataFrames de los parámetors para exportar ################################################

# Crear las columnas para todas las combinaciones
columns_alpha = []
columns_kcool = []
columns_kdis = []

for nameCF3 in name_CF3:
    for nameArdbleStar in name_Ar_dbleStar:
        col = f"{nameCF3}_{nameArdbleStar}"
        columns_alpha.append(col)
        columns_kcool.append(col)
        columns_kdis.append(col)

# Inicializar dataframes con una sola fila (que luego irás rellenando)
df_alpha = pd.DataFrame(columns=columns_alpha, index=[0])
df_kcool = pd.DataFrame(columns=columns_kcool, index=[0])
df_kdis  = pd.DataFrame(columns=columns_kdis,  index=[0])

index = -1

###################### Funciones para minimizar ################################################


def minimize_vis(x,PCF3, PAr_dbleStar, y_vis, values0, values0_yield):
    alpha = x[0]
    A = Pgamma_vis_Cociente(fCF4_real/100, PCF3, PAr_dbleStar, values0, alpha)
    chi2 = 0
    
    
    for col in y_vis.columns:
        if not("fCF4" in col) and not("Err" in col):

            
            values0_yield = y_vis[col].to_numpy()[index]
            err =  y_vis["Err " + col].to_numpy()[index]
            
            B = y_vis[col].to_numpy()
            sB = np.sqrt((y_vis["Err " + col] / values0_yield)**2 +(err * y_vis[col].to_numpy() / (values0_yield**2))**2            )
            #sB = 0.1 * B        
                
            if values0_yield!=0:
                for a,b,sb in zip(A,B,sB):
                    if b!=0 and sb!=0:
                        
                        chi2 += (((a - b/values0_yield)**2)/(sb**2))
                    
    return (chi2)

def minimize_uv(x, Pgamma_CF4_plus_star_dir, P_Ar_3rd, y_uv, values0, values0_yield):
    kcool, kdis = x

    # Modelo teórico para todos los fCF4_real
    chi2 = 0.0
    for col in y_uv.columns:
        if not ("fCF4" in col) and not ("Err" in col):
            
            n = float(col.replace("bar", ""))
            
            A = Pgamma_UV_Cociente(
                fCF4_real/100,
                Pgamma_CF4_plus_star_dir,
                P_Ar_3rd,
                n,
                values0,
                kcool,
                kdis
            )
            
            B  = y_uv[col].to_numpy()
                    
            values0_yield = y_uv[col].to_numpy()[index]
            err =  y_uv["Err " + col].to_numpy()[index]
            
            B = y_uv[col].to_numpy()
            sB = np.sqrt((y_uv["Err " + col] / values0_yield)**2 +(err * y_uv[col].to_numpy() / (values0_yield**2))**2            )
            #sB = 0.1 * B    
            
            if values0_yield!=0:
                for a, b, sb in zip(A, B, sB):
                    if b != 0 and sb != 0:
                        # Misma lógica que en el visible: normalizas por el valor a 1 bar
                        chi2 += ((a - b/values0_yield)**2)/(sb**2)

    return chi2


###################### Minimización del visible ################################################

for nameCF3 in name_CF3:
    for nameArdbleStar in ["Ar** all"]:
            
            PCF3 = np.array([])
            PAr_dbleStar = np.array([])
            
            for i in fCF4_real/100:
                # Encuentra el índice j tal que fCF4[j] <= i < fCF4[j+1]
                for j in range(len(fCF4)-1):
                    if i == fCF4[j] or i == fCF4[j+1] or (fCF4[j] < i < fCF4[j+1]) or (0 <= i < fCF4[j]):
                        
                        # Valores en j y j+1
                        y1_CF3 = pCF3[nameCF3].loc[j]
                        y2_CF3 = pCF3[nameCF3].loc[j+1]
                        y1_Ar  = pArDbleStar[nameArdbleStar].loc[j]
                        y2_Ar  = pArDbleStar[nameArdbleStar].loc[j+1]

                        # Interpolación lineal (vale también para los casos "i == fCF4[j]")
                        frac = (i - fCF4[j]) / (fCF4[j+1] - fCF4[j])

                        PCF3         = np.append(PCF3, y1_CF3 + frac * (y2_CF3 - y1_CF3))
                        PAr_dbleStar = np.append(PAr_dbleStar, y1_Ar  + frac * (y2_Ar  - y1_Ar))

                        break


            print(fCF4)      
            print(fCF4_real/100)           
            # Seleccionando -1 seleccionamos 100% CF4
            values0_vis = (fCF4_real[index]/100,PCF3[index],PAr_dbleStar[index])
            values0_yield = yield_vis["1.0bar"].to_numpy()[index]
            

            args = (PCF3, PAr_dbleStar, yield_vis, values0_vis, values0_yield)
            result = opt.differential_evolution(minimize_vis, 
                                            [(0,1)], 
                                            args=args,
                                            maxiter=100,        # 5x más rápido
                                            popsize=100)        # población más pequeña)
            alpha=(result.x)
            chi2=(result.fun)

            print("="*60)
            print("alpha=",alpha)
            print("chi²=",chi2)
            print("="*60) 
            
            
           
            colname = f"{nameCF3}_{nameArdbleStar}"
            df_alpha.loc[0, colname] = float(alpha)

###################### Minimización del ultravioleta ################################################


for nameCF4 in name_CF4:
    for nameAr3rd in name_Ar_3rd:

            
        PCF4 = np.array([])
        PAr_3rd = np.array([])
            
        for i in fCF4_real/100:
            for j in range(len(fCF4)-1):

                if (
                    i == fCF4[j] or
                    i == fCF4[j+1] or
                    (fCF4[j] < i < fCF4[j+1]) or
                    (0 <= i < fCF4[j])
                ):
                    # Valores en j y j+1
                    y1_CF4 = pCF4[nameCF4].loc[j]
                    y2_CF4 = pCF4[nameCF4].loc[j+1]

                    y1_Ar  = pAr3rd[nameAr3rd].loc[j]
                    y2_Ar  = pAr3rd[nameAr3rd].loc[j+1]

                    # Interpolación lineal (cubriendo todos los casos)
                    frac = (i - fCF4[j]) / (fCF4[j+1] - fCF4[j])

                    PCF4      = np.append(PCF4,      y1_CF4 + frac * (y2_CF4 - y1_CF4))
                    PAr_3rd   = np.append(PAr_3rd,   y1_Ar  + frac * (y2_Ar  - y1_Ar))

                    break

                        
        # Valores de referencia para el cociente (punto 100% CF4, mismo criterio que en el visible)
        
        n0    =  1
        values0_uv     = (fCF4_real[index]/100, PCF4[index], PAr_3rd[index],n0)
        values0_yield  = yield_uv["%.1fbar"%n0].to_numpy()[index]

        # Minimización respecto a kcool y kdis
        args_uv = (PCF4, PAr_3rd, yield_uv, values0_uv, values0_yield)
        result_uv = opt.differential_evolution(
            minimize_uv,
            bounds=[(0.0, 1000.0),   # límites para kcool (ajusta)
                    (0.0, 1000.0)],  # límites para kdis  (ajusta)
            args=args_uv,
            maxiter=200,
            popsize=200
        )

        kcool_fit, kdis_fit = result_uv.x

        print("="*60)
        print("UV fit ->", nameCF4, "/", nameAr3rd)
        print("kcool =", kcool_fit)
        print("kdis  =", kdis_fit)
        print("chi2  =", result_uv.fun)
        print("="*60)
        
        
        colname = f"{nameCF4}_{nameAr3rd}"
        df_kcool.loc[0, colname] = float(kcool_fit)
        df_kdis.loc[0, colname]  = float(kdis_fit)

###################### Guardamos la información en pickle y csv ################################################

        
df_alpha.to_pickle("pickle_data/alpha_results.pkl")
df_kcool.to_pickle("pickle_data/kcool_results.pkl")
df_kdis.to_pickle("pickle_data/kdis_results.pkl")

df_alpha.to_csv("csv_data/alpha_results.csv")
df_kcool.to_csv("csv_data/kcool_results.csv")
df_kdis.to_csv("csv_data/kdis_results.csv")
