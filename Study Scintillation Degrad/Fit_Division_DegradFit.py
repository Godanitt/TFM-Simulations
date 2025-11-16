


from Amoedo_Model_DivisionFit import Pgamma_UV_Cociente,Pgamma_vis_Cociente
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.optimize as opt
import os 
import dill

""" 
Ajuste del cociente con datos de Degrad con ajuste lineal
"""

#############################################################################################################
######################## Lectura de las poblaciones de Degrad ################################################

# === Carpeta donde están los pickles ===
DATA_DIR = "pickle_data"

with open(os.path.join(DATA_DIR, "linealFun_poblations_Ar_3rd.pkl"), "rb") as f:
    lineal_pAr_3rd = dill.load(f)

with open(os.path.join(DATA_DIR, "linealFun_poblations_Ar_dbleStar.pkl"), "rb") as f:
    lineal_pAr_dbleStar = dill.load(f)

with open(os.path.join(DATA_DIR, "linealFun_poblations_CF4.pkl"), "rb") as f:
    lineal_CF4 = dill.load(f)

with open(os.path.join(DATA_DIR, "linealFun_poblations_CF3.pkl"), "rb") as f:
    lineal_CF3 = dill.load(f)

yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

######################## Lectura de las poblaciones de Degrad ################################################

name_CF4            =  lineal_CF4.columns.to_numpy()[::2]
name_CF3            =  lineal_CF3.columns.to_numpy()[::2]
name_Ar_dbleStar    =  lineal_pAr_dbleStar.columns.to_numpy()[::2]
name_Ar_3rd         =  lineal_pAr_3rd.columns.to_numpy()[::2]

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

###################### Funciones para minimizar ################################################


def minimize_vis(x,PCF3, PAr_dbleStar, y_vis, values0, values0_yield):
    alpha = x[0]
    A = Pgamma_vis_Cociente(fCF4_real/100, PCF3, PAr_dbleStar, values0, alpha)
    chi2 = 0
    
    
    for col in y_vis.columns:
        if not("fCF4" in col) and not("Err" in col):

            index=-1
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

            index=-1
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
            f = lineal_CF3[nameCF3].to_numpy()[0]
            g = lineal_pAr_dbleStar[nameArdbleStar].to_numpy()[0]
            
            PCF3 = f(fCF4_real/100)
            PAr_dbleStar = g(1-fCF4_real/100)
            
            
            # Seleccionando -1 seleccionamos 100% CF4
            index = -1
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
            
            fCF4_range=np.linspace(min(fCF4_real),max((fCF4_real)),1000)
            PCF3 = f(fCF4_range/100)
            PAr_dbleStar = g(1-fCF4_range/100)
           
            colname = f"{nameCF3}_{nameArdbleStar}"
            df_alpha.loc[0, colname] = float(alpha)

###################### Minimización del ultravioleta ################################################


for nameCF4 in name_CF4:
    for nameAr3rd in name_Ar_3rd:
        f = lineal_CF4[nameCF4].to_numpy()[0]
        g = lineal_pAr_3rd[nameAr3rd].to_numpy()[0]
        
        # Poblaciones / probabilidades en los puntos experimentales
        PCF4   = f(fCF4_real/100)
        PAr_3rd = g(1-fCF4_real/100)

        # Valores de referencia para el cociente (punto 100% CF4, mismo criterio que en el visible)
        index =  -1
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
            maxiter=100,
            popsize=100
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
