


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

name_CF4            =  lineal_CF4.columns.to_numpy()[:]
name_CF3            =  lineal_CF3.columns.to_numpy()[:]
name_Ar_dbleStar    =  lineal_pAr_dbleStar.columns.to_numpy()[:]
name_Ar_3rd         =  lineal_pAr_3rd.columns.to_numpy()[:]

name_yield_uv       =  yield_uv.columns.to_numpy()
name_yield_vis      =  yield_vis.columns.to_numpy()

print(name_CF3)

###################### Funciones para minimizar ################################################
