
import numpy as np 
import pandas as pd 
import dill
import scipy.special
import importlib
import matplotlib.pyplot as plt 

"""
Script que nos permite visualizar los ajustes lineales de poblacion-%mezcla
"""

#############################################################################################################
######################## Lectura de las poblaciones de Degrad ################################################

DATA_DIR = "pickle_data"
poblations_CF4 = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
poblations_CF3 = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_CF3.pkl"))
poblations_Ar_dbleStar  = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_Ar_dbleStar.pkl"))
poblations_Ar_3rd  = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_Ar_3rd.pkl"))

lineal_poblations_CF4 = pd.read_pickle(os.path.join(DATA_DIR, "linealFun_poblations_CF4.pkl"))
lineal_poblations_CF3 = pd.read_pickle(os.path.join(DATA_DIR,  "linealFun_poblations_CF3.pkl"))
lineal_poblations_Ar_dbleStar  = pd.read_pickle(os.path.join(DATA_DIR,  "linealFun_poblations_Ar_dbleStar.pkl"))
lineal_poblations_Ar_3rd  = pd.read_pickle(os.path.join(DATA_DIR,  "linealFun_poblations_Ar_3rd.pkl"))


