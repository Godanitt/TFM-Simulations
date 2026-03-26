import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns

models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../models'))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))

sys.path.append(models_dir)
sys.path.append(data_dir)

from ArN2 import *
from read_Degrad import read_degrad
from read_experimental import read_experimental
from read_Root import export_hlevels_to_csv,read_data_per_primary_electron
from read_secondary import read_garfield_csv_folder

# Carpeta donde están los ROOT
folder_path = "../data/Secondary_GarfieldData/ArN2/root"
table_path = "../data/Secondary_GarfieldData/levels/ArN2_level_data.csv"

# 1) Exportar hLevels a CSV
export_hlevels_to_csv(
    folder_path,
    table_path,
    object_name="hLevels",
    argon_update=True)

# 2) Leer ne y ni, sacar estadísticas y guardar histogramas
summary = read_data_per_primary_electron(folder_path)
print(summary)

##########


config = pd.DataFrame({
    "Ar_meta": {
        "name principal": ["EXC"],
        "gas": "Ar",
        "energy low": 11.0,
        "energy up": 11.6,
        "name output": "Ar_meta",
        "type": "excitation"
    },
    "Ar_res": {
        "name principal": ["EXC"],
        "gas": "Ar",
        "energy low": 11.6,
        "energy up": 11.7,
        "name output": "Ar_res",
        "type": "excitation"
    },
    "Ar**": {
        "name principal": ["EXC"],
        "gas": "Ar",
        "energy low": 11.7,
        "energy up": 100,
        "name output": "Ar_dbleStar",
        "type": "excitation"
    },
    "N2*": {
        "name principal": ["C 3PI"],
        "gas": "N2",
        "energy low": 11.0,
        "energy up": 100,
        "name output": "N2_star",
        "type": "excitation"
    }
})

garfield = read_garfield_csv_folder(
    folder_path="../data/Secondary_GarfieldData/ArN2/csv",
    dataframe=config,
    output_dir="../data/Secondary_GarfieldData/ArN2/populations",
    output_general_name="../data/Secondary_GarfieldData/ArN2/populations/ArN2_secondary",
    concentration_gas="n2"
)



##########


DATA_DIR_EXP = "../data/Experimental/ArN2/"
DATA_DIR_DEGRAD = "../data/Secondary_GarfieldData/ArN2/populations"
DATA_DIR_PAR = "../data/Parameters"

yield_N2_uv  = pd.read_csv(os.path.join(DATA_DIR_EXP, "yield_N2.csv"))

degrad_data = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArN2.csv"))

parameter_data = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArN2_primary.csv"))["parameter"].to_numpy()
print(parameter_data)
parameter_data[0] = 1
print(parameter_data)