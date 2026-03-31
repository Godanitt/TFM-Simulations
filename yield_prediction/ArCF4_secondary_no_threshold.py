import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from scipy.interpolate import PchipInterpolator
import scienceplots
plt.style.use('default')

models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../models'))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))

sys.path.append(models_dir)
sys.path.append(data_dir)

from ArCF4 import *
from read_Degrad import read_degrad
from read_experimental import read_experimental
from read_Root import export_hlevels_to_csv,read_data_per_primary_electron
from read_secondary import read_garfield_csv_folder

# Carpeta donde están los ROOT
folder_path = "../data/Secondary_GarfieldData/ArCF4/root"
table_path = "../data/Secondary_GarfieldData/levels/ArCF4_level_data.csv"

# 1) Exportar hLevels a CSV
export_hlevels_to_csv(
    folder_path,
    table_path,
    object_name="hLevels",
    argon_update=True)

# 2) Leer ne y ni, sacar estadísticas y guardar histogramas
summary = read_data_per_primary_electron(folder_path)

##########



config = pd.DataFrame({
    "CF4": {
        "name principal": ["ION"],
        "gas": "CF4",
        "energy low": 15.5,
        "energy up": 16,
        "name output": "CF4",
        "type": "ionisation"
    },
    "Ar**": {
        "name principal": ["EXC"],
        "gas": "Ar",
        "energy low": 11.8,
        "energy up": 100,
        "name output": "Ar_dbleStar",
        "type": "excitation"
    },
    "CF3": {
        "name principal": ["NEUTRAL DISS"],
        "gas": "CF4",
        "energy low": 15.6,
        "energy up": 100,
        "name output": "CF3",
        "type": "inelastic"
    },
    "Ar3rd": {
        "name principal": ["IONISATION"],
        "gas": "Ar",
        "energy low": 40,
        "energy up": 120,
        "name output": "Ar_3rd",
        "type": "ionisation"
    }
})

garfield = read_garfield_csv_folder(
    folder_path="../data/Secondary_GarfieldData/ArCF4/csv",
    dataframe=config,
    output_dir="../data/Secondary_GarfieldData/ArCF4/populations",
    output_general_name="../data/Secondary_GarfieldData/ArCF4/populations/ArCF4_secondary",
    concentration_gas="cf4"
)



##########


DATA_DIR_EXP = "../data/Experimental/ArCF4/"
DATA_DIR_GARFIELD = "../data/Secondary_GarfieldData/ArCF4/populations"
DATA_DIR_PAR = "../data/Parameters"

#yield_N2_uv  = pd.read_csv(os.path.join(DATA_DIR_EXP, "yield_ArCF4.csv"))

garfield_data = pd.read_csv(os.path.join(DATA_DIR_GARFIELD, "ArCF4_secondary.csv"))
garfield_data["concentration"] = garfield_data["concentration"]/100


parameter_data = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArCF4_primary.csv"))["parameter"].to_numpy()
print(parameter_data)
parameter_data[0] = 1
print(parameter_data[-1])
parameter_data[-1] = parameter_data[-1] * 0.8
parameter_data[1] = 0.35
parameter_data[2]  = 0.35
print(parameter_data)


fN2 = np.logspace(-3,0,1000)
electric_fields = [40,60,70,80,90,110]

####

ne = 94    

####

print(summary)

plt.figure(figsize=(6,4))
plt.style.use('science')

cmap = "viridis"
cmap_obj = plt.get_cmap(cmap)
colors = cmap_obj(np.linspace(0.15, 0.85, len(electric_fields)))

for i,electric_field in enumerate(electric_fields):

    yield_teo = (theory_yield_vis(parameter_data, garfield_data[garfield_data["electric_field"]==electric_field], fN2, 1)) / (ne * 500) * ion_potential(fN2)

    plt.plot(fN2 * 100,
        yield_teo,
        color=colors[i],
        label=f"{electric_field} kV/cm prediction"
    )

plt.errorbar([100,67,10,5],[0.1,0.3,0.36,0.33],yerr = np.array([0.1,0.3,0.36,0.33])*0.2,fmt=".")
plt.title("1 bar 100 gain secondary visible yield prediction for Ar/CF4 mixture")
plt.xscale("log")
#plt.yscale("log")
plt.ylabel("ph/e$^-$")
plt.xlabel("CF$_4$ concetration [%]")
plt.xlim(1,110)
plt.legend()
plt.savefig("plots/ArCF4_vis_secondary_no_threshold.pdf")

####

plt.figure(figsize=(6,4))
plt.style.use('science')

cmap = "viridis"
cmap_obj = plt.get_cmap(cmap)
colors = cmap_obj(np.linspace(0.15, 0.85, len(electric_fields)))

for i,electric_field in enumerate(electric_fields):

    yield_teo_uv = theory_yield_uv(parameter_data, garfield_data[garfield_data["electric_field"]==electric_field], fN2, 1)  / (ne * 500) * ion_potential(fN2)  
    yield_teo_vis = (theory_yield_vis(parameter_data, garfield_data[garfield_data["electric_field"]==electric_field], fN2, 1)) / (ne * 500) * ion_potential(fN2)
    yield_teo = yield_teo_uv/15 + yield_teo_vis/5

    #yield_teo = (1/7.2) * theory_yield_uv(parameter_data, garfield_data[garfield_data["electric_field"]==electric_field], fN2, 1)  / (ne * 500) * ion_potential(fN2)  

    plt.plot(fN2 * 100,
        yield_teo,
        color=colors[i],
        label=f"{electric_field} kV/cm prediction"
    )

plt.errorbar([100,67,10,5],[0.04,0.045,0.085,0.07],yerr = np.array([0.04,0.045,0.085,0.07])*0.2,fmt=".")
plt.title("1 bar 100 gain secondary UV yield prediction for Ar/CF4 mixture")
plt.xscale("log")
#plt.yscale("log")
plt.ylabel("ph/e$^-$")
plt.xlim(1,110)
plt.ylim(0.01,0.12)
plt.xlabel("CF$_4$ concetration [\%]")
plt.legend()
plt.savefig("plots/ArCF4_uv_secondary_no_threshold.pdf")