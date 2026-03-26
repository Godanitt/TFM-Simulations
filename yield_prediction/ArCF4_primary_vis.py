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

from ArCF4 import *


######################################33

DATA_DIR_EXP = "../data/Experimental/ArCF4/"
DATA_DIR_DEGRAD = "../data/Primary_DegradData"
DATA_DIR_PAR = "../data/Parameters"

yield_N2_uv  = pd.read_csv(os.path.join(DATA_DIR_EXP, "vis.csv"))

degrad_data = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArCF4.csv"))

parameter_data = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArCF4_primary.csv"))["parameter"].to_numpy()
print(parameter_data)
parameter_data[0] = 1

######################################33


cf4_red_E0 = [0.2, 0.4, 0.7, 1.0, 2.0, 7.0, 10.0]
y_red_E0   = [500, 700, 1050, 1450, 1950, 2400, 2550]
yerr_red_E0= [70, 70, 80, 100, 120, 160, 170]

# [400–750] nm, E = 100 V/cm (rojo abierto)
cf4_red_E100 = [0.2, 0.4, 0.7, 1.0, 2.0, 7.0, 10.0]
y_red_E100   = [450, 500, 600, 1150, 1300, 1850, 1900]
yerr_red_E100= [60, 60, 60, 90, 100, 120, 120]



cf4_pct = np.array([0, 1.0, 2.0, 5.0, 10, 20, 30, 50, 75, 100])/100

# Potencial de ionización (según la columna Ar/CF4)
ion_pot = np.array([26.4, 26.7, 26.9, 27.4, 28.1, 29.4, 30.2, 31.7, 33.0, 34.3])

def ion_potential(f):
    f_cf4 = np.asarray(f, dtype=float)
    W=np.interp(f_cf4,cf4_pct,ion_pot)
    return W


######################################33

fCF4 = np.logspace(-3,0,1000)
pressure = [1,2,3,4,5]

plt.figure()
print(parameter_data)

cmap = "viridis"
cmap_obj = plt.get_cmap(cmap)
colors = cmap_obj(np.linspace(0.15, 0.85, len(pressure)))

for i in range(len(pressure)):
    yield_teo = (theory_yield_vis(parameter_data,degrad_data,fCF4,pressure[i])) * (1/0.015) * ion_potential(fCF4)
    plt.plot(
        fCF4 * 100,
        yield_teo,
        color=colors[i],
        label=f"{pressure[i]} bar prediction"
    )


plt.errorbar(cf4_red_E100,
             y_red_E100,
             yerr=yerr_red_E100,
             marker="o",
             linestyle="none",
             label="$\\alpha$ 100 V/cm data")

plt.xscale("log")
#plt.yscale("log")
plt.ylabel("phe/MeV")
plt.xlabel("CF4 concetration [%]")
plt.legend()
plt.savefig("plots/ArCF4_primary.pdf")

# x_completed = par_completed[par_completed["type"].str.contains("value")].drop(columns=["type"]).to_numpy()


# N0 = x_completed[0,0] 
# x_completed[0,0] = 1



# # Mezclas

# fCF4 = np.logspace(-5,0,100)

# # Mejoramos los valores:
# pob_Ar3rd = np.zeros_like(fCF4)
# pob_ArdbleStar = np.zeros_like(fCF4)
# pob_CF3 = np.zeros_like(fCF4)
# pob_CF4 = np.zeros_like(fCF4)

# for i in range(len(fCF4)): 
#     pob_Ar3rd[i] = np.interp(fCF4[i],f,pAr3rd)
#     pob_ArdbleStar[i] = np.interp(fCF4[i],f,pArdbleStar) 
#     pob_CF3[i] = np.interp(fCF4[i],f,pCF3)
#     pob_CF4[i] = np.interp(fCF4[i],f,pCF4)
    
# n = 1
# phMev =  1/0.015
# n1_completed = theory_yield_vis(x_completed[0,:],fCF4,n,pob_CF3,pob_ArdbleStar,pob_CF4,pob_Ar3rd) * phMev


# fig,ax = plt.subplots(ncols=1,figsize=(8,5))
# ax.plot(fCF4*100,n1_completed,color="blue",label="New Model Prediction (1-5 bar)")
# ax.errorbar(visible["fCF4 real"],
#              visible["1.0bar"]*phMev*ion_potential(visible["fCF4 real"]/100)/N0,
#              yerr=visible["Err 1.0bar"]*phMev*ion_potential(visible["fCF4 real"]/100)/N0,
#              marker="o",
#              linestyle="none",
#              label="X-Ray data")



# ax.legend(loc="upper left")
# ax.set_xscale("log")
# ax.set_xlabel("$f_{CF_4}$ (%)")
# ax.set_yscale("log")
# ax.set_ylabel("Yield [phe/MeV]")
# ax.set_xlim(0.07,120)
# ax.set_ylim(20,2300)



# fig.tight_layout()
# fig.savefig("yield_prediction.pdf",bbox_inches="tight")
