import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
from ArCF4_Completed import * 

visible = pd.read_csv("yield_vis_cal.csv")
ultravioleta = pd.read_csv("yield_uv_cal.csv")

par_completed = pd.read_csv("Parametros_Globales.csv")

f = pd.read_csv("poblations_Ar_3rd.csv")["fCF4"].to_numpy()
pAr3rd = pd.read_csv("poblations_Ar_3rd.csv")["Ar3rd all"].to_numpy()
pArdbleStar = pd.read_csv("poblations_Ar_dbleStar.csv")["Ar** all"].to_numpy()
pCF3 = pd.read_csv("poblations_CF3.csv")["CF3 >11.4"].to_numpy()
pCF4 = pd.read_csv("poblations_CF4.csv")["CF4 all"].to_numpy()


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

# Cogemos los parámetros

x_completed = par_completed[par_completed["type"].str.contains("value")].drop(columns=["type"]).to_numpy()


N0 = x_completed[0,0] 
x_completed[0,0] = 1



# Mezclas

fCF4 = np.logspace(-5,0,100)

# Mejoramos los valores:
pob_Ar3rd = np.zeros_like(fCF4)
pob_ArdbleStar = np.zeros_like(fCF4)
pob_CF3 = np.zeros_like(fCF4)
pob_CF4 = np.zeros_like(fCF4)

for i in range(len(fCF4)): 
    pob_Ar3rd[i] = np.interp(fCF4[i],f,pAr3rd)
    pob_ArdbleStar[i] = np.interp(fCF4[i],f,pArdbleStar) 
    pob_CF3[i] = np.interp(fCF4[i],f,pCF3)
    pob_CF4[i] = np.interp(fCF4[i],f,pCF4)
    
n = 1
phMev =  1/0.015
n1_completed = theory_yield_vis(x_completed[0,:],fCF4,n,pob_CF3,pob_ArdbleStar,pob_CF4,pob_Ar3rd) * phMev


fig,ax = plt.subplots(ncols=1,figsize=(8,5))
ax.plot(fCF4*100,n1_completed,color="blue",label="New Model Prediction (1-5 bar)")
ax.errorbar(visible["fCF4 real"],
             visible["1.0bar"]*phMev*ion_potential(visible["fCF4 real"]/100)/N0,
             yerr=visible["Err 1.0bar"]*phMev*ion_potential(visible["fCF4 real"]/100)/N0,
             marker="o",
             linestyle="none",
             label="X-Ray data")
ax.errorbar(cf4_red_E100,
             y_red_E100,
             yerr=yerr_red_E100,
             marker="o",
             linestyle="none",
             label="$\\alpha$ 100 V/cm data")


ax.legend(loc="upper left")
ax.set_xscale("log")
ax.set_xlabel("$f_{CF_4}$ (%)")
#ax.set_yscale("log")
ax.set_ylabel("Yield [phe/MeV]")
ax.set_xlim(0.07,120)
ax.set_ylim(20,2300)



fig.tight_layout()
fig.savefig("yield_prediction.pdf",bbox_inches="tight")
