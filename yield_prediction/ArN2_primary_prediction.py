import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from matplotlib import colors as mcolors

models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../models'))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))

sys.path.append(models_dir)
sys.path.append(data_dir)

from ArN2 import *


#####################################################
###### Traemos los datos anteriormente generados 

DATA_DIR_EXP = "../data/Experimental/ArN2/"
DATA_DIR_DEGRAD = "../data/Primary_DegradData"
DATA_DIR_PAR = "../data/Parameters"

yield_N2_uv  = pd.read_csv(os.path.join(DATA_DIR_EXP, "yield_N2.csv"))

degrad_data = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArN2.csv"))

parameter_data = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArN2_primary.csv"))["parameter"].to_numpy()
print(parameter_data)
parameter_data[0] = 1
print(parameter_data)



######################################33

x_data = [100]
y_data = [141.1] # A mi me da 136 ahora mismo!!
sy_data = [2.1]

fN2 = np.logspace(-3,0,1000)
pressure = [1,2,3,4,5]

plt.figure()

cmap = "viridis"
cmap_obj = plt.get_cmap(cmap)
colors = cmap_obj(np.linspace(0.15, 0.85, len(pressure)))

for i in range(len(pressure)):
    yield_teo = theory_yield_N2_uv(parameter_data, degrad_data, fN2, pressure[i]) / 0.012 * W_ArN2(fN2)
    plt.plot(
        fN2 * 100,
        yield_teo,
        color=colors[i],
        label=f"{pressure[i]} bar prediction"
    )

plt.errorbar(x_data,
             y_data,
             yerr = sy_data,
             marker="o",
             linestyle="none",
             label="1 bar")

plt.xscale("log")
#plt.yscale("log")
plt.ylabel("ph/e$^-$")
plt.xlabel("N$_2$ concetration [\%]")
plt.legend()
plt.savefig("plots/ArN2_primary.pdf")