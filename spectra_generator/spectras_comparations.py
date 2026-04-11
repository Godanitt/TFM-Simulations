import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import dill 
import scienceplots

plt.style.use(['science', 'grid'])
models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../models'))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))

sys.path.append(models_dir)
sys.path.append(data_dir)

from ArN2 import *
from ArN2_infrarred import *
from ArCF4 import * 
from ArCF4_infrarred import * 


######################################33

DATA_DIR_DEGRAD = "../data/Primary_DegradData"
DATA_DIR_PAR = "../data/Parameters"

degrad_data_cf4 = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArCF4.csv"))
degrad_data_cf4_IR = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArCF4_IR.csv"))
degrad_data_n2 = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArN2.csv"))
degrad_data_n2_IR = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArN2_IR.csv"))


parameter_data_cf4 = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArCF4_primary.csv"))["parameter"].to_numpy()
parameter_data_cf4_IR = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArCF4_IR_primary.csv"))["parameter"].to_numpy()
parameter_data_n2= pd.read_csv(os.path.join(DATA_DIR_PAR, "ArN2_primary.csv"))["parameter"].to_numpy()
parameter_data_n2_IR = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArN2_IR_primary.csv"))["parameter"].to_numpy()


norm_cf4 = parameter_data_cf4[0].copy()
norm_n2  = parameter_data_n2[0].copy()

norm = (norm_cf4+norm_n2)/2

######################################

def gaussiana(x,mu,sigma):
    return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma)**2)

######################################33

with open("../data/Experimental/ArN2/N2_primary_data_final.pkl", "rb") as f:
        df = dill.load(f)


spectrum = df.loc[0,"mean_spectrum"]


cmap = "viridis"
cmap_obj = plt.get_cmap(cmap)
pressure = [1]
colors = cmap_obj(np.linspace(0.15, 0.85, len(pressure)))
concentrations = [0.1,1,10,100]

wavelength = np.linspace(200,800,2000)

equations_n2 = {
    "696": theory_yield_ArN2_Ir_696,
    "727": theory_yield_ArN2_Ir_727,
    "750": theory_yield_ArN2_Ir_750,
    "763": theory_yield_ArN2_Ir_763,
    "772": theory_yield_ArN2_Ir_772,
}


equations_cf4 = {
    "696": theory_yield_ArN2_Ir_696,
    "727": theory_yield_ArN2_Ir_727,
    "750": theory_yield_ArN2_Ir_750,
    "763": theory_yield_ArN2_Ir_763,
    "772": theory_yield_ArN2_Ir_772,
    "794": theory_yield_ArN2_Ir_794,
}




# =========================================================
# PRIMERA PASADA: calcular todos los espectros y el ymax global
# =========================================================
all_spectra_n2 = []
global_ymax = 0

for con in concentrations:
    spectra_con = []
    factor = (1/0.012) * W_ArN2(con/100) / norm


    for pres in pressure:
        yield_N2 = theory_yield_N2_uv(
            parameter_data_n2, degrad_data_n2, np.array([con/100]), pres
        )
   
        yield_total = 0.13 * factor * yield_N2[0] * gaussiana(wavelength, 310, 3)
        yield_total += 0.42 * factor * yield_N2[0] * gaussiana(wavelength, 335, 2.5)
        yield_total += 0.3 * factor * yield_N2[0] * gaussiana(wavelength, 355, 2.5)
        yield_total += 0.1 * factor * yield_N2[0] * gaussiana(wavelength, 378, 2.5)
        yield_total += 0.05* factor * yield_N2[0] * gaussiana(wavelength, 403, 2.5)

        for name, yield_IR in equations_n2.items():
            yield_ir = yield_IR(
                parameter_data_n2_IR, degrad_data_n2_IR, np.array([con/100]), pres
            )
            yield_total += (factor) * yield_ir[0] * gaussiana(wavelength, float(name), 2.8)

        spectra_con.append((pres, yield_total))
        global_ymax = max(global_ymax, np.max(yield_total))

    all_spectra_n2.append((con, spectra_con))

# Guardamos todo primero
all_spectra_cf4 = []
global_ymax = 0

for con in concentrations:
    spectra_con = []

    for pres in pressure:

        factor = (1/0.015) * ion_potential(con/100) / norm

        yield_vis = theory_yield_vis(
            parameter_data_cf4, degrad_data_cf4, np.array([con/100]), pres
        ) * factor

        yield_uv, yield_cf4, yield_ArDbleStar, yield_cf3_uv = theory_yield_uv(
            parameter_data_cf4, degrad_data_cf4, np.array([con/100]), pres, activate_components=True
        )


        yield_uv *= factor
        yield_cf4 *= factor
        yield_ArDbleStar *= factor
        yield_cf3_uv *= factor

        yield_vis_spec = yield_vis[0] * gaussiana(wavelength, 630, 40)

        yield_cf4_230 = (0.8/1.85) * yield_cf4[0] * gaussiana(wavelength, 230, 20)
        yield_cf4_290 = (0.95/1.85) * yield_cf4[0] * gaussiana(wavelength, 290, 20)
        yield_cf4_364 = (0.10/1.85) * yield_cf4[0] * gaussiana(wavelength, 364, 40)
        yield_cf4_spec = yield_cf4_230 + yield_cf4_290 + yield_cf4_364

        yield_arDbleStar_spec = yield_ArDbleStar[0] * gaussiana(wavelength, 220, 60)
        yield_CF3_spec = yield_cf3_uv[0] * gaussiana(wavelength, 245, 60)


        yield_total = (
            yield_vis_spec
            + yield_cf4_spec
            + yield_arDbleStar_spec
            + yield_CF3_spec
        )


        for name, yield_IR in equations_cf4.items():
            yield_ir = yield_IR(
                parameter_data_cf4_IR, degrad_data_cf4_IR, np.array([con/100]), pres
            )
            yield_total += (factor) * yield_ir[0] * gaussiana(wavelength, float(name), 2.7)


        spectra_con.append((pres, yield_total))
        global_ymax = max(global_ymax, np.max(yield_total))

    all_spectra_cf4.append((con, spectra_con))



# =========================================================
# SEGUNDA PASADA: dibujar una única figura con 4 paneles
# =========================================================
fig, axs = plt.subplots(2, 2, figsize=(9, 6), sharex=True, sharey=True)
axs = axs.ravel()

for ax, (con_cf4, spectra_cf4_con), (con_n2,spectra_n2_con) in zip(axs, all_spectra_cf4,all_spectra_n2):

    for k, (pres, yield_total_cf4,) in enumerate(spectra_cf4_con):
        ax.plot(
            wavelength,
            yield_total_cf4,
            color="blue",
            label=f"CF$_4$ {pres:.1f} bar"
        )

    for k, (pres, yield_total_n2,) in enumerate(spectra_n2_con):
        ax.plot(
            wavelength,
            yield_total_n2,
            color="red",
            label=f"N$_2$ {pres:.1f} bar"
        )

    ax.set_title(f"{con_cf4:.1f} $\%$ Aditivo")
    ax.set_xlabel(r"$\lambda$ [nm]")
    ax.set_ylabel("ph/MeV/nm")
    ax.grid(True, which='major', alpha=0.3)
    ax.grid(True, which='minor', alpha=0.08)
    ax.set_ylim(0, 1.5 * global_ymax)
    ax.legend(ncol=2,loc="upper right")

fig.suptitle(r"Primary Ar-N$_2$ \& Ar-CF$_4$ Spectra Prediction", fontsize=14)
fig.tight_layout()
fig.savefig("Comparation.pdf", dpi=300, bbox_inches="tight")
plt.show()