import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)
from ScintillationClass import Scintillation
from ArCF4_Completed import *
from ArCF4_PabloModel import *


#######################################################################
# ======================= 1) LECTURA DE DATOS =========================
#######################################################################

DATA_DIR = "../pickle_data"

yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

PCF3        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF3.pkl"))
PCF4        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
PArDbleStar = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_dbleStar.pkl"))
PAr3rd      = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_3rd.pkl"))

# ------------------------------------------------------------
# Preprocesado
# ------------------------------------------------------------
fCF4_real   = yield_uv["fCF4 real"].to_numpy()
fCF4        = PCF3["fCF4"]

PCF3        = PCF3[["CF3 >11.5","Err CF3 >11.5"]]
PCF4        = PCF4[["CF4 all","Err CF4 all"]]
PArDbleStar = PArDbleStar[["Ar** all","Err Ar** all"]]
PAr3rd      = PAr3rd[["Ar3rd all","Err Ar3rd all"]]

# Diccionario Yields
yields = {
    "fCF4": fCF4_real,
    "sCF4": yield_uv["Err fCF4 real"].to_numpy(),
    "vis": yield_vis.drop(columns=["fCF4", "fCF4 real", "Err fCF4 real"]),
    "uv":  yield_uv.drop(columns=["fCF4", "fCF4 real", "Err fCF4 real"]),
}


# Diccionario poblaciones de degradación
poblation_degrad_data = {
    "fCF4": fCF4.to_numpy(),
    "CF3": PCF3,
    "Ar dbleStar": PArDbleStar,
    "CF4": PCF4,
    "Ar 3rd": PAr3rd,
}

# Diccionario modelos físicos
scintillation_theory_models = {
    "CF3 dir":          Pgamma_CF3dir,
    "CF3 Ar dbleStar":  Pgamma_CF3ArDbleStar,
    "CF4 dir":          Pgamma_CF4dir,
    "CF4 Ar 3rd":       Pgamma_CF4Ar3rd,
    "Ar 3rd":           Pgamma_Ar3rd
}

#######################################################################
# =========== 2) CONSTRUCCIÓN DEL OBJETO PRINCIPAL ====================
#######################################################################

ArCF4 = Scintillation(
    yields=yields,
    poblation_degrad=poblation_degrad_data,
    scintillation_models=scintillation_theory_models
)

# Se puede comentar
ArCF4.plotPoblationInterpolation("CF3", savefig="InterpolacionPoblationCF3.pdf")


#######################################################################
# ======================== 3) DEFINO TEORÍA ===========================
#######################################################################

scintillation = {
    "vis": theory_yield_vis,
    "uv": theory_yield_uv,
}

ArCF4.buildYieldFunctionsFromRaw(scintillation)

#######################################################################
# ======================== 4) AJUSTAMOS ===========================
#######################################################################


x0 = np.array([0.01, 0.16, 0.29, 1/30,    # parámetros VIS
               1, 1, 0.99, 49, 4.0])  # parámetros UV

lower = [0.0, 0.0, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0]
upper = [1.0, 1.0, 1.0, 10.0, 10.0, 10.0, 1.0, 100.0, 10.0]

#lower = [    0.001,    0.15,    0.15,    1/31,    0.9,    0.9,    0.99,    48,    3.9,]
#upper = [    0.1,    0.2,    0.3,    1/29,    1.1,   1.1,    1.0,    49.1,    4.1,]




bounds=(lower, upper)

popt = ArCF4.fitParametersGlobalRaw_residuals(bands=["vis", "uv"], x0=x0, bounds=bounds)

print("Parámetros globales:", popt)

#######################################################################
# ======================== 5) GRAFICAMOS ===========================
#######################################################################

ArCF4.choosePlotNormalization("vis", mode="handle_global")
ArCF4.choosePlotNormalization("uv", mode="handle_global")

ArCF4.enableExperimentalData("vis", 1.0)
ArCF4.enableTeoCurve("vis", 1.0)

ArCF4.enableExperimentalData("uv", 1.0)
ArCF4.enableTeoCurve("uv", 1.0)


ArCF4.enableExperimentalData("vis", 3.0)
ArCF4.enableTeoCurve("vis", 3.0)

ArCF4.enableExperimentalData("uv", 3.0)
ArCF4.enableTeoCurve("uv", 3.0)


ArCF4.enableExperimentalData("vis", 5.0)
ArCF4.enableTeoCurve("vis", 5.0)

ArCF4.enableExperimentalData("uv", 5.0)
ArCF4.enableTeoCurve("uv", 5.0)

ArCF4.plotTeoCurve("vis", savefig="Ajuste_VIS_GlobalNorm.pdf")
ArCF4.plotTeoCurve("uv",  savefig="Ajuste_UV_GlobalNorm.pdf")
