import pandas as pd 
import os 
import matplotlib.pyplot as plt

from ScintillationClass import Scintillation
from ArCF4_model import *


DATA_DIR    = os.path.join("pickle_data")

yield_uv    = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis   = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

PCF3        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF3.pkl"))
PCF4        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
PArDbleStar = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_dbleStar.pkl"))
PAr3rd      = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_3rd.pkl"))

fCF4_real   = yield_uv["fCF4 real"].to_numpy()
fCF4        = PCF3["fCF4"]
PCF3        = PCF3[["CF3 >11.5","Err CF3 >11.5"]]
PCF4        = PCF4[["CF4 all","Err CF4 all"]]
PArDbleStar = PArDbleStar[["Ar** all","Err Ar** all"]]
PAr3rd      = PAr3rd[["Ar3rd all","Err Ar3rd all"]]

# Diccionarios iniciales
yields = {
    "fCF4": yield_uv["fCF4 real"].to_numpy(),
    "sCF4": yield_uv["Err fCF4 real"].to_numpy(),
    "vis": yield_vis.drop(columns=["fCF4", "fCF4 real", "Err fCF4 real"]),  
    "uv": yield_uv.drop(columns=["fCF4", "fCF4 real", "Err fCF4 real"])  
}

poblation_degrad_data = {
    "fCF4": fCF4.to_numpy(),
    "CF3": PCF3,
    "Ar dbleStar": PArDbleStar,
    "CF4": PCF4,
    "Ar 3rd": PAr3rd,
}

scintillation_theory_models = {
    "CF3 dir": Pgamma_CF3dir,
    "CF3 Ar dbleStar": Pgamma_CF3ArDbleStar,
    "CF4 dir": Pgamma_CF4dir,
    "CF4 Ar 3rd": Pgamma_CF4Ar3rd,
    "Ar 3rd": Pgamma_Ar3rd
}

# Inicializamos el contenedor con los datos
ArCF4 = Scintillation(
    yields = yields,
    poblation_degrad=poblation_degrad_data,
    scintillation_models=scintillation_theory_models
)

#################################################################
#################################################################
#################################################################

ArCF4.plotPoblationInterpolation("CF3",savefig="plot.pdf")

# Hacemos el ajuste: 


scintillation_vis = {
    "CF3 dir": ["Probabilidad"],
    "CF3 Ar dbleStar": ["Probabilidad"]
}

scintillation_uv = {
    "CF4 dir": ["Relajacion","Centelleo"],
    "CF4 Ar 3rd": [""],
    "Ar 3rd": [""],
}

scintillation= {
    "vis": scintillation_vis,
    "uv": scintillation_uv,
}

# 1) Construyo teoría
ArCF4.build_theory_functions(scintillation)

# 2) Ajusto visible
x0_vis = [0.1,0.2]
res = ArCF4.fit_parameters_chooseNorma("vis", x0_vis, n0=1.0, idx_ref=-1)

# 3) Configuración de gráficos
ArCF4.choosePlotNormalization("vis", mode="N0", value=1.0)
ArCF4.EnableExperimentalData("vis", 1.0)
ArCF4.EnableTeoCurve("vis", 1.0)
ArCF4.EnableTeoCurve("vis", 2.0)

# 4) Plot
ArCF4.plot_teoCurve("vis")
