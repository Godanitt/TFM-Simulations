import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

plt.style.use("default")

models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../models"))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))

sys.path.append(models_dir)
sys.path.append(data_dir)

from ArCF4 import *
from read_Degrad import read_degrad
from read_experimental import read_experimental
from read_Root import export_hlevels_to_csv,read_data_per_primary_electron
from read_secondary import read_garfield_csv_folder


# ============================================================
# RUTAS
# ============================================================
folder_path = "../data/Secondary_GarfieldData/ArCF4/root"
table_path = "../data/Secondary_GarfieldData/levels/ArCF4_level_data.csv"

csv_folder = "../data/Secondary_GarfieldData/ArCF4/csv"
populations_dir = "../data/Secondary_GarfieldData/ArCF4/populations"
plots_dir = "plots"

os.makedirs(populations_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

export_hlevels_to_csv(
    folder_path,
    table_path,
    object_name="hLevels",
    argon_update=True
)

# ============================================================
# 2) LEER GANANCIAS ne y ni
#    IMPORTANTE: usar el mismo gas_concentration que luego en
#    read_garfield_csv_folder para que el merge sea consistente
# ============================================================
summary = read_data_per_primary_electron(
    folder_path,
    gas_concentration="cf4"
)


# ============================================================
# 3) CONFIGURACIÓN DE POBLACIONES
# ============================================================
config = pd.DataFrame({
    "CF4": {
        "name principal": "ION",
        "gas": "CF4",
        "energy low": 15.5,
        "energy up": 16,
        "name output": "CF4",
        "type": "ionisation"
    },
    "Ar**": {
        "name principal": "EXC",
        "gas": "Ar",
        "energy low": 11.8,
        "energy up": 100,
        "name output": "Ar_dbleStar",
        "type": "excitation"
    },
    "CF3": {
        "name principal": "NEUTRAL DISS",
        "gas": "CF4",
        "energy low": 15.6,
        "energy up": 100,
        "name output": "CF3",
        "type": "inelastic"
    },
    "Ar3rd": {
        "name principal": "IONISATION",
        "gas": "Ar",
        "energy low": 40,
        "energy up": 120,
        "name output": "Ar_3rd",
        "type": "ionisation"
    }
})


# ============================================================
# 5) POBLACIONES NORMALIZADAS POR ne
# ============================================================
garfield_norm_ne = read_garfield_csv_folder(
    folder_path=csv_folder,
    dataframe=config,
    output_dir=populations_dir,
    output_general_name=os.path.join(populations_dir, "ArCF4_secondary"),
    gas_concentration="cf4",
    gain_summary=summary,
    normalized="ni"
)


# ============================================================
# 6) CARGA DE DATOS PARA EL MODELO
# ============================================================
DATA_DIR_EXP = "../data/Experimental/ArCF4/"
DATA_DIR_GARFIELD = populations_dir
DATA_DIR_PAR = "../data/Parameters"

garfield_data = pd.read_csv(os.path.join(DATA_DIR_GARFIELD, "ArCF4_secondary.csv"))
garfield_data["concentration"] = garfield_data["concentration"] / 100.0

parameter_data = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArCF4_primary.csv"))["parameter"].to_numpy()

print("parameter_data original:")
print(parameter_data)

parameter_data[0] = 1
parameter_data[1] = 0.38
parameter_data[2] = 0.38
parameter_data[-1] = parameter_data[-1] * 0.8

print("parameter_data modificado:")
print(parameter_data)


# ============================================================
# 7) MALLA DE CONCENTRACIONES Y CAMPOS
# ============================================================
fN2 = np.logspace(-3, 0, 1000)
electric_fields = [40, 60, 70, 80, 90, 110]


# ============================================================
# 8) YIELD VISIBLE
# ============================================================
plt.figure(figsize=(6, 4))
plt.style.use("science")

cmap = "viridis"
cmap_obj = plt.get_cmap(cmap)
colors = cmap_obj(np.linspace(0.15, 0.85, len(electric_fields)))

for i, electric_field in enumerate(electric_fields):
    subset = garfield_data[garfield_data["electric_field"] == electric_field].copy()

    yield_teo = (
        theory_yield_vis(parameter_data, subset, fN2, 1)
        / 500.0
        * ion_potential(fN2)
    )

    plt.plot(
        fN2 * 100,
        yield_teo,
        color=colors[i],
        label=f"{electric_field} kV/cm prediction"
    )

plt.errorbar(
    [100, 67, 10, 5],
    [0.1, 0.3, 0.39, 0.38],
    yerr=np.array([0.1, 0.3, 0.39, 0.38]) * 0.25,
    fmt="."
)

plt.title("1 bar 100 gain secondary visible yield prediction for Ar/CF4 mixture")
plt.xscale("log")
plt.ylabel("ph/e$^-$")
plt.xlabel("CF$_4$ concentration [\%]")
plt.xlim(1, 110)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "ArCF4_vis_secondary_no_threshold.pdf"))



# ============================================================
# 9) YIELD UV
# ============================================================
plt.figure(figsize=(6, 4))
plt.style.use("science")

cmap = "viridis"
cmap_obj = plt.get_cmap(cmap)
colors = cmap_obj(np.linspace(0.15, 0.85, len(electric_fields)))

for i, electric_field in enumerate(electric_fields):
    subset = garfield_data[garfield_data["electric_field"] == electric_field].copy()

    yield_teo_uv = (
        theory_yield_uv(parameter_data, subset, fN2, 1)
        / 500.0
        * ion_potential(fN2)
    )
    yield_teo_vis = (
        theory_yield_vis(parameter_data, subset, fN2, 1)
        / 500.0
        * ion_potential(fN2)
    )

    yield_teo = yield_teo_uv / 50 + yield_teo_vis / 4.2

    plt.plot(
        fN2 * 100,
        yield_teo,
        color=colors[i],
        label=f"{electric_field} kV/cm prediction"
    )

plt.errorbar(
    [100, 67, 10, 5],
    [0.04, 0.045, 0.085, 0.068],
    yerr=np.array([0.04, 0.045, 0.085, 0.068]) * 0.25,
    fmt="."
)

plt.title("1 bar 100 gain secondary UV yield prediction for Ar/CF4 mixture")
plt.xscale("log")
plt.ylabel("ph/e$^-$")
plt.xlim(1, 110)
plt.ylim(0.01, 0.12)
plt.xlabel("CF$_4$ concentration [\%]")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "ArCF4_uv_secondary_no_threshold.pdf"))
