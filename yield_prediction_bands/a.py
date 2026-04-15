import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(["science", "grid"])

# =========================================================
# PATHS
# =========================================================
BASE_DIR = os.path.dirname(__file__)

models_dir = os.path.abspath(os.path.join(BASE_DIR, "../models"))
data_dir   = os.path.abspath(os.path.join(BASE_DIR, "../data"))
fit_dir    = os.path.abspath(os.path.join(BASE_DIR, "../primary_fits"))

sys.path.append(models_dir)
sys.path.append(data_dir)
sys.path.append(fit_dir)

from ArCF4 import *
from read_experimental import read_experimental
from fiting import fitParameters

# =========================================================
# CONFIG
# =========================================================
archivo_entrada = "../data/Experimental/ArCF4/CF4_primary_data_final.pkl"
yields = ["vis", "UV"]
presiones = [1, 2, 2.5, 3, 4, 5]
output_dir = "../data/Experimental/ArCF4/"

pressure_plot = 1.0
SCALE_TO_PHOTONS = True   # False -> normalized yield ; True -> ph/MeV

DATA_DIR_DEGRAD = "../data/Primary_DegradData"
degrad_data = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArCF4.csv"))

cmap = plt.get_cmap("viridis")
colors = cmap(np.linspace(0.12, 0.85, 6))

# =========================================================
# PARÁMETROS DEL AJUSTE
# =========================================================
x0 = np.array([
    0.14,
    0.10, 0.99, 3, 0.037 * 3,
    1.0, 0.065, 0.48, 50.10, 0.37,
    0.00001
], dtype=float)

lower = [
    0.0,
    0.0, 0.0, 0.0, 0.0,
    0.00, 0.065, 0.0, 50, 0.0,
    0.0
]

upper = [
    0.99,
    1.0, 1.0, 10000.0, 10000.0,
    10000.0, 0.066, 1.0, 50.2, 1.0,
    0.0001
]

bounds = (lower, upper)

fixed_idx = [6, 8, 10]
fixed_values = [0.065, 50.05, 0.0001]
fixed_error = 0.01

equations = {
    "vis": theory_yield_vis,
    "uv": theory_yield_uv
}

y_cols = ["1.0bar", "2.0bar", "2.5bar", "3.0bar", "4.0bar", "5.0bar"]
err_cols = ["Err 1.0bar", "Err 2.0bar", "Err 2.5bar", "Err 3.0bar", "Err 4.0bar", "Err 5.0bar"]

# =========================================================
# HELPERS
# =========================================================
def load_nominal_data():
    read_experimental(
        archivo_entrada,
        yields,
        presiones,
        output_dir,
        uncertainty_mode="all"
    )

    data_dir_exp = "../data/Experimental/ArCF4/"
    df_uv = pd.read_csv(os.path.join(data_dir_exp, "UV.csv"))
    df_vis = pd.read_csv(os.path.join(data_dir_exp, "vis.csv"))

    # Igual que en tu script
    df_uv.loc[0, "fCF4"] = 0.001
    df_vis = df_vis.fillna(0)

    return {
        "uv": df_uv,
        "vis": df_vis
    }


def load_systematic_triplet():
    read_experimental(
        archivo_entrada,
        yields,
        presiones,
        output_dir,
        uncertainty_mode="systematic"
    )

    data_dir_exp = "../data/Experimental/ArCF4/"
    df_uv = pd.read_csv(os.path.join(data_dir_exp, "UV.csv"))
    df_vis = pd.read_csv(os.path.join(data_dir_exp, "vis.csv"))

    df_uv.loc[0, "fCF4"] = 0.001
    df_vis = df_vis.fillna(0)

    # central sistemático
    df_uv_mid = df_uv.copy(deep=True)
    df_vis_mid = df_vis.copy(deep=True)

    # low/up
    err_uv = df_uv[err_cols].to_numpy(dtype=float)
    err_vis = df_vis[err_cols].to_numpy(dtype=float)

    df_uv_low = df_uv.copy(deep=True)
    df_uv_up = df_uv.copy(deep=True)
    df_vis_low = df_vis.copy(deep=True)
    df_vis_up = df_vis.copy(deep=True)

    df_uv_low.loc[:, y_cols] = df_uv_low.loc[:, y_cols].to_numpy(dtype=float) - err_uv
    df_uv_up.loc[:, y_cols]  = df_uv_up.loc[:, y_cols].to_numpy(dtype=float) + err_uv

    df_vis_low.loc[:, y_cols] = df_vis_low.loc[:, y_cols].to_numpy(dtype=float) - err_vis
    df_vis_up.loc[:, y_cols]  = df_vis_up.loc[:, y_cols].to_numpy(dtype=float) + err_vis

    mid = {"uv": df_uv_mid, "vis": df_vis_mid}
    low = {"uv": df_uv_low, "vis": df_vis_low}
    up  = {"uv": df_uv_up,  "vis": df_vis_up}

    return mid, low, up


def fit_global_model(experimental_data):
    popt = fitParameters(
        equations,
        experimental_data,
        degrad_data,
        x0=x0,
        bounds=bounds,
        fixed_idx=fixed_idx,
        fixed_values=fixed_values,
        fixed_error=fixed_error
    )
    return popt.x.copy(), popt


def get_scale_factor_curve(fcf4_fraction, norm):
    return ion_potential(fcf4_fraction) / (0.015 * norm)


def get_scale_factor_points(fcf4_percent, norm):
    return ion_potential(np.asarray(fcf4_percent, dtype=float) / 100.0) / (0.015 * norm)


# =========================================================
# CARGA DE DATOS Y AJUSTES
# =========================================================
nominal_data = load_nominal_data()
sys_mid_data, sys_low_data, sys_up_data = load_systematic_triplet()

# Ajuste nominal "mid"
par_mid, popt_mid = fit_global_model(nominal_data)

# Ajustes low/up con los datasets desplazados
par_low, popt_low = fit_global_model(sys_low_data)
par_up,  popt_up  = fit_global_model(sys_up_data)

print("=" * 60)
print("Parámetros MID:\n", par_mid)
print("Parámetros LOW:\n", par_low)
print("Parámetros UP:\n", par_up)
print("=" * 60)

# =========================================================
# CURVAS UV 1 bar
# =========================================================
fCF4 = np.logspace(-5, 0, 400)

y_mid = theory_yield_uv(par_mid, degrad_data, fCF4, pressure_plot)
y_low = theory_yield_uv(par_low, degrad_data, fCF4, pressure_plot)
y_up  = theory_yield_uv(par_up,  degrad_data, fCF4, pressure_plot)

# =========================================================
# DATOS EXPERIMENTALES UV 1 bar
# =========================================================
df_mid_nom = nominal_data["uv"]      # datos usados para el ajuste nominal
df_low_sys = sys_low_data["uv"]      # datos usados para low
df_up_sys  = sys_up_data["uv"]       # datos usados para up

x_mid = df_mid_nom["fCF4"].to_numpy(dtype=float)
y_mid_exp = df_mid_nom["1.0bar"].to_numpy(dtype=float)
e_mid_exp = df_mid_nom["Err 1.0bar"].to_numpy(dtype=float)

x_low = df_low_sys["fCF4"].to_numpy(dtype=float)
y_low_exp = df_low_sys["1.0bar"].to_numpy(dtype=float)
e_low_exp = df_low_sys["Err 1.0bar"].to_numpy(dtype=float)

x_up = df_up_sys["fCF4"].to_numpy(dtype=float)
y_up_exp = df_up_sys["1.0bar"].to_numpy(dtype=float)
e_up_exp = df_up_sys["Err 1.0bar"].to_numpy(dtype=float)

# =========================================================
# ESCALADO OPCIONAL A ph/MeV
# =========================================================
if SCALE_TO_PHOTONS:
    # Igual que en tu script: se usa la normalización nominal
    norm = par_mid[0]

    factor_curve = get_scale_factor_curve(fCF4, norm)
    y_mid = y_mid * factor_curve
    y_low = y_low * factor_curve
    y_up  = y_up  * factor_curve

    factor_mid_pts = get_scale_factor_points(x_mid, norm)
    factor_low_pts = get_scale_factor_points(x_low, norm)
    factor_up_pts  = get_scale_factor_points(x_up, norm)

    y_mid_exp = y_mid_exp * factor_mid_pts
    e_mid_exp = e_mid_exp * factor_mid_pts

    y_low_exp = y_low_exp * factor_low_pts
    e_low_exp = e_low_exp * factor_low_pts

    y_up_exp = y_up_exp * factor_up_pts
    e_up_exp = e_up_exp * factor_up_pts

    ylabel = "ph/MeV"
    outname = "plots/ArCF4_UV_mid_low_up_1bar_phMeV.pdf"
    title = "Ar-CF$_4$ UV 1 bar: mid / low / up"
else:
    ylabel = "Normalized yield"
    outname = "plots/ArCF4_UV_mid_low_up_1bar_normalized.pdf"
    title = "Ar-CF$_4$ UV 1 bar: mid / low / up"

# =========================================================
# PLOT
# =========================================================
os.makedirs("plots", exist_ok=True)

plt.figure(figsize=(6.4, 4.4))

# Curvas
plt.plot(fCF4 * 100, y_mid, lw=2.2, color=colors[2], label="mid fit")
plt.plot(fCF4 * 100, y_low, lw=2.0, ls="--", color=colors[0], label="low fit")
plt.plot(fCF4 * 100, y_up,  lw=2.0, ls="-.", color=colors[4], label="up fit")

# Datos usados en cada ajuste
plt.errorbar(
    x_mid, y_mid_exp, yerr=e_mid_exp,
    marker="o", linestyle="none", ms=4,
    color=colors[2], ecolor=colors[2],
    elinewidth=1, capsize=2,
    label="mid data"
)

plt.errorbar(
    x_low, y_low_exp, yerr=e_low_exp,
    marker="s", linestyle="none", ms=4,
    color=colors[0], ecolor=colors[0],
    elinewidth=1, capsize=2,
    label="low data"
)

plt.errorbar(
    x_up, y_up_exp, yerr=e_up_exp,
    marker="^", linestyle="none", ms=4,
    color=colors[4], ecolor=colors[4],
    elinewidth=1, capsize=2,
    label="up data"
)

plt.xscale("log")
if not SCALE_TO_PHOTONS:
    plt.yscale("log")

plt.xlabel("CF$_4$ concentration [%]")
plt.ylabel(ylabel)
plt.title(title)
plt.grid(True, which="major", alpha=0.3)
plt.grid(True, which="minor", alpha=0.08)
plt.legend()
plt.tight_layout()
plt.yscale("log")
plt.savefig(outname, dpi=300, bbox_inches="tight")
plt.show()