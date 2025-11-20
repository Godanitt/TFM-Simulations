

import numpy as np
import pandas as pd
import dill
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import sys, os

PARENT_DIR = os.path.abspath(os.path.join(".."))    # proyecto/
sys.path.append(PARENT_DIR)
from Amoedo_Model_DivisionFit import Pgamma_UV_Cociente, Pgamma_vis_Cociente


"""
Script que leyendo:
    - ajustes lineales de poblaciones de Degrad (con esto obtendremos las poblaciones)
    - yields experimentales
    - parámetros ajustados (alpha, kcool, kdis)
vuelve a generar los mismos plots del ajuste.
"""


# ================================
# RUTAS A LOS PICKLES
# ================================
DATA_DIR = "pickle_data"
FIT_DIR  = "pickle_data"


# === Leer poblaciones ===
with open(os.path.join(DATA_DIR, "linealFun_poblations_Ar_3rd.pkl"), "rb") as f:
    lineal_pAr_3rd = dill.load(f)

with open(os.path.join(DATA_DIR, "linealFun_poblations_Ar_dbleStar.pkl"), "rb") as f:
    lineal_pAr_dbleStar = dill.load(f)

with open(os.path.join(DATA_DIR, "linealFun_poblations_CF4.pkl"), "rb") as f:
    lineal_CF4 = dill.load(f)

with open(os.path.join(DATA_DIR, "linealFun_poblations_CF3.pkl"), "rb") as f:
    lineal_CF3 = dill.load(f)


# === Leer yields experimentales ===
DATA_DIR = os.path.join("..", "pickle_data")
yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))


# === Leer parámetros ajustados ===
df_alpha = pd.read_pickle(os.path.join(FIT_DIR, "alpha_results.pkl"))
df_beta = pd.read_pickle(os.path.join(FIT_DIR, "beta_results.pkl"))
df_kcool = pd.read_pickle(os.path.join(FIT_DIR, "kcool_results.pkl"))
df_kdis  = pd.read_pickle(os.path.join(FIT_DIR, "kdis_results.pkl"))


# ============================
# Variables útiles
# ============================
name_CF4         = lineal_CF4.columns.to_numpy()[::2]
name_CF3         = lineal_CF3.columns.to_numpy()[::2]
name_Ar_dbleStar = lineal_pAr_dbleStar.columns.to_numpy()[::2]
name_Ar_3rd      = lineal_pAr_3rd.columns.to_numpy()[::2]

fCF4_real = yield_uv["fCF4 real"].to_numpy()
index_ref = -1   # 100% CF4


# ============================
# Colormaps para las figuras
# ============================
cmap_exp = cm.get_cmap("viridis")
cmap_fit = cm.get_cmap("hot")

i_exp = 0
i_fit = 1


# ===========================================================
#   PLOTEAR VISIBLE USANDO α GUARDADO
# ===========================================================

plt.figure()

for nameCF3 in name_CF3:
    for nameAr in name_Ar_dbleStar:

        colname = f"{nameCF3}_{nameAr}"
        if colname not in df_alpha.columns:
            continue

        alpha = df_alpha.loc[0, colname]
        beta = df_beta.loc[0, colname]
        if pd.isna(alpha):
            continue
        if pd.isna(beta):
            continue
        
        # === Recuperar funciones originales ===
        f = lineal_CF3[nameCF3].to_numpy()[0]
        g = lineal_pAr_dbleStar[nameAr].to_numpy()[0]

        PCF3         = f(fCF4_real / 100)
        PAr_dbleStar = g(1 - fCF4_real / 100)

        values0_vis   = (fCF4_real[index_ref]/100, PCF3[index_ref], PAr_dbleStar[index_ref],1.0)
        values0_yield = yield_vis["1.0bar"].to_numpy()[index_ref]

        # === Curva suave ===
        fCF4_range = np.linspace(min(fCF4_real), max(fCF4_real), 1000)
        PCF3_range = f(fCF4_range / 100)
        PAr_range  = g(1 - fCF4_range / 100)

        model = Pgamma_vis_Cociente(
            fCF4_range / 100,
            PCF3_range,
            PAr_range,
            1.0,
            values0_vis,
            alpha,
            beta
        )

        color_fit = cmap_fit(i_fit / 10)
        i_fit += 1

        plt.plot(
            fCF4_range,
            model,
            color=color_fit,
            label=f"{colname}"
        )

        model = Pgamma_vis_Cociente(
            fCF4_range / 100,
            PCF3_range,
            PAr_range,
            1.0,
            values0_vis,
            1,
            1
        )

        color_fit = cmap_fit(i_fit / 10)
        i_fit += 1

        plt.plot(
            fCF4_range,
            model,
            color=color_fit,
            label=f"{colname} sin ajuste"
        )

# === Añadir datos experimentales VIS ===
#bars = ["1.0bar", "2.0bar", "3.0bar", "4.0bar", "5.0bar"]
bars = ["1.0bar","4.0bar"]

for b in bars:
    """
    yerr = np.sqrt(
            (yield_vis["Err " + b] / yield_vis["1.0bar"])**2 +
            (yield_vis["Err 1.0bar"] * yield_vis[b].to_numpy() /
            yield_vis["1.0bar"].to_numpy()[index_ref]**2)**2
        )
    """
    yerr = np.sqrt(
            (yield_vis["Err " + b] / yield_vis[b].to_numpy()[index_ref])**2 +
            (yield_vis["Err "+  b] * yield_vis[b].to_numpy() /
           (yield_vis[b].to_numpy()[index_ref])**2)**2
        )
    
    color_exp = cmap_exp(i_exp / len(bars))
    i_exp += 1

    if yield_vis[b].to_numpy()[index_ref] == 0.0:
        plt.errorbar(
            fCF4_real[:-1],
            yield_vis[b].to_numpy()[:-1] / yield_vis["1.0bar"].to_numpy()[index_ref],
            yerr=yerr[:-1],
            fmt="o",
            capsize=5,
            label=b,
            color=color_exp
        )
    else:
        plt.errorbar(
            fCF4_real,
            yield_vis[b].to_numpy() / yield_vis[b].to_numpy()[index_ref],
            yerr=yerr,
            fmt="o",
            capsize=5,
            label=b,
            color=color_exp
        )


plt.xscale("log")
plt.yscale("log")
plt.ylim(0.01, 5)
plt.legend(fontsize=7, markerscale=0.8, framealpha=0.9)
plt.title("Visible — Ajuste con α guardado")
plt.savefig("output/VIS-Normalizado_100%CF4_1bar.pdf")



# ===========================================================
#   PLOTEAR ULTRAVIOLETA USANDO kcool, kdis GUARDADOS
# ===========================================================

plt.figure()

i_fit = 0
i_exp = 0

for nameCF4 in name_CF4:
    for nameAr in name_Ar_3rd:

        colname = f"{nameCF4}_{nameAr}"
        if colname not in df_kcool.columns:
            continue

        kcool = df_kcool.loc[0, colname]
        kdis  = df_kdis.loc[0, colname]
        if pd.isna(kcool) or pd.isna(kdis):
            continue

        #kcool,kdis=0.5,1
        # === Recuperar funciones originales ===
        f = lineal_CF4[nameCF4].to_numpy()[0]
        g = lineal_pAr_3rd[nameAr].to_numpy()[0]

        PCF4    = f(fCF4_real / 100)
        PAr_3rd = g(1 - fCF4_real / 100)

        values0_uv    = (fCF4_real[index_ref]/100, PCF4[index_ref], PAr_3rd[index_ref], 1)
        values0_yield = yield_uv["1.0bar"].to_numpy()[index_ref]

        fCF4_range = np.linspace(min(fCF4_real), max(fCF4_real), 1000)
        PCF4_range = f(fCF4_range/100)
        PAr3_range = g(1 - fCF4_range/100)

        model_1bar = Pgamma_UV_Cociente(
            fCF4_range/100, PCF4_range, PAr3_range,
            1, values0_uv, kcool, kdis
        )

        model_1bar_nofit = Pgamma_UV_Cociente(
            fCF4_range/100, PCF4_range, PAr3_range,
            1, values0_uv, 0, 0
        )

        model_4bar = Pgamma_UV_Cociente(
            fCF4_range/100, PCF4_range, PAr3_range,
            4, values0_uv, kcool, kdis
        )

        i_fit += 1
        color_fit = cmap_fit(i_fit / 6)
        plt.plot(fCF4_range, model_1bar, label=f"1bar {colname}", color=color_fit)

        i_fit += 1
        color_fit = cmap_fit(i_fit / 6)
        plt.plot(fCF4_range, model_1bar_nofit, label=f"1bar {colname} - sin ajuste", color=color_fit)

        i_fit += 1
        color_fit = cmap_fit(i_fit / 6)
        plt.plot(fCF4_range, model_4bar, label=f"4bar {colname}", color=color_fit)


# === Añadir datos experimentales UV ===
for b in bars:
    """
    yerr = np.sqrt(
        (yield_uv["Err " + b] / yield_uv["1.0bar"])**2 +
        (yield_uv["Err 1.0bar"] * yield_uv[b].to_numpy() /
         yield_uv["1.0bar"].to_numpy()[index_ref]**2)**2
    )
    """
    b0 = "1.0bar"
    b0 = b
    yerr = np.sqrt(
        (yield_uv["Err " + b].to_numpy() / yield_uv[b0].to_numpy()[index_ref])**2 +
        (yield_uv["Err " + b0].to_numpy() * yield_uv[b].to_numpy() /
         yield_uv[b0].to_numpy()[index_ref]**2)**2
    )

    color_exp = cmap_exp(i_exp / len(bars))
    i_exp += 1

    if yield_uv[b].to_numpy()[index_ref] == 0.0:
        plt.errorbar(
            fCF4_real[:-1],
            yield_uv[b].to_numpy()[:-1] / yield_uv[b0].to_numpy()[index_ref],
            yerr=yerr[:-1],
            fmt="o",
            capsize=5,
            label=b,
            color=color_exp
        )
    else:
        plt.errorbar(
            fCF4_real,
            yield_uv[b].to_numpy() / yield_uv[b0].to_numpy()[index_ref],
            yerr=yerr,
            fmt="o",
            capsize=5,
            label=b,
            color=color_exp
        )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("fCF4 real [%]")
plt.ylabel("Cociente UV normalizado")
plt.title("UV — Ajuste con kcool, kdis guardados")
plt.legend(fontsize=7, markerscale=0.8, framealpha=0.9)
plt.savefig("output/UV-Normalizado_100%CF4_1bar.pdf", dpi=300)

print("Listo. Gráficas generadas en output/")
