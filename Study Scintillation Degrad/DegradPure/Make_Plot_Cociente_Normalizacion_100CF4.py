import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys

PARENT_DIR = os.path.abspath(os.path.join(".."))    # proyecto/
sys.path.append(PARENT_DIR)
from Amoedo_Model_DivisionFit import (
    Pgamma_UV_Cociente,
    Pgamma_vis_Cociente,
)

# ===========================================================
#  LECTURA DE DATOS
# ===========================================================

# === Poblaciones Degrad ===
DATA_DIR = os.path.join("..", "pickle_data")
pCF4        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
pCF3        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF3.pkl"))
pArDbleStar = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_dbleStar.pkl"))
pAr3rd      = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_3rd.pkl"))

# === Yields experimentales ===
yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

# === Parámetros ajustados ===
FIT_DIR  = "pickle_data"
df_alpha = pd.read_pickle(os.path.join(FIT_DIR, "alpha_results.pkl"))
df_kcool = pd.read_pickle(os.path.join(FIT_DIR, "kcool_results.pkl"))
df_kdis  = pd.read_pickle(os.path.join(FIT_DIR, "kdis_results.pkl"))

# ===========================================================
#  VARIABLES ÚTILES
# ===========================================================

# malla de Degrad (0–1, fracción CF4)
fCF4_grid = pCF4["fCF4"].to_numpy()

name_CF4         = pCF4.columns.to_numpy()[1::2]
name_CF3         = pCF3.columns.to_numpy()[1::2]
name_Ar_dbleStar = pArDbleStar.columns.to_numpy()[1::2]
name_Ar_3rd      = pAr3rd.columns.to_numpy()[1::2]

fCF4_real = yield_uv["fCF4 real"].to_numpy()
index_ref = -1   # 100% CF4

# Rango suave para curvas teóricas
fCF4_range = np.linspace(np.min(fCF4_real), np.max(fCF4_real), 1000)

# Colormaps
cmap_exp = cm.get_cmap("viridis")
cmap_fit = cm.get_cmap("hot")

# Asegurar carpeta de salida
os.makedirs("output", exist_ok=True)

# ===========================================================
#  PLOTEAR VISIBLE (α guardado)
# ===========================================================

plt.figure()
i_fit = 0
i_exp = 0

for nameCF3 in name_CF3:
    # en tu ajuste sólo usas "Ar** all", pero aquí dejamos general
    for nameAr in name_Ar_dbleStar:

        colname = f"{nameCF3}_{nameAr}"
        if colname not in df_alpha.columns:
            continue

        alpha = df_alpha.loc[0, colname]
        if pd.isna(alpha):
            continue

        # --- Poblaciones en puntos experimentales ---
        pop_CF3_grid = pCF3[nameCF3].to_numpy()
        pop_Ar_grid  = pArDbleStar[nameAr].to_numpy()

        PCF3_exp         = np.interp(fCF4_real/100.0, fCF4_grid, pop_CF3_grid)
        PAr_dbleStar_exp = np.interp(fCF4_real/100.0, fCF4_grid, pop_Ar_grid)

        values0_vis = (
            fCF4_real[index_ref] / 100.0,
            PCF3_exp[index_ref],
            PAr_dbleStar_exp[index_ref],
        )

        # --- Poblaciones en rango suave ---
        PCF3_range = np.interp(fCF4_range/100.0, fCF4_grid, pop_CF3_grid)
        PAr_range  = np.interp(fCF4_range/100.0, fCF4_grid, pop_Ar_grid)

        # Modelo con α ajustado
        model_fit = Pgamma_vis_Cociente(
            fCF4_range/100.0,
            PCF3_range,
            PAr_range,
            values0_vis,
            alpha
        )

        color_fit = cmap_fit(i_fit / 10)
        i_fit += 1

        plt.plot(
            fCF4_range,
            model_fit,
            color=color_fit,
            label=f"{colname}"
        )

        # Modelo sin ajuste (α = 1)
        model_nofit = Pgamma_vis_Cociente(
            fCF4_range/100.0,
            PCF3_range,
            PAr_range,
            values0_vis,
            1.0
        )

        color_fit = cmap_fit(i_fit / 10)
        i_fit += 1

        plt.plot(
            fCF4_range,
            model_nofit,
            color=color_fit,
            linestyle="--",
            label=f"{colname} sin ajuste"
        )

# --- Datos experimentales VIS ---
bars = ["1.0bar", "4.0bar"]

for b in bars:
    # misma estructura de error que en tu script de plotting anterior
    yerr = np.sqrt(
        (yield_vis["Err " + b] / yield_vis[b])**2 +
        (yield_vis["Err " + b] * yield_vis[b].to_numpy() /
         yield_vis[b].to_numpy()[index_ref]**2)**2
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
plt.xlabel("fCF4 real [%]")
plt.ylabel("Cociente VIS normalizado")
plt.legend(fontsize=7, markerscale=0.8, framealpha=0.9)
plt.title("Visible — Ajuste con α guardado")
plt.tight_layout()
plt.savefig("output/VIS-Normalizado_Degrad.pdf", dpi=300)


# ===========================================================
#  PLOTEAR ULTRAVIOLETA (kcool, kdis guardados)
# ===========================================================

plt.figure()
i_fit = 0
i_exp = 0

for nameCF4 in name_CF4:
    for nameAr3 in name_Ar_3rd:

        colname = f"{nameCF4}_{nameAr3}"
        if colname not in df_kcool.columns:
            continue

        kcool = df_kcool.loc[0, colname]
        kdis  = df_kdis.loc[0, colname]
        if pd.isna(kcool) or pd.isna(kdis):
            continue

        # --- Poblaciones en puntos experimentales ---
        pop_CF4_grid = pCF4[nameCF4].to_numpy()
        pop_Ar3_grid = pAr3rd[nameAr3].to_numpy()

        PCF4_exp  = np.interp(fCF4_real/100.0, fCF4_grid, pop_CF4_grid)
        PAr3_exp  = np.interp(fCF4_real/100.0, fCF4_grid, pop_Ar3_grid)

        n0 = 1.0
        values0_uv = (
            fCF4_real[index_ref] / 100.0,
            PCF4_exp[index_ref],
            PAr3_exp[index_ref],
            n0
        )

        # --- Poblaciones en rango suave ---
        PCF4_range = np.interp(fCF4_range/100.0, fCF4_grid, pop_CF4_grid)
        PAr3_range = np.interp(fCF4_range/100.0, fCF4_grid, pop_Ar3_grid)

        # 1 bar con ajuste
        model_1bar = Pgamma_UV_Cociente(
            fCF4_range/100.0,
            PCF4_range,
            PAr3_range,
            1.0,
            values0_uv,
            kcool,
            kdis
        )

        # 1 bar sin ajuste (kcool=kdis=0)
        model_1bar_nofit = Pgamma_UV_Cociente(
            fCF4_range/100.0,
            PCF4_range,
            PAr3_range,
            1.0,
            values0_uv,
            0.0,
            0.0
        )

        # 4 bar con ajuste
        model_4bar = Pgamma_UV_Cociente(
            fCF4_range/100.0,
            PCF4_range,
            PAr3_range,
            4.0,
            values0_uv,
            kcool,
            kdis
        )

        i_fit += 1
        color_fit = cmap_fit(i_fit / 6)
        plt.plot(fCF4_range, model_1bar, label=f"1bar {colname}", color=color_fit)

        i_fit += 1
        color_fit = cmap_fit(i_fit / 6)
        plt.plot(
            fCF4_range,
            model_1bar_nofit,
            label=f"1bar {colname} - sin ajuste",
            color=color_fit,
            linestyle="--"
        )

        i_fit += 1
        color_fit = cmap_fit(i_fit / 6)
        plt.plot(fCF4_range, model_4bar, label=f"4bar {colname}", color=color_fit)


# --- Datos experimentales UV ---
for b in bars:
    yerr = np.sqrt(
        (yield_uv["Err " + b] / yield_uv[b])**2 +
        (yield_uv["Err " + b] * yield_uv[b].to_numpy() /
         yield_uv[b].to_numpy()[index_ref]**2)**2
    )

    color_exp = cmap_exp(i_exp / len(bars))
    i_exp += 1

    if yield_uv[b].to_numpy()[index_ref] == 0.0:
        plt.errorbar(
            fCF4_real[:-1],
            yield_uv[b].to_numpy()[:-1] / yield_uv[b].to_numpy()[index_ref],
            yerr=yerr[:-1],
            fmt="o",
            capsize=5,
            label=b,
            color=color_exp
        )
    else:
        plt.errorbar(
            fCF4_real,
            yield_uv[b].to_numpy() / yield_uv[b].to_numpy()[index_ref],
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
plt.title("UV — Ajuste con kcool, kdis guardados (Degrad)")
plt.legend(fontsize=7, markerscale=0.8, framealpha=0.9)
plt.tight_layout()
plt.savefig("output/UV-Normalizado_Degrad.pdf", dpi=300)

print("Listo. Gráficas generadas en output/ (VIS-Normalizado_Degrad.pdf, UV-Normalizado_Degrad.pdf)")
