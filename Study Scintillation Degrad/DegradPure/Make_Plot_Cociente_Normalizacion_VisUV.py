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
    Pgamma_CF3_refined,
    Pgamma_CF4_refined,
    Pgamma_Ar3rd_refined,
)

"""
Script que leyendo:
    - poblaciones de Degrad (tabuladas, NO funciones lineales)
    - yields experimentales
    - parámetros ajustados (alpha, kcool, kdis)
vuelve a generar los mismos plots "refined" del ajuste,
normalizando a (Yvis + Yuv)@1bar, fCF4 = 100%.
"""

# ================================
# RUTAS A LOS PICKLES
# ================================
DATA_DIR = os.path.join("..", "pickle_data")
FIT_DIR  = "pickle_data"

# === Leer poblaciones Degrad (tabuladas) ===
pCF4        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
pCF3        = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF3.pkl"))
pArDbleStar = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_dbleStar.pkl"))
pAr3rd      = pd.read_pickle(os.path.join(DATA_DIR, "poblations_Ar_3rd.pkl"))

# === Leer yields experimentales ===
yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

# === Leer parámetros ajustados ===
df_alpha = pd.read_pickle(os.path.join(FIT_DIR, "alpha_results.pkl"))
df_kcool = pd.read_pickle(os.path.join(FIT_DIR, "kcool_results.pkl"))
df_kdis  = pd.read_pickle(os.path.join(FIT_DIR, "kdis_results.pkl"))

# ============================
# Variables útiles
# ============================
# malla de Degrad en fracción (0–1)
fCF4_grid = pCF4["fCF4"].to_numpy()

# nombres de columnas de poblaciones (como en el script de ajuste)
name_CF4         = pCF4.columns.to_numpy()[1::2]
name_CF3         = pCF3.columns.to_numpy()[1::2]
name_Ar_dbleStar = pArDbleStar.columns.to_numpy()[1::2]
name_Ar_3rd      = pAr3rd.columns.to_numpy()[1::2]

fCF4_real = yield_uv["fCF4 real"].to_numpy()  # en %
index_ref = -1   # 100% CF4

# ============================
# Normalización experimental
# ============================
bars = ["1.0bar", "2.0bar", "3.0bar", "4.0bar", "5.0bar"]

norm_den_exp = (
    yield_vis["1.0bar"].to_numpy()[index_ref] +
    yield_uv ["1.0bar"].to_numpy()[index_ref]
)

print("=="*40)
print("VIS 1bar 100%CF4:", yield_vis["1.0bar"].to_numpy()[index_ref])
print("UV  1bar 100%CF4:", yield_uv ["1.0bar"].to_numpy()[index_ref])
print("Denominador experimental:", norm_den_exp)
print("=="*40)

# Colormaps
cmap_fit_ref = cm.get_cmap("hot")
cmap_exp_ref = cm.get_cmap("viridis")

# Rango suave para fCF4 (en %)
f_range = np.linspace(np.min(fCF4_real), np.max(fCF4_real), 5000)

# Asegurar carpeta de salida
os.makedirs("output", exist_ok=True)

# ===========================================================
#  FIGURA A — VIS REFINED NORMALIZADO
# ===========================================================
print("Generando VIS refined normalizado a (Yvis+Yuv)@1bar...")

plt.figure(figsize=(7, 5))

i_fit = 0
i_exp = 0

for nameCF3 in name_CF3:
    for nameAr in name_Ar_dbleStar:

        colname = f"{nameCF3}_{nameAr}"
        if colname not in df_alpha.columns:
            continue

        alpha = df_alpha.loc[0, colname]
        if pd.isna(alpha):
            continue

        # --- Poblaciones para esta combinación (tabuladas) ---
        pop_CF3_grid   = pCF3[nameCF3].to_numpy()
        pop_Ardbl_grid = pArDbleStar[nameAr].to_numpy()

        # Interpolamos a f_range/100 (pasamos de % a fracción)
        PCF3_range   = np.interp(f_range/100.0, fCF4_grid, pop_CF3_grid)
        PArdbl_range = np.interp(f_range/100.0, fCF4_grid, pop_Ardbl_grid)

        # --- Modelo VIS refined ---
        VIS_model_ref = Pgamma_CF3_refined(
            f_range/100.0,
            PCF3_range,
            PArdbl_range,
            alpha
        )

        # === Construimos UV_model_ref "a la opción A"
        #    (primer CF4_Ar3rd con kcool,kdis válidos) ===
        UV_model_ref = None
        for nameCF4 in name_CF4:
            for nameAr3rd in name_Ar_3rd:
                colname2 = f"{nameCF4}_{nameAr3rd}"
                if colname2 not in df_kcool.columns:
                    continue

                kcool_tmp = df_kcool.loc[0, colname2]
                kdis_tmp  = df_kdis.loc[0, colname2]
                if pd.isna(kcool_tmp) or pd.isna(kdis_tmp):
                    continue

                pop_CF4_grid = pCF4[nameCF4].to_numpy()
                pop_Ar3_grid = pAr3rd[nameAr3rd].to_numpy()

                PCF4_range = np.interp(f_range/100.0, fCF4_grid, pop_CF4_grid)
                PAr3_range = np.interp(f_range/100.0, fCF4_grid, pop_Ar3_grid)

                UV_model_ref = Pgamma_CF4_refined(
                    f_range/100.0,
                    PCF4_range,
                    PAr3_range,
                    1.0,
                    kcool_tmp,
                    kdis_tmp
                )
                break
            if UV_model_ref is not None:
                break

        if UV_model_ref is None:
            # no hay ningún UV válido -> saltamos
            continue

        # --- Normalización teórica (opción A: VIS_ref + UV_ref en el último punto) ---
        norm_den_fit = VIS_model_ref[-1] + UV_model_ref[-1]
        VIS_model_norm = VIS_model_ref / norm_den_fit

        color_fit = cmap_fit_ref(i_fit / 6)
        plt.plot(
            f_range,
            VIS_model_norm,
            color=color_fit,
            label=f"modelo VIS {colname}"
        )
        i_fit += 1

# --- Datos VIS experimentales, misma normalización experimental ---
for b in bars:
    y_vis_norm     = yield_vis[b].to_numpy()      / norm_den_exp
    yerr_vis_norm  = yield_vis["Err " + b].to_numpy() / norm_den_exp

    color_exp = cmap_exp_ref(i_exp / len(bars))
    i_exp += 1

    plt.errorbar(
        fCF4_real,
        y_vis_norm,
        yerr=yerr_vis_norm,
        fmt="o",
        capsize=4,
        color=color_exp,
        label=b
    )

plt.xscale("log")
plt.yscale("log")
plt.ylim(2e-3, 1)
plt.xlabel("fCF4 real [%]")
plt.ylabel("Yvis / (Yvis + Yuv @1bar)")
plt.title("VIS — Modelo REFINED normalizado (opción A)")
plt.legend(fontsize=6, ncol=2, framealpha=0.9)
plt.tight_layout()
plt.savefig("output/VIS-Normalizado_VisUV_Degrad.pdf", dpi=300)


# ===========================================================
#  FIGURA B — UV REFINED NORMALIZADO
# ===========================================================
print("Generando UV refined normalizado a (Yvis+Yuv)@1bar...")

plt.figure(figsize=(7, 5))

i_fit = 0
i_exp = 0

for nameCF4 in name_CF4:
    for nameAr3rd in name_Ar_3rd:

        colname = f"{nameCF4}_{nameAr3rd}"
        if colname not in df_kcool.columns:
            continue

        kcool = df_kcool.loc[0, colname]
        kdis  = df_kdis.loc[0, colname]
        if pd.isna(kcool) or pd.isna(kdis):
            continue

        pop_CF4_grid = pCF4[nameCF4].to_numpy()
        pop_Ar3_grid = pAr3rd[nameAr3rd].to_numpy()

        PCF4_range = np.interp(f_range/100.0, fCF4_grid, pop_CF4_grid)
        PAr3_range = np.interp(f_range/100.0, fCF4_grid, pop_Ar3_grid)

        for n in [1.0, 4.0]:
            # Modelo UV refined (CF4 + Ar3rd)
            UV_model_ref = (
                Pgamma_CF4_refined(
                    f_range/100.0,
                    PCF4_range,
                    PAr3_range,
                    n,
                    kcool,
                    kdis
                )
                + Pgamma_Ar3rd_refined(
                    f_range/100.0,
                    PAr3_range,
                    n
                )
            )

            # --- VIS_model_ref "global" como en opción A:
            VIS_model_ref = None
            for nameCF3 in name_CF3:
                for nameAr in name_Ar_dbleStar:
                    colname_vis = f"{nameCF3}_{nameAr}"
                    if colname_vis not in df_alpha.columns:
                        continue

                    alpha_tmp = df_alpha.loc[0, colname_vis]
                    if pd.isna(alpha_tmp):
                        continue

                    pop_CF3_grid   = pCF3[nameCF3].to_numpy()
                    pop_Ardbl_grid = pArDbleStar[nameAr].to_numpy()

                    PCF3_range = np.interp(f_range/100.0, fCF4_grid, pop_CF3_grid)
                    PArdbl_range = np.interp(f_range/100.0, fCF4_grid, pop_Ardbl_grid)

                    VIS_model_ref = Pgamma_CF3_refined(
                        f_range/100.0,
                        PCF3_range,
                        PArdbl_range,
                        alpha_tmp
                    )
                    break
                if VIS_model_ref is not None:
                    break

            if VIS_model_ref is None:
                continue

            # Normalización teórica (opción A)
            norm_den_fit = VIS_model_ref[-1] + UV_model_ref[-1]
            UV_model_norm = UV_model_ref / norm_den_fit

            color_fit = cmap_fit_ref(i_fit / 6)
            plt.plot(
                f_range,
                UV_model_norm,
                color=color_fit,
                label=f"modelo UV {colname} {n:.1f} bar"
            )
            i_fit += 1

# --- Datos UV experimentales, misma normalización experimental ---
for b in bars:
    y_uv_norm    = yield_uv[b].to_numpy()        / norm_den_exp
    yerr_uv_norm = yield_uv["Err " + b].to_numpy() / norm_den_exp

    color_exp = cmap_exp_ref(i_exp / len(bars))
    i_exp += 1

    plt.errorbar(
        fCF4_real,
        y_uv_norm,
        yerr=yerr_uv_norm,
        fmt="o",
        capsize=4,
        color=color_exp,
        label=b
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("fCF4 real [%]")
plt.ylabel("Yuv / (Yvis + Yuv @1bar)")
plt.title("UV — Modelo REFINED normalizado (opción A)")
plt.legend(fontsize=6, ncol=2, framealpha=0.9)
plt.tight_layout()
plt.savefig("output/UV-Normalizado_VisUV_Degrad.pdf", dpi=300)

print("Listo. Gráficas refined generadas en output/:")
print("  - VIS-Normalizado_VisUV_Degrad.pdf")
print("  - UV-Normalizado_VisUV_Degrad.pdf")
