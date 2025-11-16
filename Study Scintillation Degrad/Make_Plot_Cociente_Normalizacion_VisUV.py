

import numpy as np
import pandas as pd
import dill
import matplotlib.pyplot as plt
import os
from Amoedo_Model_DivisionFit import Pgamma_UV_Cociente, Pgamma_vis_Cociente,Pgamma_CF3_refined,Pgamma_CF4_refined,Pgamma_Ar3rd_refined
import matplotlib.cm as cm


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
yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))


# === Leer parámetros ajustados ===
df_alpha = pd.read_pickle(os.path.join(FIT_DIR, "alpha_results.pkl"))
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
#   NUEVOS PLOTS — AJUSTES usando funciones REFINED
#   normalizando a (Yvis + Yuv) @ 1bar, fCF4 = 100%
# ===========================================================

print("Generando gráficos con ajustes (refined) normalizados a (Yvis+Yuv)@1bar...")

bars = ["1.0bar", "2.0bar", "3.0bar", "4.0bar", "5.0bar"]
index_ref = -1

# Normalización experimental (fijada y conocida)
norm_den_exp = (
    yield_vis["1.0bar"].to_numpy()[index_ref] +
    yield_uv ["1.0bar"].to_numpy()[index_ref]
)

print("=="*40)
print("VIS 1bar 100%CF4:", yield_vis["1.0bar"].to_numpy()[index_ref])
print("UV  1bar 100%CF4:", yield_uv ["1.0bar"].to_numpy()[index_ref])
print("Denominador experimental:", norm_den_exp)
print("=="*40)


# ===========================================================
#  FIGURA A — VIS REFINED NORMALIZADO
# ===========================================================

plt.figure(figsize=(7, 5))

cmap_fit_ref = plt.colormaps["hot"]
cmap_exp_ref = plt.colormaps["viridis"]
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
        # --- poblaciones ---
        f_CF3   = lineal_CF3[nameCF3].to_numpy()[0]
        f_Ardbl = lineal_pAr_dbleStar[nameAr].to_numpy()[0]

        f_range = np.linspace(min(fCF4_real), max(fCF4_real), 5000)
        PCF3_range   = f_CF3(f_range / 100)
        PArdbl_range = f_Ardbl(1 - f_range / 100)

        # --- Modelo VIS REFINED ---
        VIS_model_ref = Pgamma_CF3_refined(
            f_range/100,
            PCF3_range,
            PArdbl_range,
            alpha
        )

        # --- Modelo UV para misma combinación? NO existe ---
        # Aquí está el error original:
        # necesitas normalizar el VIS con **un UV_model del mismo bin de ajuste**
        # pero CF3-α y CF4-kcool no están emparejados.
        # El ajuste correcto es:
        # normalizar solo por VIS_model_ref[-1]
        # si quieres normalizar VIS y UV juntos, necesitas un pairing.
        # asumimos que quieres usar SOLO VIS_model para normalizar VIS:
        
        # pero tú querías VIS/(VIS+UV teórico)
        # Para eso elegimos el UV_model correspondiente a *misma presión 1 bar* y CF4→100%
        # pero como no hay pairing directo, usamos los valores del ajuste UV global:
        # calculamos UV_model_ref aquí:
        UV_model_ref = None
        
        # buscamos la columna en kcool df que coincida con nameCF3 versión CF4
        # pero no existe relación directa -> usamos modelo UV "all"
        # pongo un UV_model_ref dummy calculado con mejores parámetros globales:
        for nameCF4 in name_CF4:
            for nameAr3rd in name_Ar_3rd:

                colname2 = f"{nameCF4}_{nameAr3rd}"
                if colname2 not in df_kcool.columns:
                    continue

                kcool_tmp = df_kcool.loc[0, colname2]
                kdis_tmp  = df_kdis.loc[0, colname2]
                if pd.isna(kcool_tmp) or pd.isna(kdis_tmp):
                    continue
                
                #kcool_tmp,kdis_tmp=0,0

                f_CF4   = lineal_CF4[nameCF4].to_numpy()[0]
                f_Ar3rd = lineal_pAr_3rd[nameAr3rd].to_numpy()[0]

                PCF4_range  = f_CF4(f_range/100)
                PAr3_range  = f_Ar3rd(1 - f_range/100)

                UV_model_ref = Pgamma_CF4_refined(
                    f_range/100, PCF4_range, PAr3_range,
                    1.0, kcool_tmp, kdis_tmp
                ) 
                break
            if UV_model_ref is not None:
                break

        if UV_model_ref is None:
            continue

        # --- normalización correcta por modelo propio ---
        norm_den_fit = VIS_model_ref[-1] + UV_model_ref[-1]

        VIS_model_norm = VIS_model_ref / norm_den_fit

        # --- Plot model ---
        color_fit = cmap_fit_ref(i_fit / 6)
        plt.plot(
            f_range,
            VIS_model_norm,
            color=color_fit,
            label=f"modelo VIS {colname}"
        )
        i_fit += 1


# --- Datos VIS ---
for b in bars:

    y_vis_norm = yield_vis[b].to_numpy() / norm_den_exp
    yerr_vis_norm = yield_vis["Err " + b].to_numpy() / norm_den_exp

    color_exp = cmap_exp_ref(i_exp / len(bars))
    i_exp += 1

    plt.errorbar(
        fCF4_real,
        y_vis_norm,
        yerr=yerr_vis_norm,
        fmt="o", capsize=4,
        color=color_exp,
        label=b
    )

plt.xscale("log")
plt.yscale("log")
plt.ylim(2e-3,1)
plt.xlabel("fCF4 real [%]")
plt.ylabel("Yvis / (Yvis + Yuv @1bar)")
plt.title("VIS — Modelo REFINED normalizado (cada curva con su propio denom teórico)")
plt.legend(fontsize=6, ncol=2, framealpha=0.9)
plt.savefig("output/VIS-Normalizado_VisUV.pdf", dpi=300)


# ===========================================================
#  FIGURA B — UV REFINED NORMALIZADO
# ===========================================================

plt.figure(figsize=(7,5))

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
        
        #kcool,kdis=0,0

        f_CF4   = lineal_CF4[nameCF4].to_numpy()[0]
        f_Ar3rd = lineal_pAr_3rd[nameAr3rd].to_numpy()[0]

        f_range = np.linspace(min(fCF4_real), max(fCF4_real), 5000)
        PCF4_range = f_CF4(f_range/100)
        PAr3_range = f_Ar3rd(1 - f_range/100)

        for n in [1,5]:
            # Modelo UV refined
            UV_model_ref = Pgamma_CF4_refined(
                f_range/100, PCF4_range, PAr3_range,
                n, kcool, kdis
            ) + Pgamma_Ar3rd_refined( f_range/100, PAr3_range,n)
            

            # Necesitamos VIS_model para normalizar — tomamos el VIS global:
            VIS_model_ref = None
            for nameCF3 in name_CF3:
                for nameAr in name_Ar_dbleStar:
                    if f"{nameCF3}_{nameAr}" in df_alpha.columns:
                        alpha_tmp = df_alpha.loc[0, f"{nameCF3}_{nameAr}"]
                        #alpha_tmp = 1
                        f_CF3   = lineal_CF3[nameCF3].to_numpy()[0]
                        f_Ardbl = lineal_pAr_dbleStar[nameAr].to_numpy()[0]

                        PCF3_range = f_CF3(f_range/100)
                        PArdbl_range = f_Ardbl(1 - f_range/100)

                        VIS_model_ref = Pgamma_CF3_refined(
                            f_range/100, PCF3_range, PArdbl_range,
                            alpha_tmp
                        )
                        break
                if VIS_model_ref is not None:
                    break

            if VIS_model_ref is None:
                continue

            # Normalización por modelo correspondiente
            norm_den_fit = VIS_model_ref[-1] + UV_model_ref[-1]
            UV_model_norm = UV_model_ref / norm_den_fit

            color_fit = cmap_fit_ref(i_fit / 6)
            plt.plot(f_range, UV_model_norm,
                    color=color_fit,
                    label=f"modelo UV {colname} {n:.1f} bar")
            i_fit += 1


# --- datos UV ---
for b in bars:

    y_uv_norm = yield_uv[b].to_numpy() / norm_den_exp
    yerr_uv_norm = yield_uv["Err " + b].to_numpy() / norm_den_exp

    color_exp = cmap_exp_ref(i_exp / len(bars))
    i_exp += 1

    plt.errorbar(
        fCF4_real,
        y_uv_norm,
        yerr=yerr_uv_norm,
        fmt="o", capsize=4,
        color=color_exp,
        label=b
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("fCF4 real [%]")
plt.ylabel("Yuv / (Yvis + Yuv @1bar)")
plt.title("UV — Modelo REFINED normalizado (cada curva con su propio denom teórico)")
plt.legend(fontsize=6, ncol=2)
plt.savefig("output/UV-Normalizado_VisUV.pdf", dpi=300)

print("Gráficas refined generadas con normalización individual por fit.")
