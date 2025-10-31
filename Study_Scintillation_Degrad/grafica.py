import pandas as pd
import matplotlib.pyplot as plt
import os

# === Carpeta donde están los pickles ===
DATA_DIR = "pickle_data"

# === 1. Cargar los pickles ===
scCF3 = pd.read_pickle(os.path.join(DATA_DIR, "scintillationProbability_PAModel_CF3.pkl"))
scCF4 = pd.read_pickle(os.path.join(DATA_DIR, "scintillationProbability_PAModel_CF4.pkl"))
scAr  = pd.read_pickle(os.path.join(DATA_DIR, "scintillationProbability_PAModel_Ar3rd.pkl"))
yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

print("✅ Pickles cargados correctamente desde 'pickle_data/'.\n")

# === 2. Verificaciones ===
if "fCF4" not in scCF3.columns:
    raise ValueError("⚠️ Falta columna 'fCF4' en CF3.")
if "all" not in scCF4.columns or "all" not in scAr.columns:
    raise ValueError("⚠️ CF4 o Ar3rd no contienen la columna 'all'.")
if "fCF4" not in yield_uv.columns or "fCF4" not in yield_vis.columns:
    raise ValueError("⚠️ Falta columna 'fCF4' en yield_uv o yield_vis.")

# === 3. Preparar figura única ===
plt.figure(figsize=(8, 6))
fCF4 = scCF3["fCF4"]

# === 4. (CF4 + Ar) / CF3 para cada modo ===
for col in scCF3.columns:
    if col == "fCF4":
        continue
    ratio_scint = (scCF4["all"] + scAr["all"]) / scCF3[col]
    plt.plot(fCF4*100, ratio_scint, marker=',', linestyle='-', label=f"Scint {col}")

# === 5. Y_uv / Y_vis para cada columna común ===
fCF4_yield = yield_uv["fCF4"]
for col in yield_uv.columns:
    if col != "fCF4" and col in yield_vis.columns:
        ratio_yield = yield_uv[col] / yield_vis[col]
        plt.plot(fCF4_yield, ratio_yield, marker='s', linestyle='--',linewidth=2, label=f"Yield {col}")
        print(ratio_yield)

# === 6. Formato del gráfico ===
plt.title("Comparación: (CF4 + Ar) / CF3  y  Y_uv / Y_vis", fontsize=13)
plt.xlabel("fCF4", fontsize=12)
plt.ylabel("Razón (adimensional)", fontsize=12)
plt.xscale("log")  # útil si fCF4 varía logarítmicamente
plt.yscale("log")  # útil si fCF4 varía logarítmicamente
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.ylim(0,100)
plt.legend(title="Curvas", fontsize=8)
plt.tight_layout()
plt.savefig("raw_data.pdf")
