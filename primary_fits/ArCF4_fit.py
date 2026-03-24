import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns

models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../models'))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))


sys.path.append(models_dir)
from ArCF4 import *

sys.path.append(data_dir)
from read_Degrad import read_degrad
from read_experimental import read_experimental
from fiting import fitParameters
from parameter_export import export_fit_table_latex, export_to_csv
from ploting import plot_fit_vs_experiment_by_pressure

#########################################################
####### CREAMOS LOS ARCHIVOS + LOS CARGAMOS 


archivo_entrada=np.array(["/output_99.9Ar_0.1CF4.txt",
                   "/output_99.8Ar_0.2CF4.txt",
                   "/output_99.5Ar_0.5CF4.txt",
                   "/output_99Ar_1CF4.txt",
                   "/output_98Ar_2CF4.txt",
                   "/output_95Ar_5CF4.txt",
                   "/output_90Ar_10CF4.txt",
                   "/output_80Ar_20CF4.txt",
                   "/output_50Ar_50CF4.txt",
                   "/output_PureCF4.txt"])

archivo_salida_1=np.array(["/ar_degrad_output_99.9Ar_0.1CF4.csv",
                   "/ar_degrad_output_99.8Ar_0.2CF4.csv",
                   "/ar_degrad_output_99.5Ar_0.5CF4.csv",
                   "/ar_degrad_output_99Ar_1CF4.csv",
                   "/ar_degrad_output_98Ar_2CF4.csv",
                   "/ar_degrad_output_95Ar_5CF4.csv",
                   "/ar_degrad_output_90Ar_10CF4.csv",
                   "/ar_degrad_output_80Ar_20CF4.csv",
                   "/ar_degrad_output_50Ar_50CF4.csv",
                   "/ar_degrad_output_PureCF4.csv"])

archivo_salida_2=np.array(["/cf4_degrad_output_99.9Ar_0.1CF4.csv",
                   "/cf4_degrad_output_99.8Ar_0.2CF4.csv",
                   "/cf4_degrad_output_99.5Ar_0.5CF4.csv",
                   "/cf4_degrad_output_99Ar_1CF4.csv",
                   "/cf4_degrad_output_98Ar_2CF4.csv",
                   "/cf4_degrad_output_95Ar_5CF4.csv",
                   "/cf4_degrad_output_90Ar_10CF4.csv",
                   "/cf4_degrad_output_80Ar_20CF4.csv",
                   "/cf4_degrad_output_50Ar_50CF4.csv",
                   "/cf4_degrad_output_PureCF4.csv"])


prefijo = "../data/Primary_DegradData/ArCF4/txt"
archivo_entrada = np.char.add(prefijo, archivo_entrada)

prefijo = "../data/Primary_DegradData/ArCF4/csv"
archivo_salida_1 = np.char.add(prefijo, archivo_salida_1)
archivo_salida_2 = np.char.add(prefijo, archivo_salida_2)

gas1 = "ARGON"
gas2 = "CF4"
concentration = np.array([0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1])   

dataframe = pd.DataFrame(
    {
        "CF4":    [["ION CF3 +"],                            "CF4",      0, 100, "CF4"],

        "Ar**":   [["EXC"],                                  "ARGON",    0, 100, "Ar_dbleStar"],

        "CF3":    [["NEUTRAL DISS"],                         "CF4",      0, 100, "CF3"],

        "Ar3rd":  [["CHARGE STATE =2","CHARGE STATE "],      "ARGON",    0, 100, "Ar_3rd"]
        
    }, 
    index=["name principal", "gas", "energy low", "energy up", "name output"]
)

output_dir = "../data/Primary_DegradData/ArCF4/"
output_general_name =  "../data/Primary_DegradData/ArCF4"

read_degrad(archivo_entrada, archivo_salida_1, archivo_salida_2, gas1, gas2, concentration, dataframe, output_dir, output_general_name)


archivo_entrada = "../data/Experimental/ArCF4/CF4_primary_data_final.pkl"
yields = ["vis","UV"]
presiones = [1,2,2.5,3,4,5]
concentraciones_reales= None
no_sistematic = True



output_dir = "../data/Experimental/ArCF4/"

read_experimental(archivo_entrada, yields, presiones, output_dir, concentraciones_reales=concentraciones_reales, no_sistematic = True)

#####################################################
###### Traemos los datos anteriormente generados 

DATA_DIR = "../data/Experimental/ArCF4/"
yield_uv  = pd.read_csv(os.path.join(DATA_DIR, "UV.csv"))
yield_vis = pd.read_csv(os.path.join(DATA_DIR, "vis.csv"))


DATA_DIR = "../data/Primary_DegradData"
degrad_data        = pd.read_csv(os.path.join(DATA_DIR, "ArCF4.csv"))


#########################################################3
####### AJUSTE


x0 = np.array([1,
               0.11352059, 0.00156568 , 999, 
               0.01793004 ,0.1, 0.31565123, 50.1, 0.99,
               999 ])
lower = [0.0,
         0.0, 0.0, 0.0,  
         0.0, 0.065, 0.0, 50, 0.0,
         0.0]

upper = [1.0, 
         1.0, 1.0, 10000.0, 
         10.0, 10.0, 1.0, 50.2, 1.0,
         10000.0]

bounds=(lower, upper)

equations = {
    "vis": theory_yield_vis,
    "uv": theory_yield_uv
}

experimental_data = {
    "vis": yield_vis,
    "uv": yield_uv
}

popt = fitParameters(equations, experimental_data, degrad_data, x0=x0, bounds=bounds)

J = popt.jac
m, p = J.shape
s2 = 2 * popt.cost / (m - p)
cov_theta =  s2 * np.linalg.inv(J.T @ J)
chi2 = 2 * popt.cost
N_res = popt.fun.size
N_par = popt.x.size
dof   = N_res - N_par
chi2_red = chi2 / dof

print("="*60)
print("Parámetros globales:", popt.x)
print(f"Chi2 (real): {chi2}")
print(f"Grados de libertad: {dof}")
print(f"Chi2 reducido: {chi2_red}")
print("="*60)


#######################################################################
# =================== LATEX, TYPST, CSV EXPORT ========================
#######################################################################

names_tex = [
    "$N_{\\text{norm}}$",
    "$P_{\\mathrm{CF_3}}|_{\\mathrm{dir}}$",
    "$P_{\\mathrm{Ar}^{**}} $",
    "${K_{\\mathrm{Ar^{**},Q(Ar)}}} [ns]$",
    "$1 / {\\tau_{\\mathrm{disocc}} K_{\\mathrm{relax}}}$",
    "$\\tau_{\mathrm{uv}} K_{\mathrm{CF_4^(+,*)Q(CF_4)}}$",
    "$P_{\\mathrm{CF_4^{+,*}}}|_{\\mathrm{dir}}$",
    "$K_{\\mathrm{Ar^{++},Q(CF_4)}}$ [ns]",
    "$P_{\\mathrm{Ar}^{++}}$",
    "${K_{\\mathrm{Ar^{**},Q(CF_4)}}} [ns]$"
]

latex_table, _, perr, rel = export_fit_table_latex(
    result=popt,
    names=names_tex,
    filename="tex_param/ArCF4_param.tex",
    caption="Parámetros obtenidos del ajuparamste global.",
    label="tab:fit_params",
    sigfigs=4
)

names_csv = [
    "Nnorm",
    "PCF3dir$",
    "PAr**",
    "KAr**QAr",
    "1/tauDiscKrelax",
    "tauUvKCF4QCF4",
    "PCF4dir",
    "KAr++QCF4",
    "PAr++",
    "KAr**QCF4"
]

export_to_csv("../data/Parameters/ArCF4_primary.csv",popt,names_csv)

#######################################################################
# =================== PLOT ========================
#######################################################################

pressure = [1,3,5]

concentrations = np.logspace(-4, 0, 1000)
yield_vis_plot = yield_vis["1.0bar"]

fig, ax, pressure_cols = plot_fit_vs_experiment_by_pressure(
    df_exp=yield_vis,
    theory_func=theory_yield_vis,
    fit_params=popt.x,
    degrad_data=degrad_data,
    concentration_grid=concentrations,
    pressures = pressure,
    x_col="fCF4",
    x_plot_factor=100,
    min_positive_x=1e-3,
    title="Visible",
    xlabel=r"Concentration of CF$_4$ [%]",
    ylabel="Normalized Yield",
    xlim=(0.1 * 0.9, 100 * 1.1),
    ylim=(0.001, 0.4),
    xscale="log",
    yscale="log",
    cmap="inferno",
    darken_factor=-0.15,
    legend=True,
    legend_kwargs={"ncol": 2, "fontsize": 9},
    output="plots/ArCF4_visible.pdf",
    show=False,
)


concentrations = np.logspace(-6, 0, 1000)

fig, ax, pressure_cols = plot_fit_vs_experiment_by_pressure(
    df_exp=yield_uv,
    theory_func=theory_yield_uv,
    fit_params=popt.x,
    degrad_data=degrad_data,
    concentration_grid=concentrations,
    pressures = pressure,
    x_col="fCF4",
    x_plot_factor=100,
    min_positive_x=1e-5,
    title="Ultraviolet",
    xlabel=r"Concentration of CF$_4$ [%]",
    ylabel="Normalized Yield",
    xlim=(0.001 * 0.9, 100 * 1.1),
    ylim=None,
    xscale="log",
    yscale="log",
    cmap="plasma",
    darken_factor=-0.15,
    legend=True,
    legend_kwargs={"ncol": 2, "fontsize": 9},
    output="plots/ArCF4_uv.pdf",
    show=False,
)

#######################################################################
# =================== CORRELATION MATRIX ========================
#######################################################################


# Construimos matriz de correlación a partir de covarianzas
diag = np.sqrt(np.diag(cov_theta))
outer = np.outer(diag, diag)
corr = cov_theta / outer

# Seguridad numérica
corr = np.clip(corr, -1, 1)

# DataFrame para seaborn
corr_df = pd.DataFrame(corr, columns=names_tex, index=names_tex)

# --- Plot estilo seaborn ---
plt.figure(figsize=(10, 8))
sns.heatmap(
    corr_df,
    cmap="coolwarm",
    vmin=-1,
    vmax=1,
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    square=True,
    cbar_kws={"label": "Correlación"}
)
plt.title("Matriz de Correlación de Parámetros Ajustados", fontsize=14)
plt.tight_layout()

plt.savefig("plots/ArCF4_CorrelationMatrix_GlobalFit.pdf", dpi=300)