import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns

models_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../models'))
data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))

sys.path.append(models_dir)
sys.path.append(data_dir)

from ArN2 import *
from read_Degrad import read_degrad
from read_experimental import read_experimental
from fiting import fitParameters
from parameter_export import export_fit_table_latex, export_to_csv
from ploting import plot_fit_vs_experiment_by_pressure

#########################################################
####### CREAMOS LOS ARCHIVOS + LOS CARGAMOS 


archivo_entrada=np.array(["/output_Argon_0.1_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_Argon_0.5_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_Argon_1.0_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_Argon_2.0_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_Argon_5.0_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_Argon_10.0_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_Argon_20.0_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_Argon_50.0_N2_E_0.0_Vcmbar_P_1_bar_12_keV.txt",
                          "/output_100.0N2_E_0.0Vcmbar_P_1bar_12keV.txt"])

archivo_salida_1=np.array(["/ar_degrad_output_99.9Ar_0.1N2.csv",
                   "/ar_degrad_output_99.5Ar_0.5N2.csv",
                   "/ar_degrad_output_99Ar_1N2.csv",
                   "/ar_degrad_output_98Ar_2N2.csv",
                   "/ar_degrad_output_95Ar_5N2.csv",
                   "/ar_degrad_output_90Ar_1N2.csv",
                   "/ar_degrad_output_80Ar_20N2.csv",
                   "/ar_degrad_output_50Ar_50N2.csv",
                   "/ar_degrad_output_PureN2.csv"])

archivo_salida_2=np.array(["/n2_degrad_output_99.9Ar_0.1N2.csv",
                   "/n2_degrad_output_99.5Ar_0.5N2.csv",
                   "/n2_degrad_output_99Ar_1N2.csv",
                   "/n2_degrad_output_98Ar_2N2.csv",
                   "/n2_degrad_output_95Ar_5N2.csv",
                   "/n2_degrad_output_90Ar_10N2.csv",
                   "/n2_degrad_output_80Ar_20N2.csv",
                   "/n2_degrad_output_50Ar_50N2.csv",
                   "/n2_degrad_output_PureN2.csv"])


prefijo = "../data/Primary_DegradData/ArN2/txt"
archivo_entrada = np.char.add(prefijo, archivo_entrada)

prefijo = "../data/Primary_DegradData/ArN2/csv"
archivo_salida_1 = np.char.add(prefijo, archivo_salida_1)
archivo_salida_2 = np.char.add(prefijo, archivo_salida_2)

gas1 = "ARGON"
gas2 = "NITROGEN"
concentration = np.array([0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1])   

dataframe = pd.DataFrame(
    {    
        "Ar**":   [["EXC"],     "ARGON",     0, 100, "Ar_dbleStar"],
        "N2*":    [[""], "NITROGEN",  13, 14.5, "N2_star"]
    }, 
    index=["name principal", "gas", "energy low", "energy up", "name output"]
)

output_dir = "../data/Primary_DegradData/ArN2/"
output_general_name =  "../data/Primary_DegradData/ArN2"

read_degrad(archivo_entrada, archivo_salida_1, archivo_salida_2, gas1, gas2, concentration, dataframe, output_dir, output_general_name)


archivo_entrada = "../data/Experimental/ArN2/N2_primary_data_final.pkl"
yields = ["yield_N2"]
presiones = [1,2,2.5,3,4,5]
concentraciones_reales= None
no_sistematic = True


output_dir = "../data/Experimental/ArN2/"

read_experimental(archivo_entrada, yields, presiones, output_dir, concentraciones_reales=concentraciones_reales, no_sistematic = no_sistematic)

#####################################################
###### Traemos los datos anteriormente generados 

DATA_DIR = "../data/Experimental/ArN2/"
yield_N2_uv  = pd.read_csv(os.path.join(DATA_DIR, "yield_N2.csv"))

"""
columns = yield_N2_uv.columns
for i, column in enumerate(columns):
    if "Err" in column:
        yield_N2_uv[column] = 0.3 * yield_N2_uv[columns[i-1]]
"""

DATA_DIR = "../data/Primary_DegradData"
degrad_data        = pd.read_csv(os.path.join(DATA_DIR, "ArN2.csv"))


#########################################################3
####### AJUSTE


lower = [0.0, 0.0 , 0.0 , 0, 0.0, 0]

x0 = np.array([0.99, 0.9, 0, 0, 0.00001, 0.00001])

upper = [1, 1, 100, 100, 100, 1000]

bounds=(lower, upper)

equations = {
    "vis": theory_yield_N2_uv,
}

experimental_data = {
    "vis": yield_N2_uv,
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
    "$P_{\\mathrm{Ar}^{**}} $",
    "$K_{\\mathrm{Ar}^{**}Q(\\mathrm{Ar}^{**})}/K_{\\mathrm{Ar}^{**}Q(\\mathrm{N}_2)} $",
    "$1 / {\\tau_{\\mathrm{Ar}^{**}} K_{\\mathrm{Ar}Q(\\mathrm{N}_2)}} $",
    "${\\tau_{\\mathrm{N}_2} K_{\\mathrm{N}_2Q(\\mathrm{Ar}^{**})}}$",
    "${\\tau_{\\mathrm{N}_2} K_{\\mathrm{N}_2Q(\\mathrm{N}_2)}}$"
]

latex_table, _, perr, rel = export_fit_table_latex(
    result=popt,
    names=names_tex,
    filename="tex_param/fit_table.tex",
    caption="Parámetros obtenidos del ajuparamste global.",
    label="tab:fit_params",
    sigfigs=4
)

names_csv = [
    "Nnorm",
    "PAr**dir$",
    "KArQAr/KArQN2",
    "1/tau KArQN2",
    "tau KN2QAr",
    "tau KN2QN2"
]

export_to_csv("../data/Parameters/ArN2_primary.csv",popt,names_csv)

#######################################################################
# =================== PLOT ========================
#######################################################################

pressure = [1,3,5]

concentrations = np.logspace(-4, 0, 1000)

fig, ax, pressure_cols = plot_fit_vs_experiment_by_pressure(
    df_exp=yield_N2_uv,
    theory_func=theory_yield_N2_uv,
    fit_params=popt.x,
    degrad_data=degrad_data,
    concentration_grid=concentrations,
    pressures = pressure,
    x_col="N2 concentration (%)",
    x_plot_factor=100,
    min_positive_x=1e-3,
    title="Emission in Ar-N$_2$",
    xlabel=r"Concentration of N$_2$ [%]",
    ylabel="Normalized Yield",
    xlim=(0.1 * 0.9, 100 * 1.1),
    ylim=(0.01, 4),
    xscale="log",
    yscale="log",
    cmap="inferno",
    darken_factor=-0.15,
    legend=True,
    legend_kwargs={"ncol": 2, "fontsize": 9},
    output="plots/ArN2_global.pdf",
    show=False,
    activate_components = False
)


pressure = [1]

fig, ax, pressure_cols = plot_fit_vs_experiment_by_pressure(
    df_exp=yield_N2_uv,
    theory_func=theory_yield_N2_uv,
    fit_params=popt.x,
    degrad_data=degrad_data,
    concentration_grid=concentrations,
    pressures = pressure,
    x_col="N2 concentration (%)",
    x_plot_factor=100,
    min_positive_x=1e-3,
    title="Emission in Ar-N$_2$",
    xlabel=r"Concentration of N$_2$ [%]",
    ylabel="Normalized Yield",
    xlim=(0.1 * 0.9, 100 * 1.1),
    ylim=(0.01, 4),
    xscale="log",
    yscale="log",
    cmap="inferno",
    darken_factor=-0.15,
    legend=True,
    legend_kwargs={"ncol": 2, "fontsize": 9},
    output="plots/ArN2_global_components.pdf",
    show=False,
    activate_components = True
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

plt.savefig("plots/ArN2_CorrelationMatrix_GlobalFit.pdf", dpi=300)