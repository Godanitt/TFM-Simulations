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

        "Ar3rd":  [["CHARGE STATE =2"],      "ARGON",    0, 100, "Ar_3rd"]
        
    }, 
    index=["name principal", "gas", "energy low", "energy up", "name output"]
)

output_dir = "../data/Primary_DegradData/ArCF4/"
output_general_name =  "../data/Primary_DegradData/ArCF4"

read_degrad(archivo_entrada, archivo_salida_1, archivo_salida_2, gas1, gas2, concentration, dataframe, output_dir, output_general_name)


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
# 6) CARGA DE DATOS PARA EL MODELO
# ============================================================

DATA_DIR_EXP = "../data/Experimental/ArCF4/"
DATA_DIR_PAR = "../data/Parameters"
DATA_DIR_LEV = "../data/Secondary_GarfieldData/levels"
DATA_DIR_DEGRAD = "../data/Primary_DegradData"


degrad_data = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArCF4.csv"))
garfield_data = pd.read_csv(os.path.join(populations_dir, "ArCF4_secondary.csv"))
garfield_data["concentration"] = garfield_data["concentration"] / 100.0

parameter_data_og = pd.read_csv(os.path.join(DATA_DIR_PAR, "ArCF4_primary.csv"))["parameter"].to_numpy()
parameter_data_og[0] = 1
parameter_data = parameter_data_og.copy()

print("="*30)
print("parameter_data original:")
print(parameter_data)
print("="*30)


levels = pd.read_csv(os.path.join(DATA_DIR_LEV, "ArCF4_level_data.csv"))

CF3_energy = levels[levels["state_name"].str.contains("NEUTRAL", na=False)].copy()
Ar_energy  = levels[levels["state_name"].str.contains("EXC", na=False)].copy()

CF3_energy_list = CF3_energy["energy_eV"].to_numpy()
fracCF3 = np.zeros_like(CF3_energy_list)
Ar_energy_list  = Ar_energy["energy_eV"].to_numpy() 
fracAr = np.zeros_like(Ar_energy_list)



# ============================================================
# 7) MALLA DE CONCENTRACIONES Y CAMPOS
# ============================================================

fN2 = np.logspace(-3, 0, 1000)
electric_fields = [40, 60, 70, 80, 90, 110]

probabilities = np.linspace(max(parameter_data_og[1],parameter_data_og[2]),1,100)
#probabilities = np.linspace(0.28,0.7,10)

chi2 = np.zeros_like(probabilities)

electric_fields = [40, 60, 70, 80, 90, 110]
#electric_fields = [60]#40, 60, 70,]# 80, 90, 110]

standard_concentration = 0.10
normalization = "ni"

npe = 500 # numero de electrones primarios, importante

cf3_ref_value = degrad_data.loc[
        degrad_data["concentration"] == standard_concentration, "CF3"
    ].iloc[0]

ar_ref_value = degrad_data.loc[
        degrad_data["concentration"] == standard_concentration, "Ar_dbleStar"
    ].iloc[0]

CF3_completed = pd.DataFrame()
Ar_completed = pd.DataFrame()
    
    
for l,E1 in enumerate(CF3_energy_list[:]):

            dataframe = pd.DataFrame(
                {

                    "CF4":    [["ION CF3 +"],                            "CF4",      0, 100, "CF4"],

                    "Ar**":   [["EXC"],                                  "ARGON",    0, 100, "Ar_dbleStar"],

                    "CF3":    [["NEUTRAL DISS"],                         "CF4",      E1, 100, "CF3"],

                    "Ar3rd":  [["CHARGE STATE =2"],      "ARGON",    0, 100, "Ar_3rd"]
                    
                }, 
                index=["name principal", "gas", "energy low", "energy up", "name output"]
            )



            read_degrad(archivo_entrada, archivo_salida_1, archivo_salida_2, gas1, gas2, concentration, dataframe, output_dir, output_general_name)
        
            degrad_data_CF4  = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArCF4.csv"))

            fracCF3[l] = degrad_data_CF4.loc[
                degrad_data_CF4["concentration"] == standard_concentration, "CF3"
            ].iloc[0] / cf3_ref_value


for l,E1 in enumerate(Ar_energy_list[:]):
            dataframe = pd.DataFrame(
                {

                    "CF4":    [["ION CF3 +"],                            "CF4",      0, 100, "CF4"],

                    "Ar**":   [["EXC"],                                  "ARGON",    E1, 100, "Ar_dbleStar"],

                    "CF3":    [["NEUTRAL DISS"],                         "CF4",      0, 100, "CF3"],

                    "Ar3rd":  [["CHARGE STATE =2"],      "ARGON",    0, 100, "Ar_3rd"]
                    
                }, 
                index=["name principal", "gas", "energy low", "energy up", "name output"]
            )


            read_degrad(archivo_entrada, archivo_salida_1, archivo_salida_2, gas1, gas2, concentration, dataframe, output_dir, output_general_name)

            degrad_data_Ar  = pd.read_csv(os.path.join(DATA_DIR_DEGRAD, "ArCF4.csv"))

            fracAr[l] = degrad_data_Ar.loc[
                degrad_data_Ar["concentration"] == standard_concentration, "Ar_dbleStar"
            ].iloc[0] / ar_ref_value

for i, electric_field in enumerate(electric_fields):
    plt.figure(figsize=(6,4))
    for j,prob in enumerate(probabilities): 
            

            parameter_data[1] = prob # max(parameter_data_og[1],parameter_data_og[2])/ # 0.35 # 
            parameter_data[2] = prob # max(parameter_data_og[1],parameter_data_og[2])/ # 0.35 # 
            fracCF3_value = parameter_data_og[1]/parameter_data[1]
            fracAR_value = parameter_data_og[2]/parameter_data[2]

            idx = np.abs(fracCF3 - fracCF3_value).argmin()
            Ecf4 = CF3_energy_list[idx] # 15.6 #    
            idx = np.abs(fracAr - fracAR_value).argmin()
            Ear =  Ar_energy_list[idx] #11.8 #  
            print("==="*10)
            print(parameter_data[1],parameter_data[2])
            print(fracCF3_value,fracAR_value)
            print(parameter_data_og[1],parameter_data_og[2])
            print(Ecf4,Ear)
            print("==="*10)


            # parameter_data[1] = 0.35 #  max(parameter_data_og[1],parameter_data_og[2])/prob # 
            # parameter_data[2] = 0.35 #  max(parameter_data_og[1],parameter_data_og[2])/prob # 

            ### garfield
            config = pd.DataFrame({
                "CF4": {
                    "name principal": "ION",
                    "gas": "CF4",
                    "energy low": 14,
                    "energy up": 20,
                    "name output": "CF4",
                    "type": "ionisation"
                },
                "Ar**": {
                    "name principal": "EXC",
                    "gas": "Ar",
                    "energy low": Ear*0.99,
                    "energy up": 100,
                    "name output": "Ar_dbleStar",
                    "type": "excitation"
                },
                "CF3": {
                    "name principal": "NEUTRAL DISS",
                    "gas": "CF4",
                    "energy low": Ecf4*0.99,
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


            garfield_norm_ne = read_garfield_csv_folder(
                folder_path=csv_folder,
                dataframe=config,
                output_dir=populations_dir,
                output_general_name=os.path.join(populations_dir, "ArCF4_secondary"),
                gas_concentration="cf4",
                gain_summary=summary,
                normalized=normalization
            )


            garfield_data = pd.read_csv(os.path.join(populations_dir, "ArCF4_secondary.csv"))
            garfield_data = garfield_data[garfield_data["electric_field"] == electric_field].copy()
            garfield_data["concentration"] = garfield_data["concentration"] / 100.0


            ### chi 2 + grafica

            press = np.array([1,1,1,1]) #bar
            con = np.array([0.05,0.10,0.67,1]) #%
            phe = np.array([0.38, 0.39, 0.28, 0.1]) #ph/e-               
    

            _fCF4 = np.logspace(-3,0,100)

            plt.xscale("log")
            plt.xlim(4,105)
            plt.plot(_fCF4*100,theory_yield_vis(parameter_data,garfield_data,_fCF4,1) * (ion_potential(_fCF4) / npe)
                     ,label=f"{parameter_data[1]:.2f};{Ecf4:.2f},{Ear:.2f}")
            
            plt.errorbar(con*100, phe, yerr = phe*0.25, fmt=".")
            #plt.legend(loc="lower left", ncol=2)

            #mask = garfield_data["concentration"].isin(con)
            #garfield_data = garfield_data.loc[mask, "concentration"]           
            

            yy = theory_yield_vis(parameter_data, garfield_data, con, 1) * (ion_potential(con) / npe)

            phe = np.asarray(phe).ravel()
            yy = np.asarray(yy).ravel()
            print(phe-yy)

            chi2[j] = np.sum(((phe - yy) / (phe * 0.22))**2)
    
    plt.show()
    plt.figure()
    plt.plot(probabilities,chi2)    
    plt.yscale("log")     
    plt.show()   