from read_Degrad_input import read_input
from amoedo_model import Pgamma_Ar3rd,Pgamma_CF3,Pgamma_CF4
import numpy as np 
import pandas as pd 
import dill
import scipy.special
import importlib


archives=np.array(["input/output_99.9Ar_0.1CF4.txt",
                   "input/output_99Ar_1CF4.txt",
                   "input/output_90Ar10CF4.txt",
                   "input/output_PureCF4.txt"])

archives_ar=np.array(["input/ar_input_90Ar10CF4.csv",
                   "input/ar_input_99.9Ar_0.1CF4.csv",
                   "input/ar_input_99Ar_1CF4.csv",
                   "input/ar_input_PureCF4.csv"])

archives_cf4=np.array(["input/cf4_input_90Ar10CF4.csv",
                   "input/cf4_input_99.9Ar_0.1CF4.csv",
                   "input/cf4_input_99Ar_1CF4.csv",
                   "input/cf4_input_PureCF4.csv"])

# ------------------------------------------------------------------------------
# Arrays para guardar la información 

#ar_eventos      = np.array([0,0,0,0],dtype=object)
#cf4_eventos     = np.array([0,0,0,0],dtype=object)

#ar_proceso      = np.array([0,0,0,0],dtype=object)
#cf4_proceso     = np.array([0,0,0,0],dtype=object)

#ar_error        = np.array([0,0,0,0],dtype=object)
#cf4_error       = np.array([0,0,0,0],dtype=object)


# ------------------------------------------------------------------------------
# Valores por archivo (concentraciones, temperaturas, presión)

fCF4            = np.array([0.001,0.01,0.1,1])   # concentraciones
n               = 1                              # P/P0
#pressure       = 760                            # torr
#temperature    = 20                             # grados centigrados

# ------------------------------------------------------------------------------
# Creamos los dataframes que contendrán toda la información

poblations_CF4 = pd.DataFrame({"fCF4":fCF4})
poblations_CF3 = pd.DataFrame({"fCF4":fCF4})
poblations_Ar_dbleStar  = pd.DataFrame({"fCF4":fCF4})
poblations_Ar_3rd  = pd.DataFrame({"fCF4":fCF4})


# PA = PabloAmoedo

scintillationProbability_PAModel_CF3 = pd.DataFrame({"fCF4":fCF4})
scintillationProbability_PAModel_CF4 = pd.DataFrame({"fCF4":fCF4})
scintillationProbability_PAModel_Ar3rd  = pd.DataFrame({"fCF4":fCF4})

yield_uv       = pd.DataFrame()
yield_vis      = pd.DataFrame()

# ------------------------------------------------------------------------------
# Creamos los nombres de las columnas 
# Cambiar el nombre aquí va a dar error, también hay que cambiarlo cuando se indexa.
# Además las indexaciones se hacen a mano (lineas 120-140)
name_CF4         =  np.array(["all"])
name_CF3         =  np.array([">11.5",">12.5",">14.5"])
name_Ar_dbleStar =  np.array(["all",">2o Resonante"])
name_Ar_3rd      =  np.array(["all"])
name_uv          =  np.array(["fCF4","1bar","2bar","2.5bar","3bar","4bar","5bar"])
name_vis         =  np.array(["fCF4","1bar","2bar","2.5bar","3bar","4bar","5bar"])

presion = np.array([1,1,2,2.5,3,4,5])
presiones = pd.DataFrame({"Presion":presion})

name_SP_CF3      =  np.array([">11.5 all",">11.5 >2o Resonante",
                              ">12.5 all",">12.5 >2o Resonante",
                              ">14.5 all",">14.5 >2o Resonante"])

name_SP_CF4      =  np.array(["all"])
name_SP_Ar3rd    =  np.array(["all"])

# ------------------------------------------------------------------------------
# Inicializamos los nombres de las columnas de los dataFrames

for col in name_CF4:
    poblations_CF4[col] = np.nan

for col in name_CF3:
    poblations_CF3[col] = np.nan
    
for col in name_Ar_dbleStar:
    poblations_Ar_dbleStar[col] = np.nan
    
for col in name_Ar_3rd:
    poblations_Ar_3rd[col] = np.nan
    
for col in name_uv:
    yield_uv[col] = [np.nan]*9
    
for col in name_vis:
    yield_vis[col] = [np.nan]*9

for col in name_SP_CF3:
    scintillationProbability_PAModel_CF3[col] = np.nan
    
for col in name_SP_CF4:
    scintillationProbability_PAModel_CF4[col] = np.nan
    
for col in name_SP_Ar3rd:
    scintillationProbability_PAModel_Ar3rd[col] = np.nan
    
# ------------------------------------------------------------------------------
# Bucle que recorra los diferentes puntos: 

for i in range(len(archives)):
    df = read_input(
        archives[i],
        archivo_salida_Ar=archives_ar[i],      
        archivo_salida_CF4=archives_cf4[i])
    
    ar  = df.loc[lambda df: df['Gas'] == "ARGON", :]
    cf4 = df.loc[lambda df: df['Gas'] == "CF4", :]
    
    poblations_CF4.loc[i, "all"] = cf4.loc[cf4['Proceso'].str.contains("ION CF3 +"), 'Eventos'].sum()
    poblations_Ar_3rd.loc[i, "all"] = ar.loc[ar['Proceso'].str.contains("CHARGE STATE =2"), 'Eventos'].sum()

    poblations_CF3.loc[i, ">11.5"] = cf4.loc[cf4['Proceso'].str.contains("NEUTRAL DISS")  & (cf4['Energia'] > 11.5), 'Eventos'].sum()
    poblations_CF3.loc[i, ">12.5"] = cf4.loc[cf4['Proceso'].str.contains("NEUTRAL DISS")  & (cf4['Energia'] > 12.5), 'Eventos'].sum()
    poblations_CF3.loc[i, ">14.5"] = cf4.loc[cf4['Proceso'].str.contains("NEUTRAL DISS")  & (cf4['Energia'] > 14.5), 'Eventos'].sum()
    
    poblations_Ar_dbleStar.loc[i, "all"]=ar.loc[(ar['Proceso'].str.contains("EXC")),'Eventos'].iloc[:].sum()
    poblations_Ar_dbleStar.loc[i, ">2o Resonante"]=ar.loc[(ar['Proceso'].str.contains("EXC")) & (ar['Energia'] > 11.8),'Eventos'].iloc[:].sum()    
    
print("-"*60)
print(poblations_CF3)
print("-"*60)
print(poblations_Ar_dbleStar)
print("-"*60)
print(poblations_CF4)
print("-"*60)
print(poblations_Ar_3rd)
print("-"*60)
# ----------------------------------------------------------------
# no se que demonios pasaba aqui, chatgpt es el puto amo
# Cargar el módulo compilado de bajo nivel
_special_ufuncs = importlib.import_module("scipy.special._special_ufuncs")

# Lista de funciones que pueden faltar
funcs = ["erf", "erfc", "erfi", "gamma", "lgamma","wofz"]

# Inyectarlas si no existen
for name in funcs:
    if not hasattr(_special_ufuncs, name) and hasattr(scipy.special, name):
        setattr(_special_ufuncs, name, getattr(scipy.special, name))
        #print(f"🔧 Añadida función faltante: {name}")

with open("input/resultados_4picos.pkl", "rb") as f:
    df = dill.load(f)
    
    
# ----------------------------------------------------------------
# concentraciones objetivo (de la primera presión, por ejemplo)
df_pressure0 = df[df["presion"] == presion[0]].copy()
concentraciones_og = df_pressure0["concentracion"].to_numpy()

# Inicializa las tablas destino alineadas por fCF4
yield_uv  = pd.DataFrame({"fCF4": concentraciones_og})
yield_vis = pd.DataFrame({"fCF4": concentraciones_og})

for i in range(1,len(name_uv)):
    # OJO: i y presion alineados 1:1
    df_pressure = df[df["presion"] == presion[i]].copy()

    # Extrae UV y vis como Series indexadas por concentracion
    s_uv = df_pressure.set_index("concentracion")["yields_zonas"].apply(lambda d: d["UV"])
    s_vis = df_pressure.set_index("concentracion")["yields_zonas"].apply(lambda d: d["vis"])

    # Si las concentraciones son floats con posibles decimales “sucios”, puedes redondear:
    # s_uv.index  = np.round(s_uv.index.astype(float), 8)
    # s_vis.index = np.round(s_vis.index.astype(float), 8)
    # conc_idx    = np.round(concentraciones_og.astype(float), 8)
    # Luego reindexar con conc_idx

    # Reindexa para alinear por concentración objetivo
    col_uv  = name_uv[i]
    col_vis = name_vis[i]

    yield_uv[col_uv]  = pd.Series(s_uv).reindex(concentraciones_og).to_numpy()
    yield_vis[col_vis]= pd.Series(s_vis).reindex(concentraciones_og).to_numpy()

# Rellena faltantes con -1 si quieres ese marcador
yield_uv  = yield_uv.fillna(0)
yield_vis = yield_vis.fillna(0)
print("=="*50)
print(yield_uv)
    
# ------------------------------------------
# Generamos los dataframes de scintillation

for i in name_SP_CF3:
    #a = np.zeros_like(name_CF3)
    #b = np.zeros_like(name_Ar_dbleStar)
    for j in range(len(name_CF3)):
        if name_CF3[j] in i: 
            a = poblations_CF3[name_CF3[j]].to_numpy()
            
        for k in range(len(name_Ar_dbleStar)):
            if name_Ar_dbleStar[k] in i: 
                b = poblations_Ar_dbleStar[name_Ar_dbleStar[k]].to_numpy()

        #print(a)
        #print(b)
        scintillationProbability_PAModel_CF3[i]=Pgamma_CF3(fCF4,a,b)

for i in name_SP_CF4:
    for j in range(len(name_CF4)):
        if name_CF4[j] in i: 
            a = poblations_CF4[name_CF4[j]].to_numpy()
        for k in range(len(name_Ar_3rd)):
            if name_Ar_3rd[k] in i: 
                b = poblations_Ar_3rd[name_Ar_3rd[k]].to_numpy()
        scintillationProbability_PAModel_CF4[i]=Pgamma_CF4(fCF4,a,b,n)
    
for i in name_SP_Ar3rd:
    for j in range(len(name_Ar_3rd)):
        if name_Ar_3rd[j] in i: 
            a = poblations_Ar_3rd[name_Ar_3rd[j]].to_numpy()
        
        scintillationProbability_PAModel_Ar3rd[i]=Pgamma_Ar3rd(fCF4,a,n)  
    
# ------------------------------------------
# Guardamos los dataframes en pickles para su posterior análisis


# --- Guardado automático en pickle ---
# Diccionario con todas las variables que quieres guardar
to_save = {
    "poblations_CF4": poblations_CF4,
    "poblations_CF3": poblations_CF3,
    "poblations_Ar_dbleStar": poblations_Ar_dbleStar,
    "poblations_Ar_3rd": poblations_Ar_3rd,
    "scintillationProbability_PAModel_CF3": scintillationProbability_PAModel_CF3,
    "scintillationProbability_PAModel_CF4": scintillationProbability_PAModel_CF4,
    "scintillationProbability_PAModel_Ar3rd": scintillationProbability_PAModel_Ar3rd,
    "yield_uv": yield_uv,
    "yield_vis": yield_vis,
    "presiones": presiones
}

# Recorremos y guardamos cada uno como .pkl
for name, df in to_save.items():
    df.to_pickle(f"pickle_data/{name}.pkl")
    print(f"✅ Guardado: {name}.pkl")