
import numpy as np 
import pandas as pd 
import dill
import scipy.special
import importlib

"""
Script que nos permite leer los datos de los yields de visible/ultravioleta, sacándolos en formato pickle y csv.
"""


#############################################################################################################
########################## FUNCION PARA LEER LOS PICKLES ##############################################
#############################################################################################################

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



def read_experimental(archivo_entrada, yields, presiones,output_dir, concentraciones_reales=None, no_sistematic = True):

    with open(archivo_entrada, "rb") as f:
        df = dill.load(f)


    name = ["fCF4","fCF4 real", "Err fCF4 real"]

    for i in presiones:
        name += ["%.1fbar"%i]
        name += ["Err %.1fbar"%i]
        
            
    df_pressure0 = df[df["presion"] == presiones[0]].copy()
    concentraciones = df_pressure0["concentracion"].to_numpy()


    if concentraciones_reales!= None: 
        concentraciones = concentraciones_reales


    for i in yields:

        yield_out = pd.DataFrame({"fCF4": concentraciones})

        for j in range(0,len(presiones)):

            # OJO: i y presion alineados 1:1
            df_pressure = df[df["presion"] == presiones[j]].copy()
    

            # Extrae UV y vis como Series indexadas por concentracion
            s = df_pressure.set_index("concentracion")["yields_zonas"].apply(lambda d: d[i])
            err_s  = df_pressure.set_index("concentracion")["u_yields_zonas"].apply(lambda d: d[i])
            if no_sistematic:
                err_s  = df_pressure.set_index("concentracion")["uyields_estadistico"].apply(lambda d: d[i])


            # Si las concentraciones son floats con posibles decimales “sucios”, puedes redondear:
            # s_uv.index  = np.round(s_uv.index.astype(float), 8)
            # s_vis.index = np.round(s_vis.index.astype(float), 8)
            # conc_idx    = np.round(concentraciones_og.astype(float), 8)
            # Luego reindexar con conc_idxz

            # Reindexa para alinear por concentración objetivo
            col  = name[j*2+3]
            err_col  = name[j*2+1+3]

            yield_out[col] = pd.Series(s).reindex(concentraciones).to_numpy()
            yield_out[err_col] = pd.Series(err_s).reindex(concentraciones).to_numpy()


        yield_out = yield_out.fillna(0)
        yield_out.to_csv(f"{output_dir}{i}.csv", index=False)
        print(f"✅ Guardado: {i}.csv")
        