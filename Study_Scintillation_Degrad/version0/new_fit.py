from amoedo_refined_model_garbagecCollector import Pgamma_UV_refined_GC,Pgamma_vis_refined_GC
from amoedo_model import Pgamma_Ar3rd,Pgamma_CF3,Pgamma_CF4
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.optimize as opt
import os 

# === Carpeta donde están los pickles ===
DATA_DIR = "pickle_data"

# ============================================
# === 1. Cargar los pickles ===
pCF4 = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
pCF3 = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_CF3.pkl"))
pAr_dbleStar  = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_Ar_dbleStar.pkl"))
pAr_3rd  = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_Ar_3rd.pkl"))

scCF3 = pd.read_pickle(os.path.join(DATA_DIR, "scintillationProbability_PAModel_CF3.pkl"))
scCF4 = pd.read_pickle(os.path.join(DATA_DIR, "scintillationProbability_PAModel_CF4.pkl"))
scAr  = pd.read_pickle(os.path.join(DATA_DIR, "scintillationProbability_PAModel_Ar3rd.pkl"))

yield_uv  = pd.read_pickle(os.path.join(DATA_DIR, "yield_uv.pkl"))
yield_vis = pd.read_pickle(os.path.join(DATA_DIR, "yield_vis.pkl"))

presiones = pd.read_pickle(os.path.join(DATA_DIR, "presiones.pkl"))


print("✅ Pickles cargados correctamente desde 'pickle_data/'.\n")

# ============================================
# === 2. Verificaciones ===
if "fCF4" not in scCF3.columns:
    raise ValueError("⚠️ Falta columna 'fCF4' en CF3.")
if "all" not in scCF4.columns or "all" not in scAr.columns:
    raise ValueError("⚠️ CF4 o Ar3rd no contienen la columna 'all'.")
if "fCF4" not in yield_uv.columns or "fCF4" not in yield_vis.columns:
    raise ValueError("⚠️ Falta columna 'fCF4' en yield_uv o yield_vis.")

# ==================================================
# === 3. Generar nombres para correr los bucles  ===

name_CF4            =  pCF4.columns.to_numpy()
name_CF3            =  pCF3.columns.to_numpy()
name_Ar_dbleStar    =  pAr_dbleStar.columns.to_numpy()
name_Ar_3rd         =  pAr_3rd.columns.to_numpy()

name_SP_CF3         =  scCF3.columns.to_numpy()
name_SP_CF4         =  scCF4.columns.to_numpy()

name_yield_uv       =  yield_uv.columns.to_numpy()
name_yield_vis      =  yield_vis.columns.to_numpy()



#least_squares(fun, x0, loss='soft_l1', f_scale=0.1,args=(t_train, y_train))

# ============================================fCF4
# === 4. Funciones para minimizar  ===

fCF4 = pCF4["fCF4"].to_numpy()*100
x0 = np.linspace(0,len(fCF4)-1,len(fCF4))
x0 = [int(x) for x in x0]

def fun_to_min_vis(x, name):     
    autoQ, C = x  # parámetros libres

    # Selección de datos según el nombre
    for j in range(len(name_CF3)):
        if name_CF3[j] in name: 
            a = pCF3[name_CF3[j]].to_numpy()
    for k in range(len(name_Ar_dbleStar)):
        if name_Ar_dbleStar[k] in name: 
            b = pAr_dbleStar[name_Ar_dbleStar[k]].to_numpy()
    
    # Extrae las filas del dataframe correspondientes a los valores de fCF4
    residuos = []
    for i in range(0,len(fCF4)):
        c = yield_vis.loc[yield_vis['fCF4'] == fCF4[i]].to_numpy()
        for j in c[0,1:]:
            if j!=0:
                val_obs = j
                val_model = Pgamma_vis_refined_GC(fCF4[i]/100, a[i], b[i], autoQ, C)
                residuos.append(val_obs - val_model)
    
                #print("Yvis")
                #print(j)
                #print("PgammaCF3")
                #print(Pgamma_CF3(fCF4[i]/100, a[i], b[i]),"\n \n")
    return np.array(residuos)


def fun_to_min_uv(x, name):     
    kcool, C = x  # parámetros libres

    # Selección de datos según el nombre
    for j in range(len(name_CF4)):
        if name_CF4[j] in name: 
            a = pCF4[name_CF4[j]].to_numpy()
    for k in range(len(name_Ar_3rd)):
        if name_Ar_3rd[k] in name: 
            b = pAr_3rd[name_Ar_3rd[k]].to_numpy()
    
    # Extrae las filas del dataframe correspondientes a los valores de fCF4
    residuos = []
    for i in range(0,len(fCF4)):
        presion = yield_uv.columns.to_numpy()
        c = yield_uv.loc[yield_uv['fCF4'] == fCF4[i]].to_numpy()
        m=0
        for j in c[0,1:]:
            m+=1
            if j!=0:
                val_obs = j
                    
                n = float(presion[m].replace("bar", ""))
                print(n)                

                val_model = Pgamma_UV_refined_GC(fCF4[i]/100, a[i], b[i], n, kcool, C)
                residuos.append(val_obs - val_model)
    
    return np.array(residuos)



# ============================================
# === 5. Minimizacion  ===


fit_vis = pd.DataFrame({"Parametro":["autoQ","C"]})
for col in name_SP_CF3:
    fit_vis[col] = np.nan
fit_vis = fit_vis.set_index("Parametro")

fit_uv = pd.DataFrame({"Parametro":["kcool","C"]})
for col in name_SP_CF4:
    fit_uv[col] = np.nan
fit_uv = fit_uv.set_index("Parametro")

for name in name_SP_CF3[1:]:
    x0 = [1.0, 10/100]
    res = opt.least_squares(fun_to_min_vis, x0, args=(name,), bounds=([0, 1/1000], [1, 50/100]))
    fit_vis.loc["autoQ", name] = res.x[0]
    fit_vis.loc["C", name]     = res.x[1]

for name in name_SP_CF4[1:]:
    x0 = [100.0, 10/100]
    res = opt.least_squares(fun_to_min_uv, x0, args=(name,), bounds=([0.5, 0], [1000,1]))
    fit_uv.loc["kcool", name] = res.x[0]
    fit_uv.loc["C", name]     = res.x[1]

print(fit_vis)
print(fit_uv)

# ============================================
# === 6. Creacion de los fits  ===

Pgamma_vis_fit = scCF3.copy()
for i in name_SP_CF3:
    #a = np.zeros_like(name_CF3)
    #b = np.zeros_like(name_Ar_dbleStar)
    for j in range(len(name_CF3)):
        if name_CF3[j] in i: 
            a = pCF3[name_CF3[j]].to_numpy()
            
        for k in range(len(name_Ar_dbleStar)):
            if name_Ar_dbleStar[k] in i: 
                b = pAr_dbleStar[name_Ar_dbleStar[k]].to_numpy()

        Pgamma_vis_fit[i]=Pgamma_vis_refined_GC(fCF4/100,a,b,fit_vis.loc["autoQ", i] ,fit_vis.loc["C", i])

Pgamma_uv_fit = yield_uv[yield_uv['fCF4'].isin(fCF4)]

presion = presiones["Presion"].to_numpy()
presion = presion[1:]          
for i in name_SP_CF4:                   
    #a = np.zeros_like(name_CF3)    
    #b = np.zeros_like(name_Ar_dbleStar)
    for j in range(len(name_CF4)):
        if name_CF4[j] in i: 
            a = pCF4[name_CF4[j]].to_numpy()
            
        for k in range(len(name_Ar_3rd)):
            if name_Ar_3rd[k] in i: 
                b = pAr_3rd[name_Ar_3rd[k]].to_numpy()

        for l in range(len(presion)):
            Pgamma_uv_fit[name_yield_uv[l+1]]=Pgamma_UV_refined_GC(fCF4/100,a,b,presion[l],fit_uv.loc["kcool", i] ,fit_uv.loc["C", i])

    
                                
print(Pgamma_uv_fit)

# ============================================
# === 7. Ploteo de VISIBLE  ===

plt.figure()

for col in Pgamma_vis_fit.columns:
    if col == "fCF4":
        continue
    plt.plot(fCF4, Pgamma_vis_fit[col], marker=',', linestyle='-', label=f"Scint {col}")

for col in yield_vis.columns:
    if col != "fCF4" and col in yield_vis.columns:
        plt.plot(yield_vis["fCF4"], yield_vis[col], marker='s', linestyle='--',linewidth=2, label=f"Yield Vis {col}")

plt.xlabel("fCF4", fontsize=12)
plt.ylabel("Yield Visible", fontsize=12)
plt.xscale("log")  # útil si fCF4 varía logarítmicamente
plt.yscale("log")  # útil si fCF4 varía logarítmicamente
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.ylim(0.5,50)
plt.legend(title="Curvas", fontsize=8)
plt.tight_layout()
plt.savefig("output/vis_fit_data.pdf")


# ============================================
# === 8. Ploteo de ULTRAVIOLETA  ===

plt.figure()

for col in Pgamma_uv_fit.columns:
    if col == "fCF4":
        continue
    plt.plot(fCF4, Pgamma_uv_fit[col], marker=',', linestyle='-', label=f"Scint {col}")

for col in yield_uv.columns:
    if col != "fCF4" and col in yield_uv.columns:
        plt.plot(yield_uv["fCF4"], yield_uv[col], marker='s', linestyle='--',linewidth=2, label=f"Yield Vis {col}")

plt.xlabel("fCF4", fontsize=12)
plt.ylabel("Yield UltraViloeta", fontsize=12)
plt.xscale("log")  # útil si fCF4 varía logarítmicamente
plt.yscale("log")  # útil si fCF4 varía logarítmicamente
plt.grid(True, which="both", ls="--", alpha=0.5)
#plt.ylim(5,200)
plt.legend(title="Curvas", fontsize=8)
plt.tight_layout()
plt.savefig("output/uv_fit_data.pdf")

