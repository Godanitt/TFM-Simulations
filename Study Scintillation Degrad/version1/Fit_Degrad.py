
import numpy as np 
import pandas as pd 
import pickle
import scipy.optimize as opt
import os 
import dill
import matplotlib.pyplot as plt
import matplotlib.cm as cm
"""
Script que ajusta las poblaciones de Degrad a una función lineal dependiente de fCF4 o (1-fCF4)    
"""

#############################################################################################################
######################## Lectura de las poblaciones de Degrad ################################################

DATA_DIR = "pickle_data"
poblations_CF4 = pd.read_pickle(os.path.join(DATA_DIR, "poblations_CF4.pkl"))
poblations_CF3 = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_CF3.pkl"))
poblations_Ar_dbleStar  = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_Ar_dbleStar.pkl"))
poblations_Ar_3rd  = pd.read_pickle(os.path.join(DATA_DIR,  "poblations_Ar_3rd.pkl"))

poblations = [poblations_CF3,poblations_CF4,poblations_Ar_dbleStar,poblations_Ar_3rd]
names = ["CF3","CF4","Ar_dbleStar","Ar_3rd"]
################# Definimos función lineal y la incertidumbre de cada punto #################################

def fun_lineal(x,a,b):
    """ f(x) = a + b*x """
    return a + b * x

def fun_lineal_inc(x,sa,sb):
    """ f(x) =  sqrt(sa + sb*x) """
    return np.sqrt(sa**2 + (sb * x)**2)

def return_lineal_fun(a,b):
    return lambda x: a+b*x

def return_lineal_inc_fun(sa,sb):
    return lambda x: np.sqrt(sa**2 + (sb * x)**2)


#############################################################################################################

linealFun_poblations_CF4 = pd.DataFrame()
linealFun_poblations_CF3 = pd.DataFrame()
linealFun_poblations_Ar_dlbeStar = pd.DataFrame()
linealFun_poblations_Ar_3rd = pd.DataFrame()

linealFun_poblations = [linealFun_poblations_CF4,
                        linealFun_poblations_CF3,
                        linealFun_poblations_Ar_dlbeStar,
                        linealFun_poblations_Ar_3rd]

for i in range(len(poblations)): 
    poblation=poblations[i]
    linealFun=linealFun_poblations[i]
    columnas = poblation.columns
    print((len(columnas)-1))
    fCF4 = poblation["fCF4"].to_numpy()
    
    fig,ax = plt.subplots(1,1)
    
    cmap = cm.get_cmap("plasma",4)       # elige el colormap que quieras
    
    for j in range(1,(len(columnas)-1),2):
        print(columnas[j])
        print(j)
        
        if "Ar" in columnas[j]:
            fCF4 = poblation["fCF4"].to_numpy()[:-1]
            colum = poblation[columnas[j]].to_numpy()[:-1]
            scolum = poblation[columnas[j+1]].to_numpy()[:-1]          
        else:
            fCF4 = poblation["fCF4"].to_numpy()
            colum = poblation[columnas[j]].to_numpy()
            scolum = poblation[columnas[j+1]].to_numpy()
        
        
        print(fCF4,"\n",colum,"\n",scolum,"\n ============= \n")
        values,svalues = opt.curve_fit(fun_lineal,fCF4,colum,p0=[0,0],sigma=scolum)
        a,b = values
        sa,sb = np.sqrt(np.diag(svalues))
        
        linealFun[columnas[j]] = [return_lineal_fun(a,b)]
        linealFun[columnas[j+1]] = [return_lineal_inc_fun(sa,sb)]
        
        color = cmap((j-1)//2)                    # color distinto por iteración
    
        x = np.linspace(min(fCF4),max(fCF4),1000)
        y = fun_lineal(x,a,b)
        sy = fun_lineal_inc(x,sa,sb)
        ax.errorbar(fCF4,colum,scolum,fmt=".",color=color,label=columnas[j])
        ax.plot(x,y,color=color)
        ax.fill_between(x,y+sy,y-sy,color=color,alpha=0.25)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("fCF4")
        ax.set_ylabel("Poblation Degrad")
        ax.legend()
        
        
    fig.savefig("output/"+names[i]+".pdf",bbox_inches="tight",dpi=300)

    
    linealFun_poblations[i] = linealFun
    

linealFun_poblations_CF4 = linealFun_poblations[0]
linealFun_poblations_CF3 = linealFun_poblations[1]
linealFun_poblations_Ar_dlbeStar = linealFun_poblations[2]
linealFun_poblations_Ar_3rd = linealFun_poblations[3]

########################## Guardado en Pickle (con DILL) ##############################################
# Diccionario con las variables a guardar
to_save = {
    "linealFun_poblations_CF4": linealFun_poblations_CF4,
    "linealFun_poblations_CF3": linealFun_poblations_CF3,
    "linealFun_poblations_Ar_dlbeStar": linealFun_poblations_Ar_dlbeStar,
    "linealFun_poblations_Ar_3rd": linealFun_poblations_Ar_3rd,
}

# Crear carpeta si no existe
os.makedirs("pickle_data", exist_ok=True)

# Guardar cada variable como .pkl usando *dill*
for name, obj in to_save.items():
    with open(f"pickle_data/{name}.pkl", "wb") as f:
        dill.dump(obj, f)
    print(f"✅ Guardado (dill): {name}.pkl")