import subprocess
import os

def run_fatGemC(args, dir_output="build"):
    # Asegúrate de que el binario esté actualizado
    # Ejecuta ./fatGem con los argumentos pasados
    subprocess.run(["./uniformE"] + args, cwd=dir_output)


##############################
#Ejemplo de parámetros
n = 8                         

# numero de electrones primarios que se lanzan
npe = [10000] * n                    

# bar
pressure = [1, 1, 0.025, 0.025,
            1, 1, 0.025, 0.025,]       
         
# mm
gap = [0.57] * n                  

# gas1 nombre -> debe ser el código para Magboltz
gas1 = ["cf4", "cf4", "cf4", "cf4", 
        "ar", "ar", "ar", "ar",  ]   
                                
# % del gas1
mixture1 = [100.0, 100.0, 100.0, 100.0, 
            80.0, 80.0, 80.0, 80.0,  ]             

# gas2 nombre -> debe ser el código para Magboltz
gas2 = ["ar", "ar", "ar", "ar", "ar", "ar",  
        "cf4", "cf4", "cf4", "cf4", "cf4", "cf4",]       

# % del gas2        
mixture2 = [0.0, 0.0, 0.0, 0.0, 
            20.0, 20.0, 20.0, 20.0, ]              

# eV/cm
fieldE = [34400, 43000, 8000, 10000,
          23200, 29000, 5280, 6600]              
 
# n=1
# npe=[100]
# pressure=[0.025]
# gap=[0.57]
# gas1=["ar"]
# mixture1=[0.0]
# gas2=["cf4"]
# mixture2=[100.0]
# fieldE=[12000]



###############################
 

subprocess.run(["rm", "-rf", "build/"])
dir_output="build"
subprocess.run(["mkdir", "build"])
subprocess.run(["cmake", ".."], cwd=dir_output)
subprocess.run("make -j$Nproc", shell=True, cwd=dir_output)
os.makedirs("rootArchives", exist_ok=True)                      # Crea La carpeta de root si no existe

for i in range(n):
    rootFileName = f"../rootArchives/{gas1[i]}{mixture2[i]:.1f}{gas2[i]}_{fieldE[i]/pressure[i]/1000:.1f}kVcmbar_{gap[i]:.2f}cm_{npe[i]}npe.root"

    # Construimos la lista de argumentos exactamente como espera fatGem:
    args = [
        rootFileName,             # name del .root
        fieldE[i],                # campo eléctrico
        "%.2f"%gap[i],            # gap en mm
        pressure[i],              # presion en bar
        npe[i],                   # numero de electrones primarios
        gas1[i],                  # gas primario
        f"{mixture1[i]:.1f}",     # porcentaje gas1
        gas2[i],                  # gas secundario
        f"{mixture2[i]:.1f}"      # porcentaje gas2
    ]
    args_str = [str(a) for a in args]
    print("-"*40)
    print(f"Ejecutando \n cmake .. \n make \n ./uniformE {' '.join(args_str)}\n")
    print("-"*40)
    run_fatGemC(args_str)
