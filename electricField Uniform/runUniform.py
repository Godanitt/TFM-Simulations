import subprocess

def run_fatGemC(args, dir_output="build"):
    # Asegúrate de que el binario esté actualizado
    # Ejecuta ./fatGem con los argumentos pasados
    subprocess.run(["./uniformE"] + args, cwd=dir_output)

###############################
# Ejemplo de parámetros
n = 1                         # Este es para limitar los electrones que se siu
npe = [10]                    # numero de electrones primarios que se lanzan
pressure = [1]                # bar
gas1 = ["ar"]                 # gas2 nombre -> debe ser el código para Magboltz
mixture1 = [99.0]             # % del gas2
gas2 = ["cf4"]                # gas2 nombre -> debe ser el código para Magboltz
mixture2 = [1.0]              # % del gas2
fieldE = [15000]              # eV/cm
 
 

subprocess.run(["rm", "-rf", "build/"])
dir_output="build"
subprocess.run(["mkdir", "build"])
subprocess.run(["cmake", ".."], cwd=dir_output)
subprocess.run("make -j$Nproc", shell=True, cwd=dir_output)
for i in range(n):
    rootFileName = f"../rootArchives/{gas1[i]}{mixture2[i]:.1f}{gas2[i]}_{pressure[i]:.1f}bar.root"

    # Construimos la lista de argumentos exactamente como espera fatGem:
    args = [
        rootFileName,             # name del .root
        fieldE[i],                # campo eléctrico
        "1",                      # pitch en mm
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
