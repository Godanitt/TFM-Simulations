import subprocess

def run_fatGemC(args, dir_output="build"):
    # Asegúrate de que el binario esté actualizado
    subprocess.run(["cmake", ".."], cwd=dir_output)
    subprocess.run(["make"], cwd=dir_output)
    # Ejecuta ./fatGem con los argumentos pasados
    subprocess.run(["./fatGem"] + args, cwd=dir_output)


# Ejemplo de parámetros
n = 1
npe = [1, 10, 10, 10]
pressure = [1, 0.1, 1, 0.1]
gas1 = ["xe", "xe", "ar", "ar"]
mixture1 = [100.0, 100.0, 100.0, 100.0]
gas2 = ["cf4", "cf4", "cf4", "cf4"]
mixture2 = [0.0, 0.0, 0.0, 0.0]

for i in range(n):
    rootFileName = f"../rootArchives/{gas1[i]}{mixture2[i]:.1f}{gas2[i]}_{pressure[i]:.1f}bar_{npe[i]:d}npe.root"

    # Construimos la lista de argumentos exactamente como espera fatGem:
    args = [
        "../fieldMaps2bar",       # campo eléctrico
        rootFileName,             # name del .root
        pressure[i],              # presion en bar
        "5",                      # pitch en mm
        npe[i],                   # numero de electrones primarios
        gas1[i],                  # gas primario
        f"{mixture1[i]:.1f}",     # porcentaje gas1
        gas2[i],                  # gas secundario
        f"{mixture2[i]:.1f}"      # porcentaje gas2
    ]
    args_str = [str(a) for a in args]
    print("-"*40)
    print(f"Ejecutando \n cmake .. \n make \n ./fatGem {' '.join(args_str)}\n")
    print("-"*40)
    run_fatGemC(args_str)
