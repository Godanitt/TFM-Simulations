#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 13:32:27 2025

@author: pablo
"""
import subprocess
import os
from multiprocessing import Pool, cpu_count

def run_fatGemC(args):
    """Ejecuta ./uniformE con la lista de argumentos en el directorio build."""
    dir_output = "build"
    print(f"--> Lanzando simulación:\n    ./uniformE {' '.join(args)}")
    subprocess.run(["./uniformE"] + args, cwd=dir_output)
    print(f"--> Finalizó simulación:\n    ./uniformE {' '.join(args)}")


###################################
# Parámetros del usuario
n = 15

npe       = [10000] * n

pressure  = [1, 1, 1,
             0.025, 0.025, 0.025,
             1, 1, 1,
             0.025, 0.025, 0.025,
             1, 1, 1,]

gap       = [0.57] * n

gas1      = ["cf4","cf4","cf4",
             "cf4","cf4","cf4",
             "ar","ar","ar",
             "ar","ar","ar",
             "he","he","he"]

mixture1  = [100.0]*6 + [80.0]*9

gas2      = ["ar","ar","ar",
             "ar","ar","ar",
             "cf4","cf4","cf4",
             "cf4","cf4","cf4",
             "cf4","cf4","cf4",]

mixture2  = [0.0]*6 + [20.0]*9

fieldE    = [38700, 40850, 43000,
             9500, 9000, 10000,
             26100, 27550, 29000,
             6270, 5940, 6600,
             18700, 16830, 17765,]
###################################


# ---------------------------
# COMPILACIÓN ÚNICA
# ---------------------------

subprocess.run(["rm", "-rf", "build/"])
os.makedirs("build", exist_ok=True)
subprocess.run(["cmake", ".."], cwd="build")
subprocess.run("make -j$Nproc", shell=True, cwd="build")

os.makedirs("rootArchives", exist_ok=True)

# ---------------------------
# PREPARACIÓN DE ARGUMENTOS
# ---------------------------

all_jobs = []

for i in range(n):
    rootFileName = (
        f"../rootArchives/"
        f"{gas1[i]}{mixture2[i]:.1f}{gas2[i]}_"
        f"{fieldE[i]/1000:.1f}kVcm_"
        f"{pressure[i]}bar_"
        f"{gap[i]:.2f}cm_{npe[i]}npe.root"
    )

    args = [
        rootFileName,
        str(fieldE[i]),
        f"{gap[i]:.2f}",
        str(pressure[i]),
        str(npe[i]),
        gas1[i],
        f"{mixture1[i]:.1f}",
        gas2[i],
        f"{mixture2[i]:.1f}",
    ]

    all_jobs.append(args)

# ---------------------------
# MULTIPROCESSING
# ---------------------------

num_cores = min(cpu_count(), n)
print(f"\nUsando {num_cores} núcleos para las simulaciones...\n")

with Pool(processes=num_cores) as pool:
    pool.map(run_fatGemC, all_jobs)

print("\n✔ Todas las simulaciones han terminado.\n")
