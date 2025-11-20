#!/bin/bash

# Activar modo estricto: si algún comando falla, se detiene el script
set -e

# Muestra el directorio actual
echo "Ejecutando scripts en $(pwd)"
echo "-----------------------------------"

# Aquí defines el orden de ejecución manualmente
# (puedes comentar o descomentar según necesites)
cd DegradLineal/
#python3 Fit_Degrad.py
python3 Fit_Division_DegradLineal.py
python3 Make_Plot_Cociente_Normalizacion_VisUV.py
python3 Make_Plot_Cociente_Normalizacion_100CF4.py
cd ..
echo "-----------------------------------"
echo "Ejecución completada correctamente ✅"
