#!/bin/bash

# Activar modo estricto: si algún comando falla, se detiene el script
set -e

# Muestra el directorio actual
echo "Ejecutando scripts en $(pwd)"
echo "-----------------------------------"

# Aquí defines el orden de ejecución manualmente
# (puedes comentar o descomentar según necesites)
python3 Read_Degrad.py
python3 Read_Experimental_Yield.py
python3 Fit_Degrad.py
python3 Make_Poblations_with_fit.py
python3 Fit_Division_DegradData.py
python3 Fit_Division_DegradFit.py
python3 Make_fits_plots_DegradData.py
python3 Make_fits_plots_DegradFit.py
echo "-----------------------------------"
echo "Ejecución completada correctamente ✅"
