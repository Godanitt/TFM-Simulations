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
bash lineal.sh
echo "-----------------------------------"
echo "Ejecución completada correctamente ✅"
