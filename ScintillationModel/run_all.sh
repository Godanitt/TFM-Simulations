#!/bin/bash
set -e

echo "==============================="
echo " Ejecutando Read_Degrad.py"
echo "==============================="
python3 Read_Degrad.py

echo "==============================="
echo " Ejecutando Read_Experimental_Yield.py"
echo "==============================="
python3 Read_Experimental_Yield.py

echo "==============================="
echo " Ejecutando main.py"
echo "==============================="
cd NewModel/
python3 main.py
cd ..
echo "==============================="
echo "   ✔ Todos los programas ejecutados correctamente"
echo "==============================="
