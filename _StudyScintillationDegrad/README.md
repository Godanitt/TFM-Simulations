# README — Descripción rápida de archivos

Este repositorio contiene los scripts utilizados para procesar simulaciones DEGRAD, leer poblaciones, ajustar modelos teóricos y comparar con datos experimentales de centelleo.

## Archivos principales
### Amoedo_Model_Degrad.py

Implementación del modelo teórico principal aplicado a los resultados de DEGRAD (poblaciones, rutas de desexcitación, términos de probabilidad).

### Amoedo_Model_DivisionFit.py

Funciones del modelo usadas específicamente para calcular cocientes visibles y UV y utilizarlos en los procesos de ajuste.

### Amoedo_Model.py

Módulos comunes del modelo de scintilación: constantes, funciones auxiliares y estructura general del cálculo.

### Fit_Degrad.py

Script que ajusta directamente los datos de DEGRAD a los modelos teóricos (poblaciones, rutas de excitación, CF₄, Ar, etc.).

### Fit_Division_DegradFit.py

Script que realiza el ajuste del cociente experimental vs. modelo, obteniendo parámetros como α, kcool y kdis.

### main.sh

Pequeño script de shell para ejecutar la cadena completa de procesamiento/ajuste de forma automática.

### Read_Degrad.py

Lectura y tratamiento de los ficheros de salida de DEGRAD: parsing, extracción de poblaciones y conversión a formatos manejables.

### Read_Experimental_Yield.py

Carga, limpieza y estructuración de los datos experimentales de rendimiento UV y visible usados para los ajustes.