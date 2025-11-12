# Proyecto de Simulación y Análisis

Este repositorio contiene varios módulos y scripts relacionados con la simulación y análisis de campos eléctricos, atenuación de partículas y estudios de centelleo.  

A continuación se muestra una descripción general de la estructura del proyecto.

En general todo se ejecuta en linux.

---

## Índice

1. [Atenuacion alpha](#atenuacion-alpha)
2. [electricField Uniform](#electricfield-uniform)
3. [fatGem](#fatgem)
4. [Study_Scintillation_Degrad](#study_scintillation_degrad)
   - [version1](#version1)
5. [.vscode](#vscode)
6. [.gitignore](#gitignore)
7. [README.md](#readmemd)

---

## Atenuacion alpha
Carpeta que contiene simulaciones o datos relacionados con la **atenuación de partículas alfa** para hallar el grosor de la atenuación que hace falta darle en el experimento.

---

## electricField Uniform

Simulación de electrones derivando en **campo eléctrico uniforme** para cualquier mezcla de gases, presión, gap y campo eléctrico. 

---

## fatGem

Simulación de electrones derivando en **campo eléctrico uniforme** para cualquier mezcla de gases, presión, gap y campo eléctrico. 

---

## Study_Scintillation_Degrad

Estudio del centelleo de Ar-CF4 según el modelo de Pablo Amedo (https://doi.org/10.3389/fdest.2023.1282854). Nos centramos en su refinamiento y mejora, comparándo resultados con poblaciones de Degrad.

### version1

Versión en la que actualmente se está trabajando. Incluye barras de errores, más datos de Degrad, estudio de poblaciones de Degrad vs concentración de CF4 y su fit lineal para la obtención de parámetros autoquenching y relajación, así como tasa de emisión del CF3 directo. 

