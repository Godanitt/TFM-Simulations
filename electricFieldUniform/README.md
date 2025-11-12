# Simulación de Campo Eléctrico Uniforme

Este módulo realiza la simulación de avalanchas electrónicas en un **campo eléctrico uniforme** mediante el uso de **Garfield++** y **Magboltz**.  
Incluye un programa principal en C++ (`uniformE.C`) que efectúa la simulación física, y un script auxiliar en Python (`runUniform.py`) que automatiza la compilación y ejecución del código. Muy recomendable hacerlo así. Si no, hay que ejecutar en orden, desde la carpeta principal electricFieldUniform: 

```bash
rm -f build
mkdir build
cd build
cmake ..
make 
./uniformE ../rootArchives/ar1.0cf4_15.0kVcmbar_20npe.root 15000 1.00 1 20 ar 99.0 cf4 1.0
```

En el caso de ejecucción con python simplemente:

python3 runUniform.py

Las condiciones de la simulación se cambian dentro del propio archivo. 


---

## 📂 Estructura del directorio



electricFieldUniform/
│
├── build/ # Carpeta de compilación generada con CMake
├── rootArchives/ # Archivos .root con los resultados de las simulaciones
│
├── CMakeLists.txt # Configuración de compilación del código C++
├── runUniform.py # Script en Python para compilar y ejecutar las simulaciones
└── uniformE.C # Código fuente en C++ que realiza la simulación principal

## ⚙️ Dependencias

- **Garfield++**  
- **ROOT**  
- **CMake**  
- **Python 3.x**

Asegúrate de que las bibliotecas de Garfield++ y ROOT estén correctamente instaladas y configuradas en tu entorno antes de compilar.

---

## 🧩 Descripción de los archivos principales

### `uniformE.C`
Código principal en C++ que:
- Configura el gas o mezcla de gases mediante **MediumMagboltz**.  
- Define un **campo eléctrico uniforme** usando `ComponentUser`.  
- Simula avalanchas microscópicas de electrones mediante `AvalancheMicroscopic`.  
- Registra los resultados en un archivo **ROOT**, con información de:
  - Niveles de excitación,
  - Energías electrónicas,
  - Número de electrones e iones producidos.  

### `runUniform.py`

Script en Python que:
- Limpia y recompila el proyecto mediante **CMake**.  
- Define los parámetros de simulación (gases, presión, campo, etc.).  
- Ejecuta automáticamente el binario `uniformE` con los argumentos apropiados.  

Se puede seleccionar en el archivo de python la presión (bar), gap (mm), campo eléctrico (V/cm), mezcla de gases (2 gases actualmente) y número de eventos/electrones primarios: 

```python
pressure = [1]      # Presión en bar
gap = [1]           # Distancia del campo (mm)
fieldE = [15000]    # Campo eléctrico (V/cm)
gas1 = ["ar"]       # Gas1 (código Magboltz)
mixture1 = [99.0]   # % mezcla Gas1
gas2 = ["cf4"]      # Gas2 (código Magboltz)
mixture2 = [1.0]    # % mezcla Gas2
npe = [20]          # Nº de electrones primarios
```
