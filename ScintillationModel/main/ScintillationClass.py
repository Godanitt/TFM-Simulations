import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from typing import Dict, List, Any, Optional
import scipy.optimize as opt
import inspect
import re


cmap = plt.get_cmap("viridis")

def normalize_tokens(s):
    """Convierte 'Ar dble Star' → ['ar','dble','star']."""
    s = s.lower().replace("_", " ").strip()
    return [tok for tok in re.split(r"[ \-]+", s) if tok]


def match_param_to_species(param, species):
    """
    Empareja P_CF3 <-> CF3
              P_Ar_dbleStar <-> Ar dble Star
              P_Ar3rd <-> Ar 3rd

    Regla:
        TODOS los tokens del parámetro deben existir en la especie.
    """

    # quitar prefijo P_
    if param.lower().startswith("p_"):
        param_clean = param[2:]
    else:
        param_clean = param

    param_tokens   = normalize_tokens(param_clean)
    species_tokens = normalize_tokens(species)

    # Condición correcta: todos los tokens de param están en species
    return all(tok in species_tokens for tok in param_tokens)

def darken(color, factor=0.9):
    """
    Mezcla el color con negro.
    factor=0 → negro
    factor=1 → color original
    """
    r, g, b, a = color
    return (factor*r, factor*g, factor*b, a)


###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
############# CORRECCIÓN DE LAS POBLACIONES DE DEGRAD ##############################

class Scintillation:
    """
    Clase para almacenar distintos diccionarios o DataFrames relacionados con:
    - Poblaciones de yields (vis, uv)
    - Poblaciones de degradación
    - Funciones de centelleo
    Además hace:
    - Comparaión entre modelos de centelleo con poblaciones de Degrad y datos experimentales
    - Ajuste a nuevos modelos de poblaciones de Degrad
    - Expansión de poblaciones de Degrad a diferentes concentraciones a través de interpolación.
    - Grafica de modelos de Centelleo a parámetros dados o ajustados, a presión elegida, con datos experimentales. 
    - Diferentes opciones a la hora de graficar.
    """

    def __init__(
        self,
        yields: Optional[Dict[str, List[Any]]] = None,
        poblation_degrad: Optional[Dict[str, List[Any]]] = None,
        scintillation_models: Optional[Dict[str, List[Any]]] = None,
    ):
        # Entrada, con lo que se define: 
        self.yields = yields
        self.poblation_degrad = poblation_degrad
        self.scintillation_models = scintillation_models
        
        # Importantes
        self.fCF4 = self.yields["fCF4"]/100
        self.min_fCF4_10log = np.log10(np.min(self.fCF4))
        self.max_fCF4_10log = np.log10(np.max(self.fCF4))
        self.fCF4_orig = self.poblation_degrad["fCF4"]
        
        # Corrección de Degrad con interpolacion
        self.poblation_degrad_corr=self._compute_poblation_degrad_corr() 
        
        # Inicializamos donde se guardan los resultados de los fits
        self.fit_results = {}

        # Por si queremos meter parámetros a mano
        self.manual_parameters = {}  
        
        # Flags para elegir entre parámetros del fit o manuales
        self.use_manual = {}   # band -> True/False

        # Nos permite dibujar las diferentes contribuciones
        self.enabled_components = {}   # p.ej: { "vis": ["Direct","Transf"] }


        # Para las gráficas de los fits
        self.plot_settings = {
            "normalization": {},   # band → config
            "show_exp": {},        # band → list of pressures
            "show_teo": {},        # band → list of pressures
            "global_norm": {}      # NUEVO
        }

    ###################################################################################
    ############# CORRECCIÓN DE LAS POBLACIONES DE DEGRAD ##############################
    
    def _compute_poblation_degrad_corr(self) -> Dict[str, Any]:
        """
        Interpola cada DataFrame de poblation_degrad usando la malla nueva self.fCF4.
        El eje original es self.poblation_degrad["fCF4"].
        """
        result = {}

        x_orig = self.fCF4_orig   # eje original
        x_new  = self.fCF4        # eje nuevo

        for name, obj in self.poblation_degrad.items():

            # saltar la clave del eje
            if name == "fCF4":
                continue

            # solo DataFrames
            if isinstance(obj, pd.DataFrame):

                corr_df = pd.DataFrame(index=x_new)

                for col in obj.columns:

                    if "Err" in col:
                        continue   # saltamos errores

                    y_old = obj[col].to_numpy(dtype=float)
                    y_new = np.interp(x_new, x_orig, y_old)

                    corr_df[col] = y_new

                result[name] = corr_df

        return result

    ###################################################################################
    #############  GRAFICO POBLACIONES DE DEGRAD INTERPOLADAS ##############################
    
    def plotPoblationInterpolation(self, name,savefig=""):

        original = self.poblation_degrad[name]

        x_orig = self.fCF4_orig   # ESTE es el eje original real

        plt.figure(figsize=(7,5))

        for col in original.columns:

            if "Err" in col:
                continue

            y_orig = original[col].to_numpy()

            # puntos rojos: datos originales
            plt.scatter(x_orig, y_orig, color="red", s=40)

            # línea azul: recta trozo a trozo entre datos originales
            plt.plot(x_orig, y_orig, color="blue", linewidth=2)

        plt.xlabel("fCF4")
        plt.ylabel("Población")
        plt.title(f"Interpolación lineal trozo a trozo de {name}")
        plt.grid(True, alpha=0.3)
        plt.xscale("log")
        plt.yscale("log")
        plt.tight_layout()
        if not(savefig==""):
            plt.savefig(savefig,dpi=300)
        else:
            plt.show()

    #############################################################################
    ############## FUNCIONES QUE AÑADEN INFORMACION ##############################

    def add_yields(self, key: str, values: List[Any]):
        if self.yields is None:
            self.yields = pd.DataFrame({key: pd.Series(values)})
        else:
            self.yields[key] = pd.Series(values)

    def add_poblation_degrad(self, key: str, values: List[Any]):
        if self.poblation_degrad is None:
            self.poblation_degrad = pd.DataFrame({key: pd.Series(values)})
        else:
            self.poblation_degrad[key] = pd.Series(values)

    def add_scintillation_models(self, key: str, values: List[Any]):
        if self.scintillation_models is None:
            self.scintillation_models = pd.DataFrame({key: pd.Series(values)})
        else:
            self.scintillation_models[key] = pd.Series(values)

    ###########################################################################################
    ############# FUNCIONES QUE CREAN LAS COMBINACIONES TEORICAS ##############################
    ###########################################################################################
        
       
    def build_theory_functions(self, scintillation_definition):
        """
        Construye funciones teóricas combinando modelos físicos y pesos.
        Guarda las funciones resultantes en self.theory_functions.
        Sirven para obtener los parámetros del ajuste.
        """
        self.theory_functions = {}
        contribs = {}   
        
        for band_name, components in scintillation_definition.items():

            comp = list(components.items())  # para cerrar bien valores

            def make_theory_func(comp):

                def theory_func(x, fCF4, n, **kwargs):
                    total = 0.0
                    idx = 0

                    for model_name, properties in comp:

                        model_func = self.scintillation_models[model_name]
                        modes = properties   # puede haber varias (["Relajacion","Centelleo"])
                        weight = 1.0         # empezamos con factor 1

                        for mode in modes:

                            if mode == "Probabilidad":
                                weight *= x[idx]
                                idx += 1

                            elif mode == "Relajacion":
                                k = x[idx]
                                weight *= n * fCF4 / (k + n * fCF4)
                                idx += 1

                            elif mode == "Centelleo":
                                k = x[idx]
                                weight *= k / (k + n * fCF4)
                                idx += 1

                            else:
                                weight *= 1.0


                        # ===== 2. Parámetros del modelo físico =====
                        sig = inspect.signature(model_func)
                        
                        params = {}

                        for par in sig.parameters.values():
                            pname = par.name
                            

                            if pname == "f_cf4":
                                params[pname] = fCF4
                                continue

                            if pname == "n":
                                params[pname] = n
                                continue

                            # --- BÚSQUEDA GENERAL DE LA POBLACIÓN ---
                            found = False
                            
                            for species, df in self.poblation_degrad.items():

                                # ⛔ Saltar la malla fCF4 porque NO es un DataFrame
                                if species == "fCF4":
                                    continue
                                
                                flag = match_param_to_species(pname,species)
                                #flag=True
                                if flag: 
                                    valid_cols = [c for c in df.columns
                                                if "err" not in c.lower()]

                                    if not valid_cols:
                                        continue
                                    
                                    # primera columna válida (general para cualquier gas)
                                    col = valid_cols[0]
                                    y_old = df[col].to_numpy()
                                    y_new = np.interp(fCF4, self.fCF4_orig, y_old)

                                    
                                    params[pname] = y_new
                                    
                                    found = True
                                    
                            
                            #print(params)
                            if not found:
                                raise ValueError(
                                    f"No se encontró población para parámetro '{pname}' "
                                    f"en modelo '{model_name}'."
                                )

                        # ===== 3. Suma ponderada =====
                        
                        total += weight * model_func(**params)
                        contribs[model_name] = weight * model_func(**params)

                    #return total
                    if kwargs.get("return_components", False):
                        return total, contribs
                    return total

                return theory_func

            self.theory_functions[band_name] = make_theory_func(comp)

    ##########################################################################
    ######################### AJUSTE DE PARÁMETROS ###########################
    ##########################################################################

    def fit_parameters_chooseNorma(
        self,
        band: str,
        x0: np.ndarray,
        n0: float = 1.0,
        idx_ref: int = -1,
        method: str = "BFGS",
    ):
        """
        Parametros: 
            - **band**: elijes a cual ajusta de los datos experimentales ("vis","uv"...)
            - **x0**: valores iniciales de los parámetros
        Ajusta los parámetros x para el canal `band` ('vis', 'uv', ...)
        usando una normalización a un punto de referencia:
            - presión de referencia n0 (ej. 1.0 bar)
            - índice idx_ref en la malla fCF4 (ej. -1 → último punto)

        Chi² se calcula sobre cantidades normalizadas:
            y_exp_norm  = y_exp / y_exp_ref
            y_th_norm   = y_th  / y_th_ref
        """

        if not hasattr(self, "theory_functions") or band not in self.theory_functions:
            raise ValueError(
                f"No existe función teórica para el canal '{band}'. "
                "Llama antes a build_theory_functions()."
            )

        if band not in self.yields:
            raise ValueError(
                f"No hay datos de yield para el canal '{band}' en self.yields."
            )

        theory_func = self.theory_functions[band]
        df_yield = self.yields[band]
        fCF4 = self.fCF4

        # Columnas físicas (sin errores)
        cols_phys = [c for c in df_yield.columns if "err" not in c.lower()]

        # Buscamos la columna cuya presión coincida con n0 (por ejemplo "1.0bar")
        col_ref = None
        for col in cols_phys:
            try:
                n_col = float(col.replace("bar", ""))
            except Exception:
                continue
            if abs(n_col - n0) < 1e-6:
                col_ref = col
                break

        if col_ref is None:
            raise ValueError(
                f"No se encontró ninguna columna de presión ~={n0} bar para el canal '{band}'."
            )

        # Esta función construye el chi² con normalización a n0, idx_ref
        def chi2(x):
            chi2_val = 0.0

            # Primero obtenemos la referencia experimental y teórica en n0
            y_ref_exp = df_yield[col_ref].to_numpy()
            y0_exp = y_ref_exp[idx_ref]

            # si el valor de referencia experimental es 0, evitar problemas
            if y0_exp == 0:
                return 1e30  # penalizamos esta elección de x

            # error asociado a la referencia (si existe)
            err_ref = None
            for ec in [f"Err {col_ref}", f"Err_{col_ref}", f"{col_ref} Err", f"{col_ref}_Err"]:
                if ec in df_yield.columns:
                    err_ref = df_yield[ec].to_numpy()
                    break

            # Recorremos todas las columnas físicas -> Datos Experimentales 
            for col in cols_phys:
                y_exp = df_yield[col].to_numpy()

                # errores
                sigma = None
                for ec in [f"Err {col}", f"Err_{col}", f"{col} Err", f"{col}_Err"]:
                    if ec in df_yield.columns:
                        sigma = df_yield[ec].to_numpy()
                        break
                if sigma is None:
                    sigma = np.ones_like(y_exp)

                # presión de esta columna
                try:
                    n_val = float(col.replace("bar", ""))
                except Exception:
                    n_val = 1.0

                # modelo teórico
                y_th = theory_func(x=x, fCF4=fCF4, n=n_val)

                # referencia teórica en n0 (misma idx_ref)
                y_ref_th = theory_func(x=x, fCF4=fCF4, n=n0)
                y0_th = y_ref_th[idx_ref]

                if y0_th == 0:
                    return 1e30

                # normalizamos datos y modelo
                y_exp_norm = y_exp / y0_exp
                y_th_norm = y_th / y0_th

                # normalizamos errores (simplificado)
                sigma_norm = sigma / abs(y0_exp)

                mask = sigma_norm > 0
                chi2_val += np.sum(((y_th_norm[mask] - y_exp_norm[mask]) / sigma_norm[mask]) ** 2)

            return chi2_val

        result = opt.minimize(chi2, x0, method=method)

        if not hasattr(self, "fit_results"):
            self.fit_results = {}
    
        popt = result.x
        self.fit_results[band] = popt
        return popt
    
    def set_manual_parameters(self, band, params):
        """
        Establece parámetros manuales para una banda.
        """
        if band not in self.theory_functions:
            raise ValueError(f"No existe función teórica para banda '{band}'.")

        # guardar parámetros manuales
        self.fit_results[band] = np.array(params, dtype=float)

        # activar uso de parámetros manuales
        if not hasattr(self, "use_manual"):
            self.use_manual = {}
        self.use_manual[band] = True


    def use_fit_parameters(self, band):
        """
        Hace que la banda use los parámetros del fit.
        """
        if not hasattr(self, "use_manual"):
            self.use_manual = {}
        self.use_manual[band] = False

    
    def ensure_manual_flag(self, band):
        """
        Asegura que existe self.use_manual[band].
        Si no existe, lo inicializa a False.
        """
        if band not in self.use_manual:
            self.use_manual[band] = False


    ##########################################################################
    ######################### GRAFICA ###########################
    ##########################################################################
    
    def _ensure_band_in_settings(self, band):
        if band not in self.plot_settings["normalization"]:
            self.plot_settings["normalization"][band] = ("none", None)

        if band not in self.plot_settings["show_exp"]:
            self.plot_settings["show_exp"][band] = []

        if band not in self.plot_settings["show_teo"]:
            self.plot_settings["show_teo"][band] = []
            
    def choosePlotNormalization(self, band, mode="index", value=1.0, idx_ref=-1, global_bands=None):
        """
        Configura la normalización de un canal específico.

        Existen los siguientes modos de normalización:

        - **mode = "none"**
            → No se aplica ningún tipo de normalización.

        - **mode = "N0"**
            → Se normaliza usando el valor de la curva teórica evaluada 
              a la presión indicada en `value`.
              (Por ejemplo: value = 1.0 → normaliza a 1 bar)

        - **mode = "index"**
            → Se normaliza usando el valor del array en la posición `idx_ref`.
              (Por defecto idx_ref = -1, que normalmente corresponde a la 
               concentración fCF4 = 100%)

        - **mode = "global"**
            → Normalización conjunta entre varias bandas. Debes pasar 
              `global_bands=[...]` con las bandas involucradas. 
              El valor de referencia se calcula sumando los valores de estas bandas 
              en el índice `idx_ref`.

        Parámetros
        ----------
        band : str
            Nombre del canal (banda) a normalizar.

        mode : {"none", "N0", "index", "global"}
            Tipo de normalización deseada.

        value : float
            Presión (en bar) usada únicamente cuando mode="N0".

        idx_ref : int
            Índice de referencia usado en mode="index" o mode="global".

        global_bands : list[str], opcional
            Lista de bandas para normalización conjunta cuando mode="global".
        """

        # --- Validación general correcta ----
        if band not in self.theory_functions:
            raise ValueError(f"No existe banda '{band}' en theory_functions")

        # --- Asegurar la estructura interna ---
        if not hasattr(self, "plot_settings"):
            self.plot_settings = {
                "normalization": {},
                "show_exp": {},
                "show_teo": {},
                "global_norm": {},     # NUEVO
            }

        # Guardar siempre algo
        if band not in self.plot_settings["normalization"]:
            self.plot_settings["normalization"][band] = ("none", None)

        # -----------------------------
        #  MODO 1: normalización normal
        # -----------------------------
        if mode == "none":
            self.plot_settings["normalization"][band] = ("none", None)

        elif mode == "N0":
            self.plot_settings["normalization"][band] = ("N0", value)

        elif mode == "index":
            self.plot_settings["normalization"][band] = ("index", idx_ref)

        # -----------------------------
        #  MODO 2: normalización GLOBAL
        # -----------------------------
        elif mode == "global":
            if global_bands is None:
                raise ValueError("Debes proporcionar global_bands=[...] en modo global")

            # Guardamos SOLO UNA vez — solo se usa cuando se grafique
            self.plot_settings["global_norm"][band] = {
                "bands": global_bands,
                "mode": "index",   # siempre por índice (fácil y lógico)
                "idx": idx_ref
            }

            # Y marcamos al canal como normalización global
            self.plot_settings["normalization"][band] = ("global", None)

        else:
            raise ValueError("mode debe ser: 'none', 'N0', 'index', 'global'")  

        
    def EnableExperimentalData(self, band, n):
        self._ensure_band_in_settings(band)
        if n not in self.plot_settings["show_exp"][band]:
            self.plot_settings["show_exp"][band].append(n)
            
    def EnableTeoCurve(self, band, n):
        self._ensure_band_in_settings(band)
        if n not in self.plot_settings["show_teo"][band]:
            self.plot_settings["show_teo"][band].append(n)
            
    def EnableComponent(self, band, name):
        if band not in self.enabled_components:
            self.enabled_components[band] = []
        if name not in self.enabled_components[band]:
            self.enabled_components[band].append(name)

                
    def plot_teoCurve(self, band, n=1.0, figsize=(7,5), savefig=None):
        # Asegurar que el canal existe
        if band not in self.theory_functions:
            raise ValueError(f"No existe función teórica para el canal '{band}'.")

        self._ensure_band_in_settings(band)

        f_th = self.theory_functions[band]
        fCF4 = self.fCF4

        # Normalización seleccionada
        norm_mode, norm_val = self.plot_settings["normalization"][band]

        # ---- FUNCIÓN PARA NORMALIZACIÓN GLOBAL  ----
        def compute_global_reference(idx, mode="teo"):
            """
            mode = 'exp' → usa datos experimentales
            mode = 'teo' → usa curvas teóricas
            """

            info = self.plot_settings["global_norm"][band]
            bands = info["bands"]
            idx_ref = info["idx"]

            ref_sum = 0.0

            for b in bands:

                if mode == "exp":   # NORMALIZA CONTRA EXPERIMENTAL
                    dfY = self.yields[b]
                    col = f"{n}bar"

                    if col not in dfY:
                        raise ValueError(
                            f"No hay datos experimentales para banda '{b}' a {n} bar"
                        )

                    ref_sum += dfY[col].to_numpy()[idx_ref]

                else:               # NORMALIZA CONTRA TEÓRICA
                    f_b = self.theory_functions[b]
                    x_b = self.fit_results[b]
                    y_b = f_b(x=x_b, fCF4=fCF4, n=n)
                    ref_sum += y_b[idx_ref]

            return ref_sum

        # ---- NORMALIZACIÓN TEÓRICA ----
        def normalize(arr,mode="teo"):
            if norm_mode == "none":
                return arr

            elif norm_mode == "index":
                return arr / arr[norm_val]

            elif norm_mode == "global":
                ref_sum = compute_global_reference(norm_val,mode=mode)
                return arr / ref_sum

            return arr

        # ---- NORMALIZACIÓN EXP + ERRORES ----
        def normalize_pair(y, sy,mode="exp"):
            if norm_mode == "none":
                return y, sy

            elif norm_mode == "index":
                ref = y[norm_val]
                sref = sy[norm_val]
                sy_new = np.sqrt((sy/ref)**2 + (sref*y/ref**2)**2)
                return y/ref, sy_new

            elif norm_mode == "global":
                ref_sum = compute_global_reference(norm_val,mode=mode)
                return y / ref_sum, sy / ref_sum

            return y, sy

        # ---------- PLOT ----------
        plt.figure(figsize=figsize)

        # Elegimos colores
        n_theo = len(self.plot_settings["show_teo"][band])
        n_exp  = len(self.plot_settings["show_exp"][band])
        maxima = max(n_theo, n_exp)
        colors = cmap(np.linspace(0.2, 0.8, maxima))

        # --- TEORÍA ---
        k = 0
        for n_plot in self.plot_settings["show_teo"][band]:
            fCF4_array = np.logspace(self.min_fCF4_10log, self.max_fCF4_10log, num=100)
            
            y = f_th(x=self.fit_results[band], fCF4=fCF4_array, n=n_plot)
            
            y = normalize(y)
            plt.plot(
                fCF4_array*100,
                y,
                label=f"Theory {band}, {n_plot} bar",
                lw=2,
                color=darken(colors[k], factor=0.6)
            )
            
            k += 1

        # --- EXPERIMENTAL ---
        dfY = self.yields[band]
        k = 0
        for n_plot in self.plot_settings["show_exp"][band]:
            col = f"{n_plot}bar"
            scol = f"Err {n_plot}bar"

            if col not in dfY:
                print(f"[WARN] No hay datos experimentales para {col}")
                continue

            y = dfY[col].to_numpy()
            sy = dfY[scol].to_numpy()

            y, sy = normalize_pair(y, sy)

            plt.errorbar(
                fCF4*100, y, yerr=sy,
                marker="o",
                linestyle="none",
                color=colors[k],
                label=f"Exp {band}, {n_plot} bar"
            )
            k += 1

        # --- Opciones comunes ---
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("fCF4 %")
        plt.ylabel("yield normalizado")
        plt.legend()
        plt.grid(True, alpha=0.3)

        if savefig:
            plt.savefig(savefig, dpi=300)
        else:
            plt.show()
