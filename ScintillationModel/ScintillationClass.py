import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from typing import Dict, List, Any, Optional
import scipy.optimize as opt
import inspect
import re




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
    - Grafica de modelos de Centelleo a parámetros dados 
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
        self.fCF4_orig = self.poblation_degrad["fCF4"]
        
        # Corrección de Degrad con interpolacion
        self.poblation_degrad_corr=self._compute_poblation_degrad_corr() 
        
        # Inicializamos donde se guardan los resultados de los fits
        self.fit_results = {}

        # Para las gráficas de los fits
        self.plot_settings = {
            "normalization": {},   # band → config
            "show_exp": {},        # band → list of pressures
            "show_teo": {}         # band → list of pressures
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

            # Recorremos todas las columnas físicas
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
        self.fit_results[f"{band}_N0={n0}"] = result

        return result
    
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
            
    def choosePlotNormalization(self, band, mode="N0", value=1.0, idx_ref=-1):
        """
        Configura la normalización de una curva teórica o experimental.

        mode:
            "none"  → no normalizar
            "N0"    → normalizar a f(variable) en presión = value (default 1 bar)
            "index" → normalizar usando índice idx_ref
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
            }

        if band not in self.plot_settings["normalization"]:
            self.plot_settings["normalization"][band] = ("none", None)

        # --- Aplicar modo ---
        if mode == "none":
            self.plot_settings["normalization"][band] = ("none", None)

        elif mode == "N0":
            self.plot_settings["normalization"][band] = ("N0", value)

        elif mode == "index":
            self.plot_settings["normalization"][band] = ("index", idx_ref)

        else:
            raise ValueError("mode debe ser: 'none', 'N0', 'index'")

        
    def EnableExperimentalData(self, band, n):
        self._ensure_band_in_settings(band)
        if n not in self.plot_settings["show_exp"][band]:
            self.plot_settings["show_exp"][band].append(n)
            
    def EnableTeoCurve(self, band, n):
        self._ensure_band_in_settings(band)
        if n not in self.plot_settings["show_teo"][band]:
            self.plot_settings["show_teo"][band].append(n)
    def plot_teoCurve(self, band, n=1.0, figsize=(7,5), savefig=None):
        # Asegurar que el canal existe
        if band not in self.theory_functions:
            raise ValueError(f"No existe función teórica para el canal '{band}'.")

        self._ensure_band_in_settings(band)

        f_th = self.theory_functions[band]
        fCF4 = self.fCF4

        # Normalización seleccionada
        norm_mode, norm_val = self.plot_settings["normalization"][band]

        def normalize(arr):
            if norm_mode == "none":
                return arr
            elif norm_mode == "N0":
                ref = f_th(x=self.fit_results[band], fCF4=fCF4, n=norm_val)[-1]
                return arr/ref
            elif norm_mode == "index":
                return arr / arr[norm_val]
            return arr

        plt.figure(figsize=figsize)

        # -------------------------------
        # Curvas teóricas seleccionadas
        # -------------------------------
        for n_plot in self.plot_settings["show_teo"][band]:
            y = f_th(x=self.fit_results[band], fCF4=fCF4, n=n_plot)
            y = normalize(y)
            plt.plot(fCF4, y, label=f"Theory {band}, {n_plot} bar", lw=2)

        # -------------------------------
        # Datos experimentales
        # -------------------------------
        dfY = self.yields[band]
        for n_plot in self.plot_settings["show_exp"][band]:
            col = f"{n_plot}bar"
            if col not in dfY:
                print(f"[WARN] No hay datos experimentales para {col}")
                continue
            y = dfY[col].to_numpy()
            y = normalize(y)
            plt.scatter(fCF4, y, marker="o", label=f"Exp {band}, {n_plot} bar")

        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("fCF4")
        plt.ylabel("yield normalizado")
        plt.legend()
        plt.grid(True, alpha=0.3)

        if savefig:
            plt.savefig(savefig, dpi=300)
        else:
            plt.show()
