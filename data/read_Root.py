from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot


def _normalise_gas_name(name):
    """Normaliza nombres de gas para comparación."""
    if pd.isna(name):
        return None
    s = str(name).strip().lower()

    mapping = {
        "ar": "ar",
        "argon": "ar",
        "cf4": "cf4",
        "n2": "n2",
        "nitrogen": "n2",
        "co2": "co2",
        "xe": "xe",
        "ne": "ne",
        "he": "he",
        "ch4": "ch4",
        "ic4h10": "ic4h10",
        "c2h6": "c2h6",
        "c3h8": "c3h8",
    }
    return mapping.get(s, s)


def _get_active_gases_from_filename(root_file):
    """
    Extrae pares gas-fracción del nombre del archivo.
    Ejemplo:
        ar_100.0_n2_0.0_60.0kVcm_... -> ["ar"]
        ar_90.0_n2_10.0_60.0kVcm_... -> ["ar", "n2"]

    Devuelve la lista de gases con fracción > 0.
    """
    tokens = root_file.stem.lower().split("_")
    active_gases = []

    i = 0
    while i < len(tokens) - 1:
        gas = _normalise_gas_name(tokens[i])

        try:
            frac = float(tokens[i + 1])
            if frac > 0:
                active_gases.append(gas)
            i += 2
        except ValueError:
            i += 1

    return active_gases


def _build_mapping_table_for_file(table_df, active_gases, argon_update):
    """
    Construye la tabla de mapeo que se usará para un archivo concreto.

    Reglas:
    - Si es mezcla: usa la tabla global tal cual.
    - Si es gas puro: filtra ese gas y renumera desde 0.
    - Si además es Ar puro y argon_update=True:
      inserta hueco 3..7 y desplaza +5 a partir de >2.
    """
    base = table_df.copy()

    # Mezcla: no tocar numeración global
    if len(active_gases) != 1:
        return base

    # Gas puro
    pure_gas = active_gases[0]
    base = base[base["_gas_norm"] == pure_gas].copy()
    base = base.sort_values("level").reset_index(drop=True)

    # Renumerar desde 0
    base["level"] = np.arange(len(base), dtype=int)

    # Corrección específica de Ar puro
    if pure_gas == "ar" and argon_update:
        # Desplazar +5 todo lo que esté por encima del 2
        base.loc[base["level"] > 2, "level"] += 5

        # Insertar filas vacías en 3..7
        extra_cols = list(base.columns)
        gap_rows = []
        for lev in range(3, 8):
            row = {col: pd.NA for col in extra_cols}
            row["level"] = lev
            gap_rows.append(row)

        gap_df = pd.DataFrame(gap_rows, columns=extra_cols)
        base = pd.concat([base, gap_df], ignore_index=True)
        base = base.sort_values("level").reset_index(drop=True)

    return base

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot


def _normalise_gas_name(name):
    """Normaliza nombres de gas para comparación."""
    if pd.isna(name):
        return None
    s = str(name).strip().lower()

    mapping = {
        "ar": "ar",
        "argon": "ar",
        "cf4": "cf4",
        "n2": "n2",
        "nitrogen": "n2",
        "co2": "co2",
        "xe": "xe",
        "ne": "ne",
        "he": "he",
        "ch4": "ch4",
        "ic4h10": "ic4h10",
        "c2h6": "c2h6",
        "c3h8": "c3h8",
    }
    return mapping.get(s, s)


def _get_active_gases_from_filename(root_file):
    """
    Extrae pares gas-fracción del nombre del archivo.
    Ejemplo:
        ar_100.0_n2_0.0_60.0kVcm_... -> ["ar"]
        ar_90.0_n2_10.0_60.0kVcm_... -> ["ar", "n2"]

    Devuelve la lista de gases con fracción > 0.
    """
    tokens = root_file.stem.lower().split("_")
    active_gases = []

    i = 0
    while i < len(tokens) - 1:
        gas = _normalise_gas_name(tokens[i])

        try:
            frac = float(tokens[i + 1])
            if frac > 0 and gas is not None:
                active_gases.append(gas)
            i += 2
        except ValueError:
            i += 1

    return active_gases


def _apply_argon_gap_to_mapping(df):
    """
    Inserta hueco 3..7 y desplaza +5 todos los niveles > 2.
    Se usa cuando hay que corregir la ausencia de los estados de Ar.
    """
    df = df.copy()
    df.loc[df["level"] > 2, "level"] += 5

    gap_rows = []
    cols = list(df.columns)
    for lev in range(3, 8):
        row = {col: pd.NA for col in cols}
        row["level"] = lev
        gap_rows.append(row)

    gap_df = pd.DataFrame(gap_rows, columns=cols)
    df = pd.concat([df, gap_df], ignore_index=True)
    df = df.sort_values("level").reset_index(drop=True)
    return df


def _build_mapping_table_for_file(table_df, active_gases, argon_update):
    """
    Construye la tabla de mapeo que se usará para un archivo concreto.

    Reglas:
    - Gas puro:
        * filtra al gas activo
        * renumera desde 0
        * si es Ar y argon_update=True, inserta hueco 3..7
    - Mezcla:
        * filtra solo los gases activos
        * conserva numeración global de mezcla
        * si hay Ar y argon_update=True, inserta hueco 3..7 globalmente
    """
    base = table_df.copy()

    # Si se detectan gases activos, filtrar solo esos
    if active_gases:
        base = base[base["_gas_norm"].isin(active_gases)].copy()

    base = base.sort_values("level").reset_index(drop=True)

    # Caso gas puro
    if len(active_gases) == 1:
        pure_gas = active_gases[0]

        base = base[base["_gas_norm"] == pure_gas].copy()
        base = base.sort_values("level").reset_index(drop=True)

        # Renumerar desde 0 para gas puro
        base["level"] = np.arange(len(base), dtype=int)

        if pure_gas == "ar" and argon_update:
            base = _apply_argon_gap_to_mapping(base)

        return base

    # Caso mezcla
    if len(active_gases) >= 2:
        if "ar" in active_gases and argon_update:
            base = _apply_argon_gap_to_mapping(base)
        return base

    # Caso residual: si no se detectó nada útil, devolver sin tocar
    return base


def export_hlevels_to_csv(
    folder_path,
    table_path,
    object_name="hLevels",
    argon_update=True
):
    """
    Lee el histograma hLevels de todos los ROOT de una carpeta, lo cruza con una
    tabla externa con el mapeo de niveles, y exporta un CSV enriquecido.

    Parámetros
    ----------
    folder_path : str o Path
        Carpeta donde están los archivos .root
    table_path : str o Path
        Ruta al CSV/tabla con columnas al menos:
            level, gas, state_name
        Opcionalmente puede tener también:
            type, energy_eV
    object_name : str
        Nombre del histograma ROOT a leer
    argon_update : bool
        Si True:
        - en mezclas con Ar: inserta hueco global 3..7
        - en Ar puro: renumera desde 0 e inserta hueco 3..7
    """
    folder = Path(folder_path)
    table_path = Path(table_path)

    if not folder.is_dir():
        raise NotADirectoryError(f"La ruta no es una carpeta válida: {folder}")

    if not table_path.is_file():
        raise FileNotFoundError(f"No existe la tabla de niveles: {table_path}")

    table_df = pd.read_csv(table_path)

    required_cols = {"level", "gas", "state_name"}
    missing_cols = required_cols - set(table_df.columns)
    if missing_cols:
        raise ValueError(
            f"La tabla {table_path} debe contener las columnas {required_cols}. "
            f"Faltan: {missing_cols}"
        )

    table_df = table_df.copy()
    table_df["level"] = pd.to_numeric(table_df["level"], errors="coerce")
    table_df = table_df.dropna(subset=["level"]).copy()
    table_df["level"] = table_df["level"].astype(int)
    table_df["_gas_norm"] = table_df["gas"].apply(_normalise_gas_name)

    root_files = sorted(folder.glob("*.root"))
    if not root_files:
        print(f"No se encontraron archivos .root en: {folder}")
        return []

    out_dir = folder.parent / "csv"
    out_dir.mkdir(exist_ok=True)

    generated_csvs = []

    for root_file in root_files:
        try:
            with uproot.open(root_file) as f:
                keys = f.keys(cycle=False)

                if object_name not in keys:
                    print(f"[AVISO] '{object_name}' no existe en {root_file.name}")
                    continue

                h = f[object_name]
                values = h.values()

                active_gases = _get_active_gases_from_filename(root_file)

                mapping_df = _build_mapping_table_for_file(
                    table_df=table_df,
                    active_gases=active_gases,
                    argon_update=argon_update
                )

                # Los bins del histograma siempre van 0..N-1
                hist_levels = np.arange(len(values), dtype=int)

                df = pd.DataFrame({
                    "level": hist_levels,
                    "n_events": values
                })

                # Cruce con la tabla de mapeo
                df = df.merge(mapping_df, how="left", on="level")

                # Asegurar que también existan niveles del mapping aunque no aparezcan en hLevels
                if not mapping_df.empty:
                    max_level = int(max(df["level"].max(), mapping_df["level"].max()))
                else:
                    max_level = int(df["level"].max())

                df = (
                    df.set_index("level")
                      .reindex(range(0, max_level + 1))
                      .reset_index()
                      .rename(columns={"index": "level"})
                )

                # Poner 0 a los niveles insertados que no estaban en hLevels
                df["n_events"] = df["n_events"].fillna(0)

                # Reinyectar metadata del mapping
                meta_cols = [c for c in mapping_df.columns if c != "level"]
                if meta_cols:
                    meta_df = mapping_df[["level"] + meta_cols].drop_duplicates("level")
                    df = df.drop(columns=[c for c in meta_cols if c in df.columns], errors="ignore")
                    df = df.merge(meta_df, how="left", on="level")

                # Si hay gases activos, vaciar gases no activos
                if active_gases and "_gas_norm" in df.columns:
                    mask_other_gas = df["_gas_norm"].notna() & (~df["_gas_norm"].isin(active_gases))
                    for col in ["gas", "state_name", "type", "energy_eV"]:
                        if col in df.columns:
                            df.loc[mask_other_gas, col] = pd.NA

                # Eliminar columna auxiliar
                if "_gas_norm" in df.columns:
                    df = df.drop(columns=["_gas_norm"])

                preferred_order = [
                    "level",
                    "gas",
                    "state_name",
                    "type",
                    "energy_eV",
                    "n_events"
                ]
                final_cols = [c for c in preferred_order if c in df.columns] + \
                             [c for c in df.columns if c not in preferred_order]
                df = df[final_cols]

                csv_path = out_dir / f"{root_file.stem}.csv"
                df.to_csv(csv_path, index=False)
                generated_csvs.append(csv_path)

                if len(active_gases) == 1:
                    print(f"[OK] CSV generado: {csv_path.name} | gas puro: {active_gases[0]}")
                else:
                    print(f"[OK] CSV generado: {csv_path.name} | gases activos: {active_gases}")

        except Exception as e:
            print(f"[ERROR] No se pudo interpretar '{object_name}' en {root_file.name}: {e}")

    return generated_csvs


def read_data_per_primary_electron(
    folder_path,
    tree_name="dataPerPrimaryElectron",
    electron_branch="nElectrons",
    ion_branch="nIons"
):
    """
    Lee 'dataPerPrimaryElectron' en todos los ROOT de una carpeta y devuelve:
        - media y desviación estándar de nElectrons
        - media y desviación estándar de nIons

    Además guarda histogramas en:
        folder_path/gain_distribution/
    """
    folder = Path(folder_path)
    if not folder.is_dir():
        raise NotADirectoryError(f"La ruta no es una carpeta válida: {folder}")

    root_files = sorted(folder.glob("*.root"))
    if not root_files:
        print(f"No se encontraron archivos .root en: {folder}")
        return pd.DataFrame()

    plot_dir = folder / "gain_distribution"
    plot_dir.mkdir(exist_ok=True)

    summary_rows = []

    for root_file in root_files:
        try:
            with uproot.open(root_file) as f:
                keys_nocycle = f.keys(cycle=False)

                if tree_name not in keys_nocycle:
                    print(f"[AVISO] '{tree_name}' no existe en {root_file.name}")
                    continue

                tree = f[tree_name]
                branches = tree.keys()

                if electron_branch not in branches or ion_branch not in branches:
                    print(
                        f"[AVISO] En {root_file.name} faltan ramas "
                        f"'{electron_branch}' o '{ion_branch}'. "
                        f"Ramas disponibles: {branches}"
                    )
                    continue

                arrays = tree.arrays([electron_branch, ion_branch], library="np")
                ne = np.asarray(arrays[electron_branch], dtype=float)
                ni = np.asarray(arrays[ion_branch], dtype=float)

                ne_mean = np.mean(ne)
                ne_std = np.std(ne, ddof=1) if len(ne) > 1 else 0.0
                ni_mean = np.mean(ni)
                ni_std = np.std(ni, ddof=1) if len(ni) > 1 else 0.0

                summary_rows.append({
                    "file": root_file.name,
                    "ne_mean": ne_mean,
                    "ne_std": ne_std,
                    "ni_mean": ni_mean,
                    "ni_std": ni_std,
                    "n_entries": len(ne)
                })

                plt.figure(figsize=(8, 5))
                plt.hist(ne, bins="auto", edgecolor="black")
                plt.xlabel(electron_branch)
                plt.ylabel("Frecuencia")
                plt.title(f"Distribución de electrones - {root_file.stem}")
                plt.tight_layout()
                plt.savefig(plot_dir / f"{root_file.stem}_{electron_branch}.png", dpi=200)
                plt.close()

                plt.figure(figsize=(8, 5))
                plt.hist(ni, bins="auto", edgecolor="black")
                plt.xlabel(ion_branch)
                plt.ylabel("Frecuencia")
                plt.title(f"Distribución de iones - {root_file.stem}")
                plt.tight_layout()
                plt.savefig(plot_dir / f"{root_file.stem}_{ion_branch}.png", dpi=200)
                plt.close()

                print(f"[OK] Procesado: {root_file.name}")

        except Exception as e:
            print(f"[ERROR] No se pudo procesar {root_file.name}: {e}")

    summary_df = pd.DataFrame(summary_rows)

    if not summary_df.empty:
        summary_csv = folder / "dataPerPrimaryElectron_summary.csv"
        summary_df.to_csv(summary_csv, index=False)
        print(f"[OK] Resumen guardado en: {summary_csv.name}")

    return summary_df

def read_data_per_primary_electron(
    folder_path,
    tree_name="dataPerPrimaryElectron",
    electron_branch="nElectrons",
    ion_branch="nIons"
):
    """
    Lee 'dataPerPrimaryElectron' en todos los ROOT de una carpeta y devuelve:
        - media y desviación estándar de nElectrons
        - media y desviación estándar de nIons

    Además guarda histogramas en:
        folder_path/gain_distribution/
    """
    folder = Path(folder_path)
    if not folder.is_dir():
        raise NotADirectoryError(f"La ruta no es una carpeta válida: {folder}")

    root_files = sorted(folder.glob("*.root"))
    if not root_files:
        print(f"No se encontraron archivos .root en: {folder}")
        return pd.DataFrame()

    plot_dir = folder / "gain_distribution"
    plot_dir.mkdir(exist_ok=True)

    summary_rows = []

    for root_file in root_files:
        try:
            with uproot.open(root_file) as f:
                keys_nocycle = f.keys(cycle=False)

                if tree_name not in keys_nocycle:
                    print(f"[AVISO] '{tree_name}' no existe en {root_file.name}")
                    continue

                tree = f[tree_name]
                branches = tree.keys()

                if electron_branch not in branches or ion_branch not in branches:
                    print(
                        f"[AVISO] En {root_file.name} faltan ramas "
                        f"'{electron_branch}' o '{ion_branch}'. "
                        f"Ramas disponibles: {branches}"
                    )
                    continue

                arrays = tree.arrays([electron_branch, ion_branch], library="np")
                ne = np.asarray(arrays[electron_branch], dtype=float)
                ni = np.asarray(arrays[ion_branch], dtype=float)

                ne_mean = np.mean(ne)
                ne_std = np.std(ne, ddof=1) if len(ne) > 1 else 0.0
                ni_mean = np.mean(ni)
                ni_std = np.std(ni, ddof=1) if len(ni) > 1 else 0.0

                summary_rows.append({
                    "file": root_file.name,
                    "ne_mean": ne_mean,
                    "ne_std": ne_std,
                    "ni_mean": ni_mean,
                    "ni_std": ni_std,
                    "n_entries": len(ne)
                })

                plt.figure(figsize=(8, 5))
                plt.hist(ne, bins="auto", edgecolor="black")
                plt.xlabel(electron_branch)
                plt.ylabel("Frecuencia")
                plt.title(f"Distribución de electrones - {root_file.stem}")
                plt.tight_layout()
                plt.savefig(plot_dir / f"{root_file.stem}_{electron_branch}.png", dpi=200)
                plt.close()

                plt.figure(figsize=(8, 5))
                plt.hist(ni, bins="auto", edgecolor="black")
                plt.xlabel(ion_branch)
                plt.ylabel("Frecuencia")
                plt.title(f"Distribución de iones - {root_file.stem}")
                plt.tight_layout()
                plt.savefig(plot_dir / f"{root_file.stem}_{ion_branch}.png", dpi=200)
                plt.close()

                print(f"[OK] Procesado: {root_file.name}")

        except Exception as e:
            print(f"[ERROR] No se pudo procesar {root_file.name}: {e}")

    summary_df = pd.DataFrame(summary_rows)

    if not summary_df.empty:
        summary_csv = folder / "dataPerPrimaryElectron_summary.csv"
        summary_df.to_csv(summary_csv, index=False)
        print(f"[OK] Resumen guardado en: {summary_csv.name}")

    return summary_df

######################3
def read_data_per_primary_electron(folder_path, tree_name="dataPerPrimaryElectron"):
    """
    Lee el árbol 'dataPerPrimaryElectron' de todos los ROOT de una carpeta,
    extrae las ramas 'ne' y 'ni', calcula:
        - media y desviación estándar de ne
        - media y desviación estándar de ni

    Además genera las gráficas de distribución en:
        folder_path/gain_distribution/

    Parámetros
    ----------
    folder_path : str o Path
        Ruta a la carpeta con los archivos .root
    tree_name : str
        Nombre del árbol a leer (por defecto 'dataPerPrimaryElectron')

    Devuelve
    --------
    pandas.DataFrame
        Tabla resumen con estadísticas por archivo
    """
    folder = Path(folder_path)
    if not folder.is_dir():
        raise NotADirectoryError(f"La ruta no es una carpeta válida: {folder}")

    root_files = sorted(folder.glob("*.root"))
    if not root_files:
        print(f"No se encontraron archivos .root en: {folder}")
        return pd.DataFrame()

    plot_dir = folder / ".." / "gain_distribution"
    plot_dir.mkdir(exist_ok=True)

    summary_rows = []

    for root_file in root_files:
        try:
            with uproot.open(root_file) as f:
                if tree_name not in f:
                    print(f"[AVISO] '{tree_name}' no existe en {root_file.name}")
                    continue

                tree = f[tree_name]
                available_branches = tree.keys()

                if "nElectrons" not in available_branches or "nIons" not in available_branches:
                    print(
                        f"[AVISO] En {root_file.name} faltan ramas. "
                        f"Disponibles: {available_branches}"
                    )
                    continue

                arrays = tree.arrays(["nElectrons", "nIons"], library="np")
                ne = np.asarray(arrays["nElectrons"], dtype=float)
                ni = np.asarray(arrays["nIons"], dtype=float)

                ne_mean = np.mean(ne)
                ne_std = np.std(ne, ddof=1) if len(ne) > 1 else 0.0

                ni_mean = np.mean(ni)
                ni_std = np.std(ni, ddof=1) if len(ni) > 1 else 0.0

                summary_rows.append({
                    "file": root_file.name,
                    "ne_mean": ne_mean,
                    "ne_std": ne_std,
                    "ni_mean": ni_mean,
                    "ni_std": ni_std,
                    "n_entries": len(ne)
                })

                # Gráfica de ne
                plt.figure(figsize=(8, 5))
                plt.hist(ne, bins="auto", edgecolor="white")
                plt.xlabel("ne")
                plt.ylabel("Frecuencia")
                plt.title(f"Distribución de ne\n{root_file.stem}")
                plt.tight_layout()
                plt.savefig(plot_dir / f"{root_file.stem}_ne.pdf", dpi=200)
                plt.close()

                # Gráfica de ni
                # plt.figure(figsize=(8, 5))
                # plt.hist(ni, bins="auto", edgecolor="white")
                # plt.xlabel("ni")
                # plt.ylabel("Frecuencia")
                # plt.title(f"Distribución de ni\n{root_file.stem}")
                # plt.tight_layout()
                # plt.savefig(plot_dir / f"{root_file.stem}_ni.pdf", dpi=200)
                # plt.close()

                print(f"[OK] Procesado: {root_file.name}")

        except Exception as e:
            print(f"[ERROR] No se pudo procesar {root_file.name}: {e}")

    summary_df = pd.DataFrame(summary_rows)

    if not summary_df.empty:
        summary_csv = folder / ".." / "dataPerPrimaryElectron_summary.csv"
        summary_df.to_csv(summary_csv, index=False)
        print(f"[OK] Resumen guardado en: {summary_csv}")

    return summary_df


###########S