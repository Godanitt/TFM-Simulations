import os
from pathlib import Path
import numpy as np
import pandas as pd


def _parse_csv_metadata_from_name(filename, concentration_gas=None):
    """
    Extrae metadata del nombre del archivo CSV.

    Ejemplo esperado:
    ar_90.0_n2_10.0_60.0kVcm_1.000bar_0.1071mm_100npe.csv

    Devuelve:
        concentration, electric_field, pressure, gap_mm, gas_fractions
    """
    stem = Path(filename).stem
    tokens = stem.split("_")

    gas_fractions = {}
    electric_field = np.nan
    pressure = np.nan
    gap_mm = np.nan

    i = 0
    while i < len(tokens) - 1:
        # Intentar leer pares gas-fracción
        try:
            frac = float(tokens[i + 1])
            gas = tokens[i].lower()
            gas_fractions[gas] = frac
            i += 2
            continue
        except ValueError:
            pass

        token = tokens[i].lower()

        if token.endswith("kvcm"):
            try:
                electric_field = float(token.replace("kvcm", ""))
            except ValueError:
                pass

        elif token.endswith("bar"):
            try:
                pressure = float(token.replace("bar", ""))
            except ValueError:
                pass

        elif token.endswith("mm"):
            try:
                gap_mm = float(token.replace("mm", ""))
            except ValueError:
                pass

        i += 1

    concentration = np.nan
    if concentration_gas is not None:
        concentration = gas_fractions.get(concentration_gas.lower(), np.nan)
    else:
        # Si no se especifica, toma la segunda concentración si existe
        if len(gas_fractions) >= 2:
            concentration = list(gas_fractions.values())[1]
        elif len(gas_fractions) == 1:
            concentration = list(gas_fractions.values())[0]

    return {
        "concentration": concentration,
        "electric_field": electric_field,
        "pressure": pressure,
        "gap_mm": gap_mm,
        "gas_fractions": gas_fractions
    }


def read_garfield_csv_folder(
    folder_path,
    dataframe,
    output_dir,
    output_general_name,
    concentration_gas=None,
    use_poisson_error=True
):
    """
    Lee todos los CSV de una carpeta generados a partir de Garfield++ y crea
    tablas resumen con lógica similar a read_degrad.

    Se asume que cada CSV tiene columnas como:
        level, gas, state_name, type, energy_eV, n_events

    Parámetros
    ----------
    folder_path : str
        Carpeta donde están los CSV
    dataframe : pandas.DataFrame
        DataFrame de configuración con filas tipo:
            - "name principal"
            - "gas"
            - "energy up"
            - "energy low"
            - "name output"
        Opcional:
            - "type"
    output_dir : str
        Carpeta para los CSV individuales
    output_general_name : str
        Nombre base del CSV general
    concentration_gas : str or None
        Gas cuya concentración se usará como columna "concentration"
    use_poisson_error : bool
        Si True, usa sqrt(N) como error

    Devuelve
    --------
    population_gen : pandas.DataFrame
        DataFrame general con todas las poblaciones
    """
    folder = Path(folder_path)
    output_dir = Path(output_dir)

    if not folder.is_dir():
        raise NotADirectoryError(f"La ruta no es una carpeta válida: {folder}")

    output_dir.mkdir(exist_ok=True, parents=True)

    csv_files = sorted(folder.glob("*.csv"))
    if not csv_files:
        print(f"No se encontraron CSVs en: {folder}")
        return pd.DataFrame()

    # Tabla base con metadata extraída del nombre del archivo
    rows_meta = []
    for csv_file in csv_files:
        meta = _parse_csv_metadata_from_name(csv_file.name, concentration_gas=concentration_gas)
        rows_meta.append({
            "file": str(csv_file),
            "concentration": meta["concentration"],
            "electric_field": meta["electric_field"],
            "gap_mm": meta["gap_mm"],
            "pressure": meta["pressure"]
        })

    population_gen = pd.DataFrame(rows_meta)

    # Construcción de poblaciones
    for nombre_col in dataframe.columns:

        population = population_gen[
            ["concentration", "electric_field", "gap_mm", "pressure"]
        ].copy()

        name_of_state = dataframe.loc["name principal", nombre_col]
        gas = dataframe.loc["gas", nombre_col]
        energy_upper_limit = dataframe.loc["energy up", nombre_col]
        energy_lower_limit = dataframe.loc["energy low", nombre_col]
        name_of_output = dataframe.loc["name output", nombre_col]

        process_type = None
        if "type" in dataframe.index:
            process_type = dataframe.loc["type", nombre_col]

        # Forzar lista
        if isinstance(name_of_state, str):
            name_of_state = [name_of_state]
        elif pd.isna(name_of_state):
            name_of_state = []

        for i, csv_file in enumerate(population_gen["file"]):

            df = pd.read_csv(csv_file)

            expected_cols = {"gas", "state_name", "energy_eV", "n_events"}
            missing_cols = expected_cols - set(df.columns)
            if missing_cols:
                raise ValueError(
                    f"El archivo {csv_file} no contiene las columnas necesarias: {missing_cols}"
                )

            df_main_gas = df.loc[df["gas"] == gas, :].copy()

            mask = pd.Series(True, index=df_main_gas.index)

            # Filtrado por palabras clave en state_name
            for token in name_of_state:
                mask &= df_main_gas["state_name"].fillna("").str.contains(
                    str(token), case=False, na=False
                )

            # Filtrado energético
            if pd.notna(energy_lower_limit):
                mask &= df_main_gas["energy_eV"].fillna(-np.inf) > energy_lower_limit

            if pd.notna(energy_upper_limit):
                mask &= df_main_gas["energy_eV"].fillna(np.inf) < energy_upper_limit

            # Filtrado opcional por tipo
            if process_type is not None and pd.notna(process_type) and "type" in df_main_gas.columns:
                mask &= df_main_gas["type"].fillna("").str.lower() == str(process_type).lower()

            total_events = df_main_gas.loc[mask, "n_events"].fillna(0).sum()

            population.loc[i, name_of_output] = total_events
            population_gen.loc[i, name_of_output] = total_events

            if use_poisson_error:
                err = np.sqrt(total_events)
                population.loc[i, "Err" + name_of_output] = err
                population_gen.loc[i, "Err" + name_of_output] = err

        population = population.fillna(0)
        population.to_csv(output_dir / f"{name_of_output}.csv", index=False)
        print(f"✅ Guardado: {name_of_output}.csv")

    population_gen = population_gen.drop(columns=["file"]).fillna(0)

    # Orden opcional
    sort_cols = [c for c in ["concentration", "electric_field", "gap_mm", "pressure"] if c in population_gen.columns]
    if sort_cols:
        population_gen = population_gen.sort_values(sort_cols).reset_index(drop=True)

    if not str(output_general_name).endswith(".csv"):
        output_general_name = f"{output_general_name}.csv"

    population_gen.to_csv(output_general_name, index=False)
    print(f"✅ Guardado: {output_general_name}")

    return population_gen