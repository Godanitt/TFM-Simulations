import re
import pandas as pd
from pathlib import Path

#############################################################################################################
########################## FUNCION QUE LE UN ARCHIVO DE DEGRAD ##############################################
#############################################################################################################

def read_input(archivo_entrada, archivo_salida_Ar="collisions_Ar.csv", archivo_salida_CF4="collisions_CF4.csv"):
    # hecho por chatgpt, coregido por dvl
    def split_by_gas(block_text: str):
        # Encabezados tipo "ARGON ANISOTROPIC ...\n--------"
        headers = list(re.finditer(
            r"^\s*(?P<gas>[A-Z0-9]+)\s+ANISOTROPIC[^\n]*\n[-]{8,}\s*",
            block_text, flags=re.M
        ))
        parts = []
        for i, h in enumerate(headers):
            gas = h.group("gas").strip()
            start = h.end()
            end = headers[i + 1].start() if i + 1 < len(headers) else len(block_text)
            parts.append((gas, block_text[start:end]))
        return parts

    def parse_lines(gas: str, gas_block: str) -> pd.DataFrame:
        # Proceso [ELOSS/ELEVEL= ...] valor +- error %
        pattern = re.compile(
            r"^\s*(?P<proc>.+?)"
            r"(?:\s+(?:E(?:LEVEL|LOSS)=\s*(?P<energy>-?\d*\.?\d+(?:D[+-]?\d+)?)))?"
            r"\s+(?P<value>-?\d*\.?\d+)\s*\+\-\s*(?P<err>-?\d*\.?\d+)\s*%",
            flags=re.M
        )
        rows = []
        for m in pattern.finditer(gas_block):
            proc = re.sub(r"\s{2,}", " ", m.group("proc").strip())
            energy = m.group("energy")
            energy = float(energy.replace("D", "E")) if energy else None
            rows.append({
                "Gas": gas,
                "Proceso": proc,
                "Energia": energy,
                "Eventos": float(m.group("value")),
                "Error%": float(m.group("err")),
            })
        return pd.DataFrame(rows)

    def parse_file(path: str) -> pd.DataFrame:
        text = Path(path).read_text(encoding="utf-8", errors="ignore")
        title = "NUMBER OF COLLISIONS PER EVENT FOR EACH GAS"
        idx = text.find(title)
        if idx == -1:
            return pd.DataFrame()
        # Toma desde el título hasta el final (evita cortar en el primer “----”)
        sub = text[idx:]
        frames = [parse_lines(g, b) for g, b in split_by_gas(sub)]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    # --- EJECUCIÓN ---
    df_all = parse_file(archivo_entrada)

    if df_all.empty:
        print("No se encontró el bloque de colisiones en el archivo.")
    else:
        df_all[df_all["Gas"] == "ARGON"].to_csv(archivo_salida_Ar, index=False)
        df_all[df_all["Gas"] == "CF4"].to_csv(archivo_salida_CF4, index=False)
        print(f"Guardado: {archivo_salida_Ar} y {archivo_salida_CF4}")

    return df_all

#############################################################################################################
########################## PROGRAMA PRINCIPAL ##############################################
#############################################################################################################

########################## Inputs de Degrad ##############################################

archives=np.array(["input/output_99.9Ar_0.1CF4.txt",
                   "input/output_99.8Ar_0.2CF4.txt",
                   "input/output_99.5Ar_0.5CF4.txt",
                   "input/output_99Ar_1CF4.txt",
                   "input/output_98Ar_2CF4.txt",
                   "input/output_95Ar_5CF4.txt",
                   "input/output_90Ar_10CF4.txt",
                   "input/output_80Ar_20CF4.txt",
                   "input/output_50Ar_50CF4.txt",
                   "input/output_PureCF4.txt"])

archives_ar=np.array(["input/ar_degrad_output_99.9Ar_0.1CF4.txt",
                   "input/ar_degrad_output_99.8Ar_0.2CF4.txt",
                   "input/ar_degrad_output_99.5Ar_0.5CF4.txt",
                   "input/ar_degrad_output_99Ar_1CF4.txt",
                   "input/ar_degrad_output_98Ar_2CF4.txt",
                   "input/ar_degrad_output_95Ar_5CF4.txt",
                   "input/ar_degrad_output_90Ar_10CF4.txt",
                   "input/ar_degrad_output_80Ar_20CF4.txt",
                   "input/ar_degrad_output_50Ar_50CF4.txt",
                   "input/ar_degrad_output_PureCF4.txt"])

archives_cf4=np.array(["input/cf4_degrad_output_99.9Ar_0.1CF4.txt",
                   "input/cf4_degrad_output_99.8Ar_0.2CF4.txt",
                   "input/cf4_degrad_output_99.5Ar_0.5CF4.txt",
                   "input/cf4_degrad_output_99Ar_1CF4.txt",
                   "input/cf4_degrad_output_98Ar_2CF4.txt",
                   "input/cf4_degrad_output_95Ar_5CF4.txt",
                   "input/cf4_degrad_output_90Ar_10CF4.txt",
                   "input/cf4_degrad_output_80Ar_20CF4.txt",
                   "input/cf4_degrad_output_50Ar_50CF4.txt",
                   "input/cf4_degrad_output_PureCF4.txt"])


########################################################################

fCF4            = np.array([0.001,0.01,0.1,1])   # concentraciones

poblations_CF4 = pd.DataFrame({"fCF4":fCF4})
poblations_CF3 = pd.DataFrame({"fCF4":fCF4})
poblations_Ar_dbleStar  = pd.DataFrame({"fCF4":fCF4})
poblations_Ar_3rd  = pd.DataFrame({"fCF4":fCF4})


name_CF4         =  np.array(["all"])
name_CF3         =  np.array([">11.5",">12.5",">14.5"])
name_Ar_dbleStar =  np.array(["Ar** all", "Ar** 2exc"])
name_Ar_3rd      =  np.array(["all"])