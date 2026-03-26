import numpy as np
from pathlib import Path
import pandas as pd
import numpy as np
from pathlib import Path
import pandas as pd


def _extract_fit_info(result):
    """
    Extrae de forma robusta:
      - popt : parámetros completos
      - perr : errores completos
      - pcov : covarianza completa

    Prioridad:
      1) usar result.pcov / result.perr si existen
      2) si no, reconstruir desde result.jac (modo antiguo)
    """
    if not hasattr(result, "x"):
        raise ValueError("El objeto 'result' no tiene atributo 'x'.")

    popt = np.asarray(result.x, dtype=float)
    n_total = len(popt)

    # =========================================================
    # Caso nuevo: fitParameters ya devolvió todo expandido
    # =========================================================
    has_full_pcov = hasattr(result, "pcov")
    has_full_perr = hasattr(result, "perr")

    if has_full_pcov or has_full_perr:
        pcov = None
        perr = None

        if has_full_pcov:
            pcov = np.asarray(result.pcov, dtype=float)
            if pcov.shape != (n_total, n_total):
                raise ValueError(
                    f"result.pcov tiene forma {pcov.shape}, pero debería ser "
                    f"({n_total}, {n_total})."
                )

        if has_full_perr:
            perr = np.asarray(result.perr, dtype=float)
            if perr.shape != (n_total,):
                raise ValueError(
                    f"result.perr tiene longitud {len(perr)}, pero debería ser {n_total}."
                )

        if (perr is None) and (pcov is not None):
            perr = np.sqrt(np.clip(np.diag(pcov), 0.0, None))

        if (pcov is None) and (perr is not None):
            pcov = np.full((n_total, n_total), np.nan, dtype=float)
            np.fill_diagonal(pcov, perr**2)

        return popt, perr, pcov

    # =========================================================
    # Caso antiguo: reconstruir desde jac
    # =========================================================
    if not hasattr(result, "jac"):
        raise ValueError(
            "El objeto 'result' no tiene atributos 'pcov'/'perr' ni 'jac'."
        )

    J = np.asarray(result.jac, dtype=float)
    m, n = J.shape

    U, svals, VT = np.linalg.svd(J, full_matrices=False)

    if svals.size == 0:
        pcov = np.full((n, n), np.nan, dtype=float)
    else:
        threshold = np.finfo(float).eps * max(J.shape) * svals[0]
        mask = svals > threshold

        if not np.any(mask):
            pcov = np.full((n, n), np.nan, dtype=float)
        else:
            VT = VT[mask]
            svals = svals[mask]
            pcov = (VT.T / (svals**2)) @ VT

            dof = m - n
            if dof > 0:
                s_sq = 2.0 * result.cost / dof
                pcov = pcov * s_sq
            else:
                pcov[:] = np.nan

    perr = np.sqrt(np.clip(np.diag(pcov), 0.0, None))

    if len(perr) != n_total:
        raise ValueError(
            f"Incompatibilidad interna: len(result.x)={n_total} pero "
            f"la covarianza inferida desde jac da {len(perr)} parámetros. "
            f"Esto suele ocurrir si hay parámetros fijados y no se ha expandido "
            f"result.pcov/result.perr dentro de fitParameters."
        )

    return popt, perr, pcov


def export_fit_table_latex(result, names, filename, caption, label, sigfigs=4):
    """
    Exporta una tabla LaTeX con:
      - parámetro
      - valor ajustado
      - incertidumbre (1 sigma)
      - incertidumbre relativa en %

    Compatible con:
      - resultados antiguos de least_squares
      - resultados nuevos donde fitParameters ya devuelve
        result.x, result.perr y result.pcov completos
    """
    popt, perr, pcov = _extract_fit_info(result)

    if len(names) != len(popt):
        raise ValueError(
            f"len(names) = {len(names)} pero len(result.x) = {len(popt)}"
        )

    rel_err_percent = np.full_like(perr, np.nan, dtype=float)
    nonzero = (popt != 0) & np.isfinite(perr)
    rel_err_percent[nonzero] = 100.0 * np.abs(perr[nonzero] / popt[nonzero])

    def fmt_num(x, sig=sigfigs):
        if x is None or not np.isfinite(x):
            return r"--"
        return rf"\num{{{x:.{sig}g}}}"

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(r"\begin{tabular}{lccc}")
    lines.append(r"\toprule")
    lines.append(r"Parámetro & Valor & Incertidumbre & Incertidumbre (\%) \\")
    lines.append(r"\midrule")

    for name, val, err, rel in zip(names, popt, perr, rel_err_percent):
        lines.append(
            f"{name} & {fmt_num(val)} & {fmt_num(err)} & {fmt_num(rel)} \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    latex_table = "\n".join(lines)

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(latex_table, encoding="utf-8")

    return latex_table, popt, perr, rel_err_percent


def export_to_csv(output, result, names):
    """
    Exporta CSV con parámetros y errores.

    Compatible con parámetros fijados siempre que fitParameters haya dejado
    result.x/result.perr/result.pcov ya expandidos al espacio completo.
    """
    popt, perr, pcov = _extract_fit_info(result)

    if len(names) != len(popt):
        raise ValueError(
            f"len(names) = {len(names)} pero len(result.x) = {len(popt)}"
        )

    df = pd.DataFrame(
        {
            "parameter": popt,
            "err parameter": perr
        },
        index=names
    )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output)

    return df

def export_fit_table_typst(
    result,
    names,
    filename,
    caption=None,
    label=None,
    sigfigs=4
):
    """
    Exporta una tabla Typst con:
      - parámetro
      - valor ajustado
      - incertidumbre (1 sigma)
      - incertidumbre relativa en %

    Compatible con:
      - resultados antiguos de least_squares
      - resultados nuevos donde fitParameters ya devuelve
        result.x, result.perr y result.pcov completos

    Parámetros
    ----------
    result : OptimizeResult
    names : list[str]
        Nombres de parámetros en el mismo orden que result.x.
        Pueden contener sintaxis Typst, por ejemplo `$alpha$`.
    filename : str
        Fichero .typ de salida.
    caption : str | None
        Caption de la tabla.
    label : str | None
        Label opcional. Se añade como `<label>`.
    sigfigs : int
        Cifras significativas.
    """
    popt, perr, pcov = _extract_fit_info(result)

    if len(names) != len(popt):
        raise ValueError(
            f"len(names) = {len(names)} pero len(result.x) = {len(popt)}"
        )

    rel_err_percent = np.full_like(perr, np.nan, dtype=float)
    nonzero = (popt != 0) & np.isfinite(perr)
    rel_err_percent[nonzero] = 100.0 * np.abs(perr[nonzero] / popt[nonzero])

    def fmt_num(x, sig=sigfigs):
        if x is None or not np.isfinite(x):
            return '"-"'
        return f'"{x:.{sig}g}"'

    def fmt_text(s):
        s = str(s)
        # Si el usuario ya mete algo como $...$, lo dejamos tal cual para math mode.
        # Si no, lo ponemos como string de Typst.
        if s.startswith("$") and s.endswith("$"):
            return s
        return f'"{s}"'

    rows = []
    for name, val, err, rel in zip(names, popt, perr, rel_err_percent):
        rows.append(
            f"  [{fmt_text(name)}], [{fmt_num(val)}], [{fmt_num(err)}], [{fmt_num(rel)}],"
        )

    lines = []

    label_suffix = f" <{label}>" if label else ""
    caption_text = caption if caption is not None else "Fit parameters"

    lines.append(f'#figure{label_suffix}(')
    lines.append("  table(")
    lines.append("    columns: 4,")
    lines.append("    align: (left, center, center, center),")
    lines.append("    stroke: none,")
    lines.append("    table.header(")
    lines.append('      [*Parámetro*],')
    lines.append('      [*Valor*],')
    lines.append('      [*Incertidumbre*],')
    lines.append('      [*Incertidumbre (%)*],')
    lines.append("    ),")
    lines.extend(rows)
    lines.append("  ),")
    lines.append(f'  caption: [{caption_text}],')
    lines.append(")")

    typst_table = "\n".join(lines)

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(typst_table, encoding="utf-8")

    return typst_table, popt, perr, rel_err_percent