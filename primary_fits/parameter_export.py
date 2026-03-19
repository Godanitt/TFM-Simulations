import numpy as np
from pathlib import Path
import pandas as pd 

def export_fit_table_latex(result, names, filename, caption, label, sigfigs=4):
    """
    Exporta una tabla LaTeX con:
      - parámetro
      - valor ajustado
      - incertidumbre (1 sigma)
      - incertidumbre relativa en %

    Parámetros
    ----------
    result : scipy.optimize.OptimizeResult
        Resultado de scipy.optimize.least_squares.
    names : list[str]
        Lista de nombres LaTeX de los parámetros, en el mismo orden que result.x.
    filename : str
        Nombre del archivo .tex de salida.
    caption : str
        Caption de la tabla.
    label : str
        Label de la tabla (sin \\label{}).
    sigfigs : int, opcional
        Cifras significativas para los números dentro de \\num{}.

    Notas
    -----
    - La incertidumbre se estima a partir de la matriz de covarianza aproximada:
          cov ~ s^2 (J^T J)^(-1)
      usando pseudoinversa vía SVD.
    - Si los residuos ya están ponderados por sigma experimental,
      esta estimación suele ser razonable como aproximación local.
    - Requiere en LaTeX:
          \\usepackage{booktabs}
          \\usepackage{siunitx}
    """

    popt = np.asarray(result.x, dtype=float)

    if len(names) != len(popt):
        raise ValueError(
            f"len(names) = {len(names)} pero len(result.x) = {len(popt)}"
        )

    if not hasattr(result, "jac"):
        raise ValueError("El objeto 'result' no tiene atributo 'jac'.")

    J = np.asarray(result.jac, dtype=float)
    m, n = J.shape

    # --- Covarianza aproximada usando pseudoinversa robusta (SVD) ---
    # Similar a lo que hace scipy en curve_fit para evitar problemas numéricos.
    U, svals, VT = np.linalg.svd(J, full_matrices=False)

    if svals.size == 0:
        pcov = np.full((n, n), np.nan)
    else:
        threshold = np.finfo(float).eps * max(J.shape) * svals[0]
        mask = svals > threshold

        if not np.any(mask):
            pcov = np.full((n, n), np.nan)
        else:
            VT = VT[mask]
            svals = svals[mask]
            pcov = (VT.T / (svals**2)) @ VT

            # Escalado por varianza residual si hay grados de libertad
            dof = m - n
            if dof > 0:
                s_sq = 2.0 * result.cost / dof
                pcov = pcov * s_sq
            else:
                pcov[:] = np.nan

    perr = np.sqrt(np.clip(np.diag(pcov), 0.0, None))

    # % de incertidumbre relativa
    rel_err_percent = np.full_like(perr, np.nan, dtype=float)
    nonzero = popt != 0
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


def export_to_csv(output,result,names):
    popt = np.asarray(result.x, dtype=float)

    if len(names) != len(popt):
        raise ValueError(
            f"len(names) = {len(names)} pero len(result.x) = {len(popt)}"
        )

    if not hasattr(result, "jac"):
        raise ValueError("El objeto 'result' no tiene atributo 'jac'.")

    J = np.asarray(result.jac, dtype=float)
    m, n = J.shape

    # --- Covarianza aproximada usando pseudoinversa robusta (SVD) ---
    # Similar a lo que hace scipy en curve_fit para evitar problemas numéricos.
    U, svals, VT = np.linalg.svd(J, full_matrices=False)

    if svals.size == 0:
        pcov = np.full((n, n), np.nan)
    else:
        threshold = np.finfo(float).eps * max(J.shape) * svals[0]
        mask = svals > threshold

        if not np.any(mask):
            pcov = np.full((n, n), np.nan)
        else:
            VT = VT[mask]
            svals = svals[mask]
            pcov = (VT.T / (svals**2)) @ VT

            # Escalado por varianza residual si hay grados de libertad
            dof = m - n
            if dof > 0:
                s_sq = 2.0 * result.cost / dof
                pcov = pcov * s_sq
            else:
                pcov[:] = np.nan

    param = result.x
    perr = np.sqrt(np.clip(np.diag(pcov), 0.0, None))

        
    df = pd.DataFrame(
        {
            "parameter": param,
            "err parameter": perr
        },
        index=names
    )

    df.to_csv(output)
