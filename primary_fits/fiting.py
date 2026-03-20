import numpy as np
import scipy.optimize as opt

def fitParameters(equations, experimental_data, degrad_data, x0, bounds):

    concentration = degrad_data["concentration"]

    def residuals(x):
        res_list = []

        for key, theory_yield in equations.items():
            exp_data = experimental_data[key]

            # Solo columnas físicas, no columnas de error
            cols_phys = [
                c for c in exp_data.columns
                if not str(c).startswith("Err")
            ]

            for col in cols_phys:
                
                y_exp = exp_data[col].to_numpy(dtype=float)

                err_col_candidates = [
                    f"Err {col}", f"Err_{col}", f"{col} Err", f"{col}_Err"
                ]

                s_exp = None
                for ec in err_col_candidates:
                    if ec in exp_data.columns:
                        s_exp = exp_data[ec].to_numpy(dtype=float)
                        break

                if s_exp is None:
                    s_exp = np.ones_like(y_exp)

                # Evitar sigma = 0
                s_exp_eff = s_exp.copy()
                mask0 = (s_exp_eff == 0)
                if np.any(mask0):
                    s_exp_eff[mask0] = 1e12

                # Extraer presión de nombres tipo "1bar", "4bar", etc.
                try:
                    n_val = float(str(col).replace("bar", ""))
                except:
                    continue   # saltar columnas que no sean físicas

                y_th = theory_yield(x, degrad_data, concentration, n_val)

                if len(y_th) > len(y_exp):
                    n = len(y_th) - len(y_exp)
                    y_th = y_th[n:]
                elif len(y_exp) > len(y_th):
                    n = len(y_exp) - len(y_th)
                    y_exp = y_exp[n:]
                    s_exp_eff = s_exp_eff[n:]

                res = (y_exp - y_th) / s_exp_eff
                res_list.append(res)

        return np.concatenate(res_list)

    result = opt.least_squares(
        residuals,
        x0,
        bounds=bounds,
        method="trf",
        verbose=2
    )

    return result