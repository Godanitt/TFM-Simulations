import numpy as np
import scipy.optimize as opt


def fitParameters(
    equations,
    experimental_data,
    degrad_data,
    x0,
    bounds,
    is_infrared=False,
    fixed_idx=None,
    fixed_values=None,
    fixed_error=np.nan,   # usa 0.0 si prefieres que los fijados salgan con error 0
):
    """
    Ajusta parámetros con posibilidad de fijar algunos índices.

    Parámetros
    ----------
    fixed_idx : list[int] o None
        Índices de los parámetros que se quieren fijar.
    fixed_values : list[float] o None
        Valores a fijar en esos índices. Si es None, se usan los de x0.
    fixed_error : float
        Error asignado a los parámetros fijados en result.perr.
        Recomendado: np.nan. Si quieres que salgan como 0, usa 0.0.

    Devuelve
    --------
    result : OptimizeResult enriquecido con:
        - result.x         -> parámetros completos
        - result.x_free    -> parámetros libres
        - result.perr      -> errores completos
        - result.perr_free -> errores solo de libres
        - result.pcov      -> covarianza completa
        - result.pcov_free -> covarianza solo de libres
        - result.free_idx
        - result.fixed_idx
        - result.n_full
    """

    concentration = degrad_data["concentration"]

    x0 = np.asarray(x0, dtype=float).copy()
    lb = np.asarray(bounds[0], dtype=float).copy()
    ub = np.asarray(bounds[1], dtype=float).copy()

    if x0.ndim != 1:
        raise ValueError("x0 debe ser un array 1D.")
    if lb.shape != x0.shape or ub.shape != x0.shape:
        raise ValueError("bounds debe tener la misma longitud que x0.")

    n_full = len(x0)

    # -----------------------------
    # Gestión de parámetros fijados
    # -----------------------------
    if fixed_idx is None:
        fixed_idx = []

    fixed_idx = np.array(sorted(set(fixed_idx)), dtype=int)

    if fixed_idx.size > 0:
        if np.any(fixed_idx < 0) or np.any(fixed_idx >= n_full):
            raise ValueError("Hay índices en fixed_idx fuera del rango de x0.")

    if fixed_values is not None:
        fixed_values = np.asarray(fixed_values, dtype=float)
        if len(fixed_values) != len(fixed_idx):
            raise ValueError("fixed_values debe tener la misma longitud que fixed_idx.")
        x0[fixed_idx] = fixed_values

    free_mask = np.ones(n_full, dtype=bool)
    free_mask[fixed_idx] = False
    free_idx = np.where(free_mask)[0]

    # Comprobar que x0 está dentro de bounds
    if np.any(x0 < lb) or np.any(x0 > ub):
        raise ValueError("Algún valor de x0 está fuera de los bounds.")

    # Subespacio libre
    x0_free = x0[free_idx]
    lb_free = lb[free_idx]
    ub_free = ub[free_idx]

    def build_full_x(x_free):
        x_full = x0.copy()
        x_full[free_idx] = x_free
        return x_full

    def residuals_from_full_x(x, is_infrared=is_infrared):
        res_list = []

        for key, theory_yield in equations.items():
            exp_data = experimental_data[key]

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

                s_exp_eff = s_exp.copy()
                mask0 = (s_exp_eff == 0)
                if np.any(mask0):
                    s_exp_eff[mask0] = 1e12

                try:
                    n_val = float(str(col).replace("bar", ""))
                except Exception:
                    continue

                y_th = theory_yield(x, degrad_data, concentration, n_val)

                if (len(y_th) > len(y_exp)) and (not is_infrared):
                    n = len(y_th) - len(y_exp)
                    y_th = y_th[n:]
                elif (len(y_th) > len(y_exp)) and is_infrared:
                    n = len(y_th) - len(y_exp)
                    y_th = y_th[:-n]
                elif len(y_exp) > len(y_th):
                    n = len(y_exp) - len(y_th)
                    y_exp = y_exp[n:]
                    s_exp_eff = s_exp_eff[n:]

                res = (y_exp - y_th) / s_exp_eff
                res_list.append(res)

        if not res_list:
            return np.array([], dtype=float)

        return np.concatenate(res_list)

    def residuals(x_free, is_infrared=is_infrared):
        x_full = build_full_x(x_free)
        return residuals_from_full_x(x_full, is_infrared=is_infrared)

    # ---------------------------------------------------------
    # Caso extremo: todos los parámetros fijados, no se ajusta
    # ---------------------------------------------------------
    if len(free_idx) == 0:
        fun = residuals_from_full_x(x0, is_infrared=is_infrared)
        cost = 0.5 * np.dot(fun, fun)

        result = opt.OptimizeResult()
        result.x = x0.copy()
        result.x_free = np.array([], dtype=float)
        result.fun = fun
        result.cost = cost
        result.success = True
        result.status = 0
        result.message = "No hay parámetros libres: todos los parámetros están fijados."
        result.nfev = 1
        result.njev = 0
        result.jac = np.empty((len(fun), 0), dtype=float)
        result.jac_free = result.jac.copy()
        result.jac_full = np.zeros((len(fun), n_full), dtype=float)

        result.pcov_free = np.empty((0, 0), dtype=float)
        result.perr_free = np.array([], dtype=float)

        result.pcov = np.full((n_full, n_full), np.nan, dtype=float)
        result.perr = np.full(n_full, fixed_error, dtype=float)

        result.free_idx = free_idx
        result.fixed_idx = fixed_idx
        result.n_full = n_full

        return result

    # -----------------
    # Ajuste de verdad
    # -----------------
    
    result = opt.least_squares(
        residuals,
        x0_free,
        bounds=(lb_free, ub_free),
        method="trf",
        verbose=2
    )

    # =========================================================
    # Reconstrucción obligatoria al espacio completo
    # =========================================================
    
    x_free_opt = np.asarray(result.x, dtype=float).copy()
    x_full_opt = build_full_x(x_free_opt)

    J_free = np.asarray(result.jac, dtype=float)
    m, n_free = J_free.shape

    # =========================================================
    # Covarianza en el subespacio libre
    # =========================================================
    try:
        U, s, VT = np.linalg.svd(J_free, full_matrices=False)

        if s.size == 0:
            pcov_free = np.full((n_free, n_free), np.nan, dtype=float)
        else:
            threshold = np.finfo(float).eps * max(J_free.shape) * s[0]
            mask = s > threshold
            s = s[mask]
            VT = VT[mask]

            if s.size == 0:
                pcov_free = np.full((n_free, n_free), np.nan, dtype=float)
            else:
                JTJ_inv = (VT.T / (s ** 2)) @ VT

                dof_fit = m - n_free
                if dof_fit > 0:
                    s2 = 2.0 * result.cost / dof_fit
                    pcov_free = s2 * JTJ_inv
                else:
                    pcov_free = np.full((n_free, n_free), np.nan, dtype=float)

    except np.linalg.LinAlgError:
        pcov_free = np.full((n_free, n_free), np.nan, dtype=float)

    perr_free = np.sqrt(np.clip(np.diag(pcov_free), 0.0, None))

    # =========================================================
    # Expandir al espacio completo
    # =========================================================
    pcov_full = np.full((n_full, n_full), np.nan, dtype=float)
    pcov_full[np.ix_(free_idx, free_idx)] = pcov_free

    perr_full = np.full(n_full, fixed_error, dtype=float)
    perr_full[free_idx] = perr_free

    jac_full = np.zeros((m, n_full), dtype=float)
    jac_full[:, free_idx] = J_free

    # =========================================================
    # Estadísticos globales
    # =========================================================
    chi2 = 2.0 * result.cost
    N_res = result.fun.size
    N_free = n_free
    N_total = n_full
    dof = N_res - N_free
    chi2_red = chi2 / dof if dof > 0 else np.nan

    # =========================================================
    # Guardar todo en result
    # =========================================================
    result.x_free = x_free_opt
    result.x = x_full_opt

    result.jac_free = J_free
    result.jac_full = jac_full

    result.pcov_free = pcov_free
    result.perr_free = perr_free

    result.pcov = pcov_full
    result.perr = perr_full

    result.free_idx = free_idx
    result.fixed_idx = fixed_idx
    result.n_free = N_free
    result.n_total = N_total

    result.chi2 = chi2
    result.dof = dof
    result.chi2_red = chi2_red

    # =========================================================
    # Comprobaciones duras: que falle AQUÍ y no luego en export
    # =========================================================
    if len(result.x) != n_full:
        raise RuntimeError(
            f"fitParameters: result.x tiene longitud {len(result.x)} "
            f"pero debería tener {n_full}."
        )

    if result.pcov.shape != (n_full, n_full):
        raise RuntimeError(
            f"fitParameters: result.pcov tiene forma {result.pcov.shape} "
            f"pero debería ser ({n_full}, {n_full})."
        )

    if len(result.perr) != n_full:
        raise RuntimeError(
            f"fitParameters: result.perr tiene longitud {len(result.perr)} "
            f"pero debería tener {n_full}."
        )

    return result