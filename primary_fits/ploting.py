import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import scienceplots
plt.style.use('default')
plt.style.use('science')
plt.rcParams.update({
    "font.family": "serif",   # specify font family here
    "font.serif": ["Times"],  # specify font here
    "font.size": 11})          # specify font size here


def darken_color(color, factor=0.65):
    """
    factor > 0  -> oscurece
    factor < 0  -> aclara
    factor = 0  -> no cambia
    """
    rgb = np.array(mcolors.to_rgb(color), dtype=float)

    if factor > 0:
        new_rgb = rgb * (1 - factor)
    elif factor < 0:
        a = -factor
        new_rgb = rgb + (1 - rgb) * a
    else:
        new_rgb = rgb

    return tuple(np.clip(new_rgb, 0, 1))


def plot_fit_vs_experiment_by_pressure(
    df_exp,
    theory_func,
    fit_params,
    degrad_data,
    concentration_grid,
    *,
    x_col="fCF4",
    pressure_cols=None,
    pressures=None,
    pressure_regex=r"^\s*([0-9]+(?:\.[0-9]+)?)\s*bar\s*$",
    err_patterns=None,
    x_plot_factor=100.0,
    min_positive_x=None,
    cmap="viridis",
    darken_factor=0.65,
    figsize=(6, 4),
    title=None,
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    xscale="log",
    yscale="log",
    line_label_fmt=["{p:g} bar fit"],
    exp_label_fmt="{p:g} bar exp",
    legend=True,
    legend_kwargs=None,
    output=None,
    show=True,
    ax=None,
    activate_components = False,
):
    if err_patterns is None:
        err_patterns = [
            "Err {col}",
            "Err_{col}",
            "{col} Err",
            "{col}_Err",
        ]

    concentration_grid = np.asarray(concentration_grid, dtype=float)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # x experimental
    x_exp = df_exp[x_col].to_numpy(dtype=float).copy()
    if xscale == "log" and min_positive_x is not None:
        x_exp[x_exp <= 0] = min_positive_x

    x_exp_plot = x_exp 
    x_grid_plot = concentration_grid * x_plot_factor

    # detectar columnas de presión automáticamente
    if pressure_cols is None:
        pressure_cols = []
        regex = re.compile(pressure_regex)

        for col in df_exp.columns:
            col_str = str(col)

            if col_str == x_col:
                continue

            m = regex.match(col_str)
            if m:
                pressure_cols.append((float(m.group(1)), col_str))

        pressure_cols.sort(key=lambda t: t[0])

    if not pressure_cols:
        raise ValueError("No se encontraron columnas de presión válidas.")

    # filtrar solo las presiones pedidas
    if pressures is not None:
        pressures_set = {float(p) for p in pressures}
        pressure_cols = [(p, col) for p, col in pressure_cols if p in pressures_set]

        if not pressure_cols:
            raise ValueError(
                f"Ninguna de las presiones pedidas {pressures} está en los datos."
            )

    # colores
    cmap_obj = plt.get_cmap(cmap)
    colors = cmap_obj(np.linspace(0.15, 0.85, len(pressure_cols)))

    for (p, col), color in zip(pressure_cols, colors):
        point_color = darken_color(color, factor=darken_factor)

        # teoría

        y_fit = np.asarray(
            theory_func(fit_params, degrad_data, concentration_grid, p, activate_components = activate_components),
            dtype=float
        )

        linestyles = [
            "--",
            "-.",
            ":",
            (0, (1, 1)),
            (0, (5, 1)),
            (0, (3, 1, 1, 1)),
            (0, (5, 2, 1, 2)),
            (0, (10, 2)),
            (0, (3, 5, 1, 5)),
        ]
        
        if activate_components:
            for i, y in enumerate(y_fit):
                if i == 0:
                    ax.plot(
                        x_grid_plot,
                        y,
                        color=darken_color(color,0.3),
                        lw=2,
                        label=line_label_fmt[0].format(p=p)
                    )
                elif i<10: 
                    ax.plot(
                        x_grid_plot,
                        y,
                        color=darken_color(color,0.3),
                        linestyle = linestyles[i],
                        lw=2,
                        label=line_label_fmt[i].format(p=p)
                    )
        else: 
            ax.plot(
                x_grid_plot,
                y_fit,
                color=darken_color(color,0.3),
                lw=2,
                label=line_label_fmt[0].format(p=p)
            )


        # error experimental
        err_col = None
        for pattern in err_patterns:
            candidate = pattern.format(col=col)
            if candidate in df_exp.columns:
                err_col = candidate
                break

        y_exp = df_exp[col].to_numpy(dtype=float)
        yerr = df_exp[err_col].to_numpy(dtype=float) if err_col is not None else None

        ax.errorbar(
            x_exp_plot,
            y_exp,
            yerr=yerr,
            fmt="o",
            ms=4,
            color=point_color,
            ecolor=point_color,
            elinewidth=1,
            capsize=2,
            label=exp_label_fmt.format(p=p)
        )

    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    if legend:
        if legend_kwargs is None:
            legend_kwargs = {}
        ax.legend(**legend_kwargs)

    if output is not None:
        fig.savefig(output, bbox_inches="tight")

    if show:
        plt.show()

    return fig, ax, pressure_cols