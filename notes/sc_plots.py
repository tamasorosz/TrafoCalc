import numpy as np

from src.superconductor_losses import perp_loss, parallel_loss, norris_equation
from bokeh.plotting import figure, output_file, show

flux_dens = np.linspace(0.01, .1, 200)

"""
This file creates plots from the superconductor losses.
"""


def plot_pp_losses():
    """Plots the parallel and the perpendicular losses, the parameters are based on:
    N. Magnussona,*, A. Wolfbrandt. AC losses in high-temperature superconducting tapes exposed
    to longitudinal magnetic fields

    Bi-2223/Ag, HTS tape produced by American Superconductor. The cross-section was 0.21 mm × 4.29 mm
    and the self-field critical current at 77 K, determined by the 1 μV/cm criterion, was 115 A.

    The losses due to perpendicular magnetic fields follow (2) with Bc=11 mT and K=1.35.
    Both the parallel and the longitudinal field losses follow (3), with B equal to B∥ or B=.
    """
    p = figure(title="Losses", y_axis_type="log", x_range=(0.01, 1), y_range=(0.01, 30))
    # background_fill_color="#fafafa")

    p.line(
        flux_dens,
        perp_loss(1.35, 50, 4.29 * 1e-3, 0.015, flux_dens),
        legend_label="perpendicular loss",
        line_color="tomato",
        line_dash="dashed",
        line_width=2.5,
    )

    f = 50.0
    c = (0.75,)
    bp = 0.0344
    Ac = 0.21 * 4.29 * 1e-6

    bpv = np.vectorize(parallel_loss)
    p.line(
        flux_dens,
        bpv(f, c, Ac, flux_dens, bp),
        legend_label="parallel losses",
        line_dash="dotted",
        line_color="indigo",
        line_width=2.5,
    )

    p.legend.location = "top_left"

    output_file("logplot.html", title="log plot example")

    show(p)


def plot_self_field_losses():
    """Plots the parallel and the perpendicular losses, the parameters are based on:
    N. Magnussona,*, A. Wolfbrandt. AC losses in high-temperature superconducting tapes exposed
    to longitudinal magnetic fields
    """

    Ic = 115.0  # A
    current = np.linspace(2.0, 114.9, 100)

    p = figure(x_range=(0.01, 1), y_range=(0.01, 0.12), x_axis_label="Ic [A]", y_axis_label="P [W/m]")

    p.line(
        current / Ic,
        norris_equation(50, current, Ic),
        legend_label="self-field losses",
        line_color="teal",
        line_dash="dashed",
        line_width=2.5,
    )

    p.legend.location = "top_left"

    # p.output_backend = "svg"
    output_file("logplot2.html", title="log plot example")
    show(p)


if __name__ == "__main__":
    plot_pp_losses()
    plot_self_field_losses()
