import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_winding_flux(fluxes: list, z_min, z_max, label='HV', dev=False):
    """
    Plots the radial and axial flux along the windings
    :param fluxes: list of tuples which contains the (radial,axial) axial fluxes in the given points
    :param z_min: minimal axial position of the winding
    :param z_max: maximal axial position of the winding
    :return:
    """

    sns.set_theme(style="whitegrid")

    z = np.linspace(z_min, z_max, len(fluxes))
    data = pd.DataFrame(fluxes, z, columns=["Axial Flux", "Radial Flux"])
    fig, axes = plt.subplots(2, 1)
    fig.suptitle('Flux distribution in {} winding'.format(label))
    avg = data["Axial Flux"].mean()
    ## axial flux
    sns.lineplot(ax=axes[1], data=data["Axial Flux"], linewidth=2.5, palette=sns.color_palette("mako_r", 10),
                 legend="full")
    axes[1].axhline(data["Axial Flux"].mean(), color='red')
    ## radial flux
    sns.lineplot(ax=axes[0], data=data["Radial Flux"], linewidth=2.5, palette=sns.color_palette("mako_r", 6),
                 legend="full")
    axes[0].axhline(data["Radial Flux"].mean(), color='red',)

    if not dev:
        plt.show()
