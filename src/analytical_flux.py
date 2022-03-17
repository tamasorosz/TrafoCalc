import typing
from math import exp, log

from scipy.constants import mu_0, pi

"""
 This file contains analytical flux density calculatoins to estimate the maximum value for the axial and the radial 
 losses of the different windings.
 
 References: - https://iopscience.iop.org/article/10.1088/1742-6596/97/1/012318/pdf
               equation (3)
             - https://www.transform.ru/sst/$articles.en/ae000002.pdf
"""


def calc_b_parallel(N: float, I: float, h: float, g: float = 1) -> typing.Any:
    """
    Approximates the peak value of the magnetic flux density at the middle of the coil system.

    :param N: the number of the turns in the coil (#)
    :param I: the nominal current in the coil (A)
    :param g: number of groups of balanced ampere-turns
    :param h: the length of the coil (m)
    :return:
    """

    return 2.0 ** 0.5 * I * N * mu_0 / (g * h)


def calc_b_perpendicular(N: float, I: float, h: float, w: float, g: float = 1) -> typing.Any:
    """
    Approximates the maximum value of the radial flux at the winding ends.

    :param N: the number of the turns in the coil (N)
    :param I: the nominal current in the coil (A)
    :param h: the length of the coil (m)
    :param w: w is the width of the tape (mm)
    :param g: number of groups of balanced ampere-turns
    :return: the maximum value of the perpendicular magnetic field at the winding ends
    """
    return mu_0 * N * I / (2.0 ** 0.5 * pi * g * h) * log(2.0 * h / w)


def rogowski(t_lv: float, t_hv: float, gap: float, ls: float) -> typing.Any:
    """
    Calculates the rogowski factor for a more accurate analytical prediction of the magnetic field.
    :param t_lv: thickness of the low voltage winding
    :param t_hv: thickness of the high voltage winding
    :param gap: main insulation distance between the low voltage and the high voltage windings
    :param ls: length of the drop channel
    :return: rogowski factor
    """
    a = t_hv + t_lv + gap

    return 1 - a / pi / ls * (1 - exp(-ls / a))


def calc_current_density(Nt: float, height: float, thickness: float, i_ph: float):
    """
    Calculates the average current density of a superconducting winding, if the filling factor considered constant.

    :param Nt: number of active turns in the winding [#]
    :param height: height of the winding [mm]
    :param thickness: winding width [mm]
    :param i_ph: phase_current [A]

    :return:
    """

    return Nt * i_ph / (height * thickness)


if __name__ == "__main__":
    # 1.25 MVA transformer data
    # HV winding
    # l = 349 / 342.5 mm
    # Ip/Is = 69.0/1804.0
    # Nr of turns 262 - 22 x 10 turns in each pancakes, the pancekes should be parallel connected

    # lv winding
    Ns = 10.0
    Is = 1804.0/2.
    hs = 0.3425

    Np = 262
    hp = 0.355
    Ip = 69.0/2.

    NI_LV = Ns * Is
    NI_HV = Np * Ip

    print('current density lv: ', calc_current_density(Np, 342.5, 8, 34.5))
    print('current density hv: ', calc_current_density(Ns, 355.0, 13.5, 902))

    print("jlv ", NI_LV / (13 * hs * 1e3))
    print("jhv ", NI_HV / (8 * hp * 1e3))
    print("rogowski", rogowski(0.013, 0.08, 0.034, hs))
    print("LV winding: ", NI_LV, "HV:", NI_HV)
    print("axial flux (LV)", calc_b_parallel(Ns, Is, hs) * 1e3 * rogowski(0.013, 0.08, 0.034, hs), "mT")
    print("axial flux (HV)", calc_b_parallel(Np, Ip, hp) * 1e3, "mT")

    print("radial flux (LV)", calc_b_perpendicular(Ns, Is, hs, 0.013) * 1e3, "mT")

    # measured maximum of axial flux density --- 630 kVA transformer
    # Ip/Is = 34.64/909.33
    # 262 / 10
    Ns = 10.0
    Is = 909.33
    hs = 0.355
    print(
        "axial flux (LV)",
        calc_b_parallel(Ns, Is, hs) * 1e3,
        "mT",
        " diff:",
        65.8 - calc_b_parallel(Ns, Is, hs) * 1e3 / 65.8 * 100,
    )
