from math import exp, log

from scipy.constants import mu_0, pi

# References: - https://iopscience.iop.org/article/10.1088/1742-6596/97/1/012318/pdf
#               equation (3)
#             - https://www.transform.ru/sst/$articles.en/ae000002.pdf


def calc_b_parallel(N: float, I: float, h: float, g: float = 1):
    """
    Approximates the peak value of the magnetic flux density at the middle of the coil system.

    :param N: the number of the turns in the coil (#)
    :param I: the nominal current in the coil (A)
    :param g: number of groups of balanced ampere-turns
    :param h: the length of the coil (m)
    :return:
    """

    return 2.0**0.5 * I * N * mu_0 / (g * h)


def calc_b_perpendicular(N, I, h, w, g=1):
    """
    Approximates the maximum value of the radial flux at the winding ends.

    :param N: the number of the turns in the coil (N)
    :param I: the nominal current in the coil (A)
    :param h: the length of the coil (m)
    :param w: w is the width of the tape (mm)
    :param g: number of groups of balanced ampere-turns
    :return:
    """
    return mu_0 * N * I / (2.0**0.5 * pi * g * h) * log(2.0 * h / w)


def rogowski(t_lv, t_hv, gap, ls):
    a = t_hv + t_lv + gap

    return 1 - a / pi / ls * (1 - exp(-ls / a))


if __name__ == "__main__":
    # 1.25 MVA transformer data
    # HV winding
    # l = 349 / 342.5 mm
    # Ip/Is = 69.0/1804.0
    # Nr of turns 262 - 22 x 10 turns in each pancakes, the pancekes should be parallel connected

    # lv winding
    Ns = 10.0
    Is = 1804.0
    hs = 0.3425

    Np = 262
    hp = 0.349
    Ip = 69.0

    NI_LV = Ns * Is
    NI_HV = Np * Ip

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
