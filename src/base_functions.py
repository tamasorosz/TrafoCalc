from scipy.constants import mu_0, pi

""" This module contains the basic calculations for a two winding transformer design and optimization."""
C_RHO = 8.9 * 1e-6  # kg/mm3
C_RHO_CU = 2.42  # resistivity constant in 75 C
C_RHO_FE = 7.65 * 1e-6  # kg/mm3
C_MU_0 = 4. * pi * 10 ** -7.  # Vs/Am
C_RHO_BSSCO = 6.4 * 1e-6  # kg/mm3  source: shorturl.at/rGIN7 -- sigmaaldrich.com

INFEASIBLE = -1
C_WIN_MIN = 10.  # minimal width for the windings

# constants
RHO_CU = 0.0216  # 0.0216        # ohm * mm2 / m
RHO_COPPER = 8960.0  # kg/m3


def winding_mass(m, r_m, t, h, ff):
    """
    Winding mass in m phase.

     m [#] is the number of phases
     r_m in [mm] is the mean radius of the winding.
     t in [mm] is the thickness of the winding.
     h in [mm] is the height of the winding.
     ff is the copper filling factor
     The typical values of the filling factors are .40<ff<<.70

    """

    return m * r_m * 2. * pi * t * h * ff * C_RHO


def winding_dc_loss(mass, j):
    """
    This function is estimate the loss of the winding from the geometry of
    the winding and the copper filling factor.

    - j in [A/mm2] is the current density.
    """

    dc_loss = C_RHO_CU * mass * j ** 2.

    return dc_loss * 1e-3


def core_loss_unit(ind, m_c, f_bf):
    """
    This function calculates the core loss in kW/ at the given induction
    based on mediumloss

    The calculation based on a posynomial formula:

    P_nll = M_c*f_bf*[a + c*B + d*B^b + e*B^3 + f*B^4 + g*B^4] * 10^-3

    x P_nll is the no-load loss
    x M_c is the core mass in [kg]
    x f_bf is the building factor
    x B is the induction in the column (ind)

    Fitted values (medium loss):

    """

    a = 1. * 0.0417580576
    b = 3.1
    c = 1.73506161432 * 0.0417580576
    d = 1.50521940274 * 0.0417580576
    e = 0.87054946894 * 0.0417580576
    f = 0.377614241733 * 0.0417580576
    g = 0.13103679517 * 0.0417580576

    return m_c * f_bf * (a + c * ind + d * ind ** b + e * ind ** 3. + f * ind ** 4. + g * ind ** 5.) * 10 ** (-3.)


def core_mass(r_c, ff_c, h, ei, s, m):
    """
    This function estimates the mass of a 3 legged transformer core.

    mass = m_column + m_yoke + m_corner

    x m_column - mass of columns in  [kg]
    x m_yoke   - mass of yokes in    [kg]
    x m_corner - mass of corners in  [kg]

    a = r_c**2.*pi*ff_c*rho_fe

    m_corner = a*(c*r_c*sigma**2.*gamma)
    m_column = a*o*(h+ei)
    m_yoke = a*(s*sn+m*mn)

    x h    - inner winding height in [mm]
    x ei   - end insulation bottom in[mm]
    x s    - window width in [mm]
    x sn   - number of s in a transformer shape [#]
    x m    - phase insulation [mm]
    x mn   - number of main insulation in the core [#]
    """
    gamma = 1.025

    a = r_c ** 2. * pi * ff_c * C_RHO_FE

    m_corner = a * (6.0 * r_c * gamma + 6.0 * r_c)
    m_column = a * 3 * (h + ei)
    m_yoke = a * (s * 8.0 + m * 4.0)

    return m_column + m_yoke + m_corner


def window_width(g_core, t_in, t_out, g, t_r, g_r):
    """
    is the calculated winding width if there

    x g_core - is the distance between the core and the inner winding in [mm]
    x t_in   - is the thickness of the inner winding  in [mm]
    x t_out  - is the thickness of the outer winding  in [mm]
    x g is   - is the main gap in [mm]
    x t_r    - is the width of the regulating winding in [mm]
    x g_r    - is the distance between the outer winding and the regulating winding [mm]

    g is considered as a phase distance at the end of the windings
    """

    return g_core + t_in + t_out + g + t_r + g_r + g


def turn_voltage(ind, r_c, ff_c, freq):
    """
    This function calculates the turn voltage from the core area, the frequency and the flux density in the core.
    [V]
    """
    area = r_c ** 2. * pi * ff_c

    return ind * area * 4.44 * 1e-6 * freq


def short_circuit_impedance(b_pow, p_num, freq, alpha, turn_v, h, s,
                            r_in, t_in, r_ou, t_ou, g):
    """
    Short-circuit impedance calculation

    b_pow - built-in power  [kVA]
    p_num - phase number [#]
    freq  - frequency  [Hz]
    alpha - ratio of the outer and inner winding
    ff_c  - core filling factor
    turn_v- turn voltage [V]
    h     - height of the inner window [mm]
    s     - width of working window [mm]
    g     . main insulation width [mm]
    """

    p_pow = b_pow / p_num
    imp_con = 4. * pi ** 2. * mu_0 * freq * p_pow / turn_v ** 2. / (h * (1 + alpha) / 2. + 0.32 * s)
    a = r_in * t_in / 3.
    b = r_ou * t_ou / 3.
    c = (r_in + t_in / 2. + g / 2.) * g

    return imp_con * (a + b + c)


def inner_winding_radius(r_c, g_core, t_in):
    """
    Calculates the inner winding radius in the following cases:

    Core || Inner Main || Outer Main
    Core || Inner Main || Outer Main || Regulating
    Core || Inner Main || Regulating || Outer Main

    """
    return r_c + g_core + t_in / 2.


def outer_winding_radius(r_in, t_in, g, t_out):
    """
    Calculates the inner winding radius in the following cases:

    Core || Inner Main || Outer Main
    Core || Inner Main || Outer Main || Regulating

    """
    return r_in + t_in / 2. + g + t_out / 2.


def winding_power(width, height, ff_w, j_, u_t):
    """
    Calculates the power of the winding

    x width - winding width in [mm]
    x height- winding height in [mm]
    x ff_w  - filling factor [-]
    """

    return width * height * u_t * ff_w * j_ * 1e-3  # kVA


def calc_inner_width(s_p, h_, ff_w, j_, u_t):
    """
    Calculates the inner width from the height and the power

    s_p   - phase power in kVA
    h_    - winding height
    ff_w  - winding filling factor
    j_    - current density [A/mm2]
    """

    return s_p / h_ / ff_w / j_ / u_t * 1e3  # mm


def calculate_turn_num(win_voltage, turn_voltage):
    """
    win_voltage in [kV]
    turn_voltage in [V]
    """

    return round(win_voltage / turn_voltage * 1e3, 1)


def homogenous_insulation_ff(ff):
    """
    If we assumes that the horizontal ff = vertical ff, the
    function give back the insulation horizontal filling factor
    """

    return (1. - ff) ** 0.5


def opt_win_eddy_loss(v_k, k):
    """
    This function approximates the optimal eddy loss in the assumed optimal winding system.
    From the insulation, and winding width parameters.

    Source:
    Elektrotechnika, Újházy Géza, 1969/10-11,
    Erőátviteli Transzformátorok tekercsrendszerének a méretezése

    """

    return v_k / (3. * v_k + 2. * k)


def sum_winding_loss(dc_loss, eddy_loss):
    """
    This function give back the sum of the estimated eddy and dc loss in the
    winding.
    """
    return dc_loss * (1. + eddy_loss)


def phase_current(sb, ub, con_fact):
    """

    :param sb: nominal power [MVA]
    :param ub: line voltage
    :param con_fact: connection factor --- 1 for delta --- sqrt(3) for star connection --- sqrt(2)/2. for zig-zag
    :return:
    """
    return sb * 1e3 / ub / 3. ** 0.5 / con_fact


def capitalized_cost(c_mass, c_material_price,
                     w_mass_in, w_c_in,
                     w_mass_ou, w_c_out,
                     ll, ll_cost,
                     nll, nll_cost):
    """
    Objective function
    Simple capitalized cost calculation, only the active part with filling factors
    """

    return c_mass * c_material_price + w_mass_in * w_c_in + w_mass_ou * w_c_out + ll * ll_cost + nll * nll_cost
