from unittest import TestCase
from src.transformer_calculations import winding_mass, winding_dc_loss, core_loss_unit, core_mass, window_width, turn_voltage, \
    short_circuit_impedance, inner_winding_radius, outer_winding_radius, winding_power, calc_inner_width, \
    calculate_turn_num, homogenous_insulation_ff, opt_win_eddy_loss, sum_winding_loss


class TestFunctions(TestCase):

    def test_winding_mass(self):
        # reference: Karsai et al, Nagytranszformátorok - Hu example page 86. (example 3.1) 6.3 MVA transformer
        m = 3.  # phase
        r_m = 437. / 2.  # [mm]
        ff = 0.5  # [-]
        h = 1202.  # [mm]
        t = 37.  # [mm]
        assert round(winding_mass(m, r_m, t, h, ff), 1) == 815.1  # kg

    def test_winding_dc_loss(self):
        # reference: Karsai et al, Nagytranszformátorok-Hu example page 86.
        m = 3.
        r_m = 437. / 2.  # [mm]
        ff = 0.5  # [-]
        h = 1202.  # [mm]
        t = 37.  # [mm]
        j = 3.02  # [A/mm2]

        assert round(winding_dc_loss(winding_mass(m, r_m, t, h, ff), j), 3) == 17.991  # [W]

    # ----------------------------------------------------------------------
    def test_core_loss_unit(self):
        ind = 1.71  # [T]
        m_c = 24900.  # [kg]
        f_bf = 1.2

        assert round(core_loss_unit(ind, m_c, f_bf), 1) == 26.7

    def test_core_mass_on_40mva_transformer(self):
        r_c = 669.5 / 2.
        ff_c = 0.91
        s = 258.5
        h = 1140.5
        ei = 170.
        m = 24.

        assert round(core_mass(r_c, ff_c, h, ei, s, m), 0) == 24906.

    def test_core_mass_10mva(self):
        # 10 MVA trafo
        r_c = 401. / 2.
        ff_c = 0.89
        s = 178.
        h = 783.
        ei = 150.
        m = 20.

        assert round(core_mass(r_c, ff_c, h, ei, s, m), 0) == 5795.

    def test_window_width(self):
        g_core = 17.
        t_in = 54.5
        t_out = 57.5
        g = 20.
        t_r = 10.5
        g_r = 20.

        assert round(window_width(g_core, t_in, t_out, g, t_r, g_r), 0) == 200

    def test_turn_voltage(self):
        r_c = 200.3
        ff_c = 0.90
        ind = 1.73
        freq = 50.

        assert round(turn_voltage(ind, r_c, ff_c, freq), 1) == 43.6

    def test_short_circuit_impedance(self):
        b_pow = 10000.
        p_num = 3.
        freq = 50.
        alpha = 1.
        t_in = 50.5
        r_in = 483. / 2.
        t_ou = 52.5
        r_ou = 312.
        g = 20.
        h = 736.
        s = 170.
        turn_v = 43.015

        assert round(short_circuit_impedance(b_pow, p_num, freq, alpha, turn_v, h, s, r_in, t_in, r_ou, t_ou, g),
                     3) == 0.085

    def test_inner_winding_radius(self):
        g_core = 17.
        r_c = 495. / 2.
        t_in = 86.

        assert inner_winding_radius(r_c, g_core, t_in) == 615. / 2

    def test_outer_winding_radius(self):
        r_in = 615 / 2.
        t_in = 86.
        g = 20.
        t_out = 98.

        assert outer_winding_radius(r_in, t_in, g, t_out) == 839. / 2.

    def test_winding_power(self):
        j_ = 1.423
        height = 586.
        width = 86.
        u_t = 66.4
        ff_w = 0.7

        assert round(winding_power(width, height, ff_w, j_, u_t), 0) == 3333.

    def test_inner_width(self):
        j_ = 1.423
        h_ = 586.
        width = 86.
        u_t = 66.4
        ff_w = 0.7
        s_p = 3333.

        assert round(calc_inner_width(s_p, h_, ff_w, j_, u_t), 0) == width

    def test_turn_num(self):
        win_voltage = 100.
        turn_voltage = 100.

        assert calculate_turn_num(win_voltage, turn_voltage) == 1000.

    def test_homogenous_insulation_ff(self):
        assert round(homogenous_insulation_ff(0.5), 3) == 0.707

    def test_opt_win_eddy_loss(self):
        ff = 0.5  # [-]
        t = 37.  # [mm]

        assert round(opt_win_eddy_loss(t * homogenous_insulation_ff(ff), t), 2) == 0.17

    def test_sum_winding_loss(self):
        # Karsai et al Nagytranszformátorok-Hu example page 86.
        m = 3.
        r_m = 437. / 2.  # [mm]
        ff = 0.5  # [-]
        h = 1202.  # [mm]
        t = 37.  # [mm]
        j = 3.02  # [A/mm2]

        assert round(sum_winding_loss(winding_dc_loss(winding_mass(m, r_m, t, h, ff), j),
                                      opt_win_eddy_loss(t * homogenous_insulation_ff(ff), t)), 3) == 21.077
