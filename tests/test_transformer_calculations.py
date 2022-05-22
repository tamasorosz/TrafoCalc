from unittest import TestCase

from src.base_functions import (
    calc_inner_width,
    calculate_turn_num,
    core_loss_unit,
    core_mass,
    homogenous_insulation_ff,
    inner_winding_radius,
    opt_win_eddy_loss,
    outer_winding_radius,
    short_circuit_impedance,
    sum_winding_loss,
    turn_voltage,
    winding_dc_loss,
    winding_mass,
    winding_power,
    window_width,
    sc_factor,
    sc_current
)


class TestFunctions(TestCase):
    def test_winding_mass(self):
        # reference: Karsai et al, Nagytranszformátorok - Hu example page 86. (example 3.1) 6.3 MVA transformer
        # the mass corrigated due to the high thickness of the coils
        m = 3.0  # phase
        r_m = 437.0 / 2.0  # [mm]
        ff = 0.5  # [-]
        h = 1202.0  # [mm]
        t = 37.0  # [mm]
        assert round(winding_mass(m, r_m, t, h, ff), 1) == 849.6  # 815.1  kg

    def test_winding_dc_loss(self):
        # reference: Karsai et al, Nagytranszformátorok-Hu example page 86.
        # the mass of the winding corrigated due to the high thickness
        m = 3.0
        r_m = 437.0 / 2.0  # [mm]
        ff = 0.5  # [-]
        h = 1202.0  # [mm]
        t = 37.0  # [mm]
        j = 3.02  # [A/mm2]

        assert round(winding_dc_loss(winding_mass(m, r_m, t, h, ff), j), 3) == 18.8  # [kW]

    # ----------------------------------------------------------------------
    def test_core_loss_unit(self):
        ind = 1.71  # [T]
        m_c = 24900.0  # [kg]
        f_bf = 1.2

        assert round(core_loss_unit(ind, m_c, f_bf), 1) == 26.7

    def test_core_mass_on_40mva_transformer(self):
        r_c = 669.5 / 2.0
        ff_c = 0.91
        s = 258.5
        h = 1140.5
        ei = 170.0
        m = 24.0

        assert round(core_mass(r_c, ff_c, h, ei, s, m), 0) == 24906.0

    def test_core_mass_10mva(self):
        # 10 MVA trafo
        r_c = 401.0 / 2.0
        ff_c = 0.89
        s = 178.0
        h = 783.0
        ei = 150.0
        m = 20.0

        assert round(core_mass(r_c, ff_c, h, ei, s, m), 0) == 5795.0

    def test_window_width(self):
        g_core = 17.0
        t_in = 54.5
        t_out = 57.5
        g = 20.0
        t_r = 10.5
        g_r = 20.0

        assert round(window_width(g_core, t_in, t_out, g, t_r, g_r), 0) == 200

    def test_turn_voltage(self):
        r_c = 200.3
        ff_c = 0.90
        ind = 1.73
        freq = 50.0

        assert round(turn_voltage(ind, r_c, ff_c, freq), 1) == 43.6

    def test_short_circuit_impedance(self):
        b_pow = 10000.0
        p_num = 3.0
        freq = 50.0
        alpha = 1.0
        t_in = 50.5
        r_in = 483.0 / 2.0
        t_ou = 52.5
        r_ou = 312.0
        g = 20.0
        h = 736.0
        s = 170.0
        turn_v = 43.015

        assert round(short_circuit_impedance(b_pow, p_num, freq, alpha, turn_v, h, s, r_in, t_in, r_ou, t_ou, g),
                     3) == 8.51

    def test_inner_winding_radius(self):
        g_core = 17.0
        r_c = 495.0 / 2.0
        t_in = 86.0

        assert inner_winding_radius(r_c, g_core, t_in) == 615.0 / 2

    def test_outer_winding_radius(self):
        r_in = 615 / 2.0
        t_in = 86.0
        g = 20.0
        t_out = 98.0

        assert outer_winding_radius(r_in, t_in, g, t_out) == 839.0 / 2.0

    def test_winding_power(self):
        j_ = 1.423
        height = 586.0
        width = 86.0
        u_t = 66.4
        ff_w = 0.7

        assert round(winding_power(width, height, ff_w, j_, u_t), 0) == 3333.0

    def test_inner_width(self):
        j_ = 1.423
        h_ = 586.0
        width = 86.0
        u_t = 66.4
        ff_w = 0.7
        s_p = 3333.0

        assert round(calc_inner_width(s_p, h_, ff_w, j_, u_t), 0) == width

    def test_turn_num(self):
        win_voltage = 100.0
        turn_voltage = 100.0

        assert calculate_turn_num(win_voltage, turn_voltage) == 1000.0

    def test_homogenous_insulation_ff(self):
        assert round(homogenous_insulation_ff(0.5), 3) == 0.707

    def test_opt_win_eddy_loss(self):
        ff = 0.5  # [-]
        t = 37.0  # [mm]

        assert round(opt_win_eddy_loss(t * homogenous_insulation_ff(ff), t), 2) == 0.09

    def test_sum_winding_loss(self):
        # Karsai et al Nagytranszformátorok-Hu example page 86.
        m = 3.0
        r_m = 437.0 / 2.0  # [mm]
        ff = 0.5  # [-]
        h = 1202.0  # [mm]
        t = 37.0  # [mm]
        j = 3.02  # [A/mm2]

        assert (
                round(
                    sum_winding_loss(
                        winding_dc_loss(winding_mass(m, r_m, t, h, ff), j),
                        opt_win_eddy_loss(t * homogenous_insulation_ff(ff), t)
                    ),
                    3,
                )
                == 19.544
        )

    def test_sc_factor(self):
        self.assertAlmostEqual(sc_factor(15, 1), 2.57, 2)

    def test_sc_current(self):
        # test example from karsai: Nagytranszformátorok p 135, bit lower than the given,
        # because the kappa is not rounded

        self.assertAlmostEqual(sc_current(120, 0.09, 0.09, 0.05), 3234)
