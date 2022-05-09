from unittest import TestCase
from src.superconductor_losses import parallel_loss, perp_loss, norris_equation, cryostat_losses, cryo_surface, \
    thermal_incomes, cooler_cost, sc_load_loss, magnusson_ac_loss

from math import pi


class TestLosses(TestCase):

    def test_losses(self):
        # BSCCO cable from Magnussons' paper
        # f = 50.0
        # c = 0.75
        # bp = 0.033
        # Ac = 0.31 * 4.29 * 1e-6

        bpar = 66.0 * 1e-3  # T

        # parallel losses
        # 0.133 W/m in 10.1109/TASC.2003.813123 but this  value seems ok
        self.assertAlmostEqual(parallel_loss(50, bpar), 0.177, 3)
        # perpendicular losses
        self.assertAlmostEqual(perp_loss(50, 68. * 1e-3), 2.009, 3)
        # norris equation
        self.assertAlmostEqual(norris_equation(50, 50, 115.), 0.0047, 4)

    def test_sc_transformer_losses(self):
        # Th example is coming from DOI: 10.1109/TASC.2003.813123
        # BSCCO cable from Magnussons' paper
        # f = 50.0
        # c = 0.75
        # bp = 0.033
        # Ac = 0.31 * 4.29 * 1e-6

        bpar = 0.066  # T
        bperp = 0.068  # T

        self.assertAlmostEqual(parallel_loss(50, bpar), 0.179, 2)
        self.assertAlmostEqual(perp_loss(50, bperp), 2.009, 2)  # 1.960

        # magnusson-formula based loss - 2 W/m
        self.assertAlmostEqual(magnusson_ac_loss(bpar, bperp, 50, 25), 2.125, 2)

    def test_45kVA_transformer_loss(self):
        """
        Source: DOI: 10.1109/TASC.2005.869712
        45 kVA transformer
        :return:
        """

        # The maximum axial and radial fields 30 mT, 23.7 mT
        bax = 30. * 1e-3
        brad = 23.7 * 1e-3
        i1 = 18.75  # A
        i2 = 281.25  # A

        ac = magnusson_ac_loss(bax, brad, 50, 18.75)

        self.assertAlmostEqual(ac, 0.279, 2)

        # hv length - 14 turns, 13 pancakes
        l_hv = 14 * 13 * (365 + 375.5) / 2 * 3.14 * 1e-3  # m
        l_lv = 210 * (286 + 304.5) / 2 * 3.14 * 1e-3  # m

        l = l_lv + l_hv

        self.assertAlmostEqual(ac * l, 113.5, 1)

    def test_losses_second_case_small_field(self):
        bpar = 20.0 * 1e-3  # mT

        # parallel losses
        self.assertAlmostEqual(parallel_loss(50, bpar), 0.0535, 3)
        # perpendicular losses
        self.assertAlmostEqual(perp_loss(50, 20. * 1e-3), 0.1625, 1)

    def test_cryostat_losses(self):
        r_in = 0.5 / pi
        r_ou = 1.0 / pi
        h = 1.0
        a_cs = cryo_surface(r_in, r_ou, h)  # [m2]
        P_cryo = cryostat_losses(a_cs)
        P_thermal = thermal_incomes(100, 100)

        sc_ll = sc_load_loss(100, 100, 100)

        self.assertAlmostEqual(a_cs, 3.238, 2)
        # self.assertAlmostEqual(P_cryo, 0.029537, 2)
        self.assertAlmostEqual(P_thermal, 54, 1)
        self.assertAlmostEqual(sc_ll, 5400, 1)
        # tco = modified_tco_evaluation(1, 1, 1, 0, 0, 1, 0)
        # self.assertAlmostEqual(tco, 3439.615, 1)

    def test_cooler_cost(self):
        p_loss_1 = 100.  # W
        p_loss_2 = 1000.  # W

        c1 = cooler_cost(p_loss_1)
        c2 = cooler_cost(p_loss_2)

        self.assertAlmostEqual(c1, 24984.95, 1)
        self.assertAlmostEqual(c2, 92827.91, 1)
