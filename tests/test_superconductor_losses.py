from unittest import TestCase
from src.superconductor_losses import parallel_loss, perp_loss, norris_equation, cryostat_losses, cryo_surface, \
    thermal_incomes, modified_tco_evaluation

from math import pi


class TestLosses(TestCase):

    def test_losses(self):
        # BSCCO cable from Magnussons' paper
        # f = 50.0
        # c = 0.75
        # bp = 0.033
        # Ac = 0.31 * 4.29 * 1e-6

        bpar = 66.0 * 1e-3  # mT

        # parallel losses
        # 0.133 W/m in 10.1109/TASC.2003.813123 but this  value seems ok
        self.assertAlmostEqual(parallel_loss(50, bpar), 0.177, 3)
        # perpendicular losses
        self.assertAlmostEqual(perp_loss(50, 68. * 1e-3), 2.009, 3)
        # norris equation
        self.assertAlmostEqual(norris_equation(50, 50, 115.), 0.0047, 4)

    def test_cryostat_losses(self):
        r_in = 0.5 / pi
        r_ou = 1.0 / pi
        h = 1.0
        a_cs = cryo_surface(r_in, r_ou, h)  # [m2]
        P_cryo = cryostat_losses(a_cs)
        P_thermal = thermal_incomes(100, 100)

        self.assertAlmostEqual(a_cs, 3.238, 2)
        self.assertAlmostEqual(P_cryo, 0.029537, 2)
        self.assertAlmostEqual(P_thermal, 54, 1)

        tco = modified_tco_evaluation(1, 1, 1, 0, 0, 1, 0)

        self.assertAlmostEqual(tco, 54, 1)
