from unittest import TestCase
from src.superconductor_losses import parallel_loss, perp_loss


class TestLosses(TestCase):

    def test_losses(self):
        # BSCCO cable from Magnussons' paper
        f = 50.0
        c = 0.75
        bp = 0.033
        Ac = 0.21 * 4.29 * 1e-6

        bpar = 50. * 1e-3  # mT
        self.assertAlmostEqual(parallel_loss(f, c, Ac, bpar, bp), 0.04968,3)
