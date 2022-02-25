from unittest import TestCase
from src.models import WindingModel


class TestWindingModel(TestCase):

    def test_calc(self):
        winding = WindingModel(r_in=230, h_in=1100, t_in=35, con='y', ul=22.0, j=2.6, ff=60.0)
        winding.calculate_phase_quantities(nominal_power=6300.)
        self.assertAlmostEqual(winding.ph_current, 95.57, 2)
        self.assertAlmostEqual(winding.ph_voltage, 12.716, 2)

    def test_calc_properties(self):
        winding = WindingModel(r_in=230, h_in=1100, t_in=35, con='y', ul=22.0, j=3.02, ff=53.5)
        winding.calc_properties(ph_num=3.)

        self.assertAlmostEqual(winding.mean_radius, 247.5,2)
        self.assertAlmostEqual(winding.outer_radius, 265, 2)
        self.assertAlmostEqual(winding.mass, 855.225, 2)
        self.assertAlmostEqual(winding.dc_loss, 18.876, 2)
        self.assertAlmostEqual(winding.ac_loss, 0.17, 2)

