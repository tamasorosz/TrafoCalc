from unittest import TestCase
from src.trafo_opt_conventional import TransformerModel
from src.models import transformer, independent_variables

"""6300 kVA Transformer from Karsai, Nagytranszformátorok """

ind = independent_variables(r_c=184.0, b_c=1.568,
                            j_in=3.02, j_ou=3.0,
                            h_in=979., m_gap=26.7)



class TestConvTransformerModel(TestCase):

    def test_phase_calc(self):
        trafo = TransformerModel(ind, trafo_6300)
        trafo.calc_phase_quantites()

        self.assertAlmostEqual(trafo.transformer_design['lv']['ph_voltage'], 4.0, 3)
        self.assertAlmostEqual(trafo.transformer_design['hv']['ph_voltage'], 19.075, 3)