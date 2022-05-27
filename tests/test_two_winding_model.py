from unittest import TestCase

from importlib_resources import files

from src.two_winding_model import TransformerDesign, TwoWindingModel

"""10 MVA Transformer from Karsai, Nagytranszformátorok """


class TestConvTransformerModel(TestCase):
    def test_transformer_from_json(self):
        path = files("data").joinpath("10MVA_example.json")

        import json

        with open(path) as json_file:
            data = json.load(json_file)

        transformer = TransformerDesign.from_dict(data)

        trafo_model = TwoWindingModel(input=transformer)
        trafo_model.calculate()

        # winding details
        self.assertAlmostEqual(trafo_model.results.turn_voltage, 46.64, 1)
        self.assertAlmostEqual(trafo_model.lv_winding.thickness, 35.0, 0)
        self.assertAlmostEqual(trafo_model.hv_winding.thickness, 43.4, 0)
        self.assertAlmostEqual(trafo_model.hv_winding.amper_turns, 71406.8, 0)

        # core window geometry
        self.assertAlmostEqual(trafo_model.results.window_width, 198.4, 0)
        self.assertAlmostEqual(trafo_model.results.core_mass, 7751.0, 0)

        # losses
        self.assertAlmostEqual(trafo_model.results.core_loss, 8.2, 1)
        self.assertAlmostEqual(trafo_model.hv_winding.mass, 1620., 0)

        self.assertAlmostEqual(trafo_model.results.load_loss, 49.18, 0)
        self.assertAlmostEqual(trafo_model.results.sci, 7.5, 0)

        self.assertAlmostEqual(trafo_model.results.capitalized_cost, 0, 0)
        print(trafo_model)

        # FEM calculation
        trafo_model.fem_simulation(detailed_output=False)

        self.assertAlmostEqual(trafo_model.results.fem_based_sci, 7.56, 1)

        del trafo_model
