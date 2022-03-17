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
        self.assertAlmostEqual(trafo_model.results.turn_voltage, 46.27, 1)
        self.assertAlmostEqual(trafo_model.lv_winding.thickness, 35.0, 0)
        self.assertAlmostEqual(trafo_model.hv_winding.thickness, 44.0, 0)

        # core window geometry
        self.assertAlmostEqual(trafo_model.results.window_width, 199, 0)
        self.assertAlmostEqual(trafo_model.results.core_mass, 7301.5, 0)

        # losses
        self.assertAlmostEqual(trafo_model.results.core_loss, 7.7, 1)
        self.assertAlmostEqual(trafo_model.hv_winding.mass, 1585.5, 0)

        self.assertAlmostEqual(trafo_model.results.load_loss, 44.8, 0)
        self.assertAlmostEqual(trafo_model.results.sci * 100, 7.5, 1)

        self.assertAlmostEqual(trafo_model.results.capitalized_cost, 0, 0)
        print(trafo_model)

        # FEM calculation
        trafo_model.fem_simulation()

        self.assertAlmostEqual(trafo_model.results.fem_based_sci, 7.7, 1)

        del trafo_model

    def test_sc_transformer(self):
        path = files("data").joinpath("630kVA_sc_transformer.json")

        import json

        with open(path) as json_file:
            data = json.load(json_file)

        transformer = TransformerDesign.from_dict(data)

        trafo_model = TwoWindingModel(input=transformer)
        trafo_model.calculate(is_sc=True)
        print(trafo_model)
        # FEM calculation
        trafo_model.fem_simulation()

        self.assertAlmostEqual(trafo_model.lv_winding.inner_radius, 244.0, 1)
        self.assertAlmostEqual(trafo_model.lv_winding.outer_radius, 252, 1)
        self.assertAlmostEqual(trafo_model.lv_winding.thickness, 8.0, 0)
        self.assertAlmostEqual(trafo_model.hv_winding.thickness, 13.5, 0)
        self.assertAlmostEqual(trafo_model.results.sci*100, 4.8,0)
        self.assertAlmostEqual(trafo_model.results.fem_based_sci, 6.3,1)

        del trafo_model
