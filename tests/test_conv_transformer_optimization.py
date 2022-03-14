from unittest import TestCase
from src.trafo_opt_conventional import TwoWindingModel, TransformerDesign
from pathlib import Path

"""6300 kVA Transformer from Karsai, Nagytranszformátorok """


class TestConvTransformerModel(TestCase):

    def test_transformer_from_json(self):
        path = Path().cwd().parent.joinpath('data').joinpath('6300_kVA_example.json')

        import json

        with open(path) as json_file:
            data = json.load(json_file)

        transformer = TransformerDesign.from_dict(data)

        trafo_model = TwoWindingModel(input=transformer)
        trafo_model.calculate()

        self.assertAlmostEqual(trafo_model.results.turn_voltage, 31.0, 1)
        print(trafo_model)
