from unittest import TestCase

from importlib_resources import files

from src.models import (
    IndependentVariables,
    MaterialCosts,
    TransformerDesign,
    TransformerRequirements,
    WindingDesign,
    WindingParams,
)


class TestWindingModel(TestCase):
    def test_calc(self):
        winding = WindingParams(connection="y", line_voltage=22.0)

        winding.calculate_phase_quantities(nominal_power=6300.0)
        self.assertAlmostEqual(winding.ph_current, 95.57, 2)
        self.assertAlmostEqual(winding.ph_voltage, 12.716, 2)

    def test_calc_properties(self):
        winding = WindingDesign(winding_height=1100, inner_radius=230, thickness=35, filling_factor=53.5, current_density=3.02)

        winding.calc_properties()

        self.assertAlmostEqual(winding.mean_radius, 247.5, 2)
        self.assertAlmostEqual(winding.outer_radius, 265, 2)
        self.assertAlmostEqual(winding.mass, 855.225, 2)
        self.assertAlmostEqual(winding.dc_loss, 18.876, 2)
        self.assertAlmostEqual(winding.ac_loss, 1.59, 2)


class TestTransformerToJSON(TestCase):
    def test_transformer_to_json(self):
        cost = MaterialCosts()
        req_params = TransformerRequirements(
            power=6300,
            freq=50,
            sci_req=7.34,
            drop_tol=5.0,
            hv=WindingParams(connection="y", line_voltage=33.0, filling_factor=56.0),
            lv=WindingParams(connection="y", line_voltage=22.0, filling_factor=53.5),
            min_main_gap=20.0,
            min_core_gap=14.0,
            ei=150.0,
            phase_distance=40.0,
            alpha=0.97,
            core_fillingf=83.8,
        )

        ind_params = IndependentVariables(rc=184, bc=1.568, j_in=3.02, j_ou=3.0, h_in=979.0, m_gap=26.7)
        transformer = TransformerDesign(
            description="6300 kVA Transformer from Karsai, Large Power Transformers  book (Hun)",
            required=req_params,
            design_params=ind_params,
            costs=cost,
        )

        print(transformer)

        json_string = transformer.to_json()

        path = files("data").joinpath("6300_kVA_example.json")

        with open(path, "w") as outfile:
            outfile.write(json_string)

    def test_transformer_from_json(self):
        path = files("data").joinpath("6300_kVA_example.json")

        import json

        with open(path) as json_file:
            data = json.load(json_file)
        transformer = TransformerDesign.from_dict(data)

        self.assertIn("6300", transformer.description)
        self.assertEqual(6300, transformer.required.power)
        self.assertEqual(50, transformer.required.freq)
        self.assertEqual(184, transformer.design_params.rc)

    def test_10_mva_transformer_data_from_json(self):
        path = files("data").joinpath("10MVA_example.json")

        import json

        with open(path) as json_file:
            data = json.load(json_file)
        transformer = TransformerDesign.from_dict(data)

        self.assertIn("10", transformer.description)
