from unittest import TestCase
from src.transformer_fem_model import FemModel


class TestGeoCreation(TestCase):

    def test_rectangle(self):
        model = FemModel()

        x0 = 0
        y0 = 0
        width = 100  # mm
        height = 100  # mm
        boundary = {"magnetic": "A = 0"}

        center = model.create_rectangle(x0, y0, width, height, boundary)

        self.assertEqual(center, (50*1e-3, 50*1e-3))

    def test_create_windings(self):

        model = FemModel()

        model.create_winding(10,10,20,100,'lv',0.7, 3.0)
