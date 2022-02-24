from unittest import TestCase
from src.basic_model import Model


class TestGeoCreation(TestCase):

    def test_rectangle(self):
        model = Model()

        x0 = 0
        y0 = 0
        width = 100  # mm
        height = 100  # mm
        boundary = {"magnetic": "A = 0"}

        center = model.create_rectangle(x0, y0, width, height, boundary)

        self.assertEqual(center, (50, 50))

    def test_create_windings(self):

        model = Model()

        model.create_winding(10,10,20,100,'lv',0.7, 3.0)
