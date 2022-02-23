from unittest import TestCase
from src.basic_model import create_rectangle, create_winding
from agrossuite import agros


class TestGeoCreation(TestCase):

    def test_rectangle(self):
        geo = self.init_agros()

        x0 = 0
        y0 = 0
        width = 100  # mm
        height = 100  # mm
        boundary = {"magnetic": "A = 0"}

        center = create_rectangle(geo, x0, y0, width, height, boundary)

        self.assertEqual(center, (50, 50))

    def test_create_windings(self):
        self.init_agros()
