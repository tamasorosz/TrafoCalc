from unittest import TestCase
from math import pi
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

        self.assertEqual(center, (50 * 1e-3, 50 * 1e-3))

    def test_create_windings(self):
        model = FemModel()

        model.create_winding(10, 10, 20, 100, "lv", 0.7, 3.0)

    def test_fem_simulation(self):
        # 31.5 MVA transformer
        # source: https://www.ee.iitb.ac.in/~fclab/FEM/FEM1.pdf
        simulation = FemModel()

        # airgap
        # creating the core window and the two windings
        simulation.create_rectangle(270, 540, 287, 1800, None)

        # label for the air/oil region in the transformer
        simulation.geo.add_label(0.280, 0.545, materials={"magnetic": "Air"})

        # core
        simulation.create_rectangle(0, 0, 1097, 2880, {"magnetic": "A = 0"})

        # label for the air/oil region in the transformer
        simulation.geo.add_label(0.01, 1e-3, materials={"magnetic": "Core"})

        # windings
        simulation.create_winding(
            293,
            620,
            52,
            1520,
            "lv",
            1,
            1.708299595,
        )

        simulation.create_winding(
            394,
            620,
            65,
            1520,
            "hv",
            1,
            -1.3666396,
        )

        computation = simulation.problem.computation()
        computation.solve()
        solution = computation.solution("magnetic")

        # the base quantites referred to the low voltage winding
        u_b = 132  # voltage --- kV
        s_b = 31.5  # nominal power  --- MVA
        z_b = u_b ** 2.0 / s_b  # base impedance
        i_b = 31500 / u_b / 3. ** 0.5

        omega = 2.0 * pi * 50.0
        L = 2 * solution.volume_integrals()["Wm"] / i_b ** 2.0
        print('Magnetic Energy', solution.volume_integrals()["Wm"])
        print('L:', L)
        print('zb, ib:', z_b, i_b)

        fem_based_sci = omega * L / z_b * 100.0  # the short-circuit impedance in [%] values
        print(fem_based_sci)

        self.assertAlmostEqual(1480.7, solution.volume_integrals()["Wm"], 1)
