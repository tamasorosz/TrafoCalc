from agrossuite import agros


class Model:
    """The goal of this class is to build a basic 2d axisymmetric model for transformer simulation in Agros Suite"""

    def __init__(self):

        self.problem = agros.problem(clear=True)
        self.geo = self.problem.geometry()

        self.problem.coordinate_type = "axisymmetric"
        self.problem.mesh_type = "triangle"

        self.magnetic = self.problem.field("magnetic")
        self.magnetic.analysis_type = "steadystate"
        self.magnetic.number_of_refinements = 2
        self.magnetic.polynomial_order = 2
        self.magnetic.adaptivity_type = "disabled"
        self.magnetic.solver = "linear"

        # boundaries
        self.magnetic.add_boundary("A = 0", "magnetic_potential", {"magnetic_potential_real": 0})

        # materials
        self.magnetic.add_material("Air",
                                   {"magnetic_permeability": 1, "magnetic_conductivity": 0, "magnetic_remanence": 0,
                                    "magnetic_remanence_angle": 0,
                                    "magnetic_velocity_x": 0, "magnetic_velocity_y": 0,
                                    "magnetic_velocity_angular": 0,
                                    "magnetic_current_density_external_real": 0,
                                    "magnetic_total_current_prescribed": 0, "magnetic_total_current_real": 0})

    def create_rectangle(self, x0, y0, width, height, boundary: dict = None):
        """
        A rectangle class to define the windings and the working window of the transformer.

        @param geo: geometry object
        @param x0: x coordinate of the bottom - left node
        @param y0: y coordinate of the bottom - left node
        @param height: height of the rectangle
        @param width: width of the rectangle
        @param boundary: boundary conditions, dictionary like that: {"magnetic":"A = 0"}

        The rectangle has the same bondary condition in all edges.
        """

        # unit conversion m -> mm

        x0 *= 1e-3
        y0 *= 1e-3
        height *= 1e-3
        width *= 1e-3

        if boundary is not None:

            self.geo.add_edge(x0, y0, x0 + width, y0, boundaries=boundary)
            self.geo.add_edge(x0 + width, y0, x0 + width, y0 + height, boundaries=boundary)
            self.geo.add_edge(x0 + width, y0 + height, x0, y0 + height, boundaries=boundary)
            self.geo.add_edge(x0, y0 + height, x0, y0, boundaries=boundary)

        else:

            self.geo.add_edge(x0, y0, x0 + width, y0)
            self.geo.add_edge(x0 + width, y0, x0 + width, y0 + height)
            self.geo.add_edge(x0 + width, y0 + height, x0, y0 + height)
            self.geo.add_edge(x0, y0 + height, x0, y0)

        return x0 + width / 2. * 1e3, y0 + height / 2. * 1e3  # gives back the center of the rectangle in mm-s

    def create_winding(mag, geo, x0, y0, width, height, name, filling_f, j):
        """
        @param geo: geometry object
        @param x0: x coordinate of the bottom - left node
        @param y0: y coordinate of the bottom - left node
        @param height: height of the rectangle
        @param width: width of the rectangle
        @param name: name of the winding
        @param filling_f: filling factor, between [0--1]
        @param j: current density in A/mm2
        """

        x_label, y_label = create_rectangle(geo, x0, y0, width, height)
        # print(x0, y0, width, height)
        mag.add_material(name, {"magnetic_permeability": 1,
                                "magnetic_conductivity": 57 * 1e6 * filling_f,
                                "magnetic_remanence": 0,
                                "magnetic_remanence_angle": 0,
                                "magnetic_velocity_x": 0,
                                "magnetic_velocity_y": 0,
                                "magnetic_velocity_angular": 0,
                                "magnetic_current_density_external_real": j * 1e6 * filling_f})  # j in A/m2

        # creates a label for the material definition
        geo.add_label(x_label, y_label, materials={"magnetic": "{0}".format(name)})

        return
