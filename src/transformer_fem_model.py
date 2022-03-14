from agrossuite import agros


class FemModel:
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

        return x0 + width / 2., y0 + height / 2.  # gives back the center of the rectangle in [m]-s

    def create_winding(self, x0, y0, width, height, name, filling_f, j):
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

        x_label, y_label = self.create_rectangle(x0, y0, width, height)
        # print(x0, y0, width, height)
        self.magnetic.add_material(name, {"magnetic_permeability": 1,
                                          "magnetic_conductivity": 57 * 1e6 * filling_f,
                                          "magnetic_remanence": 0,
                                          "magnetic_remanence_angle": 0,
                                          "magnetic_velocity_x": 0,
                                          "magnetic_velocity_y": 0,
                                          "magnetic_velocity_angular": 0,
                                          "magnetic_current_density_external_real": j * 1e6 * filling_f})  # j in A/m2

        # creates a label for the material definition
        self.geo.add_label(x_label, y_label, materials={"magnetic": "{0}".format(name)})

        return

    # def model(self):
    #
    #
    #
    #     # base impedance for the short circuit impedance calculation
    #     ub = param.u_in_line  # voltage --- kV
    #     sb = param.power / 1000.  # nominal power  --- MVA
    #     zb = ub ** 2. / sb  # impedance
    #     ib = param.power / ub / 1.73  # phase_current(sb, ub, param.con_fact_in)  # amps in the impedance base calculation
    #     Omega = 2. * pi * param.freq
    #
    #     # -- calculating the impedance from the volume integrals --
    #     L = 2 * solution.volume_integrals()['Wm'] / ib ** 2.
    #     dep.sci = 2. * pi * param.freq * L / zb * 100.
    #
    #     # calculating the average radial and axial flux density
    #     nx = 5  # number of the local inductances along the z coordinate
    #     ny = 5  # number of the local inductance along the x coordinate
    #
    #     # low voltage winding  ----------------------------------
    #     Br_avg = 0.
    #     Bz_avg = 0.
    #
    #     dx = dep.t_in / nx * 1e-3
    #     dy = indep.h_in / ny * 1e-3
    #
    #     for i in range(0, nx):
    #         xx = dep.r_in * 1e-3 + i * dx
    #         for j in range(0, ny):
    #             yy = param.ei / 2. * 1e-3 + j * dy
    #
    #             point = solution.local_values(xx, yy)
    #             Br_avg += point["Brr"]
    #             Bz_avg += point["Brz"]
    #
    #     Br_avg /= (nx * ny)
    #     Bz_avg /= (nx * ny)
    #
    #     b1 = (Br_avg ** 2. + Bz_avg ** 2.) ** 0.5
    #
    #     # ---------------------------------------------------------
    #     # print('bz_avg:', Bz_avg)
    #     # print('br_avg:', Br_avg)
    #
    #     # -- calculating the dc and ac losses from the model
    #     n_in = param.u_in / dep.turn_voltage * 1e3  # nr of the turns in the inner winding
    #     n_ou = param.u_out / dep.turn_voltage * 1e3  # nr of the turns in the inner winding
    #     Bg = b_gap(n_in, ib,
    #                indep.h_in)  # TODO <--- should be changed by the maximum of the inductance in the gap --- at least ...
    #
