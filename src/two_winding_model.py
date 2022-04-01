import typing
from dataclasses import dataclass, field
from math import pi

from dataclasses_json import dataclass_json

from src.base_functions import turn_voltage, calc_inner_width, inner_winding_radius, outer_winding_radius, \
    window_width, core_mass, core_loss_unit, short_circuit_impedance, capitalized_cost

from src.models import MainResults, TransformerDesign, WindingDesign
from src.transformer_fem_model import FemModel

C_WIN_MIN = 10.0  # [mm] technological limit for the thickness of the windings, it should be larger than 10 mm-s
INFEASIBLE = -1
CORE_BF = 1.2  # building factor of the core


@dataclass_json
@dataclass
class TwoWindingModel:
    input: TransformerDesign
    results: MainResults = field(default=MainResults())
    # winding models
    hv_winding: typing.Any = field(default=None)
    lv_winding: typing.Any = field(default=None)

    def calculate(self, is_sc=False):
        """
        Calculates the main geometrical parameters and invokes the WindingParameter class which calculates the searched
        parameters -> losses, masses
        :param is_sc: True if superconducting transformer considered
        :return: feasible (boolean)
        """
        # 1) phase power, assumes a 3 phased 3 legged transformer core
        ph_power = self.input.required.power / 3.0  # [kVA]

        # 2) turn voltage
        self.results.turn_voltage = turn_voltage(
            self.input.design_params.bc,
            self.input.design_params.rc,
            self.input.required.core_fillingf / 100.0,
            self.input.required.freq,
        )

        # 3) main parameters for the inner winding
        t_in = calc_inner_width(
            ph_power,
            self.input.design_params.h_in,
            self.input.required.lv.filling_factor / 100.0,
            self.input.design_params.j_in,
            self.results.turn_voltage,
        )

        r_in = inner_winding_radius(self.input.design_params.rc, self.input.required.min_core_gap, t_in)

        # 4/ outer winding parameters
        h_ou = self.input.design_params.h_in * self.input.required.alpha
        t_ou = calc_inner_width(
            ph_power,
            h_ou,
            self.input.required.hv.filling_factor / 100.0,
            self.input.design_params.j_ou,
            self.results.turn_voltage,
        )

        # outer winding radius
        r_ou = outer_winding_radius(r_in, t_in, self.input.design_params.m_gap, t_ou)

        # calculating the detailed parameters of the winding
        self.lv_winding = WindingDesign(
            inner_radius=r_in - t_in / 2.0,
            thickness=t_in,
            winding_height=self.input.design_params.h_in,
            current_density=self.input.design_params.j_in,
            filling_factor=self.input.required.lv.filling_factor,
        )
        self.hv_winding = WindingDesign(
            inner_radius=r_ou - t_ou / 2.0,
            thickness=t_ou,
            winding_height=h_ou,
            current_density=self.input.design_params.j_ou,
            filling_factor=self.input.required.hv.filling_factor,
        )

        if is_sc:
            self.lv_winding.calc_properties()
            self.hv_winding.calc_properties()
        else:
            # if the resulting thickness of the winding is smaller than the required minimum the
            # solution is not feasible
            if t_in < C_WIN_MIN or t_ou < C_WIN_MIN:
                self.results.feasible = INFEASIBLE
                raise ValueError("The winding thickness is too narrow.")

            self.lv_winding.calc_properties()
            self.hv_winding.calc_properties()

        # window width and window heights
        self.results.window_width = window_width(
            self.input.required.min_core_gap, t_in, t_ou, self.input.design_params.m_gap, 0, 0
        )

        self.results.wh = self.input.design_params.h_in + self.input.required.ei

        # core parameters
        self.results.core_mass = core_mass(
            self.input.design_params.rc,
            self.input.required.core_fillingf / 100.0,
            self.input.design_params.h_in,
            self.input.required.ei,
            self.results.window_width,
            self.input.required.phase_distance / 2.0,
        )
        self.results.core_loss = core_loss_unit(self.input.design_params.bc, self.results.core_mass, CORE_BF)

        self.results.load_loss = (
                self.lv_winding.ac_loss + self.lv_winding.dc_loss + self.hv_winding.ac_loss + self.hv_winding.dc_loss
        )

        # short circuit impedance calculation with analytical formulas
        self.results.sci = short_circuit_impedance(
            self.input.required.power,
            3.0,
            self.input.required.freq,
            self.input.required.alpha,
            self.results.turn_voltage,
            self.input.design_params.h_in,
            self.results.window_width,
            r_in,
            t_in,
            r_ou,
            t_ou,
            self.input.design_params.m_gap,
        )

        self.results.capitalized_cost = capitalized_cost(
            self.results.core_mass,
            self.input.costs.core_cost,
            self.lv_winding.mass,
            self.input.costs.lv_cost,
            self.hv_winding.mass,
            self.input.costs.hv_cost,
            self.results.load_loss,
            self.input.costs.ll_cost,
            self.results.core_loss,
            self.input.costs.nll_cost,
        )

        self.results.feasible = True

    def fem_simulation(self):
        if not self.results.feasible:
            raise ValueError("Invalid Transformer Geometry")

        # initializing the model
        simulation = FemModel()

        # creating the core window and the two windings
        simulation.create_rectangle(
            self.input.design_params.rc, 0, self.results.window_width, self.results.wh, {"magnetic": "A = 0"}
        )

        # label for the air/oil region in the transformer
        simulation.geo.add_label((self.input.design_params.rc + 10) * 1e-3, 10 * 1e-3, materials={"magnetic": "Air"})

        # windings
        simulation.create_winding(
            self.lv_winding.inner_radius,
            self.input.required.ei / 2.0,
            self.lv_winding.thickness,
            self.lv_winding.winding_height,
            "lv",
            self.lv_winding.filling_factor / 100.0,
            self.lv_winding.current_density * 2.0 ** 0.5,
        )
        # #
        simulation.create_winding(
            self.hv_winding.inner_radius,
            self.input.required.ei / 2.0,
            self.hv_winding.thickness,
            self.hv_winding.winding_height,
            "hv",
            self.hv_winding.filling_factor / 100.0,
            -self.hv_winding.current_density * 2.0 ** 0.5,
        )

        computation = simulation.problem.computation()
        computation.solve()
        solution = computation.solution("magnetic")

        # calculate the base quantites for the LV winding
        self.input.required.lv.calculate_phase_quantities(self.input.required.power)

        # the base quantites referred to the low voltage winding
        u_b = self.input.required.hv.line_voltage  # voltage --- kV
        s_b = self.input.required.power / 1000.0  # nominal power  --- MVA
        z_b = u_b ** 2.0 / s_b  # base impedance
        i_b = self.input.required.power / u_b / 1.73

        omega = 2.0 * pi * self.input.required.freq
        L = 2 * solution.volume_integrals()["Wm"] / i_b ** 2.0
        self.results.fem_based_sci = omega * L / z_b * 100.0  # the short-circuit impedance in [%] values

        # axial and radial components of the magnetic flux densities along the inner radius of the hv winding
        for i in range(
                int(self.input.required.ei / 2.0), int(self.hv_winding.winding_height + self.input.required.ei / 2.0),
                10
        ):
            point = solution.local_values(self.hv_winding.inner_radius * 1e-3, i * 1e-3)
            self.results.fem_bax = max(abs(point["Brz"]), self.results.fem_bax)

        # iterates over the top of the hv winding
        for i in range(int(self.hv_winding.inner_radius), int(self.hv_winding.inner_radius + self.hv_winding.thickness),
                       3):
            point = solution.local_values(i * 1e-3,
                                          (self.hv_winding.winding_height + self.input.required.ei / 2.0) * 1e-3)
            self.results.fem_brad = max(abs(point["Brr"]), self.results.fem_brad)

        print(self.results.fem_bax, self.results.fem_brad)
