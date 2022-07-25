import typing
from dataclasses import dataclass, field
from math import pi

from dataclasses_json import dataclass_json

from src.base_functions import turn_voltage, calc_inner_width, inner_winding_radius, outer_winding_radius, \
    window_width, core_mass, core_loss_unit, short_circuit_impedance, capitalized_cost

from src.models import MainResults, TransformerDesign, WindingDesign
from src.transformer_fem_model import FemModel
from src.superconductor_losses import cryostat_losses, sc_load_loss, cryo_surface, thermal_incomes, cooler_cost
from src.diagrams import plot_winding_flux

C_WIN_MIN = 10.0  # [mm] technological limit for the thickness of the windings, it should be larger than 10 mm-s
SC_WIN_MIN = 8.0  # [mm] sc_transformer winding minimum
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

        # outer winding radius (MEAN)
        r_ou = outer_winding_radius(r_in, t_in, self.input.design_params.m_gap, t_ou)

        # calculating the detailed parameters of the winding
        self.lv_winding = WindingDesign(
            inner_radius=round(r_in - t_in / 2.0, 1),
            thickness=t_in,
            winding_height=self.input.design_params.h_in,
            current_density=self.input.design_params.j_in,
            filling_factor=self.input.required.lv.filling_factor,
        )
        self.hv_winding = WindingDesign(
            inner_radius=round(r_ou - t_ou / 2.0, 1),
            thickness=t_ou,
            winding_height=h_ou,
            current_density=self.input.design_params.j_ou,
            filling_factor=self.input.required.hv.filling_factor,
        )

        # calculate the phase quantities of the windings
        self.input.required.lv.calculate_phase_quantities(self.input.required.power)
        self.input.required.hv.calculate_phase_quantities(self.input.required.power)

        if is_sc:
            if t_in < C_WIN_MIN or t_ou < C_WIN_MIN:
                ###
                self.lv_winding.calc_sc_properties(self.input.required.freq)
                self.hv_winding.calc_sc_properties(self.input.required.freq)
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

        self.results.load_loss = round(
            self.lv_winding.ac_loss + self.lv_winding.dc_loss + self.hv_winding.ac_loss + self.hv_winding.dc_loss, 2
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

        if is_sc:
            # there are significant losses generated by the cryostat and in the current leads
            # the approximate surface of the cryostat
            a_cs = cryo_surface(r_in, r_ou + self.input.design_params.m_gap, h_ou)  # [m2]
            cryo_loss = cryostat_losses(a_cs)
            thermal_loss = thermal_incomes(self.input.required.lv.ph_current, self.input.required.hv.ph_current)

            # the formula uses the default c factor for the calculations
            self.results.load_loss = sc_load_loss(self.results.load_loss, cryo_loss, thermal_loss)

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
        self.results.copper_mass = self.lv_winding.mass + self.hv_winding.mass
        self.results.feasible = True

    def fem_simulation(self, detailed_output=True):
        if not self.results.feasible:
            raise ValueError("Invalid Transformer Geometry")

        # initializing the model
        simulation = FemModel()

        # # creating the core window and the two windings
        simulation.create_rectangle(self.input.design_params.rc, self.input.design_params.rc, self.results.window_width,
                                    self.results.wh, None)

        # label for the air/oil region in the transformer
        simulation.geo.add_label((self.input.design_params.rc + 20) * 1e-3, (self.input.design_params.rc + 10) * 1e-3,
                                 materials={"magnetic": "Air"})

        # core
        simulation.create_rectangle(0, 0, self.results.window_width + 2 * self.input.design_params.rc,
                                    self.results.wh + 2 * self.input.design_params.rc,
                                    {"magnetic": "A = 0"})

        # label for the air/oil region in the transformer
        simulation.geo.add_label(0.01, 1e-3, materials={"magnetic": "Core"})

        # windings
        simulation.create_winding(
            self.lv_winding.inner_radius,
            self.input.required.ei / 2.0 + self.input.design_params.rc,
            self.lv_winding.thickness,
            self.lv_winding.winding_height,
            "lv",
            self.lv_winding.filling_factor / 100.0,
            self.lv_winding.current_density,
        )

        simulation.create_winding(
            self.hv_winding.inner_radius,
            self.input.required.ei / 2.0 + self.input.design_params.rc,
            self.hv_winding.thickness,
            self.hv_winding.winding_height,
            "hv",
            self.hv_winding.filling_factor / 100.0,
            -self.hv_winding.current_density,
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
        i_b = self.input.required.power / u_b / 3. ** 0.5

        omega = 2.0 * pi * self.input.required.freq
        L = 2 * solution.volume_integrals()["Wm"] / i_b ** 2.0
        print('Magnetic Energy', solution.volume_integrals()["Wm"])
        print('zb, ib:', round(z_b, 2), 'ohm', round(i_b, 2), 'A')

        self.results.fem_based_sci = round(omega * L / z_b * 100.0, 2)  # the short-circuit impedance in [%] values
        print('SCI:', self.results.fem_based_sci, '[%]')
        # axial and radial components of the magnetic flux densities along the inner radius of the hv winding

        # iterates over on horizontal slices in the hv winding and the lv winding to collect the required flux data
        self.results.fem_bax_hv = []
        self.results.fem_brad_hv = []

        i = int(self.input.design_params.rc + self.input.required.ei / 2.0)
        hv_top = int(self.input.design_params.rc + self.hv_winding.winding_height + self.input.required.ei / 2.0)
        dz = int(self.hv_winding.winding_height / 20)
        while i <= hv_top + dz / 2:

            max_rad = 0.
            max_ax = 0.

            # iterates over the winding in the radial direction and stores the max value
            for j in range(int(self.hv_winding.inner_radius),
                           int(self.hv_winding.inner_radius + self.hv_winding.thickness), 1):
                point = solution.local_values(j * 1e-3, i * 1e-3)

                max_rad = max(abs(point["Brr"]), max_rad)
                max_ax = max(abs(point["Brz"]), max_ax)

            # max values along the hv winding
            self.results.fem_brad_hv.append(round(max_rad * 1e3, 2))
            self.results.fem_bax_hv.append(round(max_ax * 1e3, 2))

            i += dz

        # create a common list from bax and brad values
        self.results.br_bax_hv = list(zip(self.results.fem_bax_hv, self.results.fem_brad_hv))

        # collecting the critical points from the lv winding
        self.results.fem_bax_lv = []
        self.results.fem_brad_lv = []

        i = self.input.design_params.rc + self.input.required.ei / 2.0
        lv_top = int(self.input.design_params.rc + self.lv_winding.winding_height + self.input.required.ei / 2.0)
        dz = self.lv_winding.winding_height / 20
        while i <= lv_top + dz / 2:

            max_rad = 0.
            max_ax = 0.

            # iterates over the winding in the radial direction and stores the max value
            for j in range(int(self.lv_winding.inner_radius),
                           int(self.lv_winding.inner_radius + self.lv_winding.thickness), 5):
                point = solution.local_values(j * 1e-3, i * 1e-3)
                max_rad = max(abs(point["Brr"]), max_rad)
                max_ax = max(abs(point["Brz"]), max_ax)

            # max values along the hv winding
            self.results.fem_brad_lv.append(round(max_rad * 1e3, 2))
            self.results.fem_bax_lv.append(round(max_ax * 1e3, 2))

            i += dz

        # create a common list from bax and brad values
        self.results.br_bax_lv = list(zip(self.results.fem_bax_lv, self.results.fem_brad_lv))

        # maximum radial and axial values in the hv winding
        self.results.fem_brad_hv = max(self.results.fem_brad_hv)
        self.results.fem_bax_hv = max(self.results.fem_bax_hv)
        # lv
        self.results.fem_brad_lv = max(self.results.fem_brad_lv)
        self.results.fem_bax_lv = max(self.results.fem_bax_lv)

        print('Bax  [HV] =', self.results.fem_bax_hv, '[mT]')
        print('Brad [HV] =', self.results.fem_brad_hv, '[mT]')

        print('Bax  [LV] =', self.results.fem_bax_lv, '[mT]')
        print('Brad [Lv] =', self.results.fem_brad_lv, '[mT]')

        if detailed_output:
            # print('Values along the hv winding:', list(self.results.br_bax_hv))
            # print('Values along the lv winding:', list(self.results.br_bax_lv))
            plot_winding_flux(self.results.br_bax_lv, 0, self.lv_winding.winding_height, label='LV')
            plot_winding_flux(self.results.br_bax_hv, 0, self.hv_winding.winding_height, label='HV')
