from dataclasses import dataclass, field

from dataclasses_json import dataclass_json

from src.models import TransformerDesign, MainResults, WindingDesign
from src.transformer_calculations import *

C_WIN_MIN = 10.0  # [mm] technological limit for the thickness of the windings, it should be larger than 10 mm-s
INFEASIBLE = -1
CORE_BF = 1.2  # building factor of the core


@dataclass_json
@dataclass
class TwoWindingModel:
    input: TransformerDesign
    results: MainResults = field(default=MainResults)
    # winding models
    hv_winding: WindingDesign = field(default=None)  # derived winding parameters and data
    lv_winding: WindingDesign = field(default=None)

    def calculate(self):
        """
        Calculates the main geometrical parameters and invokes the WindingParameter class which calculates the searched
        parameters -> losses, masses
        :return: feasible (boolean)
        """
        # 1) phase power, assumes a 3 phased 3 legged transformer core
        ph_power = self.input.required.power / 3.  # [kVA]

        # 2) turn voltage
        self.results.turn_voltage = turn_voltage(self.input.design_params.bc, self.input.design_params.rc,
                                                 self.input.required.core_fillingf / 100., self.input.required.freq)

        # 3) main parameters for the inner winding
        t_in = calc_inner_width(ph_power, self.input.design_params.h_in, self.input.required.lv.filling_factor / 100.,
                                self.input.design_params.j_in, self.results.turn_voltage)

        r_in = inner_winding_radius(self.input.design_params.rc, self.input.required.min_core_gap, t_in)

        # 4/ outer winding parameters
        h_ou = self.input.design_params.h_in * self.input.required.alpha
        t_ou = calc_inner_width(ph_power, h_ou, self.input.required.hv.filling_factor / 100.,
                                self.input.design_params.j_ou,
                                self.results.turn_voltage)

        # if the resulting thickness of the winding is smaller than the required minimum the solution is not feasible
        if t_in < C_WIN_MIN or t_ou < C_WIN_MIN:
            self.results.feasible = INFEASIBLE
            return

        # outer winding radius
        r_ou = outer_winding_radius(r_in, t_in, self.input.design_params.m_gap, t_ou)

        # calculating the detailed parameters of the winding
        self.lv_winding = WindingDesign(inner_radius=r_in, thickness=t_in, winding_height=self.input.design_params.h_in,
                                        current_density=self.input.design_params.j_in,
                                        filling_factor=self.input.required.lv.filling_factor)
        self.lv_winding.calc_properties()

        self.hv_winding = WindingDesign(inner_radius=r_ou, thickness=t_ou, winding_height=h_ou,
                                        current_density=self.input.design_params.j_ou,
                                        filling_factor=self.input.required.hv.filling_factor)

        self.hv_winding.calc_properties()

        # window width and window heights
        self.results.window_width = window_width(self.input.required.min_core_gap, t_in, t_ou,
                                                 self.input.design_params.m_gap, 0, 0)

        self.results.wh = self.input.design_params.h_in + self.input.required.ei

        # core parameters
        self.results.core_mass = core_mass(self.input.design_params.rc, self.input.required.core_fillingf / 100.,
                                           self.input.design_params.h_in,
                                           self.input.required.ei, self.results.window_width,
                                           self.input.required.phase_distance / 2.)
        self.results.core_loss = core_loss_unit(self.input.design_params.bc, self.results.core_mass,
                                                CORE_BF)

        self.results.load_loss = self.lv_winding.ac_loss + self.lv_winding.dc_loss + self.hv_winding.ac_loss + self.hv_winding.dc_loss

        # short circuit impedance calculation with analytical formulas
        self.results.sci = short_circuit_impedance(self.input.required.power, 3.,
                                                   self.input.required.freq, self.input.required.alpha,
                                                   self.results.turn_voltage,
                                                   self.input.design_params.h_in, self.results.window_width, r_in,
                                                   t_in, r_ou,
                                                   t_ou, self.input.design_params.m_gap)

        self.results.capitalized_cost = capitalized_cost(self.results.core_mass, self.input.costs.core_cost,
                                                         self.lv_winding.mass, self.input.costs.lv_cost,
                                                         self.hv_winding.mass, self.input.costs.hv_cost,
                                                         self.results.load_loss, self.input.costs.ll_cost,
                                                         self.results.core_loss, self.input.costs.nll_cost)
