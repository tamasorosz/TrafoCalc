from dataclasses import dataclass, field

from dataclasses_json import dataclass_json

from src.base_functions import homogenous_insulation_ff, opt_win_eddy_loss, phase_current, winding_dc_loss, winding_mass
from src.superconductor_losses import perp_loss, parallel_loss, norris_equation
from src.base_functions import C_RHO_BSSCO


@dataclass_json
@dataclass
class IndependentVariables:
    rc: float  # core radius in [mm]
    bc: float  # flux density in the core [T]
    j_in: float  # current density [A/mm2]
    j_ou: float  # current density [A/mm2]
    h_in: float  # height of the inner winding [mm]
    m_gap: float  # main gap - main insulation distance between the windings [mm]


@dataclass_json
@dataclass
class WindingParams:
    # required parameters
    connection: str  # it can be 'y' or delta
    line_voltage: float  # the line voltage in the winding [kV]
    filling_factor: float = field(default=0.0)  # the expected copper filling of the given winding
    # calculated parameters
    ph_current = 0.0
    ph_voltage = 0.0

    def calculate_phase_quantities(self, nominal_power):
        """
        Calculates the phase quantities to the given line voltages and line currents.
        :param nominal_power: in kVA
        :return:
        """
        # in this simple code the connection can be only 'y' or 'd'
        if str(self.connection).lower() == "y":
            connection_factor = 1.73

        elif str(self.connection).lower() == "d":
            connection_factor = 1.0

        else:
            raise ValueError("invalid connection type")

        self.ph_current = phase_current(nominal_power * 1e-3, self.line_voltage, connection_factor)
        self.ph_voltage = self.line_voltage / connection_factor


C_WIN_MIN = 10.0  # [mm] technological limit for the thickness of the windings, it should be larger than 10 mm-s
INFEASIBLE = -1


@dataclass_json
@dataclass
class WindingDesign:
    inner_radius: float = field(default=0.0)
    thickness: float = field(default=0.0)  # input - design variable
    winding_height: float = field(default=0.0)  # input - design variable
    filling_factor: float = field(default=0.0)  # input - engineers guess from the constructoins variable
    current_density: float = field(default=0.0)  # input - design variable
    # calculated parameters
    mass: float = field(default=0.0)
    dc_loss: float = field(default=0.0)
    ac_loss: float = field(default=0.0)
    outer_radius: float = field(default=0.0)
    cable_length: float = field(default=0.0)
    amper_turns: float = field(default=0.0)

    def calc_properties(self):
        # geometry
        self.mean_radius = self.inner_radius + self.thickness / 2.0
        self.outer_radius = self.inner_radius + self.thickness

        self.mass = winding_mass(3, self.mean_radius, self.thickness, self.winding_height, self.filling_factor / 100.0)

        self.dc_loss = winding_dc_loss(self.mass, self.current_density)
        self.ac_loss = (
                opt_win_eddy_loss(self.thickness * homogenous_insulation_ff(self.filling_factor / 100.0),
                                  self.thickness)
                * self.dc_loss
        )
        self.amper_turns = round(
            self.thickness * self.filling_factor / 100.0 * self.winding_height * self.current_density, 1)
        return True

    def calc_sc_properties(self, freq):
        # geometry
        self.mean_radius = self.inner_radius + self.thickness / 2.0
        self.outer_radius = self.inner_radius + self.thickness
        self.mass = winding_mass(3, self.mean_radius, self.thickness, self.winding_height, self.filling_factor / 100.0,
                                 material='BSSCO')

        ac = 0.31 * 4.1  # mm2
        self.cable_length = round(self.mass / C_RHO_BSSCO / ac *1e-3, 2) # m
        self.dc_loss = 0.  # superconducting 'loss' assumed in the normal state

        # the ac loss calculated for the assumed properties of
        self.ac_loss = perp_loss(1.35, freq, 4.29 * 1e-3, 0.011) * self.cable_length

        self.amper_turns = round(
            self.thickness * self.filling_factor / 100.0 * self.winding_height * self.current_density, 1)


@dataclass_json
@dataclass
class MaterialCosts:
    ll_cost: float = field(default=0.0)  # load loss cost of the given design
    nll_cost: float = field(default=0.0)  # no load loss cost of the given design
    lv_cost: float = field(default=0.0)  # cost of 1kg copper in the lv winding
    hv_cost: float = field(default=0.0)  # cost of 1kg copper in the hv winding
    core_cost: float = field(default=0.0)  # cost of the magnetic steel


@dataclass_json
@dataclass
class TransformerRequirements:
    """
    Required parameters for a 3 phase transformer design and technological parameters.
    """

    power: float  # nominal power in MVA
    freq: float  # network frequency in Hz
    sci_req: float  # the % value of the required short circuit impedance
    drop_tol: float  # the tolerance value of the transformer
    hv: WindingParams  # hv terminal parameters
    lv: WindingParams  # lv terminal parameters
    min_main_gap: float  # the minimal insulation distance of the main gap at the given voltage level (*)
    min_core_gap: float  # the minimal insulation distance between the core and the
    ei: float  # the total height of the end insulation below and above the windings
    phase_distance: float  # the minimum thickness of the phase distance between the windings
    alpha: float  # the ratio of the HV/LV windings
    core_fillingf: float  # % value of the steel sheets in the stacked core


@dataclass_json
@dataclass
class TransformerDesign:
    """This class contains the requiered parameters and the optimal values of the transformer design."""

    description: str
    required: TransformerRequirements
    costs: MaterialCosts
    design_params: IndependentVariables


@dataclass_json
@dataclass
class MainResults:
    wh: float = field(default=0.0)
    feasible: bool = field(default=False)
    core_loss: float = field(default=0.0)
    load_loss: float = field(default=0.0)
    copper_mass: float = field(default=0.0)
    sci: float = field(default=0.0)
    window_width: float = field(default=0.0)
    core_mass: float = field(default=0.0)
    turn_voltage: float = field(default=0.0)
    capitalized_cost: float = field(default=0.0)
    fem_based_sci: float = field(default=0.0)
    fem_bax: float = field(default=0.0)
    fem_brad: float = field(default=0.0)
    fem_bax_brad:list = field(default=0.0)
