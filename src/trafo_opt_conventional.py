from src.models import TransformerDesign, MainResults, WindingDesign, WindingParams
from src.transformer_calculations import *

from dataclasses import dataclass, field
from dataclasses_json import dataclass_json

C_WIN_MIN = 10.0  # [mm] technological limit for the thickness of the windings, it should be larger than 10 mm-s


@dataclass_json
@dataclass
class TwoWindingModel:
    input: TransformerDesign
    hv_winding: WindingDesign  # derived winding parameters and data
    lv_winding: WindingDesign
    results: MainResults

    def calculate(self):
        # 1) phase power, assumes a 3 phased 3 legged transformer core
        ph_voltage = self.input.required.power / 3.  # [kVA]

        # 2) turn voltage
        self.results.turn_voltage = turn_voltage(self.input.design_params.bc, self.input.design_params.rc,
                                                 self.input.required.core_fillingf / 100., self.input.required.freq)

        # 3/ winding thickness - inner
        self.lv_winding.thickness
        dep.t_in = calc_inner_width(dep.ph_pow, ind.h_in, para.ff_in, ind.j_in, dep.turn_voltage)

        # check the 'strength' of the coil if it's smaller than a technological limit, the solution is infeasible
        if dep.t_in < C_WIN_MIN:
            dep.feasible = INFEASIBLE
            return copy(dep)


def calc_dependent_variables(ind, para):
    """
    Calculates the non-magnetic/electrical parameters of the model.

    :param ind: the independent parameters of the model, which defines and 'individual' solution
    :param param: the required parameters or the different limits from the standards, technology ...
    :return: a dictionary with the dependent key-design variables
    """

    # 3/ winding thickness - inner
    dep.t_in = calc_inner_width(dep.ph_pow, ind.h_in, para.ff_in, ind.j_in, dep.turn_voltage)

    # check the 'strength' of the coil if it's smaller than a technological limit, the solution is infeasible
    if dep.t_in < C_WIN_MIN:
        dep.feasible = INFEASIBLE
        return copy(dep)

    # 4/ winding thickness -- outer winding
    dep.h_ou = ind.h_in * para.alpha
    dep.t_ou = calc_inner_width(dep.ph_pow, dep.h_ou, para.ff_ou, ind.j_ou, dep.turn_voltage)

    # check the strength of the coil
    if dep.t_ou < C_WIN_MIN:
        dep.feasible = INFEASIBLE
        return copy(dep)

    # 5/ regulating winding thickness
    dep.h_reg = ind.h_in * para.beta

    if para.reg_range > 1e-6 and para.bin_reg is False:
        dep.t_reg = calc_t_reg(para.reg_range, ind.h_in, para.beta,
                               dep.ph_pow, ind.j_reg, para.ff_reg, dep.turn_voltage)
    else:
        if para.bin_reg is True:
            dep.t_reg = 0.

    # 6/ innder winding radius
    dep.r_in = inner_winding_radius(ind.r_c, para.gap_core, dep.t_in)

    # 7/ outer winding radius
    dep.r_ou = outer_winding_radius(dep.r_in, dep.t_in, ind.m_gap, dep.t_ou)

    # 8/ regulating winding radius
    if para.reg_range > 1e-6 and para.bin_reg is False:
        dep.r_reg = rad_reg_winding_outer(dep.r_ou, dep.t_ou, para.gap, dep.t_reg)
    else:
        if para.bin_reg is True:
            dep.r_reg = 0.

        if para.bin_reg is True:
            dep.t_reg = dep.t_ou
            dep.h_ou = dep.h_ou * 1. - para.reg_range

    # 9/ calculating the window width wothout the regulating winding
    if para.reg_range > 1e-6:
        gap_r = ind.m_gap
    else:
        gap_r = 0.

    dep.s = window_width(para.gap_core, dep.t_in, dep.t_ou, ind.m_gap, dep.t_reg, gap_r,
                         para.phase_distance)  # TODO: this variable is missing somehow from the previous optimization model

    # 10/ calculating winding masses and losses
    dep.m_in = winding_mass(para.ph_num, dep.r_in, dep.t_in, ind.h_in, para.ff_in)
    dep.m_ou = winding_mass(para.ph_num, dep.r_ou, dep.t_ou, dep.h_ou, para.ff_ou)

    if para.reg_range > 1e-6:
        dep.m_reg = winding_mass(para.ph_num, dep.r_reg, dep.t_reg, dep.h_reg, para.ff_reg)
    else:
        dep.m_reg = 0.
    # 11/ calculating the core mass
    dep.c_mass = core_mass(ind.r_c, para.ff_c, ind.h_in, para.ei, dep.s, para.phase_distance / 2.)
    dep.core_loss = core_loss_unit(ind.b_c, dep.c_mass, para.f_bf)

    # 12/ window height
    dep.wh = ind.h_in + para.ei

    return copy(dep)
