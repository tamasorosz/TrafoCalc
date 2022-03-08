from src.models import design_parameters
from src.transformer_calculations import *


class TransformerModel:

    def __init__(self, independent_params, transformer_data: dict):
        self.indep = independent_params
        self.transformer_data = transformer_data
        self.transformer_design = design_parameters

    def calc_phase_quantites(self):
        """
        This function calculates the phase voltages and phase currents in the hv and the lv windings and the phase power
        of the transformer.
        """

        # LV voltage winding
        # in this simple code the connection can be only 'y' or 'd'
        connection_factor = 1.
        if str(self.transformer_data['con_inner']).lower() == 'y':
            connection_factor = 1.73
        self.transformer_design['lv']['ph_voltage'] = self.transformer_data['line_voltage_inner']/connection_factor
        self.transformer_design['lv']['ph_current'] = phase_current(self.transformer_data['nominal_power'] * 1e-3,
                                                                    self.transformer_data['line_voltage_inner'],
                                                                    connection_factor)

        # HV voltage winding
        connection_factor = 1.
        if str(self.transformer_data['con_outer']).lower() == 'y':
            connection_factor = 1.73

        self.transformer_design['hv']['ph_voltage'] = self.transformer_data['line_voltage_outer']/connection_factor
        self.transformer_design['hv']['ph_current'] = phase_current(self.transformer_data['nominal_power'] * 1e-3,
                                                                    self.transformer_data['line_voltage_outer'],
                                                                    connection_factor)

        self.transformer_design['general']['ph_pow'] = self.transformer_data['nominal_power'] / self.transformer_data['phase_number']




def calc_dependent_variables(ind, para):
    """
    Calculates the non-magnetic/electrical parameters of the model.

    :param ind: the independent parameters of the model, which defines and 'individual' solution
    :param param: the required parameters or the different limits from the standards, technology ...
    :return: a dictionary with the dependent key-design variables
    """
    dep = Dependent_Params()

    # calculating derived parameters

    ## 1/ phase power
    #dep.ph_pow = para.power / para.ph_num

    # 2/ turn voltage
    dep.turn_voltage = turn_voltage(ind.b_c, ind.r_c, para.ff_c, para.freq)

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


def run_fem(param, dep, indep):
    """
    :param param: problem parameters
    :param dep: dependent variables
    :param indep: independent variables
    :return:
    """

    # problem
    problem = a2d.problem(clear=True)
    problem.coordinate_type = "axisymmetric"
    problem.mesh_type = "triangle"

    magnetic = problem.field("magnetic")
    magnetic.analysis_type = "steadystate"
    magnetic.number_of_refinements = 2
    magnetic.polynomial_order = 1
    magnetic.adaptivity_type = "disabled"
    magnetic.solver = "linear"

    # boundaries
    magnetic.add_boundary("A = 0", "magnetic_potential", {"magnetic_potential_real": 0})

    # materials
    magnetic.add_material("Air", {"magnetic_permeability": 1, "magnetic_conductivity": 0, "magnetic_remanence": 0,
                                  "magnetic_remanence_angle": 0,
                                  "magnetic_velocity_x": 0, "magnetic_velocity_y": 0, "magnetic_velocity_angular": 0,
                                  "magnetic_current_density_external_real": 0,
                                  "magnetic_total_current_prescribed": 0, "magnetic_total_current_real": 0})

    geometry = problem.geometry()

    # define regions inside the working window
    geometry.add_label(indep.r_c * 1e-3 + 0.01, 0.05, materials={"magnetic": "Air"})

    # print(dep.s, dep.wh)

    create_rectangle(geometry, indep.r_c, 0, dep.s, dep.wh, {"magnetic": "A = 0"})
    # TODO: the distance of the winding from the bottom core yoke is different from the endins/2

    # r_in, r_ou and r_reg are represents the mean diameters of the different windings, for the fem simulation
    # we have to calculate the inner radiuses
    r_in = dep.r_in - dep.t_in / 2.
    r_ou = dep.r_ou - dep.t_ou / 2.
    r_reg = dep.r_reg - dep.t_reg / 2.

    create_winding(magnetic, geometry, r_in, param.ei / 2., dep.t_in, indep.h_in, "lv", param.ff_in,
                   indep.j_in * 2. ** 0.5)
    create_winding(magnetic, geometry, r_ou, param.ei / 2., dep.t_ou, dep.h_ou, "hv", param.ff_ou,
                   -indep.j_ou * 2. ** 0.5)

    # checks that the regulating winding is activated or not
    if param.reg_range > 1e-6:
        create_winding(magnetic, geometry, r_reg, param.ei / 2., dep.t_reg, dep.h_reg, "reg", param.ff_reg, 0.)

    computation = problem.computation()
    computation.solve()
    solution = computation.solution("magnetic")

    # base impedance for the short circuit impedance calculation
    ub = param.u_in_line  # voltage --- kV
    sb = param.power / 1000.  # nominal power  --- MVA
    zb = ub ** 2. / sb  # impedance
    ib = param.power / ub / 1.73  # phase_current(sb, ub, param.con_fact_in)  # amps in the impedance base calculation
    Omega = 2. * pi * param.freq

    # -- calculating the impedance from the volume integrals --
    L = 2 * solution.volume_integrals()['Wm'] / ib ** 2.
    dep.sci = 2. * pi * param.freq * L / zb * 100.

    # calculating the average radial and axial flux density
    nx = 5  # number of the local inductances along the z coordinate
    ny = 5  # number of the local inductance along the x coordinate

    # low voltage winding  ----------------------------------
    Br_avg = 0.
    Bz_avg = 0.

    dx = dep.t_in / nx * 1e-3
    dy = indep.h_in / ny * 1e-3

    for i in range(0, nx):
        xx = dep.r_in * 1e-3 + i * dx
        for j in range(0, ny):
            yy = param.ei / 2. * 1e-3 + j * dy

            point = solution.local_values(xx, yy)
            Br_avg += point["Brr"]
            Bz_avg += point["Brz"]

    Br_avg /= (nx * ny)
    Bz_avg /= (nx * ny)

    b1 = (Br_avg ** 2. + Bz_avg ** 2.) ** 0.5

    # ---------------------------------------------------------
    # print('bz_avg:', Bz_avg)
    # print('br_avg:', Br_avg)

    # -- calculating the dc and ac losses from the model
    n_in = param.u_in / dep.turn_voltage * 1e3  # nr of the turns in the inner winding
    n_ou = param.u_out / dep.turn_voltage * 1e3  # nr of the turns in the inner winding
    Bg = b_gap(n_in, ib,
               indep.h_in)  # TODO <--- should be changed by the maximum of the inductance in the gap --- at least ...

    # inner winding -- ctc
    status, m, objval = gp_model_winding(n_in, dep.r_in + dep.t_in / 2., dep.t_in, indep.h_in, param.in_ins_rad,
                                         param.in_ins_s, param.in_ins_ax, ib, Omega, abs(Br_avg), abs(Bz_avg),
                                         2., param.ff_in, 2.5, 8., param.ph_num)
    w_in = m[W]
    print('w inner: ', m[W])
    print('h inner: ', m[H])

    # outer winding
    i_ou = phase_current(sb, param.u_out_line, param.con_fact_ou)

    # high voltage winding
    Br_avg = 0.
    Bz_avg = 0.

    dx = dep.t_ou / nx * 1e-3
    dy = dep.h_ou / ny * 1e-3

    for i in range(0, nx):
        xx = dep.r_ou * 1e-3 + i * dx
        for j in range(0, ny):
            yy = param.ei / 2. * 1e-3 + j * dy * 0.8

            point = solution.local_values(xx, yy)
            Br_avg += point["Brr"]
            Bz_avg += point["Brz"]

    Br_avg /= (nx * ny)
    Bz_avg /= (nx * ny)

    b2 = (Br_avg ** 2. + Bz_avg ** 2.) ** 0.5

    # print('bz_avg:', Bz_avg)
    # print('br_avg:', Br_avg)

    status, m, objval = gp_model_winding(n_in, dep.r_in + dep.t_in / 2., dep.t_in,
                                         indep.h_in, param.in_ins_rad, param.ou_ins_s, param.ou_ins_ax, ib, Omega,
                                         abs(Br_avg), abs(Bz_avg), 1., param.ff_ou, 2.7, 16., param.ph_num)

    w_ou = m[W]
    print('w outer: ', m[W])
    print('h outer: ', m[H])

    # high voltage upper segment

    Br_avg = 0.
    Bz_avg = 0.

    dx = dep.t_ou / nx * 1e-3
    dy = dep.h_ou / ny * 1e-3

    for i in range(0, nx):
        xx = dep.r_ou * 1e-3 + i * dx
        for j in range(0, ny):
            yy = param.ei / 2. * 1e-3 + j * dy * 0.2 + ny * dy * 0.8  # upper segment

            point = solution.local_values(xx, yy)
            Br_avg += point["Brr"]
            Bz_avg += point["Brz"]

    Br_avg /= (nx * ny)
    Bz_avg /= (nx * ny)

    b3 = (Br_avg ** 2. + Bz_avg ** 2.) ** 0.5

    # print('bz_avg:', Bz_avg)
    # print('br_avg:', Br_avg)
    # print(b1,b2,b3)
    status, m, objval = gp_model_winding(n_in, dep.r_in + dep.t_in / 2., dep.t_in,
                                         indep.h_in, param.in_ins_rad, param.ou_ins_s, param.ou_ins_ax, ib, Omega,
                                         abs(Br_avg), abs(Bz_avg), 1., param.ff_ou, 2.7, 16., param.ph_num)

    w_ou_up = m[W]
    print('w outer: ', m[W])
    print('h outer: ', m[H])

    fem_analytics_improved_model(indep, param, dep, w_in, b1 * 1.41, w_ou, b2 * 1.41, w_ou_up, b3 * 1.41)

    toc = eval_objective(param, dep, param.drop_tol, param.drop, dep.sci)

    # TODO: the calculation speed to claculate the masses also in the future objval[Vol]
    # print('TOC: ', toc, 'sci', dep.sci, 'drop input', param.drop)
    del problem
    del geometry
    del solution

    return toc
