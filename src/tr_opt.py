"""
script for calculating the load losses and the short circuit impedance
"""
from agrossuite import agros as a2d

#import agros as a2d
import numpy as np
#from floor_planning_wa2 import gp_model_winding, b_gap, Ploss, Vol, W, H
from optimization_functions import independent_variables, Dependent_Params, Parameters, phase_current, \
    calc_dependent_variables, fem_analytics_improved_model, eval_objective, short_circuit_impedance
from artap.algorithm_genetic import NSGAII
#from artap.results import Results, GraphicalResults
from math import pi
import matplotlib.pyplot as plt


from artap.problem import Problem

import tempfile

##########

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

    #print(dep.s, dep.wh)

    create_rectangle(geometry, indep.r_c, 0, dep.s, dep.wh, {"magnetic": "A = 0"})
    # TODO: the distance of the winding from the bottom core yoke is different from the endins/2

    # r_in, r_ou and r_reg are represents the mean diameters of the different windings, for the fem simulation
    # we have to calculate the inner radiuses
    r_in = dep.r_in - dep.t_in/2.
    r_ou = dep.r_ou - dep.t_ou/2.
    r_reg = dep.r_reg - dep.t_reg/2.

    create_winding(magnetic, geometry, r_in, param.ei/2., dep.t_in, indep.h_in, "lv", param.ff_in, indep.j_in*2.**0.5)
    create_winding(magnetic, geometry, r_ou, param.ei/2., dep.t_ou, dep.h_ou, "hv", param.ff_ou, -indep.j_ou*2.**0.5)

    # checks that the regulating winding is activated or not
    if param.reg_range > 1e-6:
        create_winding(magnetic, geometry, r_reg, param.ei/2., dep.t_reg, dep.h_reg, "reg", param.ff_reg, 0.)

    computation = problem.computation()
    computation.solve()
    solution = computation.solution("magnetic")

    # base impedance for the short circuit impedance calculation
    ub = param.u_in_line    # voltage --- kV
    sb = param.power/1000.  # nominal power  --- MVA
    zb = ub ** 2. / sb      # impedance
    ib = param.power/ub/1.73  #phase_current(sb, ub, param.con_fact_in)  # amps in the impedance base calculation
    Omega = 2. * pi * param.freq

    # -- calculating the impedance from the volume integrals --
    L = 2 * solution.volume_integrals()['Wm'] / ib ** 2.
    dep.sci = 2. * pi * param.freq * L / zb * 100.

    # calculating the average radial and axial flux density
    nx = 5 # number of the local inductances along the z coordinate
    ny = 5 # number of the local inductance along the x coordinate

    # low voltage winding  ----------------------------------
    Br_avg = 0.
    Bz_avg = 0.

    dx = dep.t_in/nx*1e-3
    dy = indep.h_in/ny*1e-3

    for i in range(0, nx):
        xx = dep.r_in*1e-3 + i * dx
        for j in range(0, ny):
            yy = param.ei/2.*1e-3 + j * dy

            point = solution.local_values(xx, yy)
            Br_avg += point["Brr"]
            Bz_avg += point["Brz"]

    Br_avg /= (nx*ny)
    Bz_avg /= (nx*ny)

    b1 = (Br_avg ** 2. + Bz_avg ** 2.) ** 0.5

    # ---------------------------------------------------------
    #print('bz_avg:', Bz_avg)
    #print('br_avg:', Br_avg)

    # -- calculating the dc and ac losses from the model
    n_in = param.u_in / dep.turn_voltage *1e3  # nr of the turns in the inner winding
    n_ou = param.u_out / dep.turn_voltage*1e3  # nr of the turns in the inner winding
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

    b2 = (Br_avg**2.+Bz_avg**2.)**0.5

    #print('bz_avg:', Bz_avg)
    #print('br_avg:', Br_avg)

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
            yy = param.ei / 2. * 1e-3 + j * dy *0.2 + ny * dy * 0.8 # upper segment

            point = solution.local_values(xx, yy)
            Br_avg += point["Brr"]
            Bz_avg += point["Brz"]

    Br_avg /= (nx * ny)
    Bz_avg /= (nx * ny)

    b3 = (Br_avg ** 2. + Bz_avg ** 2.) ** 0.5

    #print('bz_avg:', Bz_avg)
    #print('br_avg:', Br_avg)
    #print(b1,b2,b3)
    status, m, objval = gp_model_winding(n_in, dep.r_in + dep.t_in / 2., dep.t_in,
                                         indep.h_in, param.in_ins_rad, param.ou_ins_s, param.ou_ins_ax, ib, Omega,
                                         abs(Br_avg), abs(Bz_avg), 1., param.ff_ou, 2.7, 16., param.ph_num)

    w_ou_up = m[W]
    print('w outer: ', m[W])
    print('h outer: ', m[H])

    fem_analytics_improved_model(indep, param, dep, w_in, b1*1.41, w_ou, b2*1.41, w_ou_up, b3*1.41)

    toc = eval_objective(param, dep, param.drop_tol, param.drop, dep.sci)

    # TODO: the calculation speed to claculate the masses also in the future objval[Vol]
    #print('TOC: ', toc, 'sci', dep.sci, 'drop input', param.drop)
    del problem
    del geometry
    del solution

    return toc


# problem parameters
class TestProblem(Problem):

    def __init__(self, name):
        parameters = {'r_c': {'initial_value': 250., 'bounds': [200., 600.]},
                      'b_c': {'initial_value': 1.6, 'bounds': [1.4, 1.7]},
                      'j_in': {'initial_value': 2.0, 'bounds': [1.5, 3.5]},
                      'j_ou': {'initial_value': 2.0, 'bounds': [1.5, 3.5]},
                      'j_reg': {'initial_value': 3., 'bounds': [2.0, 3.8]},
                      'h_in': {'initial_value': 1500., 'bounds': [800., 3000.]},
                      'm_gap': {'initial_value': 37., 'bounds': [37., 100.]}}
        costs = ['F']

        super().__init__(name, parameters, costs)

    def evaluate(self, x: list):

        ind = independent_variables(r_c=x[0], b_c=x[1],
                                    j_in=x[2], j_ou=x[3], j_reg=x[4],
                                    h_in=x[5],
                                    m_gap=x[6])

        # define parameters
        p = Parameters()
        p.power = 31500.
        p.freq = 50.
        p.ph_num = 3.
        p.drop_tol = 0.05
        p.drop = 14.5 # expected
        p.gap = 37.
        p.gap_core = 20.
        p.ei = 150.
        p.phase_distance = 37.
        p.ff_c = 0.9
        p.cc = 3.5
        p.u_in_line = 33.  # line voltage in the inner winding
        p.u_in = 19.05
        p.u_out = 132.
        p.ff_in = 0.6
        p.ff_ou = 0.60
        p.win_c = 10.
        p.wout_c = 8.5
        p.reg_range = 0.
        p.ff_reg = 0.65
        p.alpha = 0.97
        p.beta = 0.85

        p.ll_c = 2000.
        p.nll_c = 7000.


        p.con_fact_in = 1.
        p.con_fact_ou = 1.73
        p.u_out = 132. # phase voltage on the outer terminal [kV]
        p.u_out_line = 132.

        #### inner winding definition
        p.in_ins_ax = 3.5
        p.in_ins_rad = 0.45
        p.in_ins_s = 0.17  # ctx

        #### outer winding
        p.ou_ins_ax = 4.8
        p.ou_ins_rad = 0.5
        p.ou_ins_s = 0.3  # triple twins

        dep = calc_dependent_variables(ind, p)
        toc = run_fem(p, dep, ind)

        return [toc]


class TestProblem2(Problem):

    def __init__(self, name):
        parameters = [{'name':'r_c', 'bounds': [150., 300.]},
                      {'name':'b_c', 'initial_value': 1.6, 'bounds': [1.4, 1.7]},
                      {'name':'j_in','initial_value': 2.0, 'bounds': [1.5, 3.]},
                      {'name':'j_ou','initial_value': 2.0, 'bounds': [1.5, 3.]},
                      {'name':'j_reg', 'initial_value': 3., 'bounds': [2.0, 3.5]},
                      {'name':'h_in','initial_value': 1500., 'bounds': [1000., 2000.]},
                      {'name':'m_gap','initial_value': 37., 'bounds': [24., 70.]}]

        costs = [{'name':'F'}]

        super().__init__(name, parameters, costs)

    def evaluate(self, x: list):

        ind = independent_variables(r_c=x[0], b_c=x[1],
                                    j_in=x[2], j_ou=x[3], j_reg=x[4],
                                    h_in=x[5],
                                    m_gap=x[6])

        # define parameters
        # define parameters
        p = Parameters()
        p.power = 10000.
        p.freq = 50.
        p.ph_num = 3.
        p.drop_tol = 0.03
        p.drop = 7.5  # expected
        p.gap = 38.
        p.gap_core = 20.
        p.ei = 80.
        p.phase_distance = 50.
        p.ff_c = 0.88
        p.cc = 3.5
        p.u_in_line = 6.9  # line voltage in the inner winding
        p.u_in = 6.9 / 1.73
        p.ff_in = 0.7
        p.ff_ou = 0.60
        p.win_c = 10.
        p.wout_c = 8.5
        p.reg_range = 0.
        p.ff_reg = 0.65
        p.alpha = 1.0
        p.beta = 0.85

        p.ll_c = 2500.
        p.nll_c = 5000.

        p.con_fact_in = 1.
        p.con_fact_ou = 1.
        p.u_out = 33. / 1.73  # phase voltage on the outer terminal [kV]
        p.u_out_line = 33.

        #### inner winding definition
        p.in_ins_ax = 3.
        p.in_ins_rad = 0.45
        p.in_ins_s = 0.17  # ctx

        #### outer winding
        p.ou_ins_ax = 3.5
        p.ou_ins_rad = 0.5
        p.ou_ins_s = 0.3  # triple twins

        dep = calc_dependent_variables(ind, p)
        toc = run_fem(p, dep, ind)

        return [toc]


class TestProblem3(Problem):

    def set(self):
        self.name = 'TrOPT Test Problem'
        self.working_dir = "./"
        self.parameters = [{'name':'r_c', 'bounds': [150., 400.]},
                           {'name':'b_c', 'initial_value': 1.6, 'bounds': [1.4, 1.7]},
                           {'name':'j_in','initial_value': 2.0, 'bounds': [1.5, 3.5]},
                           {'name':'j_ou','initial_value': 2.0, 'bounds': [1.5, 3.5]},
                           {'name':'j_reg', 'initial_value': 3., 'bounds': [2.0, 3.8]},
                           {'name':'h_in','initial_value': 1500., 'bounds': [1000., 2500.]},
                           {'name':'m_gap','initial_value': 37., 'bounds': [38., 70.]}]

        self.costs = [{'name': 'F_1', 'criteria': 'minimize'}, {'name': 'F_2', 'criteria': 'maximize'}]

    def evaluate(self, individual):
        x = individual.vector
        ind = independent_variables(r_c=x[0], b_c=x[1],
                                    j_in=x[2], j_ou=x[3], j_reg=x[4],
                                    h_in=x[5],
                                    m_gap=x[6])
        individual.sci = 0.
        # define parameters
        # define parameters
        p = Parameters()
        p.power = 50000.
        p.freq = 50.
        p.ph_num = 3.
        p.drop_tol = 0.05 #1e6
        p.drop = 12.5  # expected
        p.gap = 38.
        p.gap_core = 20.
        p.ei = 300.
        p.phase_distance = 50.
        p.ff_c = 0.89
        p.cc = 3.5
        p.u_in_line = 20.  # line voltage in the inner winding
        p.u_in = 20.
        p.ff_in = 0.65
        p.ff_ou = 0.60
        p.win_c = 10.
        p.wout_c = 8.5
        p.reg_range = 10.
        p.ff_reg = 0.65
        p.alpha = 0.95
        p.beta = 0.85

        p.ll_c = 0.
        p.nll_c = 0.

        p.con_fact_in = 1.
        p.con_fact_ou = 1.
        p.u_out = 130. / 1.73  # phase voltage on the outer terminal [kV]
        p.u_out_line = 130.

        #### inner winding definition
        p.in_ins_ax = 3.
        p.in_ins_rad = 0.45
        p.in_ins_s = 0.17  # ctx

        #### outer winding
        p.ou_ins_ax = 3.5
        p.ou_ins_rad = 0.5
        p.ou_ins_s = 0.3  # triple twins

        dep = calc_dependent_variables(ind, p)
        toc = run_fem(p, dep, ind)
        individual.sci = dep.sci
        return [toc, x[6]]


###########

def run_monte_carlo():

    problem = TestProblem3()
    problem.data_store = FileDataStore(problem, database_name='monte_carlo_exclude_12.5.db', mode="write")

    algorithm = NSGAII(problem)
    algorithm.options['max_population_number'] = 1
    algorithm.options['max_population_size'] = 1000

    algorithm.options['verbose_level'] = 1
    algorithm.run()


def process_results():
    problem = TestProblem3()
    #database_names =  ["./monte_carlo_7_6.db", "./monte_carlo_7_2.db",  "./monte_carlo_0_0.db"]
    database_names =  ["./monte_carlo_exclude_12.5.db"]

    colors = ['bo', 'rx', 'g.']
    plt.figure(1)
    n = 60
    grid = np.linspace(38, 70, n)
    table = []
    for i in range(n):
        table.append([])

    for i in range(3):
        problem.data_store = FileDataStore(problem, database_name=database_names[i], mode="read")

        for population in problem.data_store.populations:
            for individual in population.individuals:
                if individual.costs[0] > -1:
                    table[round(individual.vector[6])-38].append(individual)

        for row in table:
            if len(row) is not 0:
                sorted(row, key=lambda x: individual.costs[0])
                plt.semilogy(row[-1].vector[6], row[-1].costs[0], colors[i])
                plt.grid()

    # plt.xlim([0, 0.25])
    # plt.ylim([0, 0.001e9])
    plt.xlabel("Main gap")
    plt.ylabel("Total Cost of Ownership")

    plt.figure(2)
    for i in range(3):
        problem.data_store = FileDataStore(problem, database_name=database_names[i], mode="read")
        for population in problem.data_store.populations:
            for individual in population.individuals:
                plt.semilogy(individual.sci, individual.costs[0], colors[i])

    plt.xlim([0, 50])
    # plt.ylim([0, 0.001e9])

    plt.grid()
    plt.xlabel("SCI")
    plt.ylabel("Total Cost of Ownership")
    plt.show()


def Test_6300kVA_transformer():

    # Test example from KKK
    indep = ind = independent_variables(r_c=184.0, b_c=1.568,
                                        j_in= 3.02 , j_ou= 3.0, j_reg=1.73,
                                        h_in=979., m_gap=26.7)
    #
    # # optimum 1 -- k1, k2 = [7000, 1000]
    # indep = ind = independent_variables(r_c=206., b_c=1.5,
    #                                    j_in= 2.29 , j_ou= 2.189, j_reg=2.175,
    #                                    h_in=1254., m_gap=28.)
    #
    # # define parameters
    p.ff_c = 0.838
    p.cc = 3.5
    p.u_in_line = 22.  # line voltage in the inner winding
    p.u_in = 22.
    p.ff_in = 0.535
    p.ff_ou = 0.56
    p.win_c = 10.
    p.wout_c = 8.5
    p.reg_range = 0.
    p.ff_reg = 0.65
    p.alpha = 1.
    p.beta = 0.85

    p.ll_c = 2000.
    p.nll_c = 6000.

    p.con_fact_in = 1.
    p.con_fact_ou = 1.
    p.u_out = 33./1.73  # phase voltage on the outer terminal [kV]
    p.u_out_line = 33.
    #
    # #### inner winding definition
    p.in_ins_ax = 3.5
    p.in_ins_rad = 0.45
    p.in_ins_s = 0.17  # ctx
    #
    # #### outer winding
    p.ou_ins_ax = 4.8
    p.ou_ins_rad = 0.5
    p.ou_ins_s = 0.3  # triple twins
    #
    dep = calc_dependent_variables(ind, p)
    toc = run_fem(p, dep, ind)
    print('load loss:    ', dep.load_loss, '\n'
          'turn voltage: ', dep.turn_voltage, '\n'
          'sci:          ', dep.sci, '\n'
          't_in: ', dep.t_in, '\n'
          't_ou: ', dep.t_ou, '\n'
          'm_in: ', dep.m_in, '\n'
          'dc_in: ', dep.dc_in, '\n'
          'ec_in: ', dep.ec_in, '\n'
          'd_m: ', 2*dep.r_in, '\n' 
          'd_m: ', 2*dep.r_ou, '\n' 
          'dc_ou: ', dep.dc_ou, '\n'
          'ec_ou: ', dep.ec_ou, '\n'
          'm_ou:', dep.m_ou, '\n'
          'core mass: ', dep.c_mass, '\n'
          'core loss: ', 1)

    print(dep.m_ou)
    print(dep.m_in)

    return

if __name__ == "__main__":

    Test_6300kVA_transformer()
    #run_monte_carlo()
