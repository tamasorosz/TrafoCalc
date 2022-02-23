from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller.objects import Rectangle

from copy import copy

"""
Sources: https://www.mdpi.com/2076-3417/10/4/1361
         https://pp.bme.hu/eecs/article/view/10207/7294
         
Input data:
 
    Design Parameters
        filling_inner   - Filling factor of the copper in the inner (mostly LV) winding
        filling_outer   - Filling factor of the copper in the outer (mostly HV) winding
        end_insulation  - sum of the end insulation
        core_insulation - the thickness of the core insulation
        
    Independent Variables for the two winding based transformer model:
         Rc, Jin, Jou, Gap
"""

ModelDir.set_base(__file__)


class PowerTransformer(BaseModel):
    """
    This class waits for two input keys: params and vars, params contains the design parameters, which values are fix
    for all of the calculated individuals during the optimization of the transformers. The vars  are the design variables,
    which used together with the design parameters to calculate the transformer data.
    """
    design_parameters = {'ff_in': 0,  # filling factor in the inner winding [%]
                         'ff_ou': 0,  # filling factor in the outer winding [%]
                         'alpha': 0,  # the ratio between the inner and the outer winding [-]
                         'freq': 0,  # network frequency
                         'end_ins': 0,  # the sum of the end insulation under and on the windings [mm]
                         'core_ins': 0}  # core insulation distance between the core and the windings [mm]

    design_variables = {'core_diam': 0,  # core diameter in [mm]
                        'gap': 0,  # main insulation distance between the main windings in [mm]
                        'hin': 0,  # the height of the inner winding [mm]
                        'tin': 0,  # the thickness of the inner winding [mm]
                        'tou': 0,  # the thickness of the outer winding [mm]
                        'jin': 0,  # the current density in the inner winding [A/mm2]
                        'jou': 0}  # the current density in the outer winding [A/mm2]

    def __init__(self, **kwargs):
        super(PowerTransformer, self).__init__(**kwargs)

        # self._init_directories()

        # The default values of the example script are coming from the validation example

        # design parameters
        self.design_parameters['ff_in'] = kwargs.get('ff_in', 0)
        self.design_parameters['ff_ou'] = kwargs.get('ff_ou', 0)
        self.design_parameters['alpha'] = kwargs.get('alpha', 0)
        self.design_parameters['freq'] = kwargs.get('freq', 0)
        self.design_parameters['end_ins'] = kwargs.get('end_ins', 0)
        self.design_parameters['core_ins'] = kwargs.get('core_ins', 0)

        # design variables
        self.design_variables['core_diam'] = kwargs.get('core_diam', 0)
        self.design_variables['gap'] = kwargs.get('gap', 0)
        self.design_variables['hin'] = kwargs.get('hin', 0)
        self.design_variables['tin'] = kwargs.get('tin', 0)
        self.design_variables['tou'] = kwargs.get('tou', 0)
        self.design_variables['jin'] = kwargs.get('jin', 0)
        self.design_variables['jou'] = kwargs.get('jou', 0)

        # calculating the geometrical parameters from the given data
        self.set_window_parameters()
        self.set_inner_winding_parameters()
        self.set_outer_winding_parameters()

    def set_window_parameters(self):
        # inner window of the simulated power transformer
        # the description of the variables can be found in the connected drawing
        self.w1 = self.design_parameters['core_ins'] + self.design_variables['tin'] + 2. * self.design_variables[
            'gap'] + \
                  self.design_variables['tou']
        self.h1 = self.design_parameters['end_ins'] + self.design_variables['hin']
        self.r1 = self.design_variables['core_diam'] / 2.
        self.z1 = 0

    def set_inner_winding_parameters(self):
        self.w2 = self.design_variables['tin']
        self.h2 = self.design_variables['hin']
        self.r2 = self.r1 + self.design_parameters['core_ins']
        self.z2 = self.design_parameters['end_ins'] / 2

        self.js = self.design_variables['jin'] * 1_000_000 * self.design_parameters['ff_in'] / 100.  # A/m2

    def set_outer_winding_parameters(self):
        # HV rectangle
        self.w3 = self.design_variables['tou']
        self.h3 = self.design_variables['hin'] * self.design_parameters['alpha']
        self.r3 = self.r2 + self.design_variables['tin'] + self.design_variables['gap']
        self.z3 = self.design_parameters['end_ins'] / 2

        # Excitation
        self.jp = self.design_variables['jou'] * 1_000_000 * self.design_parameters['ff_ou'] / 100.

    def setup_solver(self):
        # a 2-dimensional, axisymmetric FEMM model is used for the basic simulations
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material('air')
        coil = Material('coil')

        LV = copy(coil)
        LV.name = 'LV'
        LV.Je = self.js

        HV = copy(coil)
        HV.name = 'HV'
        HV.Je = -self.jp  # the current direction shuold be different and the amperturns should be balanced

        self.snapshot.add_material(air)
        self.snapshot.add_material(LV)
        self.snapshot.add_material(HV)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def add_postprocessing(self):
        points = [((self.r1 + self.r2) / 2, self.z1 + self.h1 / 2)]
        self.snapshot.add_postprocessing("integration", points, "Energy")

    def build_geometry(self):
        """
        The goemetry represents the inner window of the calculated power transformer, which contains three
        rectangles:
         - the inner side of the core window
         - the inner winding
         - the outer winding
        """

        r1 = Rectangle(self.r1, self.z1, width=self.w1, height=self.h1)
        r2 = Rectangle(self.r2, self.z2, width=self.w2, height=self.h2)
        r3 = Rectangle(self.r3, self.z3, width=self.w3, height=self.h3)

        self.geom.add_rectangle(r1)
        self.geom.add_rectangle(r2)
        self.geom.add_rectangle(r3)

        self.snapshot.add_geometry(self.geom)

        self.assign_material(*r2.cp, 'LV')
        self.assign_material(*r3.cp, 'HV')
        self.assign_material((self.r1 + self.r2) / 2, self.z1 + self.h1 / 2, 'air')

        self.assign_boundary(*r1.a.mean(r1.b), 'a0')
        self.assign_boundary(*r1.d.mean(r1.c), 'a0')
        self.assign_boundary(*r1.a.mean(r1.d), 'a0')
        self.assign_boundary(*r1.b.mean(r1.c), 'a0')


if __name__ == "__main__":
    m = PowerTransformer(ff_in=60, ff_ou=60, alpha=1.0, end_ins=180., core_ins=20., core_diam=420., gap=50.,
                         hin=1100, tin=35, tou=44, jin=2.6, jou=2.26, exportname="dev")
    #m = PowerTransformer(exportname="dev", jin=3.0, jou=3.02)
    print(m(cleanup=False, devmode=True))
