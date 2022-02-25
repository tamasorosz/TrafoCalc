from src.transformer_calculatoins import phase_current, winding_mass, winding_dc_loss, opt_win_eddy_loss, \
    homogenous_insulation_ff

# Calculated parameters
# ---------------------
winding_params = {'inner_radius': 0, 'mean_radius': 0, 'outer_radius': 0, 'thickness': 0, 'height': 0, 'mass': 0,
                  'dc_loss': 0, 'ac_loss': 0, 'current_density': 0, 'ph_current': 0, 'ph_voltage': 0}

core_params = {'diameter': 0, 'core_loss': 0, 'turn_voltage': 0, 'core_mass': 0, 'window_height': 0}

general_params = {'sci': 0, 'turn_voltage': 0, 'window_width': 0, 'ph_pow': 0, 'load_loss': 0, 'core_mass': 0,
                  'feasible': False}

# dictionary for storing the details of the calculated parameters of the transformer design
design_parameters = {'general': general_params, 'core': core_params, 'lv': winding_params, 'hv': winding_params}


class WindingModel:

    def __init__(self, r_in=0, h_in=0, t_in=0, con='y', ul=0, j=0, ff=0):
        self.connection = con
        self.line_voltage = ul
        self.current_density = j
        self.filling_factor = ff  # ff is given in %
        self.rin = r_in  # inner radius of the winding
        self.hin = h_in  # height of the winding
        self.tin = t_in  # thickness of the winding

    def calculate_phase_quantities(self, nominal_power):
        """
        :param nominal_power: in kVA
        :return:
        """
        # in this simple code the connection can be only 'y' or 'd'
        connection_factor = 1.
        if str(self.connection).lower() == 'y':
            connection_factor = 1.73

        self.ph_current = phase_current(nominal_power * 1e-3, self.line_voltage, connection_factor)
        self.ph_voltage = self.line_voltage / connection_factor

    def calc_properties(self, ph_num):
        # geometry
        self.mean_radius = self.rin + self.tin / 2.
        self.outer_radius = self.rin + self.tin

        self.mass = winding_mass(ph_num, self.mean_radius, self.tin, self.hin, self.filling_factor / 100.)
        self.dc_loss = winding_dc_loss(self.mass, self.current_density)
        # the ac loss of the winding with optimal conductor layout
        self.ac_loss = opt_win_eddy_loss(self.tin * homogenous_insulation_ff(self.filling_factor/100.), self.tin)
