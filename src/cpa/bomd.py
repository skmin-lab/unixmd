from __future__ import division
from cpa.cpa import CPA

class BOMD(CPA):
    """ Class for born-oppenheimer molecular dynamics (BOMD) sampling
        to run dynamics with Classical Path Approximation (CPA)

        :param object molecule: Molecule object
        :param object thermostat: Thermostat object
        :param integer istate: Electronic state
        :param double dt: Time interval
        :param integer nsteps: Total step of nuclear propagation
        :param string samp_dir: Sampling data directory
        :param string unit_dt: Unit of time step
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, \
            samp_dir="./", unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, False, \
            unit_dt, out_freq, verbosity)

        self.samp_dir = samp_dir
        if (not os.path.exists(self.samp_dir)):
            os.makedirs(self.samp_dir)
        else:
            error_message = "File already exists!"
            error_vars = f"samp_dir = {self.samp_dir}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error message} ( {error_vars} )")

    def run(self):
        """ Run MQC dynamics according to BOMD
        """
        pass

    def save_bin(self):
        """ Routine to save MD info of each step using pickle
        """

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        pass

    def update_energy(self):
        """ Routine to update the energy of molecules in BOMD
        """
        pass

    def print_init(self):
        """ Routine to print the initial information of dynamics
        """
        pass

    def print_step(self):
        """ Routine to print each steps information about dynamics
        """
        pass

