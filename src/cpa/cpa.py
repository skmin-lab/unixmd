from __future__ import division

class CPA(object):
    """ Class for nuclear propagator used to BOMD sampling

        :param object molecule: Molecule object
        :param object thermostat: Thermostat object
        :param integer istate: Initial adiabatic state
        :param double dt: Time interval
        :param integer nsteps: Nuclear step
        :param boolean l_print_dm: Logical to print BO population and coherence

    """
    def __init__(self, molecule, thermostat, istate, dt, nsteps, \
        l_print_dm, unit_dt, out_freq, verbosity):
        # Save name of MQC dynamics
        self.md_type = self.__class__.name__

        # Initialize Molecule object
        self.mol = molecule

        # Initialize Thermostat object
        self.thermo = thermostat

        # Initialize input values
        self.istate = istate
        self.nsteps = nsteps

        # Initialize time step
        self.istep = -1
        self.fstep = -1

        # Decide unit of time step
        self.unit_dt = unit_dt.lower()
        if (self.unit_dt == 'au'):
            self.dt = dt
        elif (self.unit_dt == 'fs'):
            self.dt = dt * fs_to_au
        else:
            error_message = "Invalid unit for time step!"
            error_vars = f"unit_dt = {unit_dt}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.l_print_dm = l_print_dm

        self.out_freq = out_freq
        self.verbosity = verbosity

    def run_init(self):
        """ Initialize MQC dynamics
        """
        pass

    def cl_update_position(self):
        """ Routine to update nuclear positions
        """
        pass

    def cl_update_velocity(self):
        """ Routine to update nuclear velocities
        """
        pass

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self):
        """ Routine to print the initial information of dynamics
        """
        pass

    def touch_file(self):
        """ Routine to write PyUNIxMD output files
        """
        pass

    def write_md_output(self):
        """ Write output files
        """
        pass

    def write_final_xyz(self):
        """ Write final positions and velocities
        """
        pass

    def check_qmmm(self):
        """ Routine to check compatibility between QM and MM objects
        """
        pass
        
