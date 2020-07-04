from __future__ import division
from misc import fs_to_au
import textwrap
import numpy as np

class MQC(object):
    """ Class for nuclear/electronic propagator used in MQC dynamics

        :param object molecule: molecule object
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string propagation: propagation scheme
        :param boolean l_adjnac: logical to adjust nonadiabatic coupling
    """
    def __init__(self, molecule, istate, dt, nsteps, nesteps, \
        propagation, l_adjnac):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize input values
        self.istate = istate
        self.dt = dt * fs_to_au
        self.nsteps = nsteps
        self.nesteps = nesteps

        self.propagation = propagation

        self.l_adjnac = l_adjnac

        self.rforce = np.zeros((molecule.nat, molecule.nsp))

        # Initialize coefficients and densities
        molecule.states[self.istate].coef = 1. + 0.j
        molecule.rho[self.istate, self.istate] = 1. + 0.j

    def cl_update_position(self, molecule):
        """ Routine to update nuclear positions

            :param object molecule: molecule object
        """
        self.calculate_force(molecule)

        molecule.vel += 0.5 * self.dt * self.rforce / np.column_stack([molecule.mass] * molecule.nsp)
        molecule.pos += self.dt * molecule.vel

    def cl_update_velocity(self, molecule):
        """ Routine to update nuclear velocities

            :param object molecule: molecule object
        """
        self.calculate_force(molecule)

        molecule.vel += 0.5 * self.dt * self.rforce / np.column_stack([molecule.mass] * molecule.nsp)
        molecule.update_kinetic()

#    def calculate_temperature(self, molecule):
#        """ Routine to calculate current temperature
#        """
#        pass
#        #self.temperature = molecule.ekin * 2 / float(molecule.dof) * au_to_K

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self, molecule, theory, field, thermostat, debug):
        """ Routine to print the initial information of dynamics

            :param object molecule: molecule object
            :param object theory: theory object containing on-the-fly calculation infomation
            :param object field: force field object containing MM calculation infomation
            :param object thermostat: thermostat type
            :param integer debug: verbosity level for standard output
        """
        # Print molecule information: coordinate, velocity
        molecule.print_init()

        # Print dynamics information
        dynamics_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Dynamics Information":>43s}
        {"-" * 68}
          QM Program               = {theory.qm_prog:>16s}
          QM Method                = {theory.qm_method:>16s}
        """)
        if (molecule.qmmm):
            dynamics_info += f"  MM Program               = {field.mm_prog:>16s}\n"
        dynamics_info += textwrap.indent(textwrap.dedent(f"""\

          MQC Method               = {self.md_type:>16s}
          Time Interval (fs)       = {self.dt / fs_to_au:16.6f}
          Initial State (0:GS)     = {self.istate:>16d}
          Nuclear Step             = {self.nsteps:>16d}
        """), "  ")
        if (self.md_type != "BOMD"):
            dynamics_info += f"  Electronic Step          = {self.nesteps:>16d}\n"
            dynamics_info += f"  Propagation Scheme       = {self.propagation:>16s}\n"
        print (dynamics_info, flush=True)

        # Print thermostat information
        thermostat.print_init()


