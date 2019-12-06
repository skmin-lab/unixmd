from __future__ import division
import numpy as np
from misc import fs_to_au
import textwrap

class MQC(object):
    """ Class for nuclear/electronic propagator used in MQC dynamics
    """
    def __init__(self, molecule, istate, dt, nsteps, nesteps, \
        propagation, l_adjnac):
        # Initialize input values
        self.istate = istate
        self.dt = dt * fs_to_au
        self.nsteps = nsteps
        self.nesteps = nesteps

        self.propagation = propagation

        self.l_adjnac = l_adjnac

        self.rforce = np.zeros((molecule.nat, molecule.nsp))

        # initialize coefficients and densities
        molecule.states[self.istate].coef = 1. + 0.j
        molecule.rho[self.istate, self.istate] = 1. + 0.j

    def cl_update_position(self, molecule):
        """ Routine to update nulcear postions
        """
        self.calculate_force(molecule)

        molecule.vel += 0.5 * self.dt * self.rforce / np.column_stack([molecule.mass] * molecule.nsp)
        molecule.pos += self.dt * molecule.vel

    def cl_update_velocity(self, molecule):
        """ Rotine to update nulcear velocities
        """
        self.calculate_force(molecule)

        molecule.vel += 0.5 * self.dt * self.rforce / np.column_stack([molecule.mass] * molecule.nsp)
        molecule.update_kinetic()

    def calculate_temperature(self, molecule):
        """ Routine to calculate current temperature
        """
        pass
    #    self.temperature = molecule.ekin * 2 / float(molecule.dof) * au_to_K
    
    def calculate_force(self):
        """ Routine to calculate the forces
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self, molecule, theory, thermostat):
        """ Routine to print the initial information of dynamics
        """
        molecule.print_init()

        qm_prog = str(theory.__class__).split('.')[1]
        qm_method = theory.__class__.__name__
        print(f"QM_prog       = {qm_prog}", flush=True)
        print(f"QM_method     = {qm_method}", flush=True)
        
        mdtype = self.__class__.__name__
        print(f"MD_type     = {mdtype}", flush=True)
        
        # print initial infomation of dynamics
        dynamics_info = textwrap.dedent(f"""\
        nsteps      = {self.nsteps} 
        nesteps     = {self.nesteps}
        dt          = {self.dt}
        propagation = {self.propagation}
        istate      = {self.istate}
        """)
        print(dynamics_info, flush=True) 
        
        thermostat.print_init()
        INIT = f"#INFO{'STEP':^10}{'State':^10}{'Max_prob':>12}{'rand':>15}{'Kinetic(H)':>20}{'Potential(H)':>17}{'Total(H)':>14}{'Temperature(K)':>22}{'norm':>15}"
        state = 0
        print(INIT, flush=True)
