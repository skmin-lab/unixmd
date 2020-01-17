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

    #def print_init(self, molecule, theory, thermostat, debug):
    def print_init(self, molecule, theory, thermostat):
        """ Routine to print the initial information of dynamics
        """
        molecule.print_init()
         
        qm_prog = str(theory.__class__).split('.')[1]
        qm_method = theory.__class__.__name__
        mdtype = self.__class__.__name__
        if(mdtype == "SH"):
            method = "FSSH"
        elif(mdtype == "Eh"):
            method = "Ehrenfest dynamics"
        elif(mdtype == "SHXF"):
            method = "SHXF"
        elif(mdtype == "BOMD"):
            mehtod = "BOMD"
        
        dynamics_info = textwrap.dedent(f"""\
        {'-'*118}
        {'Dynamics information':>69}
        {'-'*118}
        Electronic structure calculation is performed by \"{qm_method}\", which is implemented in \"{qm_prog}\"
        
        MQC method  = {method}
        
        Nuclei propation    = {self.nsteps}   Electron propagtion = {self.nesteps} 
        Time interval       = {self.dt} fs
        Initial state       = {self.istate} (0:GS)
        propagation         = {self.propagation}
        """)
        print(dynamics_info, flush=True) 
        
        print(f"{'-'*118}",flush=True)
        print(f"{'Thermostat information':>69}",flush=True)
        print(f"{'-'*118}",flush=True)
        thermostat.print_init()
        
        print(f"{'-'*118}",flush=True)
        print(f"{'Start dynamics':>69}",flush=True)
        print(f"{'-'*118}",flush=True)
        if(method == "FSSH" or method == "SHXF"):
            INIT = f"#INFO{'STEP':^10}{'State':^10}{'Max_prob':>12}{'rand':>15}{'Kinetic(H)':>20}{'Potential(H)':>17}{'Total(H)':>14}{'Temperature(K)':>22}{'norm':>15}"
            debug2 = f"#DEBUG2{'STEP':^8}{'Acc. Hopping Prob.':^10}"
        elif(method == "Ehrenfest dynamics"):
            INIT = f"#INFO{'STEP':^10}{'Kinetic(H)':>12}{'Potential(H)':>17}{'Total(H)':>14}{'Temperature(K)':>22}{'norm':>15}"
            debug2 = ""
        elif(method == "SHXF"):
            INIT = f"#INFO{'STEP':^10}{'State':^10}{'Max_prob':>12}{'rand':>15}{'Kinetic(H)':>20}{'Potential(H)':>17}{'Total(H)':>14}{'Temperature(K)':>22}{'norm':>15}"
        elif(method == "BOMD"):
            INIT = f"#INFO{'STEP':^10}{'State':^10}{'Kinetic(H)':>20}{'Potential(H)':>17}{'Total(H)':>14}{'Temperature(K)':>22}{'norm':>15}"
        print(INIT, flush=True)
    
        debug1 = f"#DEBUG1{'STEP':^8}"
        for ist in range(molecule.nst):
            debug1 += f"{'Potential(H)':>17}{ist}" 
        print(debug1, flush=True) 
        print(debug2, flush=True) 
        
