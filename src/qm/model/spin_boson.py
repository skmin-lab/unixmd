from __future__ import division
from qm.model.model import Model
import os, shutil, re
import numpy as np

class Spin_Boson(Model):
    """ Class for Spin-Boson model BO calculation
    """
    def __init__(self, molecule, qm_path="./", E=0.5, V=0.5, lamb=1.0, W=0.1, T=1.0):
        # Initialize model common variables
        super(Spin_Boson, self).__init__(qm_path)

        # reference unit is 300K.
        REF_TO_au = 0.00094996
        self.E = E * REF_TO_au
        self.V = V * REF_TO_au
        self.lamb = lamb * REF_TO_au
        self.W = W * REF_TO_au
        self.T = T * REF_TO_au

        molecule.l_nacme = False
        self.re_calc = False
    
    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        # H is diabatic hamiltonian
        H = np.zeros((2,2))
        F = np.zeros((2, molecule.nat))

        H_q = np.array([[self.E, self.V],[self.V, -self.E]])
        H += H_q
        
        # w: frequency, pot: harmonic potential, coup: coupling strength betwee spin and bosons
        w = 0.0; pot = 0.0; coup = 0.0

        # Here, mass is 1/w not one.
        for iat in range(molecule.nat):
            w = 1.0 / molecule.mass[iat]
            pot = 0.5 * w * molecule.pos[iat]**2
            coup = np.sqrt(self.lamb * w/ (2.0 * molecule.nat)) * molecule.pos[iat]
            
            H[0][0] += pot + coup
            H[1][1] += pot - coup

            dpot = w * molecule.pos[iat]
            dcoup = coup/molecule.pos[iat]
            F[0, iat] = -dpot - dcoup
            F[1, iat] = -dpot + dcoup
        
        # adiabatic varibles
        energy = np.zeros(2)
        force = np.zeros((2, molecule.nat))

        # nac[0] -> <gs|  |es>
        # nac[1] -> <es|  |gs>

        # Diagonalization
        # In adiabatic representation, 0 and 1 is ground state and excited state.
        energy[0] = ((H[0][0] + H[1][1]) - np.sqrt((H[0][0]-H[1][1])**2 + 4.0 * H[0][1]**2))/2.0
        energy[1] = ((H[0][0] + H[1][1]) + np.sqrt((H[0][0]-H[1][1])**2 + 4.0 * H[0][1]**2))/2.0

        for ist in range(molecule.nst):
            molecule.states[ist].energy = energy[ist]
        
        theta = np.arctan((energy[1]-H[0][0])/self.V)

        unitary = np.zeros((2,2))
        unitary[0][0] =  np.cos(theta)
        unitary[0][1] =  np.sin(theta)
        unitary[1][0] = -np.sin(theta)
        unitary[1][1] =  np.cos(theta)

        cos = np.cos(theta)
        sin = np.sin(theta)

        print(f"UMAT  {istep+1}  cos = {unitary[0][0]:15.8f}  sin = {unitary[0][1]:15.8f}", flush=True)
        for iat in range(molecule.nat):
            molecule.states[0].force[iat, 0] = sin**2 * F[0, iat] + cos**2 * F[1, iat]
            molecule.states[1].force[iat, 0] = cos**2 * F[0, iat] + sin**2 * F[1, iat]
            molecule.nac[0, 1, iat, 0] = (-F[0, iat] * sin * cos + F[1,iat] * sin * cos) / (energy[0] - energy[1])
            molecule.nac[1, 0, iat, 0] = -molecule.nac[0, 1, iat, 0]
        nacme = 0.0

#        for iat in range(molecule.nat):
#
#          nacme += molecule.nac[0,1,iat,0] * molecule.vel[iat]
#        print(nacme)#, molecule.vel)

