from __future__ import division
from qm.model.model import Model
import os, shutil, re
import numpy as np

class DAC(Model):
    """ Class for double avoided crossing (DAC) model BO calculation
    """
    def __init__(self, molecule, E0 = 0.05, A=0.1, B=0.28, C=0.015, D=0.06, qm_path="./"):
        # Initialize model common variables
        super(DAC, self).__init__(qm_path)

        # define parameters
        self.E0 = E0
        self.A = A
        self.B = B
        self.C = C
        self.D = D

        molecule.l_nacme = False
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """
        """
        # diabatic Hamiltonian
        H = np.zeros((2,2))
        dH = np.zeros((2,2))
        unitary = np.zeros((2,2))

        x = molecule.pos[0]

        # define Hamiltonian
        H[0][0] = 0.
        H[1][1] = self.E0 - self.A * np.exp(-self.B * x ** 2)
        H[1][0] = self.C * np.exp(-self.D * x ** 2)
        H[0][1] = H[1][0]

        # define a derivative of Hamiltonian
        dH[0][0] = 0.
        dH[1][1] = self.A * self.B * 2. * x * np.exp(-self.B * x ** 2)
        dH[1][0] = -2. * self.D * self.C * x * np.exp(-self.D * x ** 2)
        dH[0][1] = dH[1][0]

        # diagonalization
        a = 4. * H[1][0] * H[0][1] + (H[1][1] - H[0][0]) ** 2
        sqa = np.sqrt(a)
        tantheta = (H[1][1] - H[0][0] - sqa) / H[1][0]  * 0.5
        theta = np.arctan(tantheta)

        unitary[0][0] = np.cos(theta)
        unitary[1][0] = np.sin(theta)
        unitary[0][1] = -np.sin(theta)
        unitary[1][1] = np.cos(theta)

        molecule.states[0].energy = 0.5 * (H[0][0] + H[1][1]) - 0.5 * sqa
        molecule.states[1].energy = 0.5 * (H[0][0] + H[1][1]) + 0.5 * sqa

        molecule.states[0].force = np.dot(unitary[:, 0], np.matmul(dH, unitary[:, 0])) 
        molecule.states[1].force = np.dot(unitary[:, 1], np.matmul(dH, unitary[:, 1])) 

        molecule.nac[0, 1, 0, 0] = np.dot(unitary[:, 0], np.matmul(dH, unitary[:, 1])) / sqa
        molecule.nac[1, 0, 0, 0] = -molecule.nac[0, 1, 0, 0]
