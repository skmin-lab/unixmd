from __future__ import division
from qm.model.model import Model
import numpy as np

class SAC(Model):
    """ Class for simple avoided crossing (SAC) model BO calculation

        :param object molecule: molecule object
        :param double A: parameter for simple avoided crossing model
        :param double B: parameter for simple avoided crossing model
        :param double C: parameter for simple avoided crossing model
        :param double D: parameter for simple avoided crossing model
    """
    def __init__(self, molecule, A=0.01, B=1.6, C=0.005, D=1.):
        # Initialize model common variables
        super(SAC, self).__init__(None)

        # Define parameters
        self.A = A
        self.B = B
        self.C = C
        self.D = D

        # Set 'l_nacme' with respect to the computational method
        # SAC model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # SAC model can compute the gradient of several states simultaneously
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from simple avoided crossing model BO calculation

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Initialize diabatic Hamiltonian
        H = np.zeros((2, 2))
        dH = np.zeros((2, 2))
        unitary = np.zeros((2, 2))

        x = molecule.pos[0]

        # Define Hamiltonian
        H[0, 0] = np.sign(x) * self.A * (1. - np.exp(- self.B * abs(x)))
        H[1, 1] = - H[0, 0]
        H[1, 0] = self.C * np.exp(- self.D * x ** 2)
        H[0, 1] = H[1, 0]

        # Define a derivative of Hamiltonian
        dH[0, 0] = self.A * self.B * np.exp(- self.B * abs(x))
        dH[1, 1] = - dH[0, 0]
        dH[1, 0] = - 2. * self.D * self.C * x * np.exp(- self.D * x ** 2)
        dH[0, 1] = dH[1, 0]

        # Diagonalization
        a = 4. * H[1, 0] * H[0, 1] + (H[1, 1] - H[0, 0]) ** 2
        sqa = np.sqrt(a)
        tantheta = (H[1, 1] - H[0, 0] - sqa) / H[1, 0]  * 0.5
        theta = np.arctan(tantheta)

        unitary[0, 0] = np.cos(theta)
        unitary[1, 0] = np.sin(theta)
        unitary[0, 1] = - np.sin(theta)
        unitary[1, 1] = np.cos(theta)

        # Extract adiabatic quantities
        molecule.states[0].energy = 0.5 * (H[0, 0] + H[1, 1]) - 0.5 * sqa
        molecule.states[1].energy = 0.5 * (H[0, 0] + H[1, 1]) + 0.5 * sqa

        molecule.states[0].force = np.dot(unitary[:, 1], np.matmul(dH, unitary[:, 1]))
        molecule.states[1].force = np.dot(unitary[:, 0], np.matmul(dH, unitary[:, 0]))

        molecule.nac[0, 1, 0, 0] = np.dot(unitary[:, 0], np.matmul(dH, unitary[:, 1])) / sqa
        molecule.nac[1, 0, 0, 0] = - molecule.nac[0, 1, 0, 0]

