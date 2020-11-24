from __future__ import division
from qm.model.model import Model
import numpy as np

class ECR(Model):
    """ Class for extended coupling region with reflection (ECR) model BO calculation

        :param object molecule: molecule object
        :param double A: parameter for extended coupling region with reflection
        :param double B: parameter for extended coupling region with reflection
        :param double C: parameter for extended coupling region with reflection
    """
    def __init__(self, molecule, A=6E-4, B=0.1, C=0.9):
        # Initialize model common variables
        super(ECR, self).__init__(None)

        # Define parameters
        self.A = A
        self.B = B
        self.C = C

        # Set 'l_nacme' with respect to the computational method
        # ECR model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # ECR model can compute the gradient of several states simultaneously
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
        H[0, 0] = self.A
        H[1, 1] = - self.A
        if (x < 0):
            H[1, 0] = self.B * np.exp(self.C * x)
        else:
            H[1, 0] = self.B * (2. - np.exp(- self.C * x))
        H[0, 1] = H[1, 0]

        # Define a derivative of Hamiltonian
        dH[0, 0] = 0.
        dH[1, 1] = 0.
        dH[1, 0] = self.B * self.C * np.exp(- self.C * abs(x))
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

