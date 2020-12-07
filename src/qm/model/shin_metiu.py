from __future__ import division
from qm.model.model import Model
import numpy as np
from math import erf
from misc import eps

class Shin_Metiu(Model):
    """ Class for 1D Shin-Metiu model BO calculation in a real-space grid

        :param object molecule: molecule object
        :param integer nx: the number of grid points
        :param double xmin: lower bound of the 1D space
        :param double xmax: upper bound of the 1D space
        :param double L: the distance between two fixed nuclei
        :param double Rc: the parameter of a moving nucleus
        :param double Rl: the parameter of a fixed nucleus in the left side
        :param double Rr: the parameter of a fixed nucleus in the right side
    """
    def __init__(self, molecule, nx=401, xmin=-20.0, xmax=20.0, L=19.0, Rc=5.0, Rl=4.0, Rr=3.1):
        # Initialize model common variables
        super(Shin_Metiu, self).__init__(None)

        # Set the grid
        self.nx = nx
        self.xmin = xmin
        self.xmax = xmax

        # Parameters in au
        self.L = L + eps
        self.Rc = Rc
        self.Rl = Rl
        self.Rr = Rr

        self.dx = (self.xmax - self.xmin) / np.float(self.nx - 1)
        self.H = np.zeros((self.nx, self.nx))

        # Set 'l_nacme' with respect to the computational method
        # Shin-Metiu model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # Shin-Metiu model can compute the gradient of several states simultaneously
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from Shin-Metiu BO calculation

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Initialize Hamiltonian
        self.H = 0.

        # Add the kinetic-energy contribution (tridiagonal)
        self.H += - 0.5 * (np.diag([1.] * (self.nx - 1), - 1) + np.diag([- 2.] * self.nx, 0) + \
            np.diag([1.] * (self.nx - 1), 1)) / self.dx ** 2
 
        x = molecule.pos[0, 0]

        # Add the potential contribution (diagonal)
        xes = [self.xmin + ix * self.dx for ix in range(self.nx)]
        Vs = [self.get_V(x, xe) for xe in xes]
        self.H += np.diag(Vs)

        # Diagonalization
        ws, unitary = np.linalg.eig(self.H)

        # Sorting eigenvalues in the ascending order and the corresponding eigenvectors
        idx = np.argsort(ws)
        ws = ws[idx]
        unitary = unitary[:, idx]

        # Slicing eigenvalues and eigenvectors up to the given number of states
        ws = ws[0:molecule.nst]
        unitary = unitary[:, 0:molecule.nst]

        for ist in range(molecule.nst):
            molecule.states[ist].energy = ws[ist]

        # Extract adiabatic quantities
        dVs = [self.get_dV(x, xe) for xe in xes]
        dVijs = np.dot(np.transpose(unitary), np.dot(np.diag(dVs), unitary))

        Fs = - np.diag(dVijs)
        for ist in range(molecule.nst):
            molecule.states[ist].force = Fs[ist]

        for ist in range(molecule.nst):
            for jst in range(ist + 1, molecule.nst):
                molecule.nac[ist, jst, 0, 0] = dVijs[ist, jst] / (ws[jst] - ws[ist])
                molecule.nac[jst, ist, 0, 0] = - molecule.nac[ist, jst, 0, 0]

    def get_V(self, x, xe):
        """ Calculate potential elements of the BO Hamiltonian

            :param double x: the nuclear position
            :param double xe: the electronic position
        """
        RR = np.abs(x - xe)

        if (RR > eps):
            V = - erf(RR / self.Rc) / RR
        else:
            V = - 2. / (np.sqrt(np.pi) * self.Rc)

        V += - erf(np.abs(xe - 0.5 * self.L) / self.Rr) / np.abs(xe - 0.5 * self.L) - \
            erf(np.abs(xe + 0.5 * self.L) / self.Rl) / np.abs(xe + 0.5 * self.L) + \
            1. / np.abs(x - 0.5 * self.L) + 1. / np.abs(x + 0.5 * self.L)

        return V

    def get_dV(self, x, xe):
        """ Calculate del potential elements of the BO Hamiltonian

            :param double x: the nuclear position
            :param double xe: the electronic position
        """
        RR = np.abs(x - xe)
 
        if (RR > eps):
            dV = (x - xe) * erf(RR / self.Rc) / RR ** 3 - \
                2. * (x - xe) * np.exp(- RR ** 2 / self.Rc ** 2) / np.sqrt(np.pi) / self.Rc / RR ** 2
        else:
            dV = 0.

        dV -= (np.abs(x - 0.5 * self.L) ** (- 3)) * (x - 0.5 * self.L) + \
            (np.abs(x + 0.5 * self.L) ** (- 3)) * (x + 0.5 * self.L)

        return dV

