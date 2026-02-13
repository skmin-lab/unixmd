from __future__ import division
from qm.model.model import Model
import numpy as np
from scipy.special import erf
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

        self.dx = (self.xmax - self.xmin) / np.float64(self.nx - 1)
        self.H = np.zeros((self.nx, self.nx))

        # Set 'l_nacme' with respect to the computational method
        # Shin-Metiu model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # Shin-Metiu model can compute the gradient of several states simultaneously
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only, traj=None):
        """ Extract energy, gradient and nonadiabatic couplings from Shin-Metiu BO calculation

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
            :param object traj: Trajectory object containing the calculator and trajectory
        """
        # Initialize Hamiltonian
        self.H = 0.

        # Add the kinetic-energy contribution (tridiagonal)
        self.H += - 0.5 * (np.diag([1.] * (self.nx - 1), - 1) + np.diag([- 2.] * self.nx, 0) + \
            np.diag([1.] * (self.nx - 1), 1)) / self.dx ** 2
 
        x = molecule.pos[0, 0]

        # Add the potential contribution (diagonal)
        xes = np.linspace(self.xmin, self.xmax, self.nx)
        Vs = self.get_V(x, xes)
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
        dVs = self.get_dV(x, xes)
        dVijs = np.dot(np.transpose(unitary), np.dot(np.diag(dVs), unitary))

        Fs = - np.diag(dVijs)
        for ist in range(molecule.nst):
            molecule.states[ist].force = Fs[ist]

        # Vectorized NAC calculation: dVijs / (ws[j] - ws[i])
        energy_diff = ws[np.newaxis, :] - ws[:, np.newaxis]  # energy_diff[i,j] = ws[j] - ws[i]
        # Avoid division by zero on diagonal by setting it to 1 (result will be zeroed anyway)
        energy_diff_safe = np.where(energy_diff != 0, energy_diff, 1.)
        nac_full = dVijs / energy_diff_safe
        # Enforce antisymmetry: keep upper triangle and subtract transpose
        nac_upper = np.triu(nac_full, k=1)
        molecule.nac[:, :, 0, 0] = nac_upper - nac_upper.T

    def get_V(self, x, xes):
        """ Calculate potential elements of the BO Hamiltonian (vectorized)

            :param double x: the nuclear position
            :param double/ndarray xes: the electronic position(s)
        """
        RR = np.abs(x - xes)

        # Use limit value: lim_{r->0} erf(r/R)/r = 2/(sqrt(pi)*R)
        RR_safe = np.where(RR > eps, RR, 1.)
        V = np.where(RR > eps,
                     -erf(RR / self.Rc) / RR_safe,
                     -2. / (np.sqrt(np.pi) * self.Rc))

        # Fixed nuclei terms with proper limits
        xes_r = np.abs(xes - 0.5 * self.L)
        xes_l = np.abs(xes + 0.5 * self.L)
        xes_r_safe = np.where(xes_r > eps, xes_r, 1.)
        xes_l_safe = np.where(xes_l > eps, xes_l, 1.)

        # lim_{r->0} erf(r/R)/r = 2/(sqrt(pi)*R)
        V_r = np.where(xes_r > eps,
                       -erf(xes_r / self.Rr) / xes_r_safe,
                       -2. / (np.sqrt(np.pi) * self.Rr))
        V_l = np.where(xes_l > eps,
                       -erf(xes_l / self.Rl) / xes_l_safe,
                       -2. / (np.sqrt(np.pi) * self.Rl))

        # Nuclear repulsion terms (1/r diverges, use safe denominator)
        x_r = np.abs(x - 0.5 * self.L)
        x_l = np.abs(x + 0.5 * self.L)
        x_r_safe = np.where(x_r > eps, x_r, eps)
        x_l_safe = np.where(x_l > eps, x_l, eps)

        V += V_r + V_l + 1. / x_r_safe + 1. / x_l_safe

        return V

    def get_dV(self, x, xes):
        """ Calculate del potential elements of the BO Hamiltonian (vectorized)

            :param double x: the nuclear position
            :param double/ndarray xes: the electronic position(s)
        """
        RR = np.abs(x - xes)

        # Use safe denominator to avoid division by zero warnings
        RR_safe = np.where(RR > eps, RR, 1.)
        dV = np.where(RR > eps,
                      (x - xes) * erf(RR / self.Rc) / RR_safe ** 3
                      - 2. * (x - xes) * np.exp(-RR ** 2 / self.Rc ** 2) / np.sqrt(np.pi) / self.Rc / RR_safe ** 2,
                      0.)

        # Safe denominators for fixed nuclei terms
        x_r = np.abs(x - 0.5 * self.L)
        x_l = np.abs(x + 0.5 * self.L)
        x_r_safe = np.where(x_r > eps, x_r, eps)
        x_l_safe = np.where(x_l > eps, x_l, eps)

        dV -= ((x - 0.5 * self.L) / x_r_safe ** 3
               + (x + 0.5 * self.L) / x_l_safe ** 3)

        return dV


