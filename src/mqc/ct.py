from __future__ import division
import numpy as np
from mqc.mqc import MQC

class CT(MQC):
    """ Class for Coupled-Trajectory Mixed Quantum-Classical dynamics (CTMQC)
    """
    def __init__(self, molecules, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True, threshold=1.0e-4, dist_cutoff=2.0/0.529166):
        # Read total number of trajectories
        self.ntrajs = len(molecules)

        # Initialize input values: initial coef. and density of each trajectories
        for itraj in range(self.ntrajs):
            super().__init__(molecules[itraj], istate, dt, nsteps, nesteps, \
                propagation, l_adjnac)

        self.nsp = molecules[0].nsp
        self.nst = molecules[0].nst
        self.nat = molecules[0].nat

        # Initialize XF variables
        self.nst_pairs = self.nst * (self.nst - 1) 

        self.phase = np.zeros((self.ntrajs, self.nst, self.nat, self.nsp))
        self.qmom  = np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp))
        self.sigma = np.zeros((self.ntrajs, self.nat, self.nsp))
        self.sigma[:] = 0.08753727/self.ntrajs

        self.xfforce = np.zeros((self.ntrajs, self.nat, self.nsp))
        self.xfcdot  = np.zeros((self.ntrajs, self.nst), dtype=complex)

        # Trajectory-dependent temporaries to obtain qmom
        self.prod_g_ij = np.zeros((self.ntrajs, self.ntrajs))
        self.g_i = np.zeros((self.ntrajs))
        self.w_ij = np.zeros((self.ntrajs, self.ntrajs, self.nat, self.nsp))
        self.slope = np.zeros((self.ntrajs, self.nat, self.nsp))
        self.w0 = np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp))
        self.d0 = np.zeros((self.nst_pairs, self.nat, self.nsp))
        self.r0 = np.zeros((self.nst_pairs, self.nat, self.nsp))

        self.qmom_center = np.array(np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp)))
        self.qmom_center_old = np.array(np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp)))
        self.k_lk = np.array(np.zeros((self.ntrajs, self.nst, self.nst)))

    def run(self):
        """ Run MQC dynamics according to CTMQC
        """
        pass

    def calculate_qmom(self, molecules):
        """ Routine to calculate quantum momentum
        """
        pass

    def calculate_force(self, states):
        """ Routine to calculate force
        """
        pass

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in surface hopping
        """
        pass

    def el_propagator(self, molecule):
        """ Routine to propagate BO coefficients or density matrix
        """
        pass

    def print_init(self, molecule, theory, thermostat, debug):
        """ Routine to print the initial information of dynamics
        """
        pass

    def print_step(self, molecule, istep, debug):
        """ Routine to print each steps infomation about dynamics
        """
        pass

