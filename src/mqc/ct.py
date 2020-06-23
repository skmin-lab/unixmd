from __future__ import division
import numpy as np
from mqc.mqc import MQC

class CT(MQC):
    """ Class for Coupled-Trajectory Mixed Quantum-Classical dynamics (CTMQC)
    """
    def __init__(self, molecules, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True, threshold=1.0e-4, dist_cutoff=2.0/0.529166):

        # Initialize input values
        pass

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

