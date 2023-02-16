from __future__ import division
from qed.qed_calculator import QED_calculator
import numpy as np

class Jaynes_Cummings(QED_calculator):
    """ Class for Jaynes-Cummings Hamiltonian calculation

        :param object polariton: Polariton object
        :param double coupling_strength: Coupling strength for cavity-matter interaction
        :param boolean l_crt: Logical to include the counter-rotating terms
    """
    def __init__(self, polariton, coupling_strength=0.001, l_crt=False):
        # Save name of QED calculator
        super().__init__()

        # TODO : Should self.coup_str change to self.coupling_stregth for convention?
        self.coup_str = coupling_strength

        # For simplified JC model, it excludes counter-rotating terms
        # due to the use of rotating-wave approximation
        # If one want to use extended JC model, then the CRT must be included
        self.l_crt = l_crt

        # Generate an indicator describing the indices between polaritonic and uncoupled states
        self.get_d_ind = np.zeros((polariton.pst, 2), dtype=np.int32)
        self.cur_d_ind = np.zeros((polariton.pst, 2), dtype=np.int32)

        kst = 0
        for ist in range(polariton.nst):
            for jst in range(polariton.nphotons + 1):
                self.get_d_ind[kst] = [ist, jst]
                kst += 1

        # Initialize diabatic Hamiltonian matrix
        self.ham_d = np.zeros((polariton.pst, polariton.pst)) 
        self.ham_d_old = np.zeros((polariton.pst, polariton.pst)) 

        # Initialize unitray matrix obtained from diagonalization of Hamiltonian matrix
        self.unitary = np.zeros((polariton.pst, polariton.pst))
        self.unitary_old = np.zeros((polariton.pst, polariton.pst))
        self.unitary_dot = np.zeros((polariton.pst, polariton.pst))

        # Initialize permutation matrix to check current character of diabatic states
        self.permut = np.zeros((polariton.pst, polariton.pst), dtype=np.int32)


