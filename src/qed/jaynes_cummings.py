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

    def get_data(self, polariton, base_dir, calc_force_only):
        """ Construct the polaritonic and diabatic states from molecule, cavity, and interaction parts

            :param object polariton: Polariton object
            :param string base_dir: Base directory
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        super().get_data(base_dir, calc_force_only)
        # Step 1: Construct JC Hamiltonian matrix in uncoupled basis
        self.construct_Hamiltonian(polariton)
        self.move_dir(base_dir)

    def construct_Hamiltonian(self, polariton):
        """ Construct JC Hamiltonian matrix in uncoupled basis

            :param object polariton: Polariton object
        """
        self.ham_d = np.zeros((polariton.pst, polariton.pst)) 

        # Molecule part: Only diagonal elements exist
        tmp_ham = np.zeros((polariton.pst, polariton.pst)) 
        for ist in range(polariton.pst):
            ind_mol1 = self.get_d_ind[ist, 0]
            tmp_ham[ist, ist] = polariton.states[ind_mol1].energy

        self.ham_d += tmp_ham

        # Cavity part: Only diagonal elements exist
        tmp_ham = np.zeros((polariton.pst, polariton.pst)) 
        for ist in range(polariton.pst):
            ind_photon1 = self.get_d_ind[ist, 1]
            tmp_ham[ist, ist] = polariton.photon_freq * (ind_photon1 + 0.5)

        self.ham_d += tmp_ham

        # Interaction part: Only off-diagonal elements exist for JC Hamiltonian
        tmp_ham = np.zeros((polariton.pst, polariton.pst)) 
        for ist in range(polariton.pst):
            ind_mol1 = self.get_d_ind[ist, 0]
            ind_photon1 = self.get_d_ind[ist, 1]
            for jst in range(polariton.pst):
                ind_mol2 = self.get_d_ind[jst, 0]
                ind_photon2 = self.get_d_ind[jst, 1]
                off_term = 0.
                if (ind_mol2 == ind_mol1 + 1 and ind_photon2 == ind_photon1 - 1):
                    off_term = 1.
                elif (ind_mol2 == ind_mol1 - 1 and ind_photon2 == ind_photon1 + 1):
                    off_term = 1.
                # For extended JC Hamiltonian, the counter-rotating terms are included
                if (self.l_crt):
                    if (ind_mol2 == ind_mol1 + 1 and ind_photon2 == ind_photon1 + 1):
                        off_term = 1.
                    elif (ind_mol2 == ind_mol1 - 1 and ind_photon2 == ind_photon1 - 1):
                        off_term = 1.
                polarization = sum(polariton.field_pol_vec * polariton.tdp[ind_mol1, ind_mol2])
                tmp_ham[ist, jst] = polarization * self.coup_str * off_term

        self.ham_d += tmp_ham

        # Write 'ham_d.dat' file; Hamiltonian matrix elements in uncoupled basis
        ham_d_row = ""
        for ist in range(polariton.pst):
            ham_d_row += " ".join([f"{i:12.6f}" for i in self.ham_d[ist]]) + "\n"

        file_name = "ham_d.dat"
        with open(file_name, "w") as f:
            f.write(ham_d_row)


