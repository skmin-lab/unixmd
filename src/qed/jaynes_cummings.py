from __future__ import division
from qed.qed_calculator import QED_calculator
import os, shutil
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

    def get_data(self, polariton, base_dir, pol_list, dt, istep, calc_force_only):
        """ Construct the polaritonic and diabatic states from molecule, cavity, and interaction parts

            :param object polariton: Polariton object
            :param string base_dir: Base directory
            :param integer,list pol_list: List of polaritonic states for QED calculation
            :param double dt: Time interval
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        super().get_data(base_dir, calc_force_only)
        # Step 1: Construct JC Hamiltonian matrix in uncoupled basis
        self.construct_Hamiltonian(polariton)
        # Step 2: Diagonalize Hamiltonian matrix and obtain eneriges
        if (not calc_force_only):
            self.solve_polaritonic_states(polariton, dt, istep)
            self.save_output_files(base_dir, pol_list, istep)
        # Step 3: Calculate properties of polaritonic states; forces
        self.calculate_properties(polariton, pol_list)
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

    def solve_polaritonic_states(self, polariton, dt, istep):
        """ Diagonalize JC Hamiltonian matrix to get the unitary matrix
            The energy for polaritonic states are directly obtained from the eigenvalues

            :param object polariton: Polariton object
            :param double dt: Time interval
            :param integer istep: Current MD step
        """
        # Backup the permutation matrix and index for current diabatic state
        if (istep >= 0):
            permut_old = np.zeros((polariton.pst, polariton.pst), dtype=np.int32) 
            permut_old = np.copy(self.permut)
            cur_d_ind_old = np.zeros((polariton.pst, 2), dtype=np.int32)
            cur_d_ind_old = np.copy(self.cur_d_ind)

        # eigh is a method to diagonalize a complex Hermitian or a real-symmetric matrix
        w, v = np.linalg.eigh(self.ham_d)

        # The eigenvector becomes unitary matrix between diabatic and adiabatic states
        # For self.unitary, 1st and 2nd ranks represent diabatic and adiabatic indices, respectively
        self.unitary = np.copy(v)

        # Temporary matrix to be used to check the phase of unitary matrix
        if (istep >= 0):
            tmp_mat = np.zeros((polariton.pst, polariton.pst)) 
            tmp_mat = np.matmul(np.transpose(self.unitary_old), self.unitary)
            for ist in range(polariton.pst):
                for jst in range(polariton.pst):
                    tmp_val = round(tmp_mat[ist, jst] ** 2.)
                    # Maintain the original sign for large elements
                    if (tmp_mat[ist, jst] > 0.):
                        tmp_mat[ist, jst] = tmp_val
                    else:
                        tmp_mat[ist, jst] = - tmp_val

            tmp_sign = np.zeros(polariton.pst)
            tmp_sign = np.sum(tmp_mat, axis=0)
            # If the phase is reversed against the old step, change the sign
            for ist in range(polariton.pst):
                for jst in range(polariton.pst):
                    self.unitary[ist, jst] = self.unitary[ist, jst] * tmp_sign[jst]

        # Write 'unitary.dat' file; unitary matrix elements
        unitary_row = ""
        for ist in range(polariton.pst):
            unitary_row += " ".join([f"{i:12.6f}" for i in self.unitary[ist]]) + "\n"

        file_name = "unitary.dat"
        with open(file_name, "w") as f:
            f.write(unitary_row)

        # Calculate numerical derivative of the unitary matrix
        if (istep >= 0):
            self.unitary_dot = (self.unitary - self.unitary_old) / dt

            # Write 'unitary_dot.dat' file; unitary matrix elements
            unitary_dot_row = ""
            for ist in range(polariton.pst):
                unitary_dot_row += " ".join([f"{i:12.6f}" for i in self.unitary_dot[ist]]) + "\n"

            file_name = "unitary_dot.dat"
            with open(file_name, "w") as f:
                f.write(unitary_dot_row)

        # Calculate permutation matrix to be used to check diabatic characters
        for ist in range(polariton.pst):
            for jst in range(polariton.pst):
                self.permut[ist, jst] = round(self.unitary[ist, jst] ** 2.)

        # Calculate the character of diabatic states at current step
        # self.cur_d_ind is used to check the switching of adjacent two states
        self.cur_d_ind = np.transpose(np.matmul(np.transpose(self.get_d_ind), self.permut))

        # Write 'index_AD.dat' file; diabatic character for each polaritonic state
        cur_d_ind_row = ""
        for ist in range(polariton.pst):
            cur_d_ind_row += " ".join([f"{i:6d}" for i in self.cur_d_ind[ist]]) + "\n"

        file_name = "index_AD.dat"
        with open(file_name, "w") as f:
            f.write(cur_d_ind_row)

        # The eigenvalues correspond to the energies of polaritonic states
        for ist in range(polariton.pst):
            polariton.pol_states[ist].energy = w[ist]

    def calculate_properties(self, polariton, pol_list):
        """ Calculate properties of polaritonic states
            The forces for polaritonic states are calculated from the gradients
            of electronic states and transition dipole gradients

            :param object polariton: Polariton object
            :param integer,list pol_list: List of polaritonic states for QED calculation
        """
        for rst in pol_list:
            # The applied gradient on the target polaritonic state
            qed_grad = np.zeros((polariton.nat, polariton.ndim))

            # First term: Ehrenfest-like term
            # Loop for diabatic states
            for ist in range(polariton.pst):
                ind_mol1 = self.get_d_ind[ist, 0]
                ind_photon1 = self.get_d_ind[ist, 1]
                qed_grad -= polariton.states[ind_mol1].force * self.unitary[ist, rst] ** 2.

            # Second term: TDP gradients
            # Loop for diabatic states
            for ist in range(polariton.pst):
                ind_mol1 = self.get_d_ind[ist, 0]
                ind_photon1 = self.get_d_ind[ist, 1]
                # Loop for diabatic states
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

                    if (ist != jst):
                        field_dot_grad = np.zeros((polariton.nat, polariton.ndim))
                        for idim in range(polariton.ndim):
                            field_dot_grad[0:polariton.nat_qm] += polariton.field_pol_vec[idim] \
                                * polariton.tdp_grad[ind_mol1, ind_mol2, idim]
                        qed_grad += field_dot_grad * self.coup_str * self.unitary[ist, rst] \
                            * self.unitary[jst, rst] * off_term

            polariton.pol_states[rst].force = - np.copy(qed_grad)

    def save_output_files(self, base_dir, pol_list, istep):
        """ Save output files generated by QED calculation

            :param string base_dir: Base directory
            :param integer,list pol_list: List of polaritonic states for QED calculation
            :param integer istep: Current MD step
        """
        # Copy the output file to 'qed_log' directory
        tmp_dir = os.path.join(base_dir, "qed_log")
        if (os.path.exists(tmp_dir)):
            ham_d_step = f"ham_d.dat.{istep + 1}.{pol_list[0]}"
            shutil.copy("ham_d.dat", os.path.join(tmp_dir, ham_d_step))
            unitary_step = f"unitary.dat.{istep + 1}.{pol_list[0]}"
            shutil.copy("unitary.dat", os.path.join(tmp_dir, unitary_step))
            if (istep >= 0):
                unitary_dot_step = f"unitary_dot.dat.{istep + 1}.{pol_list[0]}"
                shutil.copy("unitary_dot.dat", os.path.join(tmp_dir, unitary_dot_step))
            index_step = f"index_AD.dat.{istep + 1}.{pol_list[0]}"
            shutil.copy("index_AD.dat", os.path.join(tmp_dir, index_step))

    def backup_qed(self):
        """ Backup Hamiltonian matrix and unitary matrix
        """
        self.ham_d_old = np.copy(self.ham_d)
        self.unitary_old = np.copy(self.unitary)

    def transform(self, polariton, mode):
        """ Transform the coefficients using unitary operation

            :param object polariton: Polariton object
            :param string mode: Transformation mode for coefficients (a2d = to diabatic, d2a = to adiabatic)
        """
        # For self.unitary, 1st and 2nd ranks represent diabatic and adiabatic indices, respectively
        if (mode == "a2d"):

            # D = U * C
            # Index for diabatic states
            for ist in range(polariton.pst):
                tmp_real = 0.
                tmp_imag = 0.
                # Index for adiabatic states
                for jst in range(polariton.pst):
                    tmp_real += self.unitary[ist, jst] * polariton.pol_states[jst].coef_a.real
                    tmp_imag += self.unitary[ist, jst] * polariton.pol_states[jst].coef_a.imag
                polariton.pol_states[ist].coef_d = complex(tmp_real, tmp_imag)

            for ist in range(polariton.pst):
                for jst in range(polariton.pst):
                    polariton.rho_d[ist, jst] = polariton.pol_states[ist].coef_d.conjugate() \
                        * polariton.pol_states[jst].coef_d

        elif (mode == "d2a"):

            # C = U^T * D; U^-1 = U^T
            # Index for adiabatic states
            for ist in range(polariton.pst):
                tmp_real = 0.
                tmp_imag = 0.
                # Index for diabatic states
                for jst in range(polariton.pst):
                    tmp_real += self.unitary[jst, ist] * polariton.pol_states[jst].coef_d.real
                    tmp_imag += self.unitary[jst, ist] * polariton.pol_states[jst].coef_d.imag
                polariton.pol_states[ist].coef_a = complex(tmp_real, tmp_imag)

            for ist in range(polariton.pst):
                for jst in range(polariton.pst):
                    polariton.rho_a[ist, jst] = polariton.pol_states[ist].coef_a.conjugate() \
                        * polariton.pol_states[jst].coef_a


