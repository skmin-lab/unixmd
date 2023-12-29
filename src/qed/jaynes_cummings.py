from __future__ import division
from qed.qed_calculator import QED_calculator
from misc import call_name, eps
import textwrap
import os, shutil
import numpy as np

class Jaynes_Cummings(QED_calculator):
    """ Class for Jaynes-Cummings Hamiltonian calculation

        :param object polariton: Polariton object
        :param double coupling_strength: Coupling strength for cavity-matter interaction
        :param string force_level: Level for calculation of force of target polaritonic state
        :param boolean l_check_crossing: Logical to check diabatic character of polaritonic states
        :param boolean l_crt: Logical to include the counter-rotating terms
    """
    def __init__(self, polariton, coupling_strength=0.001, force_level="nac", \
        l_check_crossing=False, l_crt=False):
        # Save name of QED calculator
        super().__init__()

        # TODO : Should self.coup_str change to self.coupling_stregth for convention?
        self.coup_str = coupling_strength

        self.force_level = force_level.lower()
        if not (self.force_level in ["energy", "nac", "tdp"]):
            error_message = "Invalid method for calculation of polaritonic state gradient!"
            error_vars = f"force_level = {self.force_level}"
            raise ValueError (f"( {self.qed_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (polariton.l_nacme and self.force_level == "nac"):
            error_message = "Chosen force computation needs evaluation of NACVs, check your QM object or QED force level!"
            error_vars = f"(QM) qm_prog.qm_method = {qm.qm_prog}.{qm.qm_method}, force_level = {self.force_level}"
            raise ValueError (f"( {self.qed_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (self.force_level in ["energy", "nac"]):
            print ("\n\n WARNING: Chosen force computation may not be accurate for strong coupling case, carefully choose force computation! \n\n", flush=True)

        # Check the diabatic character of adjacent two states
        self.l_check_crossing = l_check_crossing

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

        self.l_trivial = False

        # Set 'l_pnacme' with respect to the Hamiltonian model
        # pNACs can be computed within JC model
        polariton.l_pnacme = True
        if (self.force_level == "tdp" and not polariton.l_nacme):
            polariton.l_pnacme = False

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
            self.solve_polaritonic_states(polariton, pol_list, dt, istep)
            self.save_output_files(base_dir, pol_list, istep)
        # Step 3: Calculate properties of polaritonic states; forces, NAC
        self.calculate_properties(polariton, pol_list, calc_force_only)
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

    def solve_polaritonic_states(self, polariton, pol_list, dt, istep):
        """ Diagonalize JC Hamiltonian matrix to get the unitary matrix
            The energy for polaritonic states are directly obtained from the eigenvalues

            :param object polariton: Polariton object
            :param integer,list pol_list: List of polaritonic states for QED calculation
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

        # Check the diabatic character of adjacent two states
        if (istep >= 0 and self.l_check_crossing):
            self.check_trivial_crossing(polariton, pol_list, permut_old, cur_d_ind_old)

    def check_trivial_crossing(self, polariton, pol_list, permut_old, cur_d_ind_old):
        """ Check the diabatic character of adjacent two states using permutation matrix

            :param object polariton: Polariton object
            :param integer,list pol_list: List of polaritonic states for QED calculation
            :param permut_old: Permutation matrix at old step
            :type permut_old: integer, 2D list
            :param cur_d_ind_old: Indicator for diabatic character at old step
            :type cur_d_ind_old: integer, 2D list
        """
        # Check whether the swap of two states occurs
        l_swap = False
        permut_diff_norm = np.linalg.norm(self.permut - permut_old, ord='fro')
        if (permut_diff_norm > 1.):
            l_swap = True

        # Check whether the trivial crossing of two states occurs
        self.l_trivial = False
        if (l_swap):
            # Diabatic character of running state at old step
            ind_mol1 = cur_d_ind_old[pol_list[0], 0]
            ind_photon1 = cur_d_ind_old[pol_list[0], 1]
            # Diabatic character of running state at current step
            ind_mol2 = self.cur_d_ind[pol_list[0], 0]
            ind_photon2 = self.cur_d_ind[pol_list[0], 1]

            # To decide the trivial crossing, we need the Hamiltonian matrix element
            # For this, we must know the corresponding indices for the swap
            # The order of Hamiltonian matrix is matched with self.get_d_ind
            for ist in range(polariton.pst):
                ind_mol = self.get_d_ind[ist, 0]
                ind_photon = self.get_d_ind[ist, 1]
                if (ind_mol == ind_mol1 and ind_photon == ind_photon1):
                    ind_d1 = ist
                if (ind_mol == ind_mol2 and ind_photon == ind_photon2):
                    ind_d2 = ist

            # This condition checks the coupling between two swapped states
            if (ind_d2 != ind_d1 and abs(self.ham_d[ind_d2, ind_d1]) < eps):
                self.l_trivial = True
                for ist in range(polariton.pst):
                    ind_mol = self.cur_d_ind[ist, 0]
                    ind_photon = self.cur_d_ind[ist, 1]
                    if (ind_mol == ind_mol1 and ind_photon == ind_photon1):
                        # This index is used for trivial crossing hopping
                        self.trivial_state = ist

    def calculate_properties(self, polariton, pol_list, calc_force_only):
        """ Calculate properties of polaritonic states
            The forces for polaritonic states are calculated from the gradients
            of electronic states and transition dipole gradients
            The NACs between polaritonic states are calculated from the gradients,
            NACs, and transition dipole gradients

            :param object polariton: Polariton object
            :param integer,list pol_list: List of polaritonic states for QED calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        for rst in pol_list:
            # The applied gradient on the target polaritonic state
            qed_grad = np.zeros((polariton.nat, polariton.ndim))

            # For trivial crossing, the index to be changed must be used
            # for conservation of the total energy
            if (self.l_trivial):
                rst_new = self.trivial_state
            else:
                rst_new = rst

            # First term: Ehrenfest-like term
            # Loop for diabatic states
            for ist in range(polariton.pst):
                ind_mol1 = self.get_d_ind[ist, 0]
                ind_photon1 = self.get_d_ind[ist, 1]
                qed_grad -= polariton.states[ind_mol1].force * self.unitary[ist, rst_new] ** 2.

            if (self.force_level == "nac"):
                # Second term: NACVs multiplied by energy difference; Hellmann-Feynman force contribution
                # Loop for diabatic states
                for ist in range(polariton.pst):
                    ind_mol1 = self.get_d_ind[ist, 0]
                    ind_photon1 = self.get_d_ind[ist, 1]
                    # Loop for diabatic states
                    for jst in range(polariton.pst):
                        ind_mol2 = self.get_d_ind[jst, 0]
                        ind_photon2 = self.get_d_ind[jst, 1]
                        if (ist != jst):
                            qed_grad += polariton.nac[ind_mol2, ind_mol1] * self.unitary[ist, rst_new] \
                                * (self.ham_d[ist, ist] - self.ham_d[jst, jst]) * self.unitary[jst, rst_new]

            elif (self.force_level == "tdp"):
                # Second term: TDP gradients; Pulay force contribution
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
                            qed_grad += field_dot_grad * self.coup_str * self.unitary[ist, rst_new] \
                                * self.unitary[jst, rst_new] * off_term

            polariton.pol_states[rst].force = - np.copy(qed_grad)

        # The nonadiabatic coupling vectors between polaritonic states
        if (not calc_force_only and not polariton.l_pnacme):

            # Loop for adiabatic states
            for ist in range(polariton.pst):
                # Loop for adiabatic states
                for jst in range(ist + 1, polariton.pst):

                    if (ist != jst):
                        qed_nac = np.zeros((polariton.nat, polariton.ndim))

                        # First term: CI term (diagonal contribution)
                        # Loop for diabatic states
                        for ast in range(polariton.pst):
                            ind_mol1 = self.get_d_ind[ast, 0]
                            qed_nac -= polariton.states[ind_mol1].force * self.unitary[ast, ist] * self.unitary[ast, jst]

                        # First term: CI term (off-diagonal contribution)
                        # Loop for diabatic states
                        for ast in range(polariton.pst):
                            ind_mol1 = self.get_d_ind[ast, 0]
                            ind_photon1 = self.get_d_ind[ast, 1]
                            # Loop for diabatic states
                            for bst in range(polariton.pst):
                                ind_mol2 = self.get_d_ind[bst, 0]
                                ind_photon2 = self.get_d_ind[bst, 1]

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

                                if (ast != bst):
                                    field_dot_grad = np.zeros((polariton.nat, polariton.ndim))
                                    for idim in range(polariton.ndim):
                                        field_dot_grad[0:polariton.nat_qm] += polariton.field_pol_vec[idim] \
                                            * polariton.tdp_grad[ind_mol1, ind_mol2, idim]
                                    qed_nac += field_dot_grad * self.coup_str * self.unitary[ast, ist] \
                                        * self.unitary[bst, jst] * off_term

                        qed_nac /= (polariton.pol_states[jst].energy - polariton.pol_states[ist].energy)

                        # Second term: CSF term
                        # Loop for diabatic states
                        for ast in range(polariton.pst):
                            ind_mol1 = self.get_d_ind[ast, 0]
                            ind_photon1 = self.get_d_ind[ast, 1]
                            # Loop for diabatic states
                            for bst in range(polariton.pst):
                                ind_mol2 = self.get_d_ind[bst, 0]
                                ind_photon2 = self.get_d_ind[bst, 1]
                                if (ind_mol1 != ind_mol2 and ind_photon1 == ind_photon2):
                                    qed_nac += polariton.nac[ind_mol1, ind_mol2] * self.unitary[ast, ist] \
                                        * self.unitary[bst, jst]

                        polariton.pnac[ist, jst] = qed_nac
                        polariton.pnac[jst, ist] = - qed_nac

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

    def calculate_pnacme(self, polariton):
        """ Calculate NACME between polaritonic states
            It consists of NACMEs and unitary matrix derivatives

            :param object polariton: Polariton object
        """
        tmp_nacme = np.zeros((polariton.pst, polariton.pst)) 

        for ist in range(polariton.pst):
            ind_mol1 = self.get_d_ind[ist, 0]
            ind_photon1 = self.get_d_ind[ist, 1]
            for jst in range(ist + 1, polariton.pst):
                ind_mol2 = self.get_d_ind[jst, 0]
                ind_photon2 = self.get_d_ind[jst, 1]

                if (ind_photon1 == ind_photon2):
                    tmp_nacme[ist, jst] = polariton.nacme[ind_mol1, ind_mol2]
                    tmp_nacme[jst, ist] = - tmp_nacme[ist, jst]

        polariton.pnacme += np.matmul(np.transpose(self.unitary), np.matmul(tmp_nacme, self.unitary))
        polariton.pnacme += np.matmul(np.transpose(self.unitary), self.unitary_dot)

    def backup_qed(self, polariton):
        """ Backup Hamiltonian matrix, unitary matrix and polaritonic state energies

            :param object polariton: Polariton object
        """
        for states in polariton.pol_states:
            states.energy_old = states.energy

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

    def print_init(self):
        """ Print initial information about JC model
        """
        qed_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"QED Information":>40s}
        {"-" * 68}
          QED Method               = {self.qed_method:>16s}
          Coupling Strength (au)   = {self.coup_str:>16.6f}
          Force Computation Method = {self.force_level:>16s}
        """)

        # Print RWA information
        if (not self.l_crt):
            qed_info += f"  Rotating Wave Approx.    = {'Yes':>16s}\n"
        else:
            qed_info += f"  Rotating Wave Approx.    = {'No':>16s}\n"
        print (qed_info, flush=True)


