from __future__ import division
from build.el_propagator_ct import el_run
from mqc.mqc import MQC
from misc import eps, au_to_K, au_to_A, call_name, typewriter, gaussian1d
import os, shutil, textwrap
import numpy as np
import pickle

class CT(MQC):
    """ Class for coupled-trajectory mixed quantum-classical (CTMQC) dynamics

        :param object,list molecules: List for molecule objects
        :param object thermostat: Thermostat object
        :param integer,list istates: List for initial state
        :param double dt: Time interval
        :param integer nsteps: Total step of nuclear propagation
        :param integer nesteps: Total step of electronic propagation
        :param string elec_object: Electronic equation of motions
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: Logical to print BO population and coherence
        :param boolean l_adj_nac: Adjust nonadiabatic coupling to align the phases
        :param double rho_threshold: Electronic density threshold for decoherence term calculation
        :param double sigma_threshold: Sigma threshold for quantum momentum calculation
        :param double dist_cutoff: Distance cutoff for quantum momentum calculation
        :param double dist_parameter: Distance parameter to determine quantum momentum center
        :param double sigma: Sigma to determine quantum momentum center
        :param init_coefs: Initial BO coefficient
        :type init_coefs: double, 2D list or complex, 2D list
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecules, thermostat=None, istates=None, dt=0.5, nsteps=1000, nesteps=20, \
        elec_object="coefficient", propagator="rk4", l_print_dm=True, l_adj_nac=True, \
        rho_threshold=0.01, sigma_threshold=0.25, dist_cutoff=0.5, dist_parameter=10., sigma=0.3, \
        init_coefs=None, unit_dt="fs", out_freq=1, verbosity=0):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize input values
        self.mols = molecules
        self.ntrajs = len(self.mols)
        self.digit = len(str(self.ntrajs))

        self.nst = self.mols[0].nst
        self.nat_qm = self.mols[0].nat_qm
        self.ndim = self.mols[0].ndim

        # Check compatibility between istates and init_coefs
        self.istates = istates
        self.init_coefs = init_coefs
        self.check_istates()

        # Initialize input values and coefficient for first trajectory
        super().__init__(self.mols[0], thermostat, self.istates[0], dt, nsteps, nesteps, \
            elec_object, propagator, l_print_dm, l_adj_nac, self.init_coefs[0], unit_dt, out_freq, verbosity)

        if (self.elec_object != "coefficient"):
            error_message = "Electronic equation motion in CTMQC is only solved with respect to coefficient!"
            error_vars = f"elec_object = {self.elec_object}"
            raise NotImplementedError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Initialize coefficient for other trajectories
        for itraj in range(1, self.ntrajs):
            self.mols[itraj].get_coefficient(self.init_coefs[itraj], self.istates[itraj])

        # Initialize variables for CTMQC
        self.phase = np.zeros((self.ntrajs, self.nst, self.nat_qm, self.ndim))
        self.nst_pair = int(self.nst * (self.nst - 1) / 2)
        self.qmom = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim))
        self.K_lk = np.zeros((self.ntrajs, self.nst, self.nst))

        # Initialize variables to calculate quantum momentum 
        self.count_ntrajs = np.zeros((self.ntrajs, self.nat_qm))
        self.sigma_lk = np.ones((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim))
        self.slope_i = np.zeros((self.ntrajs, self.nat_qm, self.ndim))
        self.center_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim))

        # Determine parameters to calculate decoherenece effect
        self.small = 1.0E-08

        self.rho_threshold = rho_threshold

        self.upper_th = 1. - self.rho_threshold
        self.lower_th = self.rho_threshold

        self.sigma_threshold = sigma_threshold
        self.dist_cutoff = dist_cutoff

        self.dist_parameter = dist_parameter
        self.sigma = sigma

        self.dotpopnac = np.zeros((self.ntrajs, self.nst))
        self.dotpopdec = np.zeros((self.ntrajs, self.nst))

    def run(self, qm, mm=None, output_dir="./", l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to CTMQC dynamics

            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string output_dir: Name of directory where outputs to be saved.
            :param boolean l_save_qm_log: Logical for saving QM calculation log
            :param boolean l_save_mm_log: Logical for saving MM calculation log
            :param boolean l_save_scr: Logical for saving scratch directory
            :param string restart: Option for controlling dynamics restarting
        """
        # Initialize PyUNIxMD
        base_dirs, unixmd_dirs, qm_log_dirs, mm_log_dirs =\
             self.run_init(qm, mm, output_dir, l_save_qm_log, l_save_mm_log, l_save_scr, restart)

        bo_list = [ist for ist in range(self.nst)]
        qm.calc_coupling = True

        self.print_init(qm, mm, restart)

        if (restart == None):
            # Calculate initial input geometry for all trajectories at t = 0.0 s
            self.istep = -1
            for itraj in range(self.ntrajs):
                self.mol = self.mols[itraj]

                self.mol.reset_bo(qm.calc_coupling)
                qm.get_data(self.mol, base_dirs[itraj], bo_list, self.dt, self.istep, calc_force_only=False)

                # TODO: QM/MM
                self.mol.get_nacme()

                self.update_energy()

                self.get_phase(itraj)

            self.calculate_qmom(self.istep)

            for itraj in range(self.ntrajs):

                self.mol = self.mols[itraj]

                self.write_md_output(itraj, unixmd_dirs[itraj], self.istep)

                self.print_step(self.istep, itraj)

        #TODO: restart
        else: 
            error_message = "Restart option with CTMQC not implemented!"
            error_vars = f"restart = {restart}"
            raise NotImplementedError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.istep += 1

        # Main MD loop
        for istep in range(self.istep, self.nsteps):
            for itraj in range(self.ntrajs):
                self.mol = self.mols[itraj]

                self.calculate_force(itraj)
                self.cl_update_position()

                self.mol.backup_bo()
                self.mol.reset_bo(qm.calc_coupling)
                qm.get_data(self.mol, base_dirs[itraj], bo_list, self.dt, istep, calc_force_only=False)
                #TODO: QM/MM

                if (not self.mol.l_nacme and self.l_adj_nac):
                    self.mol.adjust_nac()

                self.calculate_force(itraj)
                self.cl_update_velocity()

                self.mol.get_nacme()

                el_run(self, itraj)

                #TODO: thermostat
                #if (self.thermo != None):
                #    self.thermo.run(self)

                self.update_energy()

                self.get_phase(itraj)

                #TODO: restart
                #self.fstep = istep
                #restart_file = os.path.join(base_dir, "RESTART.bin")
                #with open(restart_file, 'wb') as f:
                #    pickle.dump({'qm':qm, 'md':self}, f)

            self.calculate_qmom(istep)

            for itraj in range(self.ntrajs):
                self.mol = self.mols[itraj]

                if ((istep + 1) % self.out_freq == 0):
                    self.write_md_output(itraj, unixmd_dirs[itraj], istep)
                    self.print_step(istep, itraj)
                if (istep == self.nsteps - 1):
                    self.write_final_xyz(unixmd_dirs[itraj], istep)

        # Delete scratch directory
        if (not l_save_scr):
            for itraj in range(self.ntrajs):
                tmp_dir = os.path.join(unixmd_dirs[itraj], "scr_qm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def calculate_force(self, itrajectory):
        """ Routine to calculate force

            :param integer itrajectory: Index for trajectories
        """
        self.rforce = np.zeros((self.nat_qm, self.ndim))

        # Derivatives of energy
        for ist, istate in enumerate(self.mols[itrajectory].states):
            self.rforce += istate.force * self.mol.rho.real[ist, ist]

        # Non-adiabatic forces 
        for ist in range(self.nst):
            for jst in range(ist + 1, self.nst):
                self.rforce += 2. * self.mol.nac[ist, jst] * self.mol.rho.real[ist, jst] \
                    * (self.mol.states[ist].energy - self.mol.states[jst].energy)

        # CT forces
        ctforce = np.zeros((self.nat_qm, self.ndim))
        for ist in range(self.nst):
            for jst in range(self.nst):
                ctforce += 0.5 * self.K_lk[itrajectory, ist, jst] * \
                    (self.phase[itrajectory, jst] - self.phase[itrajectory, ist]) * \
                    self.mol.rho.real[ist, ist]  * self.mol.rho.real[jst, jst]

        # Finally, force is Ehrenfest force + CT force
        self.rforce += ctforce

    def update_energy(self):
        """ Routine to update the energy of molecules in CTMQC dynamics
        """
        # Update kinetic energy
        self.mol.update_kinetic()
        self.mol.epot = 0.
        for ist, istate in enumerate(self.mol.states):
            self.mol.epot += self.mol.rho.real[ist, ist] * istate.energy
        self.mol.etot = self.mol.epot + self.mol.ekin

    def get_phase(self, itrajectory):
        """ Routine to calculate phase

            :param integer itrajectory: Index for trajectories
        """
        for ist in range(self.nst):
            rho = self.mol.rho[ist, ist].real
            if (rho > self.upper_th or rho < self.lower_th):
                self.phase[itrajectory, ist] = np.zeros((self.nat_qm, self.ndim))
            else:
                self.phase[itrajectory, ist] += self.mol.states[ist].force * self.dt

    def calculate_qmom(self, istep):
        """ Routine to calculate quantum momentum

            :param integer istep: Current MD step
        """
        # _lk means state_pair dependency.
        # i and j are trajectory index.
        # -------------------------------------------------------------------
        # 1. Calculate variances for each trajectory
        # TODO: method to calculate sigma
        self.sigma_lk = np.ones((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim)) # TODO: state-pair
        for itraj in range(self.ntrajs):
            # Variable to determine how many trajecories are in cutoff.
            self.count_ntrajs[itraj] = np.zeros((self.nat_qm)) 

            R2_tmp = np.zeros((self.nat_qm, self.ndim)) # Temporary variable for R**2
            R_tmp = np.zeros((self.nat_qm, self.ndim))  # Temporary variable for R

            for jtraj in range(self.ntrajs):
                pos_diff = self.mols[jtraj].pos - self.mols[itraj].pos # Dimension = (self.nat_qm, self.ndim)
                pos_diff2 = np.sum(pos_diff * pos_diff, axis=1) # Dimension = (self.nat_qm)

                for iat in range(self.nat_qm):
                    distance = np.sqrt(pos_diff2[iat]) # Distance between i-th atom in itraj and jtraj
                    if (distance <= self.dist_cutoff):
                        R_tmp[iat] += self.mols[jtraj].pos[iat] # Dimension = (self.nat_qm, self.ndim)
                        R2_tmp[iat] += self.mols[jtraj].pos[iat] * self.mols[jtraj].pos[iat] # Dimension = (self.nat_qm, self.ndim)
                        self.count_ntrajs[itraj, iat] += 1

            for iat in range(self.nat_qm):
                avg_R = R_tmp[iat] / self.count_ntrajs[itraj, iat]
                avg_R2 = R2_tmp[iat] / self.count_ntrajs[itraj, iat]
                for idim in range(self.ndim):
                    self.sigma_lk[itraj, 0, iat, idim] = np.sqrt((avg_R2[idim] - avg_R[idim] ** 2)) \
                        / np.sqrt(np.sqrt(self.count_ntrajs[itraj, iat])) # / np.sqrt(np.sqrt(count_ntrajs)) is artifact to modulate sigma.
                    if (self.sigma_lk[itraj, 0, iat, idim] <= self.sigma_threshold):
                        self.sigma_lk[itraj, 0, iat, idim] = self.dist_cutoff

        # 2. Calculate slope
        # (2-1) Calculate w_ij
        # g_i means nuclear density at the position of i-th classical trajectory.
        # prod_g_i is to multiply gaussians with respect to atoms and spaces.
        g_i = np.zeros((self.ntrajs)) 
        prod_g_i = np.ones((self.ntrajs, self.ntrajs))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat_qm):
                    for idim in range(self.ndim):
                        # gaussian1d(x, pre-factor, sigma, mean)
                        # gaussian1d(R^{itraj}, 1.0, sigma^{jtraj}, R^{jtraj})
                        prod_g_i[itraj, jtraj] *= gaussian1d(self.mols[itraj].pos[iat, idim], 1., \
                            self.sigma_lk[jtraj, 0, iat, idim], self.mols[jtraj].pos[iat, idim])
                g_i[itraj] += prod_g_i[itraj, jtraj]

        # w_ij is defined as W_IJ in SI of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        w_ij = np.zeros((self.ntrajs, self.ntrajs, self.nat_qm, self.ndim))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat_qm):
                    for idim in range(self.ndim):
                        w_ij[itraj, jtraj, iat, idim] = prod_g_i[itraj, jtraj] /\
                        (2. * self.sigma_lk[jtraj, 0, iat, idim] ** 2 * g_i[itraj])

        # (2-2) Calculate slope_i
        # the slope is calculated as a sum over j of w_ij
        self.slope_i = np.zeros((self.ntrajs, self.nat_qm, self.ndim))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                self.slope_i[itraj] -= w_ij[itraj, jtraj]

        # 3. Calculate the center of quantum momentum
        rho = np.zeros((self.ntrajs, self.nst))
        for itraj in range(self.ntrajs):
            for ist in range(self.nst):
                rho[itraj, ist] = self.mols[itraj].rho[ist, ist].real

        # (3-1) Compute denominator
        deno_lk = np.zeros((self.nst_pair, self.nat_qm, self.ndim)) # denominator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):  
                    index_lk += 1
                    for iat in range(self.nat_qm):
                        for idim in range(self.ndim):
                            deno_lk[index_lk, iat, idim] += rho[itraj, ist] * rho[itraj, jst] * \
                                (self.phase[itraj, ist, iat, idim] - self.phase[itraj, jst, iat, idim]) * self.slope_i[itraj, iat, idim]

        # (3-2) Compute numerator
        ratio_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim)) # numerator / denominator
        numer_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim)) # numerator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat_qm):
                        for idim in range(self.ndim):
                            numer_lk[itraj, index_lk, iat, idim] = rho[itraj, ist] * rho[itraj, jst] * self.mols[itraj].pos[iat, idim] * \
                                (self.phase[itraj, ist, iat, idim] - self.phase[itraj, jst, iat, idim]) * self.slope_i[itraj, iat, idim]
                            if (abs(deno_lk[index_lk, iat, idim]) <= self.small):
                                ratio_lk[itraj, index_lk, iat, idim] = 0.
                            else:
                                ratio_lk[itraj, index_lk, iat, idim] = numer_lk[itraj, index_lk, iat, idim] / \
                                    deno_lk[index_lk, iat, idim]

        # Center of quantum momentum is calculated by Eq.(S28) of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        center_old_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat_qm):
                        for idim in range(self.ndim):
                            for jtraj in range(self.ntrajs):
                                center_old_lk[itraj, index_lk, iat, idim] += ratio_lk[jtraj, index_lk, iat, idim]
                            if ((abs(self.slope_i[itraj, iat, idim]) <= self.small) or (center_old_lk[itraj, index_lk, iat, idim] == 0.)):
                                center_old_lk[itraj, index_lk, iat, idim] = self.mols[itraj].pos[iat, idim]

        # Center of quantum momentum is calculated by Eq.(S21) of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        center_new_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat_qm):
                        for idim in range(self.ndim):
                            if (abs(self.slope_i[itraj, iat, idim]) <= self.small):
                                center_new_lk[itraj, index_lk, iat, idim] = self.mols[itraj].pos[iat, idim]
                            else:
                                for jtraj in range(self.ntrajs):
                                    center_new_lk[itraj, index_lk, iat, idim] += self.mols[jtraj].pos[iat, idim] * prod_g_i[itraj, jtraj] /\
                                        (2. * self.sigma_lk[jtraj, 0, iat, idim] ** 2 * g_i[itraj] * (- self.slope_i[itraj, iat, idim]))

        # (3-3) Determine qauntum momentum center TODO: atomistic flag
        self.center_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim)) # Finally, qmom_center
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat_qm):
                        for idim in range(self.ndim):
                            # test how far calculated center of quantum momentum is from current atomic position.
                            # tmp_var is deviation between position of classical trajectory and quantum momentum center.
                            tmp_var = center_old_lk[itraj, index_lk, iat, idim] - self.mols[itraj].pos[iat, idim]
                            if (abs(tmp_var) > self.dist_parameter * self.sigma): 
                                tmp_var = center_new_lk[itraj, index_lk, iat, idim] - self.mols[itraj].pos[iat, idim]
                                if (abs(tmp_var) > self.dist_parameter * self.sigma): 
                                    self.center_lk[itraj, index_lk, iat, idim] = self.mols[itraj].pos[iat, idim]
                                else:
                                    self.center_lk[itraj, index_lk, iat, idim] = center_new_lk[itraj, index_lk, iat, idim]
                            else: 
                                self.center_lk[itraj, index_lk, iat, idim] = center_old_lk[itraj, index_lk, iat, idim]

        # 4. Compute quantum momentum
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    self.qmom[itraj, index_lk] = self.slope_i[itraj] * (self.mols[itraj].pos - self.center_lk[itraj, index_lk])

        # 5. Calculate 2 * Qmom * phase / mass
        self.K_lk = np.zeros((self.ntrajs, self.nst, self.nst))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    self.K_lk[itraj, ist, jst] += 2. * np.sum(1. / self.mol.mass[0:self.nat_qm] * \
                        np.sum(self.qmom[itraj, index_lk] * self.phase[itraj, ist], axis = 1))
                    self.K_lk[itraj, jst, ist] += 2. * np.sum(1. / self.mol.mass[0:self.nat_qm] * \
                        np.sum(self.qmom[itraj, index_lk] * self.phase[itraj, jst], axis = 1))

    def check_istates(self):
        """ Routine to check istates and init_coefs
        """
        if (self.istates != None):
            if (isinstance(self.istates, list)):
                if (len(self.istates) != self.ntrajs):
                    error_message = "Number of elements of initial states must be equal to number of trajectories!"
                    error_vars = f"len(istates) = {len(self.istates)}, ntrajs = {self.ntrajs}"
                    raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
                else:
                    self.init_coefs = [None] * self.ntrajs
            else:
                error_message = "The type of initial states must be list!"
                error_vars = f"istates = {self.istates}"
                raise TypeError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            if (self.init_coefs == None):
                error_message = "Either initial states or coefficients must be given!"
                error_vars = f"istates = {self.istates}, init_coefs = {self.init_coefs}"
                raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                if (isinstance(self.init_coefs, list)):
                    if (len(self.init_coefs) != self.ntrajs):
                        error_message = "Number of elements of initial coefficients must be equal to number of trajectories!"
                        error_vars = f"len(init_coefs) = {len(self.init_coefs)}, ntrajs = {self.ntrajs}"
                        raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
                    else:
                        self.istates = [None] * self.ntrajs
                else:
                    error_message = "Type of initial coefficients must be list!"
                    error_vars = f"init_coefs = {self.init_coefs}"
                    raise TypeError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

    def write_md_output(self, itrajectory, unixmd_dir, istep):
        """ Write output files

            :param integer itrajectory: Index for trajectories
            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, istep)

        # Write time-derivative BO population
        self.write_dotpop(itrajectory, unixmd_dir, istep)

        # Write decoherence information
        self.write_dec(itrajectory, unixmd_dir, istep)

    def write_dotpop(self, itrajectory, unixmd_dir, istep):
        """ Write time-derivative BO population

            :param integer itrajectory: Index for trajectories
            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        if (self.verbosity >= 1):
            # Write NAC term in DOTPOPNAC
            tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopnac[itrajectory]])
            typewriter(tmp, unixmd_dir, "DOTPOPNAC", "a")

            # Write decoherence term in DOTPOPDEC
            tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopdec[itrajectory]])
            typewriter(tmp, unixmd_dir, "DOTPOPDEC", "a")

    def write_dec(self, itrajectory, unixmd_dir, istep):
        """ Write CT-based decoherence information

            :param integer itrajectory: Index for trajectories
            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        if (self.verbosity >= 1):
            # Write K_lk
            for ist in range(self.nst):
                for jst in range(self.nst):
                    if (ist != jst):
                        tmp = f'{istep + 1:9d}{self.K_lk[itrajectory, ist, jst]:15.8f}'
                        typewriter(tmp, unixmd_dir, f"K_lk_{ist}_{jst}", "a")

        # Write detailed quantities related to decoherence
        if (self.verbosity >= 2):
            tmp = f'{self.nat_qm:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}sigma_x{"":5s}sigma_y{"":5s}sigma_z{"":5s}count_ntrajs' + \
                "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                "".join([f'{self.sigma_lk[itrajectory, 0, iat, idim]:15.8f}' for idim in range(self.ndim)]) + \
                f'{self.count_ntrajs[itrajectory, iat]:15.8f}' for iat in range(self.nat_qm)])
            typewriter(tmp, unixmd_dir, f"SIGMA", "a")

            tmp = f'{self.nat_qm:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}slope' + \
                "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                "".join([f'{self.slope_i[itrajectory, iat, idim]:15.8f}' for idim in range(self.ndim)]) for iat in range(self.nat_qm)])
            typewriter(tmp, unixmd_dir, f"SLOPE", "a")

            # Write quantum momenta
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    tmp = f'{self.nat_qm:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Momentum center (au)' + \
                        "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                        "".join([f'{self.center_lk[itrajectory, index_lk, iat, idim]:15.8f}' for idim in range(self.ndim)]) for iat in range(self.nat_qm)])
                    typewriter(tmp, unixmd_dir, f"CENTER_{ist}_{jst}", "a")

                    tmp = f'{self.nat_qm:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Momentum (au)' + \
                        "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                        "".join([f'{self.qmom[itrajectory, index_lk, iat, idim]:15.8f}' for idim in range(self.ndim)]) for iat in range(self.nat_qm)])
                    typewriter(tmp, unixmd_dir, f"QMOM_{ist}_{jst}", "a")

            # Write Phase
            for ist in range(self.mol.nst):
                tmp = f'{self.nat_qm:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Phase (au)' + \
                    "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                    "".join([f'{self.phase[itrajectory, ist, iat, idim]:15.8f}' for idim in range(self.ndim)]) for iat in range(self.nat_qm)])
                typewriter(tmp, unixmd_dir, f"PHASE_{ist}", "a")

    def print_init(self, qm, mm, restart):
        """ Routine to print the initial information of dynamics

            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string restart: Option for controlling dynamics restarting
        """
        # Print initial information about molecule, qm, mm and thermostat
        super().print_init(qm, mm, restart)

        # Print dynamics information for start line
        dynamics_step_info = textwrap.dedent(f"""\

        {"-" * 118}
        {"Start Dynamics":>65s}
        {"-" * 118}
        """)

        # Print INIT for each trajectory at each step
        INIT = f" #INFO_TRAJ{'STEP':>8s}{'Kinetic(H)':>15s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'norm':>8s}"
        dynamics_step_info += INIT

        print (dynamics_step_info, flush=True)

    def print_step(self, istep, itrajectory):
        """ Routine to print each trajectory infomation at each step about dynamics

            :param integer istep: Current MD step
            :param integer itrajectory: Current trajectory
        """
        ctemp = self.mol.ekin * 2. / float(self.mol.ndof) * au_to_K
        norm = 0.
        for ist in range(self.mol.nst):
            norm += self.mol.rho.real[ist, ist]

        # Print INFO for each step
        INFO = f" INFO_{itrajectory+1}{istep + 1:>9d}"
        INFO += f"{self.mol.ekin:14.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        INFO += f"{norm:11.5f}"
        print (INFO, flush=True)
