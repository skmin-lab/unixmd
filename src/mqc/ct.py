from __future__ import division
from lib.libctmqc import el_run
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
        :param init_coefs: Initial BO coefficient
        :type init_coefs: double, 2D list or complex, 2D list
        :param double dist_parameter: Distance parameter to contruct Gaussian and determine quantum momentum center
        :param double min_sigma: Minimum sigma value
        :param double const_dist_cutoff: Distance cutoff to construct Gaussian
        :param double const_center_cutoff: Distance cutoff to determine quantum momentum center
        :param string unit_dt: Unit of time interval
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecules, thermostat=None, istates=None, dt=0.5, nsteps=1000, nesteps=20, \
        elec_object="coefficient", propagator="rk4", l_print_dm=True, l_adj_nac=True, rho_threshold=0.01, \
        init_coefs=None, dist_parameter=10., min_sigma=0.3, const_dist_cutoff=None, const_center_cutoff=None, \
        l_en_cons=False, unit_dt="fs", out_freq=1, verbosity=0):
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

        # Exception for electronic propagation
        if (self.elec_object != "coefficient"):
            error_message = "Electronic equation motion in CTMQC is only solved with respect to coefficient!"
            error_vars = f"elec_object = {self.elec_object}"
            raise NotImplementedError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Exception for thermostat
        if (self.thermo != None):
            error_message = "Thermostat is not implemented yet!"
            error_vars = f"thermostat = {self.thermo}"
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
        self.count_ntrajs = np.zeros((self.ntrajs, self.nat_qm, self.ndim))
        self.sigma_lk = np.ones((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim))
        self.slope_i = np.zeros((self.ntrajs, self.nat_qm, self.ndim))
        self.g_i = np.zeros((self.ntrajs)) 
        self.prod_g_i = np.ones((self.ntrajs, self.ntrajs))
        self.center_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat_qm, self.ndim))

        # Determine parameters to calculate decoherenece effect
        self.small = 1.0E-08

        self.rho_threshold = rho_threshold
        self.upper_th = 1. - self.rho_threshold
        self.lower_th = self.rho_threshold

        self.min_sigma = min_sigma
        self.const_dist_cutoff = const_dist_cutoff
        self.dist_parameter = dist_parameter
        self.const_center_cutoff = const_center_cutoff

        self.l_en_cons = l_en_cons
        self.dotpopnac = np.zeros((self.ntrajs, self.nst))
        self.dotpopdec = np.zeros((self.ntrajs, self.nst))

        # Initialize event to print
        self.event = {"DECO": []}

    def run(self, qm, mm=None, output_dir="./", l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to CTMQC dynamics

            :param object qm: QM object containing on-the-fly calculation information
            :param object mm: MM object containing MM calculation information
            :param string output_dir: Name of directory where outputs to be saved.
            :param boolean l_save_qm_log: Logical for saving QM calculation log
            :param boolean l_save_mm_log: Logical for saving MM calculation log
            :param boolean l_save_scr: Logical for saving scratch directory
            :param string restart: Option for controlling dynamics restarting
        """
        # Initialize PyUNIxMD
        abs_path_output_dir = os.path.join(os.getcwd(), output_dir)
        base_dirs, unixmd_dirs, traj_bin_dirs, qm_log_dirs, mm_log_dirs = \
            self.run_init(qm, mm, output_dir, False, l_save_qm_log, l_save_mm_log, \
            l_save_scr, restart)

        bo_list = [ist for ist in range(self.nst)]
        qm.calc_coupling = True
        qm.calc_tdp = False
        qm.calc_tdp_grad = False

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

                self.write_md_output(itraj, unixmd_dirs[itraj], qm.calc_coupling, self.istep)

                self.print_step(self.istep, itraj)

        #TODO: restart
        elif (restart == "write"):
            # Reset initial time step to t = 0.0 s
            self.istep = -1
            for itraj in range(self.ntrajs):
                self.write_md_output(itraj, unixmd_dirs[itraj], qm.calc_coupling, self.istep)
                self.print_step(self.istep, itraj)

        elif (restart == "append"):
            # Set initial time step to last successful step of previous dynamics
            self.istep = self.fstep

        self.istep += 1

        # Main MD loop
        for istep in range(self.istep, self.nsteps):
            for itraj in range(self.ntrajs):
                self.mol = self.mols[itraj]

                self.calculate_force(itraj)
                self.cl_update_position()

                self.mol.backup_bo(qm.calc_coupling)
                self.mol.reset_bo(qm.calc_coupling)

                qm.get_data(self.mol, base_dirs[itraj], bo_list, self.dt, istep, calc_force_only=False)

                if (not self.mol.l_nacme and self.l_adj_nac):
                    self.mol.adjust_nac()

                #TODO: QM/MM

                self.calculate_force(itraj)
                self.cl_update_velocity()

                self.mol.get_nacme()

                el_run(self, itraj)

                #TODO: thermostat
                #if (self.thermo != None):
                #    self.thermo.run(self, self.mol)

                self.update_energy()

                self.get_phase(itraj)

                self.check_decoherence(itraj)

            self.calculate_qmom(istep)

            for itraj in range(self.ntrajs):
                self.mol = self.mols[itraj]

                if ((istep + 1) % self.out_freq == 0):
                    self.write_md_output(itraj, unixmd_dirs[itraj], qm.calc_coupling, istep)
                    self.print_step(istep, itraj)
                if (istep == self.nsteps - 1):
                    self.write_final_xyz(unixmd_dirs[itraj], istep)

            self.fstep = istep
            restart_file = os.path.join(abs_path_output_dir, "RESTART.bin")
            with open(restart_file, 'wb') as f:
                pickle.dump({'qm':qm, 'md':self}, f)

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
                    self.mol.rho.real[ist, ist] * self.mol.rho.real[jst, jst] * \
                    (self.mol.rho.real[ist, ist] + self.mol.rho.real[jst, jst])
        ctforce /= self.nst - 1

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
        
        if (self.l_en_cons and not (self.istep == -1)):
            alpha = (self.mol.etot - self.mol.epot)
            factor = alpha / self.mol.ekin

            self.mol.vel *= np.sqrt(factor)
            self.mol.update_kinetic()

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

    def check_decoherence(self, itrajectory):
        """ Routine to check decoherence among BO states

            :param integer itrajectory: Index for trajectories
        """
        for ist in range(self.mol.nst):
            if (np.sum(abs(self.phase[itrajectory, ist])) > eps):
                rho = self.mol.rho.real[ist, ist]
                if (rho > self.upper_th):
                    self.set_decoherence(ist)
                    self.event["DECO"].append(f"{itrajectory + 1:8d}: decohered to {ist} state")
                    return

    def set_decoherence(self, one_st):
        """ Routine to reset coefficient/density if the state is decohered

            :param integer one_st: State index that its population is one
        """
        self.mol.rho = np.zeros((self.mol.nst, self.mol.nst), dtype=np.complex64)
        self.mol.rho[one_st, one_st] = 1. + 0.j
        
        if (self.elec_object == "coefficient"):
            for ist in range(self.mol.nst):
                if (ist == one_st):
                    self.mol.states[ist].coef /= np.absolute(self.mol.states[ist].coef).real
                else:
                    self.mol.states[ist].coef = 0. + 0.j

    def calculate_qmom(self, istep):
        """ Routine to calculate quantum momentum

            :param integer istep: Current MD step
        """
        # _lk means state_pair dependency.
        # i and j are trajectory index.
        # -------------------------------------------------------------------
        # 1. Calculate variances for each trajectory
        self.calculate_sigma(istep)

        # 2. Calculate slope
        self.calculate_slope()

        # 3. Calculate the center of quantum momentum
        self.calculate_center()

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

    def calculate_sigma(self, istep):
        """ Routine to calculate variances for each trajectories

            :param integer istep: Current MD step
        """
        threshold = self.dist_parameter * self.min_sigma #/ self.ntrajs
        cutoff = np.ones((self.ntrajs, self.nat_qm, self.ndim))
        for itraj in range(self.ntrajs):
            if (self.const_dist_cutoff == None):
                if (istep == -1):
                    cutoff *= threshold
                else:
                    cutoff[itraj] = self.dist_parameter * self.sigma_lk[itraj, 0] 
            else:
                cutoff = self.const_dist_cutoff * np.ones((self.ntrajs, self.nat_qm, self.ndim))

            # Variable to determine how many trajecories are in cutoff.
            self.count_ntrajs[itraj] = np.zeros((self.nat_qm, self.ndim)) 

            R2_tmp = np.zeros((self.nat_qm, self.ndim)) # Temporary variable for R**2
            R_tmp = np.zeros((self.nat_qm, self.ndim))  # Temporary variable for R

            for jtraj in range(self.ntrajs):
                pos_diff = self.mols[jtraj].pos - self.mols[itraj].pos # Dimension = (self.nat_qm, self.ndim)
                for iat in range(self.nat_qm):
                    for idim in range(self.ndim):
                        distance = abs(pos_diff[iat, idim]) # Distance between i-th dimenstion of i-th atom in itraj and jtraj
                        if (distance <= cutoff[itraj, iat, idim]):
                            R_tmp[iat, idim] += self.mols[jtraj].pos[iat, idim] # Dimension = (self.nat_qm, self.ndim)
                            R2_tmp[iat, idim] += self.mols[jtraj].pos[iat, idim] * self.mols[jtraj].pos[iat, idim] # Dimension = (self.nat_qm, self.ndim)
                            self.count_ntrajs[itraj, iat, idim] += 1

            for iat in range(self.nat_qm):
                for idim in range(self.ndim):
                    avg_R = R_tmp[iat, idim] / self.count_ntrajs[itraj, iat, idim]
                    avg_R2 = R2_tmp[iat, idim] / self.count_ntrajs[itraj, iat, idim]

                    self.sigma_lk[itraj, 0, iat, idim] = np.sqrt((avg_R2 - avg_R ** 2)) \
                        / np.sqrt(np.sqrt(self.count_ntrajs[itraj, iat, idim])) # / np.sqrt(np.sqrt(count_ntrajs)) is artifact to modulate sigma.
                    if (self.sigma_lk[itraj, 0, iat, idim] <= self.min_sigma or self.count_ntrajs[itraj, iat, idim] == 1):
                        self.sigma_lk[itraj, 0, iat, idim] = self.min_sigma

    def calculate_slope(self):
        """ Routine to calculate slope
        """
        # (2-1) Calculate w_ij
        # g_i means nuclear density at the position of i-th classical trajectory.
        # prod_g_i is to multiply gaussians with respect to atoms and spaces.
        self.g_i = np.zeros((self.ntrajs)) 
        self.prod_g_i = np.ones((self.ntrajs, self.ntrajs))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat_qm):
                    for idim in range(self.ndim):
                        # gaussian1d(x, pre-factor, sigma, mean)
                        # gaussian1d(R^{itraj}, 1.0, sigma^{jtraj}, R^{jtraj})
                        self.prod_g_i[itraj, jtraj] *= gaussian1d(self.mols[itraj].pos[iat, idim], 1., \
                            self.sigma_lk[jtraj, 0, iat, idim], self.mols[jtraj].pos[iat, idim])
                self.g_i[itraj] += self.prod_g_i[itraj, jtraj]

        # w_ij is defined as W_IJ in SI of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        w_ij = np.zeros((self.ntrajs, self.ntrajs, self.nat_qm, self.ndim))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat_qm):
                    for idim in range(self.ndim):
                        w_ij[itraj, jtraj, iat, idim] = self.prod_g_i[itraj, jtraj] /\
                        (2. * self.sigma_lk[jtraj, 0, iat, idim] ** 2 * self.g_i[itraj])

        # Smoothing 
        self.w_k = np.zeros((self.ntrajs, self.nst))
        rho = np.zeros((self.ntrajs, self.nst))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for ist in range(self.nst):
                    self.w_k[itraj, ist] += self.prod_g_i[itraj, jtraj] * self.mols[jtraj].rho.real[ist, ist] / self.g_i[itraj]

            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1

                    l_smooth = ((self.w_k[itraj, ist] < self.lower_th) or (self.w_k[itraj, ist] > self.upper_th) 
                        or (self.w_k[itraj, jst] < self.lower_th) or (self.w_k[itraj, jst] > self.upper_th))

                    if (l_smooth):
                        self.phase[itraj, ist] = np.zeros((self.nat_qm, self.ndim))
                        self.phase[itraj, jst] = np.zeros((self.nat_qm, self.ndim))

        # (2-2) Calculate slope_i
        # the slope is calculated as a sum over j of w_ij
        self.slope_i = np.zeros((self.ntrajs, self.nat_qm, self.ndim))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                self.slope_i[itraj] -= w_ij[itraj, jtraj]

    def calculate_center(self):
        """ Routine to calculate center of quantum momentum
        """
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
                                    center_new_lk[itraj, index_lk, iat, idim] += self.mols[jtraj].pos[iat, idim] * self.prod_g_i[itraj, jtraj] /\
                                        (2. * self.sigma_lk[jtraj, 0, iat, idim] ** 2 * self.g_i[itraj] * (- self.slope_i[itraj, iat, idim]))

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
                            if (self.const_center_cutoff == None):
                                cutoff = self.dist_parameter * self.sigma_lk[itraj, 0, iat, idim]
                            else:
                                cutoff = self.const_center_cutoff

                            if (abs(tmp_var) > cutoff): 
                                tmp_var = center_new_lk[itraj, index_lk, iat, idim] - self.mols[itraj].pos[iat, idim]
                                if (abs(tmp_var) > cutoff): 
                                    self.center_lk[itraj, index_lk, iat, idim] = self.mols[itraj].pos[iat, idim]
                                else:
                                    self.center_lk[itraj, index_lk, iat, idim] = center_new_lk[itraj, index_lk, iat, idim]
                            else: 
                                self.center_lk[itraj, index_lk, iat, idim] = center_old_lk[itraj, index_lk, iat, idim]

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

    def write_md_output(self, itrajectory, unixmd_dir, calc_coupling, istep):
        """ Write output files

            :param integer itrajectory: Index for trajectories
            :param string unixmd_dir: PyUNIxMD directory
            :param boolean calc_coupling: Check whether the dynamics includes coupling calculation
            :param integer istep: Current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, calc_coupling, istep)

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
                "".join([f'{self.count_ntrajs[itrajectory, iat, idim]:15.8f}' for idim in range(self.ndim)]) for iat in range(self.nat_qm)])
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

            :param object qm: QM object containing on-the-fly calculation information
            :param object mm: MM object containing MM calculation information
            :param string restart: Option for controlling dynamics restarting
        """
        # Print initial information about molecule, qm, mm and thermostat
        super().print_init(qm, mm, False, restart)

        # Print CTMQC info.
        ct_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"CTMQC Information":>43s}
        {"-" * 68}
          rho_threshold            = {self.rho_threshold:>16f}
          dist_parameter           = {self.dist_parameter:>16f}
          min_sigma                = {self.min_sigma:>16f}
        """)

        if (self.const_dist_cutoff != None):
            ct_info += f"  const_dist_cutoff        = {self.const_dist_cutoff:>16f}\n"
        else:
            ct_info += f"  const_dist_cutoff        = {str(None):>16s}\n"

        if (self.const_center_cutoff != None):
            ct_info += f"  const_center_cutoff      = {self.const_center_cutoff:>16f}\n"
        else:
            ct_info += f"  const_center_cutoff      = {str(None):>16s}\n"
        print (ct_info, flush=True)

        # Print istate
        istate_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Initial State Information":>43s}
        {"-" * 68}
        """)
        istate_info += f"  istates (1:{self.ntrajs})             =\n"
        nlines = self.ntrajs // 6
        if (self.ntrajs % 6 == 0):
            nlines -= 1

        for iline in range(nlines + 1):
            iline1 = iline * 6
            iline2 = (iline + 1) * 6
            if (iline2 > self.ntrajs):
                iline2 = self.ntrajs
            istate_info += f"  {iline1 + 1:>4d}:{iline2:<4d};"
            istate_info += "".join([f'{str(istate):7s}' for istate in self.istates[iline1:iline2]])
            istate_info += "\n"
        print (istate_info, flush=True)

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
        """ Routine to print each trajectory information at each step about dynamics

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
        
        # Print event in CTMQC
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event["DECO"] = []
