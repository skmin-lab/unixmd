from __future__ import division
from build.el_propagator_ct import el_run
from mqc.mqc import MQC
from misc import eps, au_to_K, au_to_A, call_name, typewriter, gaussian1d
import os, shutil, textwrap
import numpy as np
import pickle

class CT(MQC):
    """ Class for coupled-trajectory mixed quantum-classical (CTMQC) dynamics
    """
    def __init__(self, molecules, thermostat=None, istates=None, dt=0.5, nsteps=1000, nesteps=20, \
        rho_threshold=0.01, elec_object="coefficient", propagator="rk4", dist_threshold=0.25, dist_cutoff=0.5, \
        dist_parameter=10., sigma=0.3, l_print_dm=True, l_adj_nac=True, init_coefs=None, unit_dt="fs", \
        out_freq=1, verbosity=2):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize input values
        self.mols = molecules
        self.ntrajs = len(self.mols)
        self.digit = len(str(self.ntrajs))

        self.nst = self.mols[0].nst
        self.nat = self.mols[0].nat
        self.ndim = self.mols[0].ndim

        if (istates == None):
            raise ValueError (f"( {self.md_type}.{call_name()} ) istates should be required! {istates}")

        if (isinstance(istates, list)):
            if (len(istates) != self.ntrajs):
                raise ValueError (f"( {self.md_type}.{call_name()} ) The length of istates should be same to total number of trajectories! {istates}")
            else:
                if (max(istates) >= self.nst):
                    raise ValueError (f"( {self.md_type}.{call_name()} ) Index for initial state must be smaller than number of states! {max(istates)}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) The type of istates should be list! {istates}")

        if (init_coefs == None):
            init_coefs = [None] * self.ntrajs
        else:
            if (len(init_coefs) != self.ntrajs):
                raise ValueError (f"( {self.md_type}.{call_name()} ) The length of init_coefs should be same to total number of trajectories! {len(init_coefs)}")

        # Initialize input values and coefficient for first trajectory
        super().__init__(self.mols[0], thermostat, istates[0], dt, nsteps, nesteps, \
            elec_object, propagator, l_print_dm, l_adj_nac, init_coefs[0], unit_dt, out_freq, verbosity)

        if (self.elec_object != "coefficient"):
            raise ValueError (f"( {self.md_type}.{call_name()} ) coefficient propagation is only valid! {self.elec_object}")

        # Initialize coefficient for other trajectories
        for itraj in range(1, self.ntrajs):
            self.mols[itraj].get_coefficient(init_coefs[itraj], istates[itraj])

        # Initialize variables for CTMQC
        self.phase = np.zeros((self.ntrajs, self.nst, self.nat, self.ndim))
        self.nst_pair = int(self.nst * (self.nst - 1) / 2)
        self.qmom = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.ndim))
        self.K_lk = np.zeros((self.ntrajs, self.nst, self.nst))

        self.upper_th = 1. - rho_threshold
        self.lower_th = rho_threshold
        
        self.dist_threshold = dist_threshold
        self.dist_cutoff = dist_cutoff

        self.dist_parameter = dist_parameter
        self.sigma = sigma

        self.dotpopd = np.zeros(self.nst)

    def run(self, qm, mm=None, output_dir="./", l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        # Initialize UNI-xMD
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

                self.write_md_output(itraj, unixmd_dirs[itraj], self.istep)
                
                self.print_traj(self.istep, itraj)

            self.calculate_qmom(self.istep, unixmd_dirs)

            self.print_step(self.istep)

        #TODO: restart
        else: 
            raise ValueError (f"( {self.md_type}.{call_name()} ) restart is not valid in CTMQC ! {restart}")

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

                if ((istep + 1) % self.out_freq == 0):
                    self.write_md_output(itraj, unixmd_dirs[itraj], istep)
                    self.print_traj(istep, itraj)
                if (istep == self.nsteps - 1):
                    self.write_final_xyz(unixmd_dirs[itraj], istep)

                #TODO: restart
                #self.fstep = istep
                #restart_file = os.path.join(base_dir, "RESTART.bin")
                #with open(restart_file, 'wb') as f:
                #    pickle.dump({'qm':qm, 'md':self}, f)

            self.calculate_qmom(istep, unixmd_dirs)

            self.print_step(istep)

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
        self.rforce = np.zeros((self.nat, self.ndim))
        
        # Derivatives of energy
        for ist, istate in enumerate(self.mols[itrajectory].states):
            self.rforce += istate.force * self.mol.rho.real[ist, ist]

        # Non-adiabatic forces 
        for ist in range(self.nst):
            for jst in range(ist + 1, self.nst):
                self.rforce += 2. * self.mol.nac[ist, jst] * self.mol.rho.real[ist, jst] \
                    * (self.mol.states[ist].energy - self.mol.states[jst].energy)

        # CT forces
        ctforce = np.zeros((self.nat, self.ndim))
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
                self.phase[itrajectory, ist] = np.zeros((self.nat, self.ndim))
            else:
                self.phase[itrajectory, ist] += self.mol.states[ist].force * self.dt

    def calculate_qmom(self, istep, dirs):
        """ Routine to calculate quantum momentum
            
            :param integer istep: Current MD step
            :param string dirs: Output directory name
        """
        #TODO: default value of parameter
        #TODO: test
        smooth_factor = 2.
        # dist_cutoff = 0.5 in sh_tddft_utils.mod.F90 
        # threshold = 0.25
        # It use when nntraj is counted.
        M_parameter = 10.
        sigma = np.ones((self.nat, self.ndim)) * 0.3
        # The value of M_parameter * sigma is used when determine quantum momentum center.

        # _lk means state_pair dependency.
        # i and j are trajectory index.
        # -------------------------------------------------------------------
        # 1. Calculate variances for each trajectory
        # TODO: method to calculate sigma
        sigma_lk = np.ones((self.ntrajs, self.nst_pair, self.nat, self.ndim)) # TODO: state-pair
        for itraj in range(self.ntrajs):
            # Variable to determine how many trajecories are in cutoff.
            nntraj = np.zeros((self.nat)) 

            R2_tmp = np.zeros((self.nat, self.ndim)) # Temporary variable for R**2
            R_tmp = np.zeros((self.nat, self.ndim))  # Temporary variable for R

            for jtraj in range(self.ntrajs):
                pos_diff = self.mols[jtraj].pos - self.mols[itraj].pos # Dimension = (self.nat, self.ndim)
                pos_diff2 = np.sum(pos_diff * pos_diff, axis=1) # Dimension = (self.nat)

                for iat in range(self.nat):
                    distance = np.sqrt(pos_diff2[iat]) # Distance between i-th atom in itraj and jtraj
                    if (distance <= smooth_factor):
                        R_tmp[iat] += self.mols[jtraj].pos[iat] # Dimension = (self.nat, self.ndim)
                        R2_tmp[iat] += self.mols[jtraj].pos[iat] * self.mols[jtraj].pos[iat] # Dimension = (self.nat, self.ndim)
                        nntraj[iat] += 1

            tmp= f'{istep+1:8d}' + \
                "".join([f'{nntraj[iat]:15.8f}' for iat in range(self.nat)])
            typewriter(tmp, dirs[itraj], f"NNTRAJ", "a")

            for iat in range(self.nat):
                avg_R = R_tmp[iat] / nntraj[iat]
                avg_R2 = R2_tmp[iat] / nntraj[iat]
                for isp in range(self.ndim):
                    sigma_lk[itraj, 0, iat, isp] = np.sqrt((avg_R2[isp] - avg_R[isp] ** 2)) / np.sqrt(np.sqrt(nntraj[iat]))
                    # / np.sqrt(np.sqrt(nntraj)) is artifact to modulate sigma.
                    if (sigma_lk[itraj, 0, iat, isp] <= 1.0E-8):
                        sigma_lk[itraj, 0, iat, isp] = smooth_factor

            tmp= f'{istep+1:8d}{sigma_lk[itraj, 0, 0, 0]:15.8f}'
            typewriter(tmp, dirs[itraj], f"SIGMA", "a")

        # 2. Calculate slope
        # (2-1) Calculate w_ij
        # g_i means nuclear density at the position of i-th classical trajectory.
        # prod_g_i is to multiply gaussians with respect to atoms and spaces.
        g_i = np.zeros((self.ntrajs)) 
        prod_g_i = np.ones((self.ntrajs, self.ntrajs))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat):
                    for isp in range(self.ndim):
                        # gaussian1d(x, pre-factor, sigma, mean)
                        # gaussian1d(R^{itraj}, 1.0, sigma^{jtraj}, R^{jtraj})
                        prod_g_i[itraj, jtraj] *= gaussian1d(self.mols[itraj].pos[iat, isp], 1., \
                            sigma_lk[jtraj, 0, iat, isp], self.mols[jtraj].pos[iat, isp])
                g_i[itraj] += prod_g_i[itraj, jtraj]

        # w_ij is defined as W_IJ in SI of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        w_ij = np.zeros((self.ntrajs, self.ntrajs, self.nat, self.ndim))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat):
                    for isp in range(self.ndim):
                        w_ij[itraj, jtraj, iat, isp] = prod_g_i[itraj, jtraj] /\
                        (2. * sigma_lk[jtraj, 0, iat, isp] ** 2 * g_i[itraj])

        # (2-2) Calculate slope_i
        # the slope is calculated as a sum over j of w_ij
        slope_i = np.zeros((self.ntrajs, self.nat, self.ndim))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                slope_i[itraj] -= w_ij[itraj, jtraj]

            tmp= f'{istep+1:8d}{slope_i[itraj, 0, 0]:15.8f}'
            typewriter(tmp, dirs[itraj], f"SLOPE", "a")

        # 3. Calculate the center of quantum momentum
        rho = np.zeros((self.ntrajs, self.nst))
        for itraj in range(self.ntrajs):
            for ist in range(self.nst):
                rho[itraj, ist] = self.mols[itraj].rho[ist, ist].real

        # (3-1) Compute denominator
        deno_lk = np.zeros((self.nst_pair, self.nat, self.ndim)) # denominator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):  
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.ndim):
                            deno_lk[index_lk, iat, isp] += rho[itraj, ist] * rho[itraj, jst] * \
                                (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * slope_i[itraj, iat, isp]

        # (3-2) Compute numerator
        ratio_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.ndim)) # numerator / denominator
        numer_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.ndim)) # numerator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.ndim):
                            numer_lk[itraj, index_lk, iat, isp] = rho[itraj, ist] * rho[itraj, jst] * self.mols[itraj].pos[iat, isp] * \
                                (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * slope_i[itraj, iat, isp]
                            if (abs(deno_lk[index_lk, iat, isp]) <= 1.0E-08):
                                ratio_lk[itraj, index_lk, iat, isp] = 0.
                            else:
                                ratio_lk[itraj, index_lk, iat, isp] = numer_lk[itraj, index_lk, iat, isp] / \
                                    deno_lk[index_lk, iat, isp]

        # Center of quantum momentum is calculated by Eq.(S28) of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        center_old_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.ndim))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.ndim):
                            for jtraj in range(self.ntrajs):
                                center_old_lk[itraj, index_lk, iat, isp] += ratio_lk[jtraj, index_lk, iat, isp]
                            if ((abs(slope_i[itraj, iat, isp]) <= 1.0E-08) or (center_old_lk[itraj, index_lk, iat, isp] == 0.)):
                                center_old_lk[itraj, index_lk, iat, isp] = self.mols[itraj].pos[iat, isp]

        # Center of quantum momentum is calculated by Eq.(S21) of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        center_new_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.ndim))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.ndim):
                            if (abs(slope_i[itraj, iat, isp]) <= 1.0E-08):
                                center_new_lk[itraj, index_lk, iat, isp] = self.mols[itraj].pos[iat, isp]
                            else:
                                for jtraj in range(self.ntrajs):
                                    center_new_lk[itraj, index_lk, iat, isp] += self.mols[jtraj].pos[iat, isp] * prod_g_i[itraj, jtraj] /\
                                        (2. * sigma_lk[jtraj, 0, iat, isp] ** 2 * g_i[itraj] * (- slope_i[itraj, iat, isp]))

        # (3-3) Determine qauntum momentum center
        center_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.ndim)) # Finally, qmom_center
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.ndim):
                            # tmp_var is deviation between position of classical trajectory and quantum momentum center.
                            tmp_var = center_old_lk[itraj, index_lk, iat, isp] - self.mols[itraj].pos[iat, isp]
                            if (abs(tmp_var) > M_parameter * sigma[iat, isp]): 
                                tmp_var = center_new_lk[itraj, index_lk, iat, isp] - self.mols[itraj].pos[iat, isp]
                                if (abs(tmp_var) > M_parameter * sigma[iat, isp]): 
                                    center_lk[itraj, index_lk, iat, isp] = self.mols[itraj].pos[iat, isp]
                                else:
                                    center_lk[itraj, index_lk, iat, isp] = center_new_lk[itraj, index_lk, iat, isp]
                            else: 
                                center_lk[itraj, index_lk, iat, isp] = center_old_lk[itraj, index_lk, iat, isp]

            tmp= f'{istep+1:8d}{center_lk[itraj, 0, 0, 0]:15.8f}'
            typewriter(tmp, dirs[itraj], f"CENTER", "a")

        # 4. Compute quantum momentum
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    self.qmom[itraj, index_lk] = slope_i[itraj] * (self.mols[itraj].pos - center_lk[itraj, index_lk])

        # 5. Calculate 2 * Qmom * phase / mass
        self.K_lk = np.zeros((self.ntrajs, self.nst, self.nst))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    self.K_lk[itraj, ist, jst] += 2. * np.sum(1. / self.mol.mass[0:self.mol.nat_qm] * \
                        np.sum(self.qmom[itraj, index_lk] * self.phase[itraj, ist], axis = 1))
                    self.K_lk[itraj, jst, ist] += 2. * np.sum(1. / self.mol.mass[0:self.mol.nat_qm] * \
                        np.sum(self.qmom[itraj, index_lk] * self.phase[itraj, jst], axis = 1))

    def write_md_output(self, itrajectory, unixmd_dir, istep):
        """ Write output files

            :param string unixmd_dir: Unixmd directory
            :param integer istep: Current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, istep)

        # Write decoherence information
        self.write_deco(itrajectory, unixmd_dir, istep)

    def write_deco(self, itrajectory, unixmd_dir, istep):
        """ Write CT-based decoherence information

            :param string unixmd_dir: Unixmd directory
            :param integer istep: Current MD step
        """
        # Write time-derivative density matrix elements in DOTPOTD
        #tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopd])
        #typewriter(tmp, unixmd_dir, "DOTPOPD", "a")

        # Write auxiliary trajectories
        if (self.verbosity >= 2):
            # Write quantum momenta
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    tmp = f'{self.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Momentum (au)' + \
                        "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                        "".join([f'{self.qmom[itrajectory, index_lk, iat, isp]:15.8f}' for isp in range(self.ndim)]) for iat in range(self.nat)])
                    typewriter(tmp, unixmd_dir, f"QMOM_{ist}_{jst}", "a")

            for ist in range(self.nst):
                for jst in range(self.nst):
                    if (ist != jst):
                        tmp = f'{istep + 1:9d}{self.K_lk[itrajectory, ist, jst]:15.8f}'
                        typewriter(tmp, unixmd_dir, f"K_lk_{ist}_{jst}", "a")

            # Write auxiliary variables
            for ist in range(self.mol.nst):
                # Write auxiliary phase
                tmp = f'{self.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Momentum (au)' + \
                    "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                    "".join([f'{self.phase[itrajectory, ist, iat, isp]:15.8f}' for isp in range(self.ndim)]) for iat in range(self.nat)])
                typewriter(tmp, unixmd_dir, f"PHASE_{ist}", "a")

    def print_init(self, qm, mm, restart):
        """ Routine to print the initial information of dynamics

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param string restart: option for controlling dynamics restarting
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

        # Print INIT for averaged quantity at each step
        DEBUG1 = f" #INFO_AVG{'STEP':>9s}"
        for ist in range(self.nst):
            DEBUG1 += f"{'BOPOP_':>13s}{ist}"

        for ist in range(self.nst):
            for jst in range(ist + 1, self.nst):
                DEBUG1 += f"{'BOCOH_':>13s}{ist}_{jst}"

        dynamics_step_info += "\n" + DEBUG1

        print (dynamics_step_info, flush=True)

    def print_traj(self, istep, itrajectory):
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

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer istep: Current MD step
        """
        rho = np.zeros((self.nst, self.nst))
        for itraj in range(self.ntrajs):
            for ist in range(self.nst):
                for jst in range(ist, self.nst):
                    if (ist == jst):
                        rho[ist, jst] += self.mols[itraj].rho[ist, jst].real
                    else:
                        rho[ist, jst] += self.mols[itraj].rho[ist, ist].real * self.mols[itraj].rho[jst, jst].real
        rho /= self.ntrajs

        DEBUG1 = f" INFO_AVG{istep + 1:9d}" + "".join([f'{rho[ist, ist]:15.8f}' for ist in range(self.nst)])
        DEBUG1 += "".join([f'{rho[ist, jst]:15.8f}' for ist in range(self.nst) for jst in range(ist + 1, self.nst)])
        print(DEBUG1, flush=True)
