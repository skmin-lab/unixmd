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
        threshold=0.01, propagation="coefficient", solver="rk4", l_pop_print=False, l_adjnac=True, \
        coefficients=None, unit_dt="fs", out_freq=1, verbosity=2):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize input values
        self.mols = molecules
        self.ntrajs = len(self.mols)
        self.digit = len(str(self.ntrajs))

        self.nst = self.mols[0].nst
        self.nat = self.mols[0].nat
        self.nsp = self.mols[0].nsp

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

        if (coefficients == None):
            coefficients = [None] * self.ntrajs
        else:
            if (self.ntrajs != len(coefficients)):
                raise ValueError (f"( {self.md_type}.{call_name()} ) The length of coefficients should be same to total number of trajectories! {coefficients}")

        # Initialize input values and coefficient for first trajectory
        super().__init__(self.mols[0], thermostat, istates[0], dt, nsteps, nesteps, \
            propagation, solver, l_pop_print, l_adjnac, coefficients[0], unit_dt, out_freq, verbosity)

        if (self.propagation != "coefficient"):
            raise ValueError (f"( {self.md_type}.{call_name()} ) coefficient propagation is only valid! {self.propagation}")

        # Initialize coefficient for other trajectories
        for itraj in range(1, self.ntrajs):
            self.mols[itraj].get_coefficient(coefficients[itraj], istates[itraj])

        # Initialize variables for CTMQC
        self.phase = np.zeros((self.ntrajs, self.nst, self.nat, self.nsp))
        self.nst_pair = int(self.nst * (self.nst - 1) / 2)
        self.qmom = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp))
        self.K_lk = np.zeros((self.ntrajs, self.nst, self.nst))

        self.upper_th = 1. - threshold
        self.lower_th = threshold

        self.dotpopd = np.zeros(self.nst)

    def run(self, qm, mm=None, input_dir="./", save_qm_log=False, save_mm_log=False, save_scr=True, restart=None):
        # Initialize UNI-xMD
        base_dirs, unixmd_dirs, qm_log_dirs, mm_log_dirs =\
             self.run_init(qm, mm, input_dir, save_qm_log, save_mm_log, save_scr, restart)

        bo_list = [ist for ist in range(self.nst)]
        qm.calc_coupling = True

        # TODO: output control
        #self.print_init()

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

            self.calculate_qmom(self.istep, unixmd_dirs)

            # TODO: output control
            #self.print_step(self.istep)

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
                    #self.print_traj(istep)
                if (istep == self.nsteps - 1):
                    self.write_final_xyz(unixmd_dirs[itraj], istep)

                #TODO: restart
                #self.fstep = istep
                #restart_file = os.path.join(base_dir, "RESTART.bin")
                #with open(restart_file, 'wb') as f:
                #    pickle.dump({'qm':qm, 'md':self}, f)

            self.calculate_qmom(istep, unixmd_dirs)

            # TODO: output control
            #self.print_step(istep)

        # Delete scratch directory
        if (not save_scr):
            for itraj in range(self.ntrajs):
                tmp_dir = os.path.join(unixmd_dirs, "scr_qm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def calculate_force(self, itrajectory):
        """ Routine to calculate force
            
            :param integer itrajectory: Index for trajectories
        """
        self.rforce = np.zeros((self.nat, self.nsp))
        
        # Derivatives of energy
        for ist, istate in enumerate(self.mols[itrajectory].states):
            self.rforce += istate.force * self.mol.rho.real[ist, ist]

        # Non-adiabatic forces 
        for ist in range(self.nst):
            for jst in range(ist + 1, self.nst):
                self.rforce += 2. * self.mol.nac[ist, jst] * self.mol.rho.real[ist, jst] \
                    * (self.mol.states[ist].energy - self.mol.states[jst].energy)

        # CT forces
        ctforce = np.zeros((self.nat, self.nsp))
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
                self.phase[itrajectory, ist] = np.zeros((self.nat, self.nsp))
            else:
                self.phase[itrajectory, ist] += self.mol.states[ist].force * self.dt

    def calculate_qmom(self, istep, dirs):
        """ Routine to calculate quantum momentum
            
            :param integer istep: Current MD step
            :param string dirs: Output directory name
        """
        #TODO: parameter
        smooth_factor = 2.
        M_parameter = 10.
        sigma = np.ones((self.nat, self.nsp)) * 0.3

        # _lk means state_pair dependency.
        # i and j are trajectory index.
        # -------------------------------------------------------------------
        # 1. Calculate sigma for each trajectory
        sigma_lk = np.ones((self.ntrajs, self.nst_pair, self.nat, self.nsp)) # TODO: state-pair
        for itraj in range(self.ntrajs):
            # Variable to determine how many trajecories are in cutoff.
            nntraj = np.zeros((self.nat)) 

            R2_tmp = np.zeros((self.nat, self.nsp)) # Temporary variable for R**2
            R_tmp = np.zeros((self.nat, self.nsp))  # Temporary variable for R

            for jtraj in range(self.ntrajs):
                pos_diff = self.mols[jtraj].pos - self.mols[itraj].pos # Dimension = (self.nat, self.nsp)
                pos_diff2 = np.sum(pos_diff * pos_diff, axis=1) # Dimension = (self.nat)

                for iat in range(self.nat):
                    distance = np.sqrt(pos_diff2[iat]) # Distance between i-th atom in itraj and jtraj
                    if (distance <= smooth_factor):
                        R_tmp += self.mols[jtraj].pos # Dimension = (self.nat, self.nsp)
                        R2_tmp += self.mols[jtraj].pos * self.mols[jtraj].pos # Dimension = (self.nat, self.nsp)
                        nntraj[iat] += 1

            tmp= f'{istep+1:8d}' + \
                "".join([f'{nntraj[iat]:15.8f}' for iat in range(self.nat)])
            typewriter(tmp, dirs[itraj], f"NNTRAJ", "a")

            for iat in range(self.nat):
                avg_R = R_tmp[iat] / nntraj[iat]
                avg_R2 = R2_tmp[iat] / nntraj[iat]
                for isp in range(self.nsp):
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
                    for isp in range(self.nsp):
                        # gaussian1d(x, pre-factor, sigma, mean)
                        # gaussian1d(R^{itraj}, 1.0, sigma^{jtraj}, R^{jtraj})
                        prod_g_i[itraj, jtraj] *= gaussian1d(self.mols[itraj].pos[iat, isp], 1., \
                            sigma_lk[jtraj, 0, iat, isp], self.mols[jtraj].pos[iat, isp])
                g_i[itraj] += prod_g_i[itraj, jtraj]

        # w_ij is defined as W_IJ in SI of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        w_ij = np.zeros((self.ntrajs, self.ntrajs, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat):
                    for isp in range(self.nsp):
                        w_ij[itraj, jtraj, iat, isp] = prod_g_i[itraj, jtraj] /\
                        (2. * sigma_lk[jtraj, 0, iat, isp] ** 2 * g_i[itraj])

        # (2-2) Calculate slope_i
        # the slope is calculated as a sum over j of w_ij
        slope_i = np.zeros((self.ntrajs, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                slope_i[itraj] -= w_ij[itraj, jtraj]

            tmp= f'{istep+1:8d}{slope_i[itraj, 0, 0]:15.8f}'
            typewriter(tmp, dirs[itraj], f"SLOPE", "a")

        # 3. Calculate the center of quantum momentum
        rho_i = np.zeros((self.ntrajs, self.nst))
        for itraj in range(self.ntrajs):
            for ist in range(self.nst):
                rho_i[itraj, ist] = self.mols[itraj].rho[ist, ist].real

        # (3-1) Compute denominator
        deno_lk = np.zeros((self.nst_pair, self.nat, self.nsp)) # denominator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):  
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            deno_lk[index_lk, iat, isp] += rho_i[itraj, ist] * rho_i[itraj, jst] * (rho_i[itraj, ist] + rho_i[itraj, jst]) * \
                                (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * slope_i[itraj, iat, isp]
                            
                            #deno_lk[index_lk, iat, isp] += rho[itraj, ist] * rho[itraj, jst] * \
                            #    (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * slope_i[itraj, iat, isp]

        # (3-2) Compute numerator
        ratio_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp)) # numerator / denominator
        numer_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp)) # numerator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            numer_lk[itraj, index_lk, iat, isp] = rho_i[itraj, ist] * rho_i[itraj, jst] * (rho_i[itraj, ist] + rho_i[itraj, jst]) * \
                                self.mols[itraj].pos[iat, isp] * (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * \
                                slope_i[itraj, iat, isp]

                            #numer_lk[itraj, index_lk, iat, isp] = rho[itraj, ist] * rho[itraj, jst] * self.mols[itraj].pos[iat, isp] * \
                            #    (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * slope_i[itraj, iat, isp]
                            if (abs(deno_lk[index_lk, iat, isp]) <= 1.0E-08):
                                ratio_lk[itraj, index_lk, iat, isp] = 0.
                            else:
                                ratio_lk[itraj, index_lk, iat, isp] = numer_lk[itraj, index_lk, iat, isp] / \
                                    deno_lk[index_lk, iat, isp]

        # Center of quantum momentum is calculated by Eq.(S28) of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        center_old_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            for jtraj in range(self.ntrajs):
                                center_old_lk[itraj, index_lk, iat, isp] += ratio_lk[jtraj, index_lk, iat, isp]
                            if ((abs(slope_i[itraj, iat, isp]) <= 1.0E-08) or (center_old_lk[itraj, index_lk, iat, isp] == 0.)):
                                center_old_lk[itraj, index_lk, iat, isp] = self.mols[itraj].pos[iat, isp]

        # Center of quantum momentum is calculated by Eq.(S21) of J. Phys. Chem. Lett., 2017, 8, 3048-3055.
        center_new_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            if (abs(slope_i[itraj, iat, isp]) <= 1.0E-08):
                                center_new_lk[itraj, index_lk, iat, isp] = self.mols[itraj].pos[iat, isp]
                            else:
                                for jtraj in range(self.ntrajs):
                                    center_new_lk[itraj, index_lk, iat, isp] += self.mols[jtraj].pos[iat, isp] * prod_g_i[itraj, jtraj] /\
                                        (2. * sigma_lk[jtraj, 0, iat, isp] ** 2 * g_i[itraj] * (- slope_i[itraj, iat, isp]))

        # (3-3) Determine qauntum momentum center
        center_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp)) # Finally, qmom_center
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
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
                        "".join([f'{self.qmom[itrajectory, index_lk, iat, isp]:15.8f}' for isp in range(self.nsp)]) for iat in range(self.nat)])
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
                    "".join([f'{self.phase[itrajectory, ist, iat, isp]:15.8f}' for isp in range(self.nsp)]) for iat in range(self.nat)])
                typewriter(tmp, unixmd_dir, f"PHASE_{ist}", "a")

    def print_init(self):
        """ Routine to print the initial information of dynamics
        """
        pass

    def print_traj(self):
        """ Routine to print each trajectory infomation at each step about dynamics
        """
        pass

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer istep: Current MD step
        """
        rho = np.zeros((self.nst, self.nst), dtype=np.complex_)
        for itraj in range(self.ntrajs):
            for ist in range(self.nst):
                for jst in range(ist, self.nst):
                    rho[ist, jst] += self.mols[itraj].rho[ist, jst]
        rho /= self.ntrajs

        print(f'RHO{istep+1:8d}{rho[0, 0].real:15.8f}{rho[1, 1].real:15.8f}', flush=True)
