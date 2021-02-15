from __future__ import division
from build.el_propagator_ct import el_run
from mqc.mqc import MQC
from misc import eps, au_to_K, au_to_A, call_name, typewriter, elapsed_time
import os, shutil, textwrap
import numpy as np
import pickle

class CT(MQC):
    """ Class for coupled-trajectory mixed quantum-classical dynamics
    """
    def __init__(self, molecules, thermostat=None, istates=None, dt=0.5, nsteps=1000, nesteps=20, \
        threshold=0.01, propagation="density", solver="rk4", l_pop_print=False, l_adjnac=True, \
        coefficients=None, unit_dt="fs", out_freq=1, verbosity=2):
        # Initialize input values
        self.mols = molecules
        self.ntrajs = len(self.mols)
        self.istates = istates
        if ((self.istates != None) and (self.ntrajs != len(self.istates))):
            raise ValueError("Error: istates!")

        if (coefficients == None):
            coefficients = [None] * self.ntrajs
        else:
            if (self.ntrajs != len(coefficients)):
                raise ValueError("Error: coefficients!")

        # Initialize input values and coefficient for first trajectory
        super().__init__(self.mols[0], thermostat, istates[0], dt, nsteps, nesteps, \
            propagation, solver, l_pop_print, l_adjnac, coefficients[0], unit_dt, out_freq, verbosity)

        # Initialize coefficient for other trajectories
        for itraj in range(1, self.ntrajs):
            self.mols[itraj].get_coefficient(coefficients[itraj], self.istates[itraj])

        # Initialize variables for CTMQC
        self.nst = self.mols[0].nst
        self.nat = self.mols[0].nat
        self.nsp = self.mols[0].nsp
        self.rho_threshold = threshold

        self.phase = np.zeros((self.ntrajs, self.nst, self.nat, self.nsp))
        self.nst_pair = int(self.nst * (self.nst - 1) / 2)
        self.qmom = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp))
        self.K_lk = np.zeros((self.ntrajs, self.nst, self.nst))

        self.dotpopd = np.zeros(self.mol.nst)

    @elapsed_time
    def run(self, qm, mm=None, input_dir="./", save_qm_log=False, save_mm_log=False, save_scr=True, restart=None):
        # Initialize UNI-xMD
        # TODO
        base_dir, unixmd_dir, qm_log_dir, mm_log_dir =\
             self.run_init(qm, mm, input_dir, save_qm_log, save_mm_log, save_scr, restart)
        unixmd_dirs = [''] * self.ntrajs
        qm_log_dirs= [''] * self.ntrajs
        for itraj in range(self.ntrajs):
            unixmd_dirs[itraj] = f'{unixmd_dir}_{itraj+1:03d}'
            if (os.path.exists(unixmd_dirs[itraj])):
                shutil.move(unixmd_dirs[itraj], unixmd_dirs[itraj] + "_old_" + str(os.getpid()))
            os.makedirs(unixmd_dirs[itraj])

        bo_list = [ist for ist in range(self.mol.nst)]
        qm.calc_coupling = True
        self.print_init()

        if (restart == None):
            # Calculate initial input geometry for all trajectories at t = 0.0 s
            self.istep = -1
            for itraj in range(self.ntrajs):
                self.mol = self.mols[itraj]

                self.mol.reset_bo(qm.calc_coupling)
                qm.get_data(self.mol, base_dir, bo_list, self.dt, self.istep, calc_force_only=False)
                # TODO: QM/MM
                self.mol.get_nacme()

                self.update_energy()

                self.get_phase(itraj)

                self.mols[itraj] = self.mol

                self.write_md_output(itraj, unixmd_dirs[itraj], self.istep)
                #self.print_step(self.istep)
            self.calculate_qmom(self.istep, unixmd_dirs)

            self.print_step(self.istep)

        else: 
            raise ValueError ("restart option is invalid in CTMQC yet.")

        self.istep += 1

        # Main MD loop
        for istep in range(self.istep, self.nsteps):
            for itraj in range(self.ntrajs):
                self.mol = self.mols[itraj]
                
                self.calculate_force(itraj)
                self.cl_update_position()

                self.mol.backup_bo()
                self.mol.reset_bo(qm.calc_coupling)
                qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=False)
                # TODO: QM/MM

                self.mol.adjust_nac()

                self.calculate_force(itraj)
                self.cl_update_velocity()

                self.mol.get_nacme()
                
                el_run(self, itraj)

                # TODO: thermostat
                #if (self.thermo != None):
                #    self.thermo.run(self)

                self.update_energy()

                self.get_phase(itraj)

                if ((istep + 1) % self.out_freq == 0):
                    self.write_md_output(itraj, unixmd_dirs[itraj], istep)
#                    self.print_step(self.mols[itraj], istep)
                if (istep == self.nsteps - 1):
                    self.write_final_xyz(unixmd_dirs[itraj], istep)

                self.mols[itraj] = self.mol
                # TODO: restart
                #self.fstep = istep
                #restart_file = os.path.join(base_dir, "RESTART.bin")
                #with open(restart_file, 'wb') as f:
                #    pickle.dump({'qm':qm, 'md':self}, f)

            self.calculate_qmom(istep, unixmd_dirs)

            self.print_step(istep)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

    def calculate_force(self, itrajectory):
        """ Routine to calculate force
        """
        self.rforce = np.zeros((self.nat, self.nsp))
        # 
        for ist, istate in enumerate(self.mols[itrajectory].states):
            self.rforce += istate.force * self.mol.rho.real[ist, ist]

        # Non-adiabatic forces from Ehrenfest force 
        for ist in range(self.nst):
            for jst in range(ist + 1, self.nst):
                self.rforce += 2. * self.mol.nac[ist, jst] * self.mol.rho.real[ist, jst] \
                    * (self.mol.states[ist].energy - self.mol.states[jst].energy)

        # CT force
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
        """
        for ist in range(self.nst):
            rho_ii = self.mol.rho[ist, ist].real
            if ((rho_ii > self.rho_threshold) and (rho_ii < (1. - self.rho_threshold))):
                self.phase[itrajectory, ist] += self.mol.states[ist].force * self.dt
            else:
                self.phase[itrajectory, ist] = np.zeros((self.nat, self.nsp))

    def calculate_qmom(self, istep, dirs):
        """ Routine to calculate quantum momentum
        """
        #TODO
        M_parameter = 10.
        sigma = np.ones((self.nat, self.nsp)) * 0.3

        # _lk means state_pair dependency.
        # i and j are trajectory index.
        # -------------------------------------------------------------------
        # 1. Calculate sigma for each trajectory
        smooth_factor = 2. #TODO: Q. How to determine cutoff ?
        sigma_lk = np.ones((self.ntrajs, self.nst_pair, self.nat, self.nsp)) * sigma_x # TODO: state-pair
        for itraj in range(self.ntrajs):
            nntraj = np.zeros((self.nat)) # Count  
            R2_tmp = np.zeros((self.nat, self.nsp)) # temporary variable for R**2
            R_tmp = np.zeros((self.nat, self.nsp))  # temporary variable for R

            for jtraj in range(self.ntrajs):
                pos_diff = self.mols[jtraj].pos - self.mols[itraj].pos # dimension = (self.nat, self.nsp)
                pos_diff2 = np.sum(pos_diff * pos_diff, axis=1) # dimension = (self.nat)

                for iat in range(self.nat):
                    distance = np.sqrt(pos_diff2[iat]) # distance between i-th atom in itraj and jtraj
                    if (distance <= smooth_factor):
                        R_tmp += self.mols[jtraj].pos # dimension = (self.nat, self.nsp)
                        R2_tmp += self.mols[jtraj].pos * self.mols[jtraj].pos # dimension = (self.nat, self.nsp)
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
        # g_i means single gaussian for i-th trajectory.
        # prod_g_i is to multiply gaussians with respect to atoms.
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

        w_ij = np.zeros((self.ntrajs, self.ntrajs, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                for iat in range(self.nat):
                    for isp in range(self.nsp):
                        w_ij[itraj, jtraj, iat, isp] = prod_g_i[itraj, jtraj] /\
                        (2. * sigma_lk[jtraj, 0, iat, isp] ** 2 * g_i[itraj])

        # (2-2) Calculate slope_i
        slope_i = np.zeros((self.ntrajs, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            for jtraj in range(self.ntrajs):
                slope_i[itraj] -= w_ij[itraj, jtraj] # TODO: minus?

            tmp= f'{istep+1:8d}{slope_i[itraj, 0, 0]:15.8f}'
            typewriter(tmp, dirs[itraj], f"SLOPE", "a")

        # 3. Calculate the center of quantum momentum
        rho = np.zeros((self.ntrajs, self.nst))
        for itraj in range(self.ntrajs):
            for ist in range(self.nst):
                rho[itraj, ist] = self.mols[itraj].rho[ist, ist].real

        # (3-1) Compute denominator
        deno_lk = np.zeros((self.nst_pair, self.nat, self.nsp)) # denominator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):  
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            deno_lk[index_lk, iat, isp] += rho[itraj, ist] * rho[itraj, jst] * (rho[itraj, ist] + rho[itraj, jst]) *\
                                (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * slope_i[itraj, iat, isp]

        # (3-2) Compute numerator
        ratio = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp))
        numer_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp)) # numerator
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            numer_lk[itraj, index_lk, iat, isp] = rho[itraj, ist] * rho[itraj, jst] * (rho[itraj, ist] + rho[itraj, jst]) * \
                                self.mols[itraj].pos[iat, isp] * (self.phase[itraj, ist, iat, isp] - self.phase[itraj, jst, iat, isp]) * \
                                slope_i[itraj, iat, isp]
                            if (abs(deno_lk[index_lk, iat, isp]) <= 1.0e-08):
                                ratio[itraj, index_lk, iat, isp] = 0.
                            else:
                                ratio[itraj, index_lk, iat, isp] = numer_lk[itraj, index_lk, iat, isp] / \
                                    deno_lk[index_lk, iat, isp]

        # Center of quantum momentum is calculated by Eq.(S28) of paper Min et al.
        center_old_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            for jtraj in range(self.ntrajs):
                                center_old_lk[itraj, index_lk, iat, isp] += ratio[jtraj, index_lk, iat, isp]
                            if ((abs(slope_i[itraj, iat, isp]) <= 1.0e-08) or (center_old_lk[itraj, index_lk, iat, isp] == 0.)):
                                center_old_lk[itraj, index_lk, iat, isp] = self.mols[itraj].pos[iat, isp]

        # Center of quantum momentum is calculated by Eq.(S21) of paper Min et al.
        center_new_lk = np.zeros((self.ntrajs, self.nst_pair, self.nat, self.nsp))
        for itraj in range(self.ntrajs):
            index_lk = -1
            for ist in range(self.nst):
                for jst in range(ist + 1, self.nst):
                    index_lk += 1
                    for iat in range(self.nat):
                        for isp in range(self.nsp):
                            if (abs(slope_i[itraj, iat, isp]) <= 1.0e-08):
                                center_new_lk[itraj, index_lk, iat, isp] = self.mols[itraj].pos[iat, isp]
                            else:
                                for jtraj in range(self.ntrajs):
                                    center_new_lk[itraj, index_lk, iat, isp] += self.mols[jtraj].pos[iat, isp] * prod_g_i[itraj, jtraj] /\
                                        (2. * sigma_lk[jtraj, 0, iat, isp] ** 2 * g_i[itraj] * (- slope_i[itraj, iat, isp])) #TODO:

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

            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, istep)

        # Write decoherence information
        self.write_deco(itrajectory, unixmd_dir, istep)

    def write_deco(self, itrajectory, unixmd_dir, istep):
        """ Write CT-based decoherence information

            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
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

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics
        """
        rho = np.zeros((self.nst, self.nst), dtype=np.complex_)
        for itraj in range(self.ntrajs):
            for ist in range(self.nst):
                for jst in range(ist, self.nst):
                    rho[ist, jst] += self.mols[itraj].rho[ist, jst]
        rho /= self.ntrajs

        print(f'RHO{istep+1:8d}{rho[0, 0].real:15.8f}{rho[1, 1].real:15.8f}', flush=True)

#TODO: move to misc.py?
def gaussian1d(x, const, sigma, x0):
    if (sigma < 0.0):
        return -1
    else:
        res = const / (sigma * np.sqrt(2. * np.pi)) * np.exp(- (x - x0) ** 2 / (2. * sigma ** 2))
        return res 
