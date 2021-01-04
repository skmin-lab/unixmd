from __future__ import division
from build.el_propagator import el_run
from mqc.mqc import MQC
from misc import eps, au_to_K, call_name, typewriter
import random, os, shutil, textwrap
import numpy as np
import pickle

class CT(MQC):
    """ Class for coupled-trajectory mixed quantum-classical dynamics
    """
    def __init__(self, molecules, thermostat=None, istates=None, dt=0.5, nsteps=1000, nesteps=20, \
        propagation="density", solver="rk4", l_pop_print=False, l_adjnac=True, \
        coefficients=None, unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        self.mols = molecules
        self.ntrajs = len(self.mols)
        self.istates = istates
        if ((self.istates != None) and (self.ntrajs != len(self.istates))):
            raise ValueError("Error: istates!")

        if (coefficients == None):
            coefficients = [None] * self.ntrajs
        else:
            if (self.ntrajs != len(self.coefficients)):
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
        self.phase = np.zeros((self.ntrajs, self.nst, self.nat, self.nsp))

        self.nstates_pair = int(self.nst * (self.nst - 1) / 2)
        self.qmom = np.zeros((self.ntrajs, self.nstates_pair, self.nat, self.nsp))

    def run(self, qm, mm=None, input_dir="./", save_qm_log=False, save_mm_log=False, save_scr=True, restart=None):
        # Initialize UNI-xMD
        base_dir, unixmd_dir, qm_log_dir, mm_log_dir =\
             self.run_init(qm, mm, input_dir, save_qm_log, save_mm_log, save_scr, restart)
        bo_list = [ist for ist in range(self.mol.nst)]
        qm.calc_coupling = True
        self.print_init()
        
        if (restart == None):
            # Calculate initial input geometry for all trajectories at t = 0.0 s
            self.istep = -1
            for itraj in range(self.ntrajs):
                self.mols[itraj].reset_bo(qm.calc_coupling)
                qm.get_data(self.mols[itraj], base_dir, bo_list, self.dt, self.istep, calc_force_only=False)
                # TODO: QM/MM
                self.mols[itraj].get_nacme()

                self.update_energy(itraj)
                
                self.get_phase(itraj)

                #self.write_md_output(unixmd_dir, self.istep)
                #self.print_step(self.istep)

        else: 
            raise ValueError ("restart option is invalid in CTMQC yet.")
           
        self.istep += 1

        # Main MD loop
        for istep in range(self.istep, self.nsteps):
            self.calculate_qmom()
            for itraj in range(self.ntrajs):
                self.cl_update_position(itraj)

                self.mols[itraj].backup_bo()
                self.mols[itraj].reset_bo(qm.calc_coupling)
                qm.get_data(self.mols[itraj], base_dir, bo_list, self.dt, istep, calc_force_only=False)
                # TODO: QM/MM

                self.mols[itraj].adjust_nac()

                self.cl_update_velocity(itraj)

                self.mols[itraj].get_nacme()

                # TODO: electronic propagation
                # el_run(self)

                # TODO: thermostat
                #if (self.thermo != None):
                #    self.thermo.run(self)

                self.update_energy(itraj)

                #if ((istep + 1) % self.out_freq == 0):
                #    self.write_md_output(unixmd_dir, istep)
                #    self.print_step(istep)
                #if (istep == self.nsteps - 1):
                #    self.write_final_xyz(unixmd_dir, istep)

                # TODO: restart
                #self.fstep = istep
                #restart_file = os.path.join(base_dir, "RESTART.bin")
                #with open(restart_file, 'wb') as f:
                #    pickle.dump({'qm':qm, 'md':self}, f)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

    def cl_update_position(self, itrajectory):
        """ Routine to update nuclear positions
        """
        self.calculate_force(itrajectory)
        
        self.mols[itrajectory].vel += 0.5 * self.dt * self.rforce / np.column_stack([self.mols[itrajectory].mass] * self.nsp)
        self.mols[itrajectory].pos += self.dt * self.mols[itrajectory].vel

    def cl_update_velocity(self, itrajectory):
        """ Routine to update nuclear velocities
        """
        self.calculate_force(itrajectory)

        self.mols[itrajectory].vel += 0.5 * self.dt * self.rforce / np.column_stack([self.mols[itrajectory].mass] * self.nsp)

    def calculate_force(self, itrajectory):
        """ Routine to calculate force
        """
        self.rforce = np.zeros((self.nat, self.nsp))

        for ist, istate in enumerate(self.mols[itrajectory].states):
            self.rforce += istate.force * self.mols[itrajectory].rho.real[ist, ist]

        for ist in range(self.nst):
            for jst in range(ist + 1, self.nst):
                self.rforce += 2. * self.mols[itrajectory].nac[ist, jst] * self.mols[itrajectory].rho.real[ist, jst] \
                    * (self.mols[itrajectory].states[ist].energy - self.mols[itrajectory].states[jst].energy)

    def update_energy(self, itrajectory):
        """ Routine to update the energy of molecules in CTMQC dynamics
        """
        # Update kinetic energy
        self.mols[itrajectory].update_kinetic()
        self.mols[itrajectory].epot = 0.
        for ist, istate in enumerate(self.mols[itrajectory].states):
            self.mols[itrajectory].epot += self.mols[itrajectory].rho.real[ist, ist] * istate.energy
        self.mols[itrajectory].etot = self.mols[itrajectory].epot + self.mols[itrajectory].ekin

    def get_phase(self, itrajectory):
        """ Routine to calculate phase
        """
        rho_threshold = 0.01
        for ist in range(self.nst):
            rho_l_real = self.mols[itrajectory].rho[ist, ist].real
        #    if ((rho_l_real < rho_threshold) or (rho_l_real > (1.0 - rho_threshold))):
        #        self.phase[itrajectory, ist] =
                

    def calculate_qmom(self):
        """ Routine to calculate quantum momentum
        """
        for ist in range(self.nst):
            for jst in range(self.nst):
                print(ist, jst)

    def print_init(self):
        """ Routine to print the initial information of dynamics
        """
        pass

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics
        """
        pass
