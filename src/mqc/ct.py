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
        self.ntrajs = len(molecules)
        self.istates = istates
        if ((self.istates != None) and (self.ntrajs != len(self.istates))):
            raise ValueError("Error: istates!")

        # Initialize initial coefficient for first trajectory
        if (coefficients == None):
            coefficients = [None] * self.ntrajs
            super().__init__(self.mols[0], thermostat, istates[0], dt, nsteps, nesteps, \
                propagation, solver, l_pop_print, l_adjnac, coefficients, unit_dt, out_freq, verbosity)
        else:
            if (self.ntrajs != len(self.coefficients)):
                raise ValueError("Error: coefficients!")
            super().__init__(self.mols[0], thermostat, istates[0], dt, nsteps, nesteps, \
                propagation, solver, l_pop_print, l_adjnac, coefficients[0], unit_dt, out_freq, verbosity)

        # Initialize initial coefficient for other trajectories
        for itraj in range(1, self.ntrajs):
            self.mols[itraj].get_coefficient(self.coefficients[itraj], self.istates[itraj])

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
            for itraj in self.mols:
                itraj.reset_bo(qm.calc_coupling)
                qm.get_data(itraj, base_dir, bo_list, self.dt, self.istep, calc_force_only=False)
                # TODO: QM/MM
                itraj.get_nacme()

                self.update_energy()

                self.write_md_output(unixmd_dir, self.istep)
                self.print_step(self.istep)
        else: 
            raise ValueError ("restart option is invalid in CTMQC yet.")
           
        self.istep += 1

        # Main MD loop
        for itraj in self.mols:
            self.calculate_qmom()
            for istep in range(self.istep, self.nsteps):
                self.cl_update_position()

                itraj.backup_bo()
                itraj.reset_bo(qm.calc_coupling)
                qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=False)
                # TODO: QM/MM

                itraj.adjust_nac()

                self.cl_update_velocity()

                itraj.get_nacme()

                # TODO: electronic propagation
                # el_run(self)

                # TODO: thermostat
                #if (self.thermo != None):
                #    self.thermo.run(self)

                self.update_energy()

                if ((istep + 1) % self.out_freq == 0):
                    self.write_md_output(unixmd_dir, istep)
                    self.print_step(istep)
                if (istep == self.nsteps - 1):
                    self.write_final_xyz(unixmd_dir, istep)

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

    def calculate_qmom(self):
        """ Routine to calculate quantum momentum
        """
        pass

    def calculate_force(self):
        """ Routine to calculate force
        """
        pass

    def update_energy(self):
        """ Routine to update the energy of molecules in CTMQC dynamics
        """
        pass

    def print_init(self):
        """ Routine to print the initial information of dynamics
        """
        pass

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics
        """
        pass
