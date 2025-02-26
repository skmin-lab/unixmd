from __future__ import division
from mqc.mqc import MQC
from misc import au_to_K, call_name
import os, shutil, textwrap
import numpy as np
import pickle

class BOMD(MQC):
    """ Class for born-oppenheimer molecular dynamics (BOMD)

        :param object molecule: Molecule object
        :param object thermostat: Thermostat object
        :param integer istate: Initial state
        :param double dt: Time interval
        :param integer nsteps: Total step of nuclear propagation
        :param boolean l_adj_nac: Adjust nonadiabatic coupling to align the phases
        :param string unit_dt: Unit of time step
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, \
        l_adj_nac=False, unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, None, None, None, \
            False, l_adj_nac, None, unit_dt, out_freq, verbosity)

    def run(self, qm, mm=None, output_dir="./", l_coupling=False, l_save_bin=False, \
        l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to BOMD

            :param object qm: QM object containing on-the-fly calculation information
            :param object mm: MM object containing MM calculation information
            :param string output_dir: Name of directory where outputs to be saved.
            :param boolean l_coupling: Logical for calculation of nonadiabatic couplings
            :param boolean l_save_bin: Logical for saving the calculator and trajectory
            :param boolean l_save_qm_log: Logical for saving QM calculation log
            :param boolean l_save_mm_log: Logical for saving MM calculation log
            :param boolean l_save_scr: Logical for saving scratch directory
            :param string restart: Option for controlling dynamics restarting
        """
        # Initialize PyUNIxMD
        base_dir, unixmd_dir, samp_bin_dir, qm_log_dir, mm_log_dir = \
            self.run_init(qm, mm, output_dir, l_coupling, l_save_bin, l_save_qm_log, \
            l_save_mm_log, l_save_scr, restart)
        bo_list = [self.istate]
        qm.calc_coupling = l_coupling
        qm.calc_tdp = False
        qm.calc_tdp_grad = False
        self.print_init(qm, mm, l_coupling, restart)

        if (restart == None):
            # Calculate initial input geometry at t = 0.0 s
            self.istep = -1
            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, self.istep, calc_force_only=False)
            if (self.mol.l_qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, self.istep, calc_force_only=False)
            if (qm.calc_coupling):
                if (not self.mol.l_nacme):
                    self.mol.get_nacme()

            self.update_energy()
            if (l_save_bin):
                self.save_bin(samp_bin_dir, self.istep)

            self.write_md_output(unixmd_dir, qm.calc_coupling, self.istep)
            self.print_step(self.istep)

        elif (restart == "write"):
            # Reset initial time step to t = 0.0 s
            self.istep = -1
            self.write_md_output(unixmd_dir, qm.calc_coupling, self.istep)
            self.print_step(self.istep)

        elif (restart == "append"):
            # Set initial time step to last successful step of previous dynamics
            self.istep = self.fstep

        self.istep += 1

        # Main MD loop
        for istep in range(self.istep, self.nsteps):

            self.calculate_force()
            self.cl_update_position()

            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=False)
            if (self.mol.l_qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep, calc_force_only=False)

            if (qm.calc_coupling):
                if (not self.mol.l_nacme and self.l_adj_nac):
                    self.mol.adjust_nac()

            self.calculate_force()
            self.cl_update_velocity()

            if (qm.calc_coupling):
                if (not self.mol.l_nacme):
                    self.mol.get_nacme()

            if (self.thermo != None):
                self.thermo.run(self, self.mol)

            self.update_energy()
            if (l_save_bin):
                self.save_bin(samp_bin_dir, istep)

            if ((istep + 1) % self.out_freq == 0):
                self.write_md_output(unixmd_dir, qm.calc_coupling, istep)
                self.print_step(istep)
            if (istep == self.nsteps - 1):
                self.write_final_xyz(unixmd_dir, istep)

            self.fstep = istep
            restart_file = os.path.join(base_dir, "RESTART.bin")
            with open(restart_file, 'wb') as f:
                pickle.dump({'qm':qm, 'md':self}, f)

        # Delete scratch directory
        if (not l_save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            if (self.mol.l_qmmm and mm != None):
                tmp_dir = os.path.join(unixmd_dir, "scr_mm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def save_bin(self, samp_bin_dir, l_coupling, istep):
        """ Routine to save the calculator and trajectory for each time step

            :param string samp_bin_dir: Name of directory where the calculator and trajectory are saved
            :param boolean l_coupling: Logical for calculation of nonadiabatic couplings
            :param integer istep: Current MD step
        """
        filename = os.path.join(samp_bin_dir, f"QM.{istep + 1}.bin")
        if (l_coupling):
            with open(filename, "wb") as f:
                pickle.dump({"energy":np.array([states.energy for states in self.mol.states]), \
                    "force":self.rforce, "nacme":self.mol.nacme}, f)
        else:
            with open(filename, "wb") as f:
                pickle.dump({"energy":np.array([states.energy for states in self.mol.states]), \
                    "force":self.rforce, f)

        filename = os.path.join(samp_bin_dir, f"RV.{istep + 1}.bin")
        with open(filename, "wb") as f:
            pickle.dump({"pos":self.mol.pos, "vel":self.mol.vel}, f)

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        self.rforce = np.copy(self.mol.states[self.istate].force)

    def update_energy(self):
        """ Routine to update the energy of molecules in BOMD
        """
        # Update kinetic energy
        self.mol.update_kinetic()
        self.mol.epot = self.mol.states[self.istate].energy
        self.mol.etot = self.mol.epot + self.mol.ekin

    def print_init(self, qm, mm, l_coupling, restart):
        """ Routine to print the initial information of dynamics

            :param object qm: QM object containing on-the-fly calculation information
            :param object mm: MM object containing MM calculation information
            :param boolean l_coupling: Logical for calculation of nonadiabatic couplings
            :param string restart: Option for controlling dynamics restarting
        """
        # Print initial information about molecule, qm, mm and thermostat
        super().print_init(qm, mm, l_coupling, restart)

        # Print dynamics information for start line
        dynamics_step_info = textwrap.dedent(f"""\

        {"-" * 118}
        {"Start Dynamics":>65s}
        {"-" * 118}
        """)

        # Print INIT for each step
        INIT = f" #INFO{'STEP':>8s}{'State':>7s}{'Kinetic(H)':>13s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}"
        dynamics_step_info += INIT

        print (dynamics_step_info, flush=True)

    def print_step(self, istep):
        """ Routine to print each steps information about dynamics

            :param integer istep: Current MD step
        """
        ctemp = self.mol.ekin * 2. / float(self.mol.ndof) * au_to_K

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d}{self.istate:>5d} "
        INFO += f"{self.mol.ekin:14.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        print (INFO, flush=True)


