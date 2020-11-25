from __future__ import division
from mqc.mqc import MQC
from misc import au_to_K, call_name
import os, shutil, textwrap
import numpy as np
import pickle

class BOMD(MQC):
    """ Class for born-oppenheimer molecular dynamics (BOMD)

        :param object molecule: molecule object
        :param object thermostat: thermostat type
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param string unit_dt: unit of time step (fs = femtosecond, au = atomic unit)
        :param integer out_freq: frequency of printing output
        :param integer verbosity: verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, None, None, None, \
            False, None, None, unit_dt, out_freq, verbosity)

    def run(self, qm, mm=None, input_dir="./", save_QMlog=False, save_MMlog=False, save_scr=True, restart=None):
        """ Run MQC dynamics according to BOMD

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param string input_dir: location of input directory
            :param boolean save_QMlog: logical for saving QM calculation log
            :param boolean save_MMlog: logical for saving MM calculation log
            :param boolean save_scr: logical for saving scratch directory
        """
        # Initialize UNI-xMD
        base_dir, unixmd_dir, QMlog_dir, MMlog_dir =\
             self.run_init(qm, mm, input_dir, save_QMlog, save_MMlog, save_scr, restart)
        bo_list = [self.istate]
        qm.calc_coupling = False
        self.print_init(qm, mm)

        if (restart == None):
            
            # Calculate initial input geometry at t = 0.0 s
            self.istep = -1
            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=self.istep, calc_force_only=False)
            if (self.mol.qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep=self.istep, calc_force_only=False)
            self.update_energy()
            self.write_md_output(unixmd_dir, istep=self.istep)
            self.print_step(istep=self.istep)
        
        elif (restart == "write"):
            self.istep = -1
            self.write_md_output(unixmd_dir, istep=self.istep)
            self.print_step(istep=self.istep)

        else:
            self.istep = self.fstep
        
        self.istep += 1

        # Main MD loop
        for istep in range(self.istep, self.nsteps):

            self.cl_update_position()

            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=istep, calc_force_only=False)
            if (self.mol.qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep=istep, calc_force_only=False)

            self.cl_update_velocity()

            if (self.thermo != None):
                self.thermo.run(self)

            self.update_energy()

            if ((istep + 1) % self.out_freq == 0):
                self.write_md_output(unixmd_dir, istep=istep)
                self.print_step(istep=istep)
            if (istep == self.nsteps - 1):
                self.write_final_xyz(unixmd_dir, istep=istep)
            
            self.fstep = istep
            restart_file = os.path.join(base_dir, "RESTART.bin")
            with open(restart_file, 'wb') as f:
                pickle.dump({'qm':qm, 'md':self}, f)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            if (self.mol.qmmm and mm != None):
                tmp_dir = os.path.join(unixmd_dir, "scr_mm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

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

    def print_init(self, qm, mm):
        """ Routine to print the initial information of dynamics

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
        """
        # Print initial information about molecule, qm, mm and thermostat
        super().print_init(qm, mm)

        # Print dynamics information for start line
        dynamics_step_info = textwrap.dedent(f"""\

        {"-" * 118}
        {"Start Dynamics":>65s}
        {"-" * 118}
        """)

        # Print INIT for each step
        INIT = f" #INFO{'STEP':>8s}{'State':>7s}{'Kinetic(H)':>13s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}"
        dynamics_step_info += INIT

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        print (dynamics_step_info, flush=True)

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer istep: current MD step
        """
        ctemp = self.mol.ekin * 2. / float(self.mol.dof) * au_to_K

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d}{self.istate:>5d} "
        INFO += f"{self.mol.ekin:14.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        print (INFO, flush=True)

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{self.mol.states[ist].energy:17.8f} "
            print (DEBUG1, flush=True)


