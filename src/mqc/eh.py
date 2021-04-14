from __future__ import division
from build.el_propagator import el_run
from mqc.mqc import MQC
from misc import au_to_K, call_name, typewriter
import os, shutil, textwrap
import numpy as np
import pickle

class Eh(MQC):
    """ Class for Ehrenfest dynamics

        :param object molecule: Molecule object
        :param object thermostat: Thermostat object
        :param integer istate: Initial state
        :param double dt: Time interval
        :param integer nsteps: Total step of nuclear propation
        :param integer nesteps: Total step of electronic propagation
        :param string elec_object: Electronic equation of motions
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: Logical to print BO population and coherence
        :param boolean l_adj_nac: Logical to adjust nonadiabatic coupling
        :param init_coef: Initial BO coefficient
        :type init_coef: Double, list or complex, list
        :param string unit_dt: Unit of time step (fs = femtosecond, au = atomic unit)
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=20, \
        elec_object="density", propagator="rk4", l_print_dm=True, l_adj_nac=True, \
        init_coef=None, unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, nesteps, \
            elec_object, propagator, l_print_dm, l_adj_nac, init_coef, unit_dt, out_freq, verbosity)

        # Debug variables
        self.dotpopnac = np.zeros(self.mol.nst)

    def run(self, qm, mm=None, output_dir="./", l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to Ehrenfest dynamics

            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string output_dir: Name of directory where outputs to be saved.
            :param boolean l_save_qm_log: Logical for saving QM calculation log
            :param boolean l_save_mm_log: Logical for saving MM calculation log
            :param boolean l_save_scr: Logical for saving scratch directory
            :param string restart: Option for controlling dynamics restarting
        """
        # Initialize PyUNIxMD
        base_dir, unixmd_dir, qm_log_dir, mm_log_dir =\
             self.run_init(qm, mm, output_dir, l_save_qm_log, l_save_mm_log, l_save_scr, restart)
        bo_list = [ist for ist in range(self.mol.nst)]
        qm.calc_coupling = True
        self.print_init(qm, mm, restart)

        if (restart == None):
            # Calculate initial input geometry at t = 0.0 s
            self.istep = -1
            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, self.istep, calc_force_only=False)
            if (self.mol.l_qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, self.istep, calc_force_only=False)
            self.mol.get_nacme()

            self.update_energy()

            self.write_md_output(unixmd_dir, self.istep)
            self.print_step(self.istep)

        elif (restart == "write"):
            # Reset initial time step to t = 0.0 s
            self.istep = -1
            self.write_md_output(unixmd_dir, self.istep)
            self.print_step(self.istep)

        elif (restart == "append"):
            # Set initial time step to last successful step of previous dynamics
            self.istep = self.fstep

        self.istep += 1

        # Main MD loop
        for istep in range(self.istep, self.nsteps):
            
            self.calculate_force()
            self.cl_update_position()

            self.mol.backup_bo()
            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=False)
            if (self.mol.l_qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep, calc_force_only=False)

            if (self.l_adj_nac):
                self.mol.adjust_nac()

            self.calculate_force()
            self.cl_update_velocity()

            self.mol.get_nacme()

            el_run(self)

            if (self.thermo != None):
                self.thermo.run(self)

            self.update_energy()

            if ((istep + 1) % self.out_freq == 0):
                self.write_md_output(unixmd_dir, istep)
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

    def calculate_force(self):
        """ Calculate the Ehrenfest force
        """
        self.rforce = np.zeros((self.mol.nat, self.mol.ndim))

        for ist, istate in enumerate(self.mol.states):
            self.rforce += istate.force * self.mol.rho.real[ist, ist]

        for ist in range(self.mol.nst):
            for jst in range(ist + 1, self.mol.nst):
                self.rforce += 2. * self.mol.nac[ist, jst] * self.mol.rho.real[ist, jst] \
                    * (self.mol.states[ist].energy - self.mol.states[jst].energy)

    def update_energy(self):
        """ Routine to update the energy of molecules in Ehrenfest dynamics
        """
        # Update kinetic energy
        self.mol.update_kinetic()
        self.mol.epot = 0.
        for ist, istate in enumerate(self.mol.states):
            self.mol.epot += self.mol.rho.real[ist, ist] * self.mol.states[ist].energy
        self.mol.etot = self.mol.epot + self.mol.ekin

    def write_md_output(self, unixmd_dir, istep):
        """ Write output files

            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, istep)

        # Write time-derivative BO population
        self.write_dotpop(unixmd_dir, istep)

    def write_dotpop(self, unixmd_dir, istep):
        """ Write time-derivative BO population

            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        # Write NAC term in DOTPOPNAC
        if (self.verbosity >= 1):
            tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopnac])
            typewriter(tmp, unixmd_dir, "DOTPOPNAC", "a")

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

        # Print INIT for each step
        INIT = f" #INFO{'STEP':>8s}{'Kinetic(H)':>15s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'norm':>8s}"
        dynamics_step_info += INIT

        print (dynamics_step_info, flush=True)

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer istep: Current MD step
        """
        ctemp = self.mol.ekin * 2. / float(self.mol.ndof) * au_to_K
        norm = 0.
        for ist in range(self.mol.nst):
            norm += self.mol.rho.real[ist, ist]

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d} "
        INFO += f"{self.mol.ekin:14.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        INFO += f"{norm:11.5f}"
        print (INFO, flush=True)


