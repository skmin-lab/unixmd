from __future__ import division
from build.el_propagator import el_run
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz
from misc import au_to_K, call_name
import os, shutil, textwrap
import numpy as np

class Eh(MQC):
    """ Class for Ehrenfest dynamics

        :param object molecule: molecule object
        :param object thermostat: thermostat type
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string propagation: propagation scheme
        :param string solver: propagation solver
        :param boolean l_pop_print: logical to print BO population and coherence
        :param boolean l_adjnac: logical to adjust nonadiabatic coupling
        :param coefficient: initial BO coefficient
        :type coefficient: double, list or complex, list
        :param string unit_dt: unit of time step (fs = femtosecond, au = atomic unit)
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", solver="rk4", l_pop_print=False, l_adjnac=True, \
        coefficient=None, unit_dt="fs"):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, nesteps, \
            propagation, solver, l_pop_print, l_adjnac, coefficient, unit_dt)

    def run(self, qm, mm=None, input_dir="./", \
        save_QMlog=False, save_MMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to Ehrenfest dynamics

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param string input_dir: location of input directory
            :param boolean save_QMlog: logical for saving QM calculation log
            :param boolean save_MMlog: logical for saving MM calculation log
            :param boolean save_scr: logical for saving scratch directory
            :param integer debug: verbosity level for standard output
        """
        # Set directory information
        input_dir = os.path.expanduser(input_dir)
        base_dir = os.path.join(os.getcwd(), input_dir)

        unixmd_dir = os.path.join(base_dir, "md")
        if (os.path.exists(unixmd_dir)):
            shutil.rmtree(unixmd_dir)
        os.makedirs(unixmd_dir)

        QMlog_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(QMlog_dir)):
            shutil.rmtree(QMlog_dir)
        if (save_QMlog):
            os.makedirs(QMlog_dir)

        if (self.mol.qmmm and mm != None):
            MMlog_dir = os.path.join(base_dir, "MMlog")
            if (os.path.exists(MMlog_dir)):
                shutil.rmtree(MMlog_dir)
            if (save_MMlog):
                os.makedirs(MMlog_dir)

        if ((self.mol.qmmm and mm == None) or (not self.mol.qmmm and mm != None)):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Both self.mol.qmmm and mm object is necessary! {self.mol.qmmm} and {mm}")

        # Check compatibility for QM and MM objects
        if (self.mol.qmmm and mm != None):
            self.check_qmmm(qm, mm)

        if (self.mol.l_nacme):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Ehrenfest dynamics needs NAC! {self.mol.l_nacme}")

        # Initialize UNI-xMD
        os.chdir(base_dir)
        bo_list = [ist for ist in range(self.mol.nst)]
        qm.calc_coupling = True

        touch_file(self.mol, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, SH_chk=False, XF_chk=False)
        self.print_init(qm, mm, debug)

        # Calculate initial input geometry at t = 0.0 s
        self.mol.reset_bo(qm.calc_coupling)
        qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=-1, calc_force_only=False)
        if (self.mol.qmmm and mm != None):
            mm.get_data(self.mol, base_dir, bo_list, istep=-1, calc_force_only=False)
        if (not self.mol.l_nacme):
            self.mol.get_nacme()

        self.update_energy()

        write_md_output(self.mol, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=-1)
        self.print_step(debug, istep=-1)

        # Main MD loop
        for istep in range(self.nsteps):

            self.cl_update_position()

            self.mol.backup_bo()
            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=istep, calc_force_only=False)
            if (self.mol.qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep=istep, calc_force_only=False)

            if (not self.mol.l_nacme):
                self.mol.adjust_nac()

            self.cl_update_velocity()

            if (not self.mol.l_nacme):
                self.mol.get_nacme()

            el_run(self)

            if (self.thermo != None):
                self.thermo.run(self)

            self.update_energy()

            write_md_output(self.mol, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=istep)
            self.print_step(debug, istep=istep)
            if (istep == self.nsteps - 1):
                write_final_xyz(self.mol, unixmd_dir, istep=istep)

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
        """ Calculate the Ehrenfest force
        """
        self.rforce = np.zeros((self.mol.nat, self.mol.nsp))

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

    def print_init(self, qm, mm, debug):
        """ Routine to print the initial information of dynamics

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param integer debug: verbosity level for standard output
        """
        # Print initial information about molecule, qm, mm and thermostat
        super().print_init(qm, mm, debug)

        # Print dynamics information for start line
        dynamics_step_info = textwrap.dedent(f"""\

        {"-" * 118}
        {"Start Dynamics":>65s}
        {"-" * 118}
        """)

        # Print INIT for each step
        INIT = f" #INFO{'STEP':>8s}{'Kinetic(H)':>15s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'norm':>8s}"
        dynamics_step_info += INIT

        # Print DEBUG1 for each step
        if (debug >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        print (dynamics_step_info, flush=True)

    def print_step(self, debug, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer debug: verbosity level for standard output
            :param integer istep: current MD step
        """
        ctemp = self.mol.ekin * 2. / float(self.mol.dof) * au_to_K
        norm = 0.
        for ist in range(self.mol.nst):
            norm += self.mol.rho.real[ist, ist]

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d} "
        INFO += f"{self.mol.ekin:14.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        INFO += f"{norm:11.5f}"
        print (INFO, flush=True)

        # Print DEBUG1 for each step
        if (debug >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{self.mol.states[ist].energy:17.8f} "
            print (DEBUG1, flush=True)


