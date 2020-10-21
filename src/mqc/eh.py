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
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", solver="rk4", l_pop_print=False, l_adjnac=True, coefficient=None):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, solver, l_pop_print, l_adjnac, coefficient)

    def run(self, molecule, qm, mm=None, thermostat=None, input_dir="./", \
        save_QMlog=False, save_MMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to Ehrenfest dynamics

            :param object molecule: molecule object
            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param object thermostat: thermostat type
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

        if (molecule.qmmm and mm != None):
            MMlog_dir = os.path.join(base_dir, "MMlog")
            if (os.path.exists(MMlog_dir)):
                shutil.rmtree(MMlog_dir)
            if (save_MMlog):
                os.makedirs(MMlog_dir)

        if ((molecule.qmmm and mm == None) or (not molecule.qmmm and mm != None)):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Both molecule.qmmm and mm object is necessary! {molecule.qmmm} and {mm}")

        # Check compatibility for QM and MM objects
        if (molecule.qmmm and mm != None):
            self.check_qmmm(qm, mm)

        if (molecule.l_nacme):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Ehrenfest dynamics needs NAC! {molecule.l_nacme}")

        # Initialize UNI-xMD
        os.chdir(base_dir)
        bo_list = [ist for ist in range(molecule.nst)]
        qm.calc_coupling = True

        touch_file(molecule, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, SH_chk=False)
        self.print_init(molecule, qm, mm, thermostat, debug)

        # Calculate initial input geometry at t = 0.0 s
        molecule.reset_bo(qm.calc_coupling)
        qm.get_data(molecule, base_dir, bo_list, self.dt, istep=-1, calc_force_only=False)
        if (molecule.qmmm and mm != None):
            mm.get_data(molecule, base_dir, bo_list, istep=-1, calc_force_only=False)
        if (not molecule.l_nacme):
            molecule.get_nacme()

        self.update_energy(molecule)

        write_md_output(molecule, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=-1)
        self.print_step(molecule, debug, istep=-1)

        # Main MD loop
        for istep in range(self.nsteps):

            self.cl_update_position(molecule)

            molecule.backup_bo()
            molecule.reset_bo(qm.calc_coupling)
            qm.get_data(molecule, base_dir, bo_list, self.dt, istep=istep, calc_force_only=False)
            if (molecule.qmmm and mm != None):
                mm.get_data(molecule, base_dir, bo_list, istep=istep, calc_force_only=False)

            if (not molecule.l_nacme):
                molecule.adjust_nac()

            self.cl_update_velocity(molecule)

            if (not molecule.l_nacme):
                molecule.get_nacme()

            el_run(self, molecule)

            if (thermostat != None):
                thermostat.run(molecule, self)

            self.update_energy(molecule)

            write_md_output(molecule, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=istep)
            self.print_step(molecule, debug, istep=istep)
            if (istep == self.nsteps - 1):
                write_final_xyz(molecule, unixmd_dir, istep=istep)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            if (molecule.qmmm and mm != None):
                tmp_dir = os.path.join(unixmd_dir, "scr_mm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def calculate_force(self, molecule):
        """ Calculate the Ehrenfest force

            :param object molecule: molecule object
        """
        self.rforce = np.zeros((molecule.nat, molecule.nsp))

        for ist, istate in enumerate(molecule.states):
            self.rforce += istate.force * molecule.rho.real[ist, ist]

        for ist in range(molecule.nst):
            for jst in range(ist + 1, molecule.nst):
                self.rforce += 2. * molecule.nac[ist, jst] * molecule.rho.real[ist, jst] \
                    * (molecule.states[ist].energy - molecule.states[jst].energy)

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in Ehrenfest dynamics

            :param object molecule: molecule object
        """
        # Update kinetic energy
        molecule.update_kinetic()
        molecule.epot = 0.
        for ist, istate in enumerate(molecule.states):
            molecule.epot += molecule.rho.real[ist, ist] * molecule.states[ist].energy
        molecule.etot = molecule.epot + molecule.ekin

    def print_init(self, molecule, qm, mm, thermostat, debug):
        """ Routine to print the initial information of dynamics

            :param object molecule: molecule object
            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param object thermostat: thermostat type
            :param integer debug: verbosity level for standard output
        """
        # Print initial information about molecule, qm, mm and thermostat
        super().print_init(molecule, qm, mm, thermostat, debug)

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
            for ist in range(molecule.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        print (dynamics_step_info, flush=True)

    def print_step(self, molecule, debug, istep):
        """ Routine to print each steps infomation about dynamics

            :param object molecule: molecule object
            :param integer debug: verbosity level for standard output
            :param integer istep: current MD step
        """
        ctemp = molecule.ekin * 2. / float(molecule.dof) * au_to_K
        norm = 0.
        for ist in range(molecule.nst):
            norm += molecule.rho.real[ist, ist]

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d} "
        INFO += f"{molecule.ekin:14.8f}{molecule.epot:15.8f}{molecule.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        INFO += f"{norm:11.5f}"
        print (INFO, flush=True)

        # Print DEBUG1 for each step
        if (debug >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            for ist in range(molecule.nst):
                DEBUG1 += f"{molecule.states[ist].energy:17.8f} "
            print (DEBUG1, flush=True)


