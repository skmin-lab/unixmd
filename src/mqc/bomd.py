from __future__ import division
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz
from misc import au_to_K
import os, shutil
import numpy as np

class BOMD(MQC):
    """ Class for born-oppenheimer molecular dynamics (BOMD)

        :param object molecule: molecule object
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, None, None, None)

    def run(self, molecule, theory, thermostat, input_dir="./", \
        save_QMlog=False, save_scr=True):
        """ Run MQC dynamics according to BOMD

            :param object molecule: molecule object
            :param object theory: theory object containing on-the-fly calculation infomation
            :param object thermostat: thermostat type
            :param string input_dir: location of input directory
            :param boolean save_QMlog: logical for saving QM calculation log
            :param boolean save_scr: logical for saving scratch directory
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

        # Initialize UNI-xMD
        os.chdir(base_dir)
        bo_list = [self.istate]
        theory.calc_coupling = False

        touch_file(molecule, theory.calc_coupling, None, unixmd_dir, SH_chk=False)
        self.print_init(molecule, theory, thermostat)

        # Calculate initial input geometry at t = 0.0 s
        theory.get_bo(molecule, base_dir, -1, bo_list, calc_force_only=False)

        self.update_energy(molecule)

        write_md_output(molecule, theory.calc_coupling, -1, None, unixmd_dir)
        self.print_step(molecule, -1)

        # Main MD loop
        for istep in range(self.nsteps):

            self.cl_update_position(molecule)

            theory.get_bo(molecule, base_dir, istep, bo_list, calc_force_only=False)

            self.cl_update_velocity(molecule)

            thermostat.run(molecule, self)

            self.update_energy(molecule)

            write_md_output(molecule, theory.calc_coupling, istep, None, unixmd_dir)
            self.print_step(molecule, istep)
            if (istep == self.nsteps - 1):
                write_final_xyz(molecule, istep, unixmd_dir)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

    def calculate_force(self, molecule):
        """ Routine to calculate the forces

            :param object molecule: molecule object
        """
        self.rforce = np.copy(molecule.states[self.istate].force)

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in BOMD

            :param object molecule: molecule object
        """
        # Update kinetic energy
        molecule.update_kinetic()
        molecule.epot = molecule.states[self.istate].energy
        molecule.etot = molecule.epot + molecule.ekin

    # TODO : add argument
    #def print_step(self, molecule, istep, debug=0):
    def print_step(self, molecule, istep):
        """ Routine to print each steps infomation about dynamics

            :param object molecule: molecule object
            :param integer istep: current MD step
        """
        ctemp = molecule.ekin * 2. / float(molecule.dof) * au_to_K

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d}{self.istate:>5d} "
        INFO += f"{molecule.ekin:14.8f}{molecule.epot:15.8f}{molecule.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        print (INFO, flush=True)

        # Print DEBUG1 for each step
        # TODO : if (debug=1):
        # TODO : debug option print (add argument)
        DEBUG1 = f" DEBUG1{istep + 1:>7d}"
        for ist in range(molecule.nst):
            DEBUG1 += f"{molecule.states[ist].energy:17.8f} "
        print (DEBUG1, flush=True)



