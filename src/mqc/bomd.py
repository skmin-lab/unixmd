from __future__ import division
import numpy as np
import os, shutil
from mqc.mqc import MQC
from fileio import unixmd_init, write_md_output, write_final_xyz

class BOMD(MQC):    
    """ born-oppenheimer molecular dynamics
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, None, None, None)

    def run(self, molecule, theory, thermostat, input_dir="./", \
        save_QMlog=False, save_scr=True):
        """ run MQC dynamics according to BOMD
        """
        # set directory information
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

        # initialize unixmd
        os.chdir(base_dir)
        bo_list = [self.istate]
        theory.calc_coupling = False
        unixmd_init(molecule, theory.calc_coupling, None, unixmd_dir, SH_chk=False)

        # calculate initial input geometry at t = 0.0 s
        theory.get_bo(molecule, base_dir, -1, bo_list, calc_force_only=False)

        self.update_energy(molecule)

        write_md_output(molecule, theory.calc_coupling, -1, None, unixmd_dir)

        # main MD loop
        for istep in range(self.nsteps):

            self.cl_update_position(molecule)

            theory.get_bo(molecule, base_dir, istep, bo_list, calc_force_only=False)

            self.cl_update_velocity(molecule)

            thermostat.run(molecule, self)

            self.update_energy(molecule)

            write_md_output(molecule, theory.calc_coupling, istep, None, unixmd_dir)
            if (istep == self.nsteps - 1):
                write_final_xyz(molecule, istep, unixmd_dir)

        # delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(base_dir, "md/scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

    def calculate_force(self, molecule):
        """ Routine to calculate the forces
        """
        self.rforce = np.copy(molecule.states[self.istate].force)

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in BOMD
        """
        # update kinetic energy
        molecule.update_kinetic()
        molecule.epot = molecule.states[self.istate].energy
        molecule.etot = molecule.epot + molecule.ekin


