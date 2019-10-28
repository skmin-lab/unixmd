from __future__ import division
import numpy as np
import os, shutil
from mqc.mqc import MQC
from fileio import unixmd_init, write_md_output, write_final_xyz
from mqc.el_prop.el_propagator import *

class Eh(MQC):
    """ ehrenfest dynamics
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, l_adjnac)

    def run(self, molecule, theory, thermostat, input_dir="./", \
        save_QMlog=False, save_scr=True):
        """ run MQC dynamics according to Ehrenfest dynamics
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
        bo_list = [ist for ist in range(molecule.nst)]
        theory.calc_coupling = True
        unixmd_init(molecule, theory.calc_coupling, self.propagation, \
            unixmd_dir, SH_chk=False)
        if (molecule.l_nacme):
            raise ValueError ("Ehrenfest requries NAC calculation")

        # calculate initial input geometry at t = 0.0 s
        theory.get_bo(molecule, base_dir, -1, bo_list, calc_force_only=False)
        if (not molecule.l_nacme):
            molecule.get_nacme()

        self.update_energy(molecule)

        write_md_output(molecule, theory.calc_coupling, -1, \
            self.propagation, unixmd_dir)

        # main MD loop
        for istep in range(self.nsteps):

            self.cl_update_position(molecule)

            molecule.backup_bo()
            theory.get_bo(molecule, base_dir, istep, bo_list, calc_force_only=False)

            if (not molecule.l_nacme):
                molecule.adjust_nac()

            self.cl_update_velocity(molecule)
            
            if (not molecule.l_nacme):
                molecule.get_nacme()

            self.el_propagator(molecule)

            thermostat.run(molecule, self)

            self.update_energy(molecule)

            write_md_output(molecule, theory.calc_coupling, istep, \
                self.propagation, unixmd_dir)
            if (istep == self.nsteps - 1):
                write_final_xyz(molecule, istep, unixmd_dir)

        # delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(base_dir, "md/scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

    def calculate_force(self, molecule):
        """ Calculate the Ehrenfest forces
        """
        self.rforce = np.zeros((molecule.nat, molecule.nsp))

        for ist, istate in enumerate(molecule.states):
            self.rforce += istate.force * molecule.rho.real[ist, ist]

        for ist in range(molecule.nst):
            for jst in range(ist + 1, molecule.nst):
                self.rforce += 2. * molecule.nac[ist, jst] * molecule.rho.real[ist, jst] \
                    * (molecule.states[ist].energy - molecule.states[jst].energy)

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in ehrenfest
        """
        # update kinetic energy
        molecule.update_kinetic()
        molecule.epot = 0.
        for ist, istate in enumerate(molecule.states):
            molecule.epot += molecule.rho.real[ist, ist] * molecule.states[ist].energy
        molecule.etot = molecule.epot + molecule.ekin

    def el_propagator(self, molecule):
        """ Routine to propagate BO coefficients or density matrix
        """
        if (self.propagation == "coefficient"):
            el_coef(self.nesteps, self.dt, molecule)
        elif (self.propagation == "density"):
            el_rho(self.nesteps, self.dt, molecule)
        else:
            raise ValueError ("Other propagators Not Implemented")


