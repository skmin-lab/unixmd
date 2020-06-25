from __future__ import division
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz
import os, shutil, textwrap
import numpy as np

class CT(MQC):
    """ Class for Coupled-Trajectory Mixed Quantum-Classical dynamics (CTMQC)
    """
    def __init__(self, molecules, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True, threshold=1.0e-4, dist_cutoff=2.0/0.529166):
        # Read total number of trajectories
        self.ntrajs = len(molecules)

        # Initialize input values: initial coef. and density of each trajectories
        for itraj in range(self.ntrajs):
            super().__init__(molecules[itraj], istate, dt, nsteps, nesteps, \
                propagation, l_adjnac)

        self.nsp = molecules[0].nsp
        self.nst = molecules[0].nst
        self.nat = molecules[0].nat

        # Initialize XF variables
        self.nst_pairs = self.nst * (self.nst - 1) 

        self.phase = np.zeros((self.ntrajs, self.nst, self.nat, self.nsp))
        self.qmom = np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp))
        self.sigma = np.zeros((self.ntrajs, self.nat, self.nsp))
        self.sigma[:] = 0.08753727/self.ntrajs

        self.xfforce = np.zeros((self.ntrajs, self.nat, self.nsp))
        self.xfcdot = np.zeros((self.ntrajs, self.nst), dtype=complex)

        # Trajectory-dependent temporaries to obtain qmom
        self.prod_g_ij = np.zeros((self.ntrajs, self.ntrajs))
        self.g_i = np.zeros((self.ntrajs))
        self.w_ij = np.zeros((self.ntrajs, self.ntrajs, self.nat, self.nsp))
        self.slope = np.zeros((self.ntrajs, self.nat, self.nsp))
        self.w0 = np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp))
        self.d0 = np.zeros((self.nst_pairs, self.nat, self.nsp))
        self.r0 = np.zeros((self.nst_pairs, self.nat, self.nsp))

        self.qmom_center = np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp))
        self.qmom_center_old = np.zeros((self.ntrajs, self.nst_pairs, self.nat, self.nsp))
        self.k_lk = np.zeros((self.ntrajs, self.nst, self.nst))

    def run(self, molecules, theory, thermostat, input_dir='./', \
        save_QMlog=False, save_scr=True):
        """ Run MQC dynamics according to CTMQC
        """
        # Set directory information
        input_dir = os.path.expanduser(input_dir)
        base_dir = os.path.join(os.getcwd(), input_dir)

        unixmd_dir = os.path.join(base_dir, "md")
        if (os.path.exists(unixmd_dir)):
            shutil.rmtree(unixmd_dir)
        os.makedirs(unixmd_dir)

        traj_dirs = []
        for itraj in range(self.ntrajs):
            traj_dirs.append(os.path.join(unixmd_dir, "traj"+str(itraj)))
            os.makedirs(traj_dirs[itraj])

        QMlog_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(QMlog_dir)):
            shutil.rmtree(QMlog_dir)
        if (save_QMlog):
            os.makedirs(QMlog_dir)
            for itraj in range(self.ntrajs):
                traj_dir = os.path.join(QMlog_dir, "traj"+str(itraj))
                os.makedirs(traj_dir)

        # initialize unixmd
        os.chdir(base_dir)
        bo_list = [ist for ist in range(self.nst)]
        theory.calc_coupling = True
        for itraj in range(self.ntrajs):
            touch_file(molecules[itraj], theory.calc_coupling, self.propagation, \
                traj_dirs[itraj], SH_chk=False)
        self.print_init(molecule, theory, thermostat, debug)

        # Calculate initial input geometry at t = 0.0 s
        for itraj in range(self.ntrajs):
            theory.get_bo(molecules[itraj], base_dir, -1, bo_list, calc_force_only=False)
            if (not molecules[itraj].l_nacme):
                molecules[itraj].get_nacme()
            
            self.update_energy(molecules[itraj])
        
            write_md_output(molecules[itraj], theory.calc_coupling, -1, \
                self.propagation, traj_dirs[itraj])
        
        self.print_step(molecule, -1, debug)

        # Main MD loop
        #for istep in range(self.nsteps):
        #    # calculate qmom
        #    self.calculate_qmom(molecules)
        #    for itraj in range(self.ntrajs):

                

        
    def calculate_qmom(self, molecules):
        """ Routine to calculate quantum momentum
        """
        pass

    def calculate_force(self, states):
        """ Routine to calculate force
        """
        pass

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in surface hopping
        """
        # Update kinetic energy
        molecule.update_kinetic()
        molecule.epot = 0.
        for ist in enumerate(molecule.nst):
            molecule.epot += molecule.rho.real[ist, ist] * molecule.states[ist].energy
        molecule.etot = molecule.epot + molecule.ekin

    def el_propagator(self, molecule):
        """ Routine to propagate BO coefficients or density matrix
        """
        pass

    def print_init(self, molecule, theory, thermostat, debug):
        """ Routine to print the initial information of dynamics
        """
        pass

    def print_step(self, molecule, istep, debug):
        """ Routine to print each steps infomation about dynamics
        """
        pass

