from __future__ import division
from build.el_propagator_xf import el_run
from mqc.mqc import MQC
from misc import au_to_K, call_name, typewriter
import random, os, shutil, textwrap
import numpy as np
import pickle

class Auxiliary_Molecule(object):
    """ Class for auxiliary molecule that is used for the calculation of decoherence term

        :param object molecule: Molecule object
    """
    def __init__(self, molecule):
        # Initialize auxiliary molecule
        self.nat = molecule.nat_qm
        self.ndim = molecule.ndim
        self.symbols = np.copy(molecule.symbols[0:molecule.nat_qm])

        self.mass = np.copy(molecule.mass[0:molecule.nat_qm])

        self.pos = np.zeros((molecule.nst, self.nat, self.ndim))
        self.vel = np.zeros((molecule.nst, self.nat, self.ndim))
        self.vel_old = np.copy(self.vel)


class EhXF(MQC):
    """ Class for EhXF dynamics

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
        :param double rho_threshold: Electronic density threshold for decoherence term calculation
        :param sigma: Width of nuclear wave packet of auxiliary trajectory
        :type sigma: double or double,list
        :param init_coef: Initial BO coefficient
        :type init_coef: double, list or complex, list
        :param boolean l_xf_force: Logical to inlcude XF contribution to the total force
        :param boolean l_econs_state: Logical to use identical total energies for all auxiliary trajectories
        :param string unit_dt: Unit of time step (fs = femtosecond, au = atomic unit)
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=20, \
        elec_object="density", propagator="rk4", l_print_dm=True, l_adj_nac=True, \
        rho_threshold=0.01, sigma=None, init_coef=None, l_xf_force=True,  l_econs_state=True, \
        unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, nesteps, \
            elec_object, propagator, l_print_dm, l_adj_nac, init_coef, unit_dt, out_freq, verbosity)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.rand = 0.
        self.prob = np.zeros(self.mol.nst)
        self.acc_prob = np.zeros(self.mol.nst + 1)

        self.l_hop = False
        self.l_reject = False

        ##shkim no hop_rescale~
        # Initialize XF related variables
        self.force_hop = False
        self.l_xf_force = l_xf_force 
        self.l_econs_state = l_econs_state
        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst
        # TODO : l_fix?
        #self.l_fix = [False] * self.mol.nst
        self.rho_threshold = rho_threshold
        # TODO : aux_econs_viol?
        #self.aux_econs_viol = aux_econs_viol

        self.sigma = sigma
        if (self.sigma == None):
            error_message = "Sigma for auxiliary trajectories must be set in running script!"
            error_vars = f"sigma = {self.sigma}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        if (isinstance(self.sigma, float)):
            # uniform value for sigma
            pass
        elif (isinstance(self.sigma, list)):
            # atom-resolved values for sigma
            if (len(self.sigma) != self.mol.nat_qm):
                error_message = "Number of elements for sigma must be equal to number of atoms!"
                error_vars = f"len(sigma) = {len(self.sigma)}"
                raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of sigma must be float or list consisting of float!"
            error_vars = f"sigma = {self.sigma}"
            raise TypeError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.upper_th = 1. - self.rho_threshold
        self.lower_th = self.rho_threshold

        # Initialize auxiliary molecule object
        self.aux = Auxiliary_Molecule(self.mol)
        self.pos_0 = np.zeros((self.aux.nat, self.aux.ndim))
        self.phase = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))

        # Debug variables
        self.dotpopdec = np.zeros(self.mol.nst)
        self.dotpopnac = np.zeros(self.mol.nst)
        self.qmom = np.zeros((self.aux.nat, self.aux.ndim))

    def run(self, qm, mm=None, output_dir="./", l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to decoherence-induced Ehrenfest dynamics

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
        bo_list = [self.rstate]
        qm.calc_coupling = True
        self.print_init(qm, mm, restart)

        if (restart == None):
            # Initialize decoherence variables
            self.append_sigma()

            # Calculate initial input geometry at t = 0.0 s
            self.istep = -1
            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, self.istep, calc_force_only=False)
            if (self.mol.l_qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, self.istep, calc_force_only=False)
            self.mol.get_nacme()

            self.hop_prob()
            self.hop_check(bo_list)
            if (self.l_hop):
                if (qm.re_calc):
                    qm.get_data(self.mol, base_dir, bo_list, self.dt, self.istep, calc_force_only=True)
                if (self.mol.l_qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, self.istep, calc_force_only=True)

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

            self.hop_prob()
            self.hop_check(bo_list)
            if (self.l_hop):
                if (qm.re_calc):
                    qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=True)
                if (self.mol.l_qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, istep, calc_force_only=True)

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

    def hop_prob(self):
        """ Routine to calculate hopping probabilities

            :param integer istep: Current MD step
        """
        # Reset surface hopping variables
        self.rstate_old = self.rstate

        self.prob = np.zeros(self.mol.nst)
        self.acc_prob = np.zeros(self.mol.nst + 1)

        self.l_hop = False
        self.force_hop = False

        accum = 0.

        if (self.mol.rho.real[self.rstate, self.rstate] < self.lower_th):
            self.force_hop = True

        for ist in range(self.mol.nst):
            if (ist != self.rstate):
                if (self.force_hop):
                    self.prob[ist] = self.mol.rho.real[ist, ist] / self.upper_th
                else:
                    self.prob[ist] = - 2. * self.mol.rho.real[ist, self.rstate] * \
                        self.mol.nacme[ist, self.rstate] * self.dt / self.mol.rho.real[self.rstate, self.rstate]

                if (self.prob[ist] < 0.):
                    self.prob[ist] = 0.
                accum += self.prob[ist]
            self.acc_prob[ist + 1] = accum
        psum = self.acc_prob[self.mol.nst]

        if (psum > 1.):
            self.prob /= psum
            self.acc_prob /= psum

    def hop_check(self, bo_list):
        """ Routine to check hopping occurs with random number

            :param integer,list bo_list: List of BO states for BO calculation
        """
        self.rand = random.random()
        for ist in range(self.mol.nst):
            if (ist == self.rstate):
                continue
            if (self.rand > self.acc_prob[ist] and self.rand <= self.acc_prob[ist + 1]):
                self.l_hop = True
                self.rstate = ist
                bo_list[0] = self.rstate

        # Record hopping event
        if (self.rstate != self.rstate_old):
            if (self.force_hop):
                self.event["HOP"].append(f"Accept hopping: force hop {self.rstate_old} -> {self.rstate}")
            else:
                self.event["HOP"].append(f"Accept hopping: hop {self.rstate_old} -> {self.rstate}")

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

    def append_sigma(self):
        """ Routine to append sigma values when single float number is provided
        """
        # Create a list from single float number
        if (isinstance(self.sigma, float)):
            sigma = self.sigma
            self.sigma = self.aux.nat * [sigma]

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


