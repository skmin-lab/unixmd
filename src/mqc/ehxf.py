from __future__ import division
from build.el_propagator_xf import el_run
from mqc.mqc import MQC
from misc import eps, au_to_K, au_to_A, call_name, typewriter
import random, os, shutil, textwrap
import numpy as np
import pickle

class Auxiliary_Molecule(object):
    """ Class for auxiliary molecule that is used for the calculation of decoherence term

        :param object molecule: Molecule object
    """
    def __init__(self, molecule, l_xf1d):
        # Initialize auxiliary molecule
        if (l_xf1d):

            self.nat = 1
            self.ndim = 1
            self.symbols = ['XX']

            self.mass = np.zeros((self.nat))
            self.mass[0] = 1. / np.sum(1. / molecule.mass[0:molecule.nat_qm])

        else:

            self.nat = molecule.nat_qm
            self.ndim = molecule.ndim
            self.symbols = np.copy(molecule.symbols[0:molecule.nat_qm])

            self.mass = np.copy(molecule.mass[0:molecule.nat_qm])

        self.pos = np.zeros((molecule.nst, self.nat, self.ndim))
        self.vel = np.zeros((molecule.nst, self.nat, self.ndim))
        self.vel_old = np.copy(self.vel)


class EhXF(MQC):
    """ Class for Ehrenfest-XF dynamics

        :param object molecule: Molecule object
        :param object thermostat: Thermostat object
        :param integer istate: Initial state
        :param double dt: Time interval
        :param integer nsteps: Total step of nuclear propagation
        :param integer nesteps: Total step of electronic propagation
        :param string elec_object: Electronic equation of motions
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: Logical to print BO population and coherence
        :param boolean l_adj_nac: Adjust nonadiabatic coupling to align the phases
        :param string hop_rescale: Velocity rescaling method after successful hop
        :param string hop_reject: Velocity rescaling method after frustrated hop
        :param double rho_threshold: Electronic density threshold for decoherence term calculation
        :param sigma: Width of nuclear wave packet of auxiliary trajectory
        :type sigma: double or double,list
        :param init_coef: Initial BO coefficient
        :type init_coef: double, list or complex, list
        :param boolean l_econs_state: Logical to use identical total energies for all auxiliary trajectories
        :param string unit_dt: Unit of time interval
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=20, \
        elec_object="density", propagator="rk4", l_print_dm=True, l_adj_nac=True, l_asymp=True, \
        rho_threshold=0.01, sigma=None, l_xf_force=True, l_xf1d=False, init_coef=None, l_afssh=False, l_td_sigma=False,\
        l_econs_state=True, l_econs_phase=True, unit_dt="fs", out_freq=1, verbosity=0):
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

        # Initialize XF related variables
        self.force_hop = False
        self.l_econs_state = l_econs_state
        self.l_econs_phase = l_econs_phase
        self.l_xf_force = l_xf_force 
        self.l_xf1d = l_xf1d
        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst
        self.rho_threshold = rho_threshold
        self.sigma = sigma
        self.aux = Auxiliary_Molecule(self.mol, self.l_xf1d)
        self.pos_0 = np.zeros((self.aux.nat, self.aux.ndim))
        self.phase = np.array(np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim)))
        self.qmom = np.zeros((self.mol.nst, self.mol.nst, self.aux.nat, self.aux.ndim))
        self.df = np.zeros((self.mol.nst, self.mol.nst, self.aux.nat, self.aux.ndim))
        self.k_lk = np.zeros((self.mol.nst, self.mol.nst))
        self.l_asymp = l_asymp        

        # Initialize AFSSH related variables
        self.l_afssh = l_afssh
        self.l_td_sigma = l_td_sigma
        self.dR = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))
        self.dP = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))
        self.dF = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))
        self.dF_old = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))

        if (self.l_afssh and self.l_xf1d):
            raise ValueError (f"( {self.md_type}.{call_name()} ) EhXF1D is not implemented for AFSSH auxiliary trajectory propation method! {self.l_xf1d}")
        if (self.l_td_sigma and self.l_xf1d):
            raise ValueError (f"( {self.md_type}.{call_name()} ) EhXF1D is not implemented for TD sigma! {self.l_xf1d}")
        
        if (self.sigma == None):
            if (not self.l_td_sigma):
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
            if (self.l_xf1d):
                error_message = "Sigma must be float, not list in XF-1D scheme!"
                error_vars = f"sigma = {self.sigma}"
                raise TypeError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            if (not self.l_td_sigma):
                error_message = "Type of sigma must be float or list consisting of float!"
                error_vars = f"sigma = {self.sigma}"
                raise TypeError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
        
        self.upper_th = 1. - self.rho_threshold
        self.lower_th = self.rho_threshold

        # Debug variables
        self.dotpopdec = np.zeros(self.mol.nst)
        self.dotpopnac = np.zeros(self.mol.nst)

        # Initialize event to print
        self.event = {"HOP": [], "DECO": []}

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
            if (not self.mol.l_nacme):
                self.mol.get_nacme()

            self.hop_prob(self.istep)
            self.hop_check(bo_list)
            if (qm.re_calc and self.l_hop):
                qm.get_data(self.mol, base_dir, bo_list, self.dt, self.istep, calc_force_only=True)
                if (self.mol.l_qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, self.istep, calc_force_only=True)

            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()
            self.calc_k_lk()
            if (self.l_xf_force and not self.l_xf1d):
                self.calc_xf_force()

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

            if (not self.mol.l_nacme and self.l_adj_nac):
                self.mol.adjust_nac()

            self.calculate_force()
            self.cl_update_velocity()

            if (not self.mol.l_nacme):
                self.mol.get_nacme()

            el_run(self)

            self.hop_prob(istep)
            self.hop_check(bo_list)
            if (qm.re_calc and self.l_hop):
                qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=True)
                if (self.mol.l_qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, istep, calc_force_only=True)

            if (self.thermo != None):
                self.thermo.run(self)

#            if (self.l_xf_force and not self.l_xf1d):
#                self.force_econs()
            
            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()
            self.calc_k_lk()
            if (self.l_xf_force and not self.l_xf1d):
                self.calc_xf_force()
            
            self.update_energy() # TODO: update_energy is always executed twice. Use flag to avoid the duplication
            
            if ((istep + 1) % self.out_freq == 0):
                self.write_md_output(unixmd_dir, istep)
            if ((istep + 1) % self.out_freq == 0 or len(self.event["HOP"]) > 0 or len(self.event["DECO"]) > 0):
                self.print_step(istep)
            if (istep == self.nsteps - 1):
                self.write_final_xyz(unixmd_dir, istep)

            self.fstep = istep
            restart_file = os.path.join(base_dir, "RESTART.bin")
            with open(restart_file, 'wb') as f:
                pickle.dump({'qm':qm, 'md':self}, f)
            
            det = self.mol.pos[0, 0] * self.mol.vel[0, 0]
            if (self.l_asymp and det > 0. and np.abs(self.mol.pos[0, 0]) > np.abs(20.)):
                print(f"Trajectory reached asymptotic region")
                break

        # Delete scratch directory
        if (not l_save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            if (self.mol.l_qmmm and mm != None):
                tmp_dir = os.path.join(unixmd_dir, "scr_mm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def hop_prob(self, istep):
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
        """ Routine to calculate the forces
        """
        self.rforce = np.zeros((self.mol.nat, self.mol.ndim))

        for ist, istate in enumerate(self.mol.states):
            self.rforce += istate.force * self.mol.rho.real[ist, ist]

        for ist in range(self.mol.nst):
            for jst in range(ist + 1, self.mol.nst):
                self.rforce += 2. * self.mol.nac[ist, jst] * self.mol.rho.real[ist, jst] \
                    * (self.mol.states[ist].energy - self.mol.states[jst].energy)
        if (self.l_xf_force):
            self.rforce += self.xf_force

    def update_energy(self):
        """ Routine to update the energy of molecules in Ehrenfest dynamics
        """
        # Update kinetic energy
        self.mol.update_kinetic()
        self.mol.epot = 0.
        for ist, istate in enumerate(self.mol.states):
            self.mol.epot += self.mol.rho.real[ist, ist] * self.mol.states[ist].energy
        self.mol.etot = self.mol.epot + self.mol.ekin

    def check_decoherence(self):
        """ Routine to check if the electronic state is decohered
        """
        if (self.l_hop):
            if (True in self.l_coh):
                self.event["DECO"].append(f"Destroy auxiliary trajectories: hopping occurs")
            self.l_coh = [False] * self.mol.nst
            self.l_first = [False] * self.mol.nst
        else:
            for ist in range(self.mol.nst):
                if (self.l_coh[ist]):
                    rho = self.mol.rho.real[ist, ist]
                    if (rho > self.upper_th):
                        self.set_decoherence(ist)
                        self.event["DECO"].append(f"Destroy auxiliary trajectories: decohered to {ist} state")
                        return

    def check_coherence(self):
        """ Routine to check coherence among BO states
        """
        count = 0
        tmp_st = ""
        for ist in range(self.mol.nst):
            rho = self.mol.rho.real[ist, ist]
            if (rho > self.upper_th or rho < self.lower_th):
                self.l_coh[ist] = False
            else:
                if (self.l_coh[ist]):
                    self.l_first[ist] = False
                else:
                    self.l_first[ist] = True
                    tmp_st += f"{ist}, "
                self.l_coh[ist] = True
                count += 1

        if (count < 2):
            self.l_coh = [False] * self.mol.nst
            self.l_first = [False] * self.mol.nst
            tmp_st = ""

        if (len(tmp_st) >= 1):
            tmp_st = tmp_st.rstrip(', ')
            self.event["DECO"].append(f"Generate auxiliary trajectory on {tmp_st} state")

    def set_decoherence(self, one_st):
        """ Routine to reset coefficient/density if the state is decohered

            :param integer one_st: State index that its population is one
        """
        self.phase = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))
        self.qmom = np.zeros((self.mol.nst, self.mol.nst, self.aux.nat, self.aux.ndim))
        self.df = np.zeros((self.mol.nst, self.mol.nst, self.aux.nat, self.aux.ndim))
        self.k_lk = np.zeros((self.mol.nst, self.mol.nst))
        self.mol.rho = np.zeros((self.mol.nst, self.mol.nst), dtype=np.complex_)
        self.mol.rho[one_st, one_st] = 1. + 0.j

        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst

        if (self.elec_object == "coefficient"):
            for ist in range(self.mol.nst):
                if (ist == one_st):
                    self.mol.states[ist].coef /= np.absolute(self.mol.states[ist].coef).real
                else:
                    self.mol.states[ist].coef = 0. + 0.j

        epot_new = 0.
        for ist, istate in enumerate(self.mol.states):
            epot_new += self.mol.rho.real[ist, ist] * self.mol.states[ist].energy
        alpha = self.mol.etot - epot_new

        if (alpha >= 0.):
            alpha /= self.mol.ekin_qm
            alpha = np.sqrt(alpha)
        else:
            alpha = 0.

        self.mol.vel *= alpha



    def aux_propagator(self):
        """ Routine to propagate auxiliary molecule
        """
        if (self.l_afssh):
        
            # Use AFSSH EOMs
           
            self.aux.vel_old = np.copy(self.aux.vel)
            rst = self.rstate
            rho_rst = self.mol.rho.real[rst, rst]
           
            for ist in range(self.mol.nst):
                if (self.l_coh[ist]):
                    
                    rho_ist = self.mol.rho.real[ist, ist]
                    self.dF_old[ist] =  np.copy(self.dF[ist])
                    self.dF[ist] =  rho_ist * (self.mol.states[ist].force - self.mol.states[rst].force)
                    
                    if (self.l_first[ist]):
                        self.dR[ist] = 0.
                        self.dP[ist] = 0.
                        self.aux.pos[ist] = self.mol.pos[0:self.aux.nat] 
                        self.aux.vel[ist] = self.mol.vel[0:self.aux.nat] 
                    else:
                        # Velocity Verlet
                        self.dP[ist] += 0.5 * self.dt * self.dF_old[ist]
                        self.dR[ist] += self.dt * self.dP[ist] / np.column_stack([self.aux.mass] * self.aux.ndim)
                        self.dP[ist] += 0.5 * self.dt * self.dF[ist]
                        self.aux.pos[ist] = self.mol.pos[0:self.aux.nat] + self.dR[ist] / rho_ist  
                        self.aux.vel[ist] = self.mol.vel[0:self.aux.nat] + self.dP[ist] / np.column_stack([self.aux.mass] * self.aux.ndim) / rho_ist 
#                        self.aux.pos[ist] = 0.0 + self.dR[ist] / rho_ist  
#                        self.aux.vel[ist] = 0.0 + self.dP[ist] / np.column_stack([self.aux.mass] * self.aux.ndim) / rho_ist 
                else:
                    self.dR[ist] = 0.
                    self.dP[ist] = 0.
                    self.dF[ist] = 0.
                    self.aux.pos[ist] = self.mol.pos[0:self.aux.nat]
                    self.aux.vel[ist] = self.mol.vel[0:self.aux.nat]

            
            self.pos_0 = np.copy(self.aux.pos[self.rstate])
            
        else:
            # Original scheme (Ha, J. Phys. Chem. Lett. 2018, 9, 1094)
            # Get auxiliary position
            for ist in range(self.mol.nst):
                if (self.l_coh[ist]):
                    if (self.l_first[ist]):
                        if (self.l_xf1d):
                            self.aux.pos[ist] = np.zeros((self.aux.nat, self.aux.ndim))
                        else:
                            self.aux.pos[ist] = self.mol.pos[0:self.aux.nat]
                    else:
                        self.aux.pos[ist] += self.aux.vel[ist] * self.dt
    
            self.pos_0 = np.copy(self.aux.pos[self.rstate])
    
            # Get auxiliary velocity
            self.aux.vel_old = np.copy(self.aux.vel)
            for ist in range(self.mol.nst):
                # Calculate propagation factor alpha
                if (self.l_coh[ist]):
                    if (self.l_first[ist]):
                        alpha = self.mol.ekin_qm
                        if (self.l_econs_state):
                            alpha += self.mol.epot - self.mol.states[ist].energy
                    else:
                        alpha = self.mol.ekin_qm + self.mol.epot - self.mol.states[ist].energy
                    if (alpha < 0.):
                        alpha = 0.
    
                    # Calculate auxiliary velocity from alpha
                    if (self.l_xf1d):
                        alpha /= 0.5 * self.aux.mass[0]
                        self.aux.vel[ist] = np.sqrt(alpha)
                    else:
                        alpha /= self.mol.ekin_qm
                        self.aux.vel[ist] = self.mol.vel[0:self.aux.nat] * np.sqrt(alpha)
        
        
        if (self.l_td_sigma):
            # Only two-state case is implemented..
            if (self.l_first[0]):
                self.sigma = np.ones((self.aux.nat, self.aux.ndim)) * 100000.
            else:
                for iat in range(self.aux.nat):
                    for isp in range(self.aux.ndim):
                        if ((np.abs(self.aux.vel[0, iat, isp] - self.aux.vel[1, iat, isp])) * self.aux.mass[iat] > eps):
                            self.sigma[iat, isp] = np.sqrt(0.5 * np.abs(\
                               (self.aux.pos[0, iat, isp] - self.aux.pos[1, iat, isp])/\
                               (self.aux.vel[0, iat, isp] - self.aux.vel[1, iat, isp]))\
                                / self.aux.mass[iat])
                        else:
                            self.sigma = np.ones((self.aux.nat, self.aux.ndim)) * 100000.

    def get_phase(self):
        """ Routine to calculate phase term
        """
        for ist in range(self.mol.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    self.phase[ist] = 0.
                else:
                    for iat in range(self.aux.nat):
                        self.phase[ist, iat] += self.aux.mass[iat] * \
                            (self.aux.vel[ist, iat] - self.aux.vel_old[ist, iat])

    def calc_k_lk(self):
        """ Routine to calculate decoherence term k_lk
        """ 
        
        self.k_lk = np.zeros((self.mol.nst, self.mol.nst))
        self.qmom = np.zeros((self.mol.nst, self.mol.nst, self.aux.nat, self.aux.ndim))

        if (self.l_econs_phase):
            self.df = np.zeros((self.mol.nst, self.mol.nst, self.aux.nat, self.aux.ndim))
            # Calculate P*v
            fac = np.sum(self.aux.mass * np.sum(self.mol.vel ** 2, axis=1))

        # Calculate state-pair quantum momentum and phase difference
        for ist in range(self.mol.nst):
            for jst in range(ist + 1, self.mol.nst):
                if (self.l_coh[ist] and self.l_coh[jst]):
                    # Calculate quantum momentum
                    for iat in range(self.aux.nat):
                        self.qmom[ist, jst, iat] = (self.pos_0[iat] - self.aux.pos[ist, iat]) * self.mol.rho.real[ist, ist]\
                                                 + (self.pos_0[iat] - self.aux.pos[jst, iat]) * self.mol.rho.real[jst, jst]
                        self.qmom[ist, jst, iat] /= self.aux.mass[iat] * self.sigma[iat] ** 2
                    
                    # Calculate phase difference
                    if (self.l_first[ist] or self.l_first[jst]):
                        self.df[ist, jst] = np.zeros((self.aux.nat, self.aux.ndim))
                    elif (self.l_econs_phase):
                        if (fac < eps):
                            alpha = 0.
                        else:
                            alpha = - (self.mol.states[ist].energy - self.mol.states[jst].energy) / fac
                        for iat in range(self.aux.nat):
                            self.df[ist, jst, iat] = alpha * self.aux.mass[iat] * self.mol.vel[iat]
                    else:    
                        for iat in range(self.aux.nat):
                            self.df[ist, jst, iat] += self.aux.mass[iat] * (self.aux.vel[ist, iat] - self.aux.vel_old[ist, iat])
                            self.df[ist, jst, iat] -= self.aux.mass[iat] * (self.aux.vel[jst, iat] - self.aux.vel_old[jst, iat])

                self.qmom[jst, ist] = np.copy(self.qmom[ist, jst])
                self.df[jst, ist] = - 1. * self.df[ist, jst]

        self.qmom /= 2. * (self.mol.nst - 1)

        # Calculate state-pair 
        for ist in range(self.mol.nst):
            for jst in range(ist + 1, self.mol.nst):
                if (self.l_coh[ist] and self.l_coh[jst]):
                    for iat in range(self.aux.nat):
                        self.k_lk[ist, jst] += np.sum(self.qmom[ist, jst, iat] * self.df[ist, jst, iat])
                self.k_lk[jst, ist] = - 1. * self.k_lk[ist, jst]
         
    def calc_xf_force(self):
        """ Routine to calculate nuclear force from XF term
        """ 
        # TODO: qmom calculation is duplicated. Address this issue.
#        if (self.verbosity < 2):
#            self.qmom = np.zeros((self.aux.nat, self.aux.ndim))
#            for ist in range(self.mol.nst):
#                if (self.l_coh[ist]):
#                    for iat in range(self.aux.nat):
#                        self.qmom[iat, :] += 0.5 * self.mol.rho.real[ist, ist] / self.aux.mass[iat] / self.sigma[iat] **2\
#                                           * (self.pos_0[iat, :] - self.aux.pos[ist, iat, :])

        self.xf_force = np.zeros((self.mol.nat, self.mol.ndim))
        for ist in range(self.mol.nst):
            for jst in range(self.mol.nst):
                if (self.l_coh[ist] and self.l_coh[jst]):
                    self.xf_force += self.mol.rho.real[ist, ist] * self.mol.rho.real[jst, jst] \
                                     * self.k_lk[ist, jst] * self.df[ist, jst]
    
    def force_econs(self):
        """ Routine to rescale velocity to force the total energy conservation
        """
        epot = 0.
        for ist, istate in enumerate(self.mol.states):
            epot += self.mol.rho.real[ist, ist] * self.mol.states[ist].energy
        ekin = np.sum(0.5 * self.mol.mass * np.sum(self.mol.vel ** 2, axis=1))
        alpha = self.mol.etot - epot
        
        if (alpha < 0.):
            alpha = 1.
            print(f" CAUTION: Negative Etot-Epot, velocity is not rescaled.")
        else:
            if (ekin < eps):
                alpha = 1.
                print(f" CAUTION: Too small kinetic before rescaling, velocity is not rescaled.")
            else:
                alpha /= ekin
                alpha = np.sqrt(alpha)
        
        self.mol.vel *= alpha

    def append_sigma(self):
        """ Routine to append sigma values when single float number is provided
        """
        # Create a list from single float number
        if (isinstance(self.sigma, float)):
            sigma = self.sigma
            self.sigma = np.array(self.aux.nat * [self.aux.ndim * [sigma]])
        if (isinstance(self.sigma, list)):
            sigma = []
            for sgm in self.sigma:
                sigma.append(self.aux.nat * [sgm])
            self.sigma = sigma[:]
        if (self.l_td_sigma):
            self.sigma = np.array(self.aux.nat * [self.aux.ndim * [0.0]])

    def write_md_output(self, unixmd_dir, istep):
        """ Write output files

            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, istep)

        # Write hopping-related quantities
        self.write_sh(unixmd_dir, istep)

        # Write time-derivative BO population
        self.write_dotpop(unixmd_dir, istep)

        # Write decoherence information
        self.write_dec(unixmd_dir, istep)

    def write_sh(self, unixmd_dir, istep):
        """ Write hopping-related quantities into files

            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        # Write SHSTATE file
        tmp = f'{istep + 1:9d}{"":14s}{self.rstate}'
        typewriter(tmp, unixmd_dir, "SHSTATE", "a")

        # Write SHPROB file
        tmp = f'{istep + 1:9d}' + "".join([f'{self.prob[ist]:15.8f}' for ist in range(self.mol.nst)])
        typewriter(tmp, unixmd_dir, "SHPROB", "a")

    def write_dotpop(self, unixmd_dir, istep):
        """ Write time-derivative BO population

            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        if (self.verbosity >= 1):
            # Write NAC term in DOTPOPNAC
            tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopnac])
            typewriter(tmp, unixmd_dir, "DOTPOPNAC", "a")

            # Write decoherence term in DOTPOPDEC
            tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopdec])
            typewriter(tmp, unixmd_dir, "DOTPOPDEC", "a")

    def write_dec(self, unixmd_dir, istep):
        """ Write XF-based decoherence information

            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        # Write auxiliary trajectories #CAUTION!! ONLY 0,1 pair is written. Only for test.
        if (self.verbosity >= 2 and True in self.l_coh):
            # Write quantum momenta
            tmp = f'{self.aux.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Momentum (au)' + \
                "".join(["\n" + f'{self.aux.symbols[iat]:5s}' + \
                "".join([f'{self.qmom[0, 1, iat, isp]:15.8f}' for isp in range(self.aux.ndim)]) for iat in range(self.aux.nat)])
            typewriter(tmp, unixmd_dir, f"QMOM", "a")

            # Write auxiliary variables
            for ist in range(1):
                if (self.l_coh[ist]):
                    # Write auxiliary phase
                    tmp = f'{self.aux.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Phase (au)' + \
                        "".join(["\n" + f'{self.aux.symbols[iat]:5s}' + \
                        "".join([f'{self.df[0, 1, iat, isp]:15.8f}' for isp in range(self.aux.ndim)]) for iat in range(self.aux.nat)])
                    typewriter(tmp, unixmd_dir, f"AUX_PHASE_{ist}", "a")

            for ist in range(self.mol.nst):
                if (self.l_coh[ist]):
                    # Write auxiliary trajectory movie files
                    tmp = f'{self.aux.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)' + \
                        "".join(["\n" + f'{self.aux.symbols[iat]:5s}' + \
                        "".join([f'{self.aux.pos[ist, iat, isp] * au_to_A:15.8f}' for isp in range(self.aux.ndim)]) + \
                        "".join([f"{self.aux.vel[ist, iat, isp]:15.8f}" for isp in range(self.aux.ndim)]) for iat in range(self.aux.nat)])
                    typewriter(tmp, unixmd_dir, f"AUX_MOVIE_{ist}.xyz", "a")

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
        INIT = f" #INFO{'STEP':>8s}{'State':>7s}{'Max. Prob.':>14s}{'Rand.':>12s}{'Kinetic(H)':>15s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'Norm.':>8s}"
        dynamics_step_info += INIT

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        # Print DEBUG2 for each step
        if (self.verbosity >= 2):
            DEBUG2 = f" #DEBUG2{'STEP':>6s}{'Acc. Hopping Prob.':>22s}"
            dynamics_step_info += "\n" + DEBUG2
        
        # Print DEBUG3 for each step
        if (self.verbosity >= 3):
            self.dE = 0.
            DEBUG3 = f" #DEBUG3{'STEP':>6s}{'dE_tot/dt':>13s}{'dE':>6s}{'Feh':>7s}{'Fxf':>7s}"
            dynamics_step_info += "\n" + DEBUG3
        
        # Print DEBUGL for each step
        if (self.l_td_sigma and self.verbosity >= 3):
            DEBUGL = f" #DEBUGL{'STEP':>6s}{'Sigma':>9s}{'dR0':>7s}{'dR1':>7s}{'dP0':>7s}{'dP1':>7s}"
            dynamics_step_info += "\n" + DEBUGL

        print (dynamics_step_info, flush=True)

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer istep: Current MD step
        """
        if (istep == -1):
            max_prob = 0.
            hstate = self.rstate
        else:
            max_prob = max(self.prob)
            hstate = np.where(self.prob == max_prob)[0][0]

        ctemp = self.mol.ekin * 2. / float(self.mol.ndof) * au_to_K
        norm = 0.
        for ist in range(self.mol.nst):
            norm += self.mol.rho.real[ist, ist]
 
        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d}{self.rstate:>5d}{max_prob:11.5f} ({self.rstate}->{hstate}){self.rand:11.5f}"
        INFO += f"{self.mol.ekin:14.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        INFO += f"{norm:11.5f}"
        print (INFO, flush=True)

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{self.mol.states[ist].energy:17.8f} "
            print (DEBUG1, flush=True)

        # Print DEBUG2 for each step
        if (self.verbosity >= 2):
            DEBUG2 = f" DEBUG2{istep + 1:>7d}"
            for ist in range(self.mol.nst):
                DEBUG2 += f"{self.acc_prob[ist]:12.5f} ({self.rstate}->{ist})"
            print (DEBUG2, flush=True)
        
        # Print DEBUG3 for each step
        if (self.verbosity >= 3):
            self.calculate_force()
            Feh = self.rforce[0, 0] - self.xf_force[0, 0]
            dotpopdec = np.zeros((self.mol.nst))
            for ist in range(self.mol.nst):
                for mst in range(self.mol.nst):
                    dotpopdec[ist] -= 2.0 * self.k_lk[mst, ist] * self.mol.rho.real[ist, ist] * self.mol.rho.real[mst, mst]
            
            self.dEdt = 0.
            for ist in range(self.mol.nst):
                self.dEdt += self.mol.states[ist].energy * dotpopdec[ist]
                self.dEdt += self.mol.rho.real[ist, ist] * np.sum(self.xf_force * self.mol.vel)
            self.dE += self.dEdt
            DEBUG3 = f" DEBUG3{istep + 1:>7d}"
            DEBUG3 += f"{self.dEdt:14.7f}{self.dE:14.7f}{Feh:14.7f}{self.xf_force[0, 0]:14.7f}" 
            print (DEBUG3, flush=True)
        
        # Print DEBUGL for each step
        if (self.l_td_sigma and self.verbosity >= 3):
            DEBUGL = f" DEBUGL{istep + 1:>7d}"
            DEBUGL += f"{self.sigma[0][0]:12.5f}{self.dR[0, 0, 0]:12.5f}{self.dR[1, 0, 0]:12.5f}{self.dP[0, 0, 0]:12.5f}{self.dP[1, 0, 0]:12.5f}"
            print (DEBUGL, flush=True)

        # Print event in EhXF
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event["HOP"] = []
        self.event["DECO"] = []
