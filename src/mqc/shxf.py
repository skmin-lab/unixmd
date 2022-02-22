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
    def __init__(self, molecule):
        # Initialize auxiliary molecule
        self.nat = molecule.nat_qm
        self.ndim = molecule.ndim
        self.symbols = np.copy(molecule.symbols[0:molecule.nat_qm])

        self.mass = np.copy(molecule.mass[0:molecule.nat_qm])

        self.pos = np.zeros((molecule.nst, self.nat, self.ndim))
        self.vel = np.zeros((molecule.nst, self.nat, self.ndim))
        self.vel_old = np.copy(self.vel)


class SHXF(MQC):
    """ Class for DISH-XF dynamics

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
        :param string aux_econs_viol: How to treat trajectories violating the total energy conservation
        :param string unit_dt: Unit of time interval
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=20, \
        elec_object="density", propagator="rk4", l_print_dm=True, l_adj_nac=True, hop_rescale="augment", \
        hop_reject="reverse", rho_threshold=0.01, sigma=None, init_coef=None, \
        l_econs_state=True, aux_econs_viol="fix", unit_dt="fs", out_freq=1, verbosity=0):
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

        self.hop_rescale = hop_rescale.lower()
        if not (self.hop_rescale in ["energy", "velocity", "momentum", "augment"]):
            error_message = "Invalid rescaling method for accepted hop!"
            error_vars = f"hop_rescale = {self.hop_rescale}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.hop_reject = hop_reject.lower()
        if not (self.hop_reject in ["keep", "reverse"]):
            error_message = "Invalid rescaling method for frustrated hop!"
            error_vars = f"hop_reject = {self.hop_reject}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Check error for incompatible cases
        if (self.mol.l_nacme):
            # No analytical nonadiabatic couplings exist
            if (self.hop_rescale in ["velocity", "momentum", "augment"]):
                error_message = "NACVs are not available with current QM object, only isotropic rescaling is possible!"
                error_vars = f"hop_rescale = {self.hop_rescale}"
                raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
            # TODO : This error will be used after adding the 'flip' option for hop_reject
#            if (self.hop_reject == "reverse"):
#                error_message = "NACVs are not available with current QM object, only keep rescaling is possible!"
#                error_vars = f"hop_reject = {self.hop_reject}"
#                raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Initialize XF related variables
        self.force_hop = False
        self.l_econs_state = l_econs_state
        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst
        self.l_fix = [False] * self.mol.nst
        self.l_collapse = False
        self.rho_threshold = rho_threshold
        self.aux_econs_viol = aux_econs_viol

        if not (self.aux_econs_viol in ["fix", "collapse"]):
            error_message = "Invalid method to treat auxiliary trajectories that violate the total energy conservation!"
            error_vars = f"aux_econs_viol = {self.aux_econs_viol}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

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

        # Initialize event to print
        self.event = {"HOP": [], "DECO": []}

    def run(self, qm, mm=None, output_dir="./", l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to decoherence-induced surface hopping dynamics

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

            self.hop_prob()
            self.hop_check(bo_list)
            self.evaluate_hop(bo_list)
            if (self.l_hop):
                if (qm.re_calc):
                    qm.get_data(self.mol, base_dir, bo_list, self.dt, self.istep, calc_force_only=True)
                if (self.mol.l_qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, self.istep, calc_force_only=True)

            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()
            if (self.l_collapse):
                self.check_decoherence()
                self.check_coherence()

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

            self.hop_prob()
            self.hop_check(bo_list)
            self.evaluate_hop(bo_list)
            if (self.l_hop):
                if (qm.re_calc):
                    qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=True)
                if (self.mol.l_qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, istep, calc_force_only=True)

            if (self.thermo != None):
                self.thermo.run(self)

            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()
            if (self.l_collapse):
                self.check_decoherence()
                self.check_coherence()

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

    def evaluate_hop(self, bo_list):
        """ Routine to evaluate hopping and velocity rescaling

            :param integer,list bo_list: List of BO states for BO calculation
            :param integer istep: Current MD step
        """
        if (self.l_hop):
            # Calculate potential difference between hopping states
            pot_diff = self.mol.states[self.rstate].energy - self.mol.states[self.rstate_old].energy

            # Solve quadratic equation for scaling factor of velocities
            a = 1.
            b = 1.
            det = 1.
            if (self.hop_rescale == "velocity"):
                a = np.sum(self.mol.mass[0:self.mol.nat_qm] * np.sum(self.mol.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                b = 2. * np.sum(self.mol.mass[0:self.mol.nat_qm] * np.sum(self.mol.nac[self.rstate_old, self.rstate] \
                    * self.mol.vel[0:self.mol.nat_qm], axis=1))
                c = 2. * pot_diff
                det = b ** 2. - 4. * a * c
            elif (self.hop_rescale == "momentum"):
                a = np.sum(1. / self.mol.mass[0:self.mol.nat_qm] * np.sum(self.mol.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                b = 2. * np.sum(np.sum(self.mol.nac[self.rstate_old, self.rstate] * self.mol.vel[0:self.mol.nat_qm], axis=1))
                c = 2. * pot_diff
                det = b ** 2. - 4. * a * c
            elif (self.hop_rescale == "augment"):
                a = np.sum(1. / self.mol.mass[0:self.mol.nat_qm] * np.sum(self.mol.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                b = 2. * np.sum(np.sum(self.mol.nac[self.rstate_old, self.rstate] * self.mol.vel[0:self.mol.nat_qm], axis=1))
                c = 2. * pot_diff
                det = b ** 2. - 4. * a * c

            # Default: hopping is allowed
            self.l_reject = False

            # Velocities cannot be adjusted when zero kinetic energy is given
            if (self.hop_rescale == "energy" and self.mol.ekin_qm < eps):
                self.l_reject = True
            # Clasically forbidden hop due to lack of kinetic energy
            if (self.mol.ekin_qm < pot_diff):
                self.l_reject = True
            # Kinetic energy is enough, but there is no solution for scaling factor
            if (det < 0.):
                self.l_reject = True
            # When kinetic energy is enough, velocities are always rescaled in 'augment' case
            if (self.hop_rescale == "augment" and self.mol.ekin_qm > pot_diff):
                self.l_reject = False

            if (self.l_reject):
                # Record event for frustrated hop
                if (self.mol.ekin_qm < pot_diff):
                    self.event["HOP"].append(f"Reject hopping: smaller kinetic energy than potential energy difference between {self.rstate} and {self.rstate_old}")
                # Set scaling constant with respect to 'hop_reject'
                if (self.hop_reject == "keep"):
                    self.event["HOP"].append("Reject hopping: no solution to find rescale factor, velocity is not changed")
                elif (self.hop_reject == "reverse"):
                    # x = - 1 when 'hop_rescale' is 'energy', otherwise x = - b / a
                    self.event["HOP"].append("Reject hopping: no solution to find rescale factor, velocity is reversed along coupling direction")
                    x = - b / a
                # Recover old running state
                self.l_hop = False

                if (self.force_hop):
                    self.event["HOP"].append(f"Collapse density: reset the density according to the current state {self.rstate_old}")
                    self.set_decoherence(self.rstate_old)

                self.force_hop = False

                self.rstate = self.rstate_old
                bo_list[0] = self.rstate
            else:
                if (self.hop_rescale == "energy" or (det < 0. and self.hop_rescale == "augment")):
                    if (det < 0.):
                        self.event["HOP"].append("Accept hopping: no solution to find rescale factor, but velocity is simply rescaled")
                    x = np.sqrt(1. - pot_diff / self.mol.ekin_qm)
                else:
                    if (b < 0.):
                        x = 0.5 * (- b - np.sqrt(det)) / a
                    else:
                        x = 0.5 * (- b + np.sqrt(det)) / a

            # Rescale velocities for QM atoms
            if (not (self.hop_reject == "keep" and self.l_reject)):
                if (self.hop_rescale == "energy"):
                    self.mol.vel[0:self.mol.nat_qm] *= x

                elif (self.hop_rescale == "velocity"):
                    self.mol.vel[0:self.mol.nat_qm] += x * self.mol.nac[self.rstate_old, self.rstate]

                elif (self.hop_rescale == "momentum"):
                    self.mol.vel[0:self.mol.nat_qm] += x * self.mol.nac[self.rstate_old, self.rstate] / \
                        self.mol.mass[0:self.mol.nat_qm].reshape((-1, 1))

                elif (self.hop_rescale == "augment"):
                    if (det > 0. or self.mol.ekin_qm < pot_diff):
                        self.mol.vel[0:self.mol.nat_qm] += x * self.mol.nac[self.rstate_old, self.rstate] / \
                            self.mol.mass[0:self.mol.nat_qm].reshape((-1, 1))
                    else:
                        self.mol.vel[0:self.mol.nat_qm] *= x

            # Update kinetic energy
            self.mol.update_kinetic()

        # Record hopping event
        if (self.rstate != self.rstate_old):
            if (self.force_hop):
                self.event["HOP"].append(f"Accept hopping: force hop {self.rstate_old} -> {self.rstate}")
            else:
                self.event["HOP"].append(f"Accept hopping: hop {self.rstate_old} -> {self.rstate}")

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        self.rforce = np.copy(self.mol.states[self.rstate].force)

    def update_energy(self):
        """ Routine to update the energy of molecules in surface hopping dynamics
        """
        # Update kinetic energy
        self.mol.update_kinetic()
        self.mol.epot = self.mol.states[self.rstate].energy
        self.mol.etot = self.mol.epot + self.mol.ekin

    def check_decoherence(self):
        """ Routine to check if the electronic state is decohered
        """
        if (self.l_hop):
            if (True in self.l_coh):
                self.event["DECO"].append(f"Destroy auxiliary trajectories: hopping occurs")
            self.l_coh = [False] * self.mol.nst
            self.l_first = [False] * self.mol.nst
            self.l_fix = [False] * self.mol.nst
        else:
            for ist in range(self.mol.nst):
                if (self.l_coh[ist]):
                    rho = self.mol.rho.real[ist, ist]
                    if (rho > self.upper_th):
                        self.set_decoherence(ist)
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
        self.mol.rho = np.zeros((self.mol.nst, self.mol.nst), dtype=np.complex128)
        self.mol.rho[one_st, one_st] = 1. + 0.j

        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst
        self.l_fix = [False] * self.mol.nst

        self.event["DECO"].append(f"Destroy auxiliary trajectories: decohered to {one_st} state")

        if (self.elec_object == "coefficient"):
            for ist in range(self.mol.nst):
                if (ist == one_st):
                    self.mol.states[ist].coef /= np.absolute(self.mol.states[ist].coef).real
                else:
                    self.mol.states[ist].coef = 0. + 0.j

    def aux_propagator(self):
        """ Routine to propagate auxiliary molecule
        """
        # Get auxiliary position
        for ist in range(self.mol.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    self.aux.pos[ist] = self.mol.pos[0:self.aux.nat]
                else:
                    if (ist == self.rstate):
                        self.aux.pos[ist] = self.mol.pos[0:self.aux.nat]
                    else:
                        self.aux.pos[ist] += self.aux.vel[ist] * self.dt

        self.pos_0 = np.copy(self.aux.pos[self.rstate])

        # Get auxiliary velocity
        self.l_collapse = False
        self.aux.vel_old = np.copy(self.aux.vel)
        for ist in range(self.mol.nst):
            # Calculate propagation factor alpha
            if (self.l_coh[ist]):
                if (self.l_fix[ist]):
                    alpha = 0.
                else:
                    if (ist == self.rstate):
                        alpha = self.mol.ekin_qm
                    else:
                        if (self.l_first[ist]):
                            alpha = self.mol.ekin_qm
                            if (self.l_econs_state):
                                alpha += self.mol.states[self.rstate].energy - self.mol.states[ist].energy
                        else:
                            ekin_old = np.sum(0.5 * self.aux.mass * np.sum(self.aux.vel_old[ist] ** 2, axis=1))
                            alpha = ekin_old + self.mol.states[ist].energy_old - self.mol.states[ist].energy
                    if (alpha < 0.):
                        alpha = 0.
                        if (self.aux_econs_viol == "fix"):
                            self.l_fix[ist] = True
                            self.event["DECO"].append(f"Energy conservation violated, the auxiliary trajectory on state {ist} is fixed.")
                        elif (self.aux_econs_viol == "collapse"):
                            self.l_collapse = True
                            self.collapse(ist)
                            self.event["DECO"].append(f"Energy conservation violated, collapse the {ist} state coefficient/density to zero.")

                # Calculate auxiliary velocity from alpha
                alpha /= self.mol.ekin_qm
                self.aux.vel[ist] = self.mol.vel[0:self.aux.nat] * np.sqrt(alpha)

    def collapse(self, cstate):
        """ Routine to collapse coefficient/density of a state to zero
        """
        fac = 1. - self.mol.rho.real[cstate, cstate]

        if (self.elec_object == "coefficient"):
            for ist in range(self.mol.nst):
                if (ist == cstate):
                    self.mol.states[ist].coef = 0. + 0.j
                else:
                    self.mol.states[ist].coef /= np.sqrt(fac)

        self.mol.rho[cstate,:] = 0. + 0.j
        self.mol.rho[:,cstate] = 0. + 0.j
        self.mol.rho /= fac
         
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
        # Write auxiliary trajectories
        if (self.verbosity >= 2 and True in self.l_coh):
            # Write quantum momenta
            tmp = f'{self.aux.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Momentum (au)' + \
                "".join(["\n" + f'{self.aux.symbols[iat]:5s}' + \
                "".join([f'{self.qmom[iat, isp]:15.8f}' for isp in range(self.aux.ndim)]) for iat in range(self.aux.nat)])
            typewriter(tmp, unixmd_dir, f"QMOM", "a")

            # Write auxiliary variables
            for ist in range(self.mol.nst):
                if (self.l_coh[ist]):
                    # Write auxiliary phase
                    tmp = f'{self.aux.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Phase (au)' + \
                        "".join(["\n" + f'{self.aux.symbols[iat]:5s}' + \
                        "".join([f'{self.phase[ist, iat, isp]:15.8f}' for isp in range(self.aux.ndim)]) for iat in range(self.aux.nat)])
                    typewriter(tmp, unixmd_dir, f"AUX_PHASE_{ist}", "a")

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
        INIT = f" #INFO{'STEP':>8s}{'State':>7s}{'Kinetic(H)':>14s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'Norm.':>8s}"
        dynamics_step_info += INIT

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}{'Rand.':>11s}{'Acc. Hopping Prob.':>28s}"
            dynamics_step_info += "\n" + DEBUG1

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
        INFO = f" INFO{istep + 1:>9d}{self.rstate:>5d}"
        INFO += f"{self.mol.ekin:16.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        INFO += f"{norm:11.5f}"
        print (INFO, flush=True)

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            DEBUG1 += f"{self.rand:11.5f}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{self.acc_prob[ist]:12.5f} ({self.rstate}->{ist})"
            print (DEBUG1, flush=True)

        # Print event in SHXF
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event["HOP"] = []
        self.event["DECO"] = []


