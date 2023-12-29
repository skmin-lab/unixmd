from __future__ import division
from build_qed.el_propagator_xf import el_run
from mqc_qed.mqc import MQC_QED
from misc import eps, au_to_K, au_to_A, call_name, typewriter
import random, os, shutil, textwrap
import numpy as np
import pickle

class Auxiliary_Molecule(object):
    """ Class for auxiliary molecule that is used for the calculation of decoherence term

        :param object polariton: Polariton object
    """
    def __init__(self, polariton):
        # Initialize auxiliary molecule
        self.nat = polariton.nat_qm
        self.ndim = polariton.ndim
        self.symbols = np.copy(polariton.symbols[0:polariton.nat_qm])

        self.mass = np.copy(polariton.mass[0:polariton.nat_qm])

        self.pos = np.zeros((polariton.pst, self.nat, self.ndim))
        self.vel = np.zeros((polariton.pst, self.nat, self.ndim))
        self.vel_old = np.copy(self.vel)


class SHXF(MQC_QED):
    """ Class for DISH-XF dynamics coupled to confined cavity mode

        :param object polariton: Polariton object
        :param object thermostat: Thermostat object
        :param integer istate: Initial state
        :param double dt: Time interval
        :param integer nsteps: Total step of nuclear propagation
        :param integer nesteps: Total step of electronic propagation
        :param string elec_object: Electronic equation of motions
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: Logical to print BO population and coherence
        :param boolean l_adj_nac: Adjust nonadiabatic coupling to align the phases
        :param boolean l_adj_tdp: Adjust transition dipole moments to align the phases
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
    def __init__(self, polariton, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=20, \
        elec_object="density", propagator="rk4", l_print_dm=True, l_adj_nac=True, l_adj_tdp=True, \
        hop_rescale="augment", hop_reject="reverse", rho_threshold=0.01, sigma=None, init_coef=None, \
        l_econs_state=True, aux_econs_viol="fix", unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        super().__init__(polariton, thermostat, istate, dt, nsteps, nesteps, elec_object, \
            propagator, l_print_dm, l_adj_nac, l_adj_tdp, init_coef, unit_dt, out_freq, verbosity)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.rand = 0.
        self.prob = np.zeros(self.pol.pst)
        self.acc_prob = np.zeros(self.pol.pst + 1)

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
        if (self.pol.l_nacme):
            # No analytical nonadiabatic couplings exist
            if (self.hop_rescale in ["velocity", "momentum", "augment"]):
                error_message = "pNACVs are not available with current QED object, only isotropic rescaling is possible!"
                error_vars = f"hop_rescale = {self.hop_rescale}"
                raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
            if (self.hop_reject == "reverse"):
                error_message = "pNACVs are not available with current QED object, only keep rescaling is possible!"
                error_vars = f"hop_reject = {self.hop_reject}"
                raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Initialize XF related variables
        self.force_hop = False
        self.l_econs_state = l_econs_state
        self.l_coh = [False] * self.pol.pst
        self.l_first = [False] * self.pol.pst
        self.l_fix = [False] * self.pol.pst
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
            if (len(self.sigma) != self.pol.nat_qm):
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
        self.aux = Auxiliary_Molecule(self.pol)
        self.pos_0 = np.zeros((self.aux.nat, self.aux.ndim))
        self.phase = np.zeros((self.pol.pst, self.aux.nat, self.aux.ndim))

        # Debug variables
        self.dotpopdec_d = np.zeros(self.pol.pst)
        self.dotpopnac_d = np.zeros(self.pol.pst)
        self.qmom = np.zeros((self.aux.nat, self.aux.ndim))

        # Initialize event to print
        self.event = {"HOP": [], "DECO": []}

    def run(self, qed, qm, mm=None, output_dir="./", l_save_qed_log=False, l_save_qm_log=False, \
        l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to decoherence-induced surface hopping dynamics

            :param object qed: QED object containing cavity-molecule interaction
            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string output_dir: Name of directory where outputs to be saved.
            :param boolean l_save_qed_log: Logical for saving QED calculation log
            :param boolean l_save_qm_log: Logical for saving QM calculation log
            :param boolean l_save_mm_log: Logical for saving MM calculation log
            :param boolean l_save_scr: Logical for saving scratch directory
            :param string restart: Option for controlling dynamics restarting
        """
        # Initialize PyUNIxMD
        base_dir, unixmd_dir, qed_log_dir, qm_log_dir, mm_log_dir = \
            self.run_init(qed, qm, mm, output_dir, l_save_qed_log, l_save_qm_log, l_save_mm_log, l_save_scr, restart)
        bo_list = [ist for ist in range(self.pol.nst)]
        pol_list = [self.rstate]
        qm.calc_coupling = True
        qm.calc_tdp = True
        qm.calc_tdp_grad = False
        if (qed.force_level == "tdp"):
            qm.calc_tdp_grad = True
        self.print_init(qed, qm, mm, restart)

        if (restart == None):
            # Initialize decoherence variables
            self.append_sigma()

            # Calculate initial input geometry at t = 0.0 s
            self.istep = -1
            self.pol.reset_bo(qm.calc_coupling, qm.calc_tdp)
            self.pol.reset_qed(qm.calc_coupling)

            qm.get_data(self.pol, base_dir, bo_list, self.dt, self.istep, calc_force_only=False)
            if (self.pol.l_qmmm and mm != None):
                mm.get_data(self.pol, base_dir, bo_list, self.istep, calc_force_only=False)
            if (not self.pol.l_nacme):
                self.pol.get_nacme()

            qed.get_data(self.pol, base_dir, pol_list, self.dt, self.istep, calc_force_only=False)
            if (not self.pol.l_pnacme):
                self.pol.get_pnacme()
            qed.transform(self.pol, mode="a2d")

            self.hop_prob(qed)
            self.hop_check(pol_list)
            self.evaluate_hop(qed, pol_list)
            if (self.l_hop):
                if (self.pol.l_qmmm and mm != None):
                    mm.get_data(self.pol, base_dir, bo_list, self.istep, calc_force_only=True)
                qed.get_data(self.pol, base_dir, pol_list, self.dt, self.istep, calc_force_only=True)

            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()
            if (self.l_collapse):
                self.check_decoherence()
                self.check_coherence()
            qed.transform(self.pol, mode="a2d")

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

            self.pol.backup_bo()
            qed.backup_qed(self.pol)
            self.pol.reset_bo(qm.calc_coupling, qm.calc_tdp)
            self.pol.reset_qed(qm.calc_coupling)

            qm.get_data(self.pol, base_dir, bo_list, self.dt, istep, calc_force_only=False)
            if (self.pol.l_qmmm and mm != None):
                mm.get_data(self.pol, base_dir, bo_list, istep, calc_force_only=False)

            if (not self.pol.l_nacme and self.l_adj_nac):
                self.pol.adjust_nac()
            if (self.l_adj_tdp):
                self.pol.adjust_tdp()
            qed.get_data(self.pol, base_dir, pol_list, self.dt, istep, calc_force_only=False)

            self.calculate_force()
            self.cl_update_velocity()

            if (not self.pol.l_nacme):
                self.pol.get_nacme()
            if (not self.pol.l_pnacme):
                self.pol.get_pnacme()
            else:
                qed.calculate_pnacme(self.pol)

            el_run(self, qed)
            qed.transform(self.pol, mode="d2a")

            self.hop_prob(qed)
            self.hop_check(pol_list)
            self.evaluate_hop(qed, pol_list)
            if (self.l_hop):
                if (self.pol.l_qmmm and mm != None):
                    mm.get_data(self.pol, base_dir, bo_list, istep, calc_force_only=True)
                qed.get_data(self.pol, base_dir, pol_list, self.dt, istep, calc_force_only=True)

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
            qed.transform(self.pol, mode="a2d")

            if ((istep + 1) % self.out_freq == 0):
                self.write_md_output(unixmd_dir, istep)
            if ((istep + 1) % self.out_freq == 0 or len(self.event["HOP"]) > 0 or len(self.event["DECO"]) > 0):
                self.print_step(istep)
            if (istep == self.nsteps - 1):
                self.write_final_xyz(unixmd_dir, istep)

            self.fstep = istep
            restart_file = os.path.join(base_dir, "RESTART.bin")
            with open(restart_file, 'wb') as f:
                pickle.dump({'qed':qed, 'qm':qm, 'md':self}, f)

        # Delete scratch directory
        if (not l_save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qed")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            if (self.pol.l_qmmm and mm != None):
                tmp_dir = os.path.join(unixmd_dir, "scr_mm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def hop_prob(self, qed):
        """ Routine to calculate hopping probabilities

            :param object qed: QED object containing cavity-molecule interaction
        """
        # Reset surface hopping variables
        self.rstate_old = self.rstate

        self.prob = np.zeros(self.pol.pst)
        self.acc_prob = np.zeros(self.pol.pst + 1)

        self.l_hop = False
        self.force_hop = False

        accum = 0.

        if (self.pol.rho_a.real[self.rstate, self.rstate] < self.lower_th):
            self.force_hop = True

        # tmp_ham = U^+ * H * U
        tmp_ham = np.zeros((self.pol.pst, self.pol.pst)) 
        tmp_ham = np.matmul(np.transpose(qed.unitary), np.matmul(qed.ham_d, qed.unitary))
        # self.pol.pnacme = U^+ * K * U + U^+ * U_dot
        # H and K are Hamiltonian and NACME in uncoupled basis

        if (not qed.l_trivial):
            for ist in range(self.pol.pst):
                if (ist != self.rstate):
                    if (self.force_hop):
                        self.prob[ist] = self.pol.rho_a.real[ist, ist] / self.upper_th
                    else:
                        self.prob[ist] = - 2. * (self.pol.rho_a.imag[self.rstate, ist] * tmp_ham[self.rstate, ist] \
                            - self.pol.rho_a.real[self.rstate, ist] * self.pol.pnacme[self.rstate, ist]) \
                            * self.dt / self.pol.rho_a.real[self.rstate, self.rstate]

                    if (self.prob[ist] < 0.):
                        self.prob[ist] = 0.
                    accum += self.prob[ist]
                self.acc_prob[ist + 1] = accum
            psum = self.acc_prob[self.pol.pst]
        else:
            for ist in range(self.pol.pst):
                if (ist != self.rstate):
                    if (ist == qed.trivial_state):
                        self.prob[ist] = 1.
                    else:
                        self.prob[ist] = 0.

                    accum += self.prob[ist]
                self.acc_prob[ist + 1] = accum
            psum = self.acc_prob[self.pol.pst]

        if (psum > 1.):
            self.prob /= psum
            self.acc_prob /= psum

    def hop_check(self, pol_list):
        """ Routine to check hopping occurs with random number

            :param integer,list pol_list: List of polaritonic states for QED calculation
        """
        self.rand = random.random()
        for ist in range(self.pol.pst):
            if (ist == self.rstate):
                continue
            if (self.rand > self.acc_prob[ist] and self.rand <= self.acc_prob[ist + 1]):
                self.l_hop = True
                self.rstate = ist
                pol_list[0] = self.rstate

    def evaluate_hop(self, qed, pol_list):
        """ Routine to evaluate hopping and velocity rescaling

            :param object qed: QED object containing cavity-molecule interaction
            :param integer,list pol_list: List of polaritonic states for QED calculation
        """
        if (self.l_hop):
            if (not qed.l_trivial):
                # Calculate potential difference between hopping states
                pot_diff = self.pol.pol_states[self.rstate].energy - self.pol.pol_states[self.rstate_old].energy

                # Solve quadratic equation for scaling factor of velocities
                a = 1.
                b = 1.
                det = 1.
                if (self.hop_rescale == "velocity"):
                    a = np.sum(self.pol.mass[0:self.pol.nat_qm] * np.sum(self.pol.pnac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(self.pol.mass[0:self.pol.nat_qm] * np.sum(self.pol.pnac[self.rstate_old, self.rstate] \
                        * self.pol.vel[0:self.pol.nat_qm], axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c
                elif (self.hop_rescale == "momentum"):
                    a = np.sum(1. / self.pol.mass[0:self.pol.nat_qm] * np.sum(self.pol.pnac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(np.sum(self.pol.pnac[self.rstate_old, self.rstate] * self.pol.vel[0:self.pol.nat_qm], axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c
                elif (self.hop_rescale == "augment"):
                    a = np.sum(1. / self.pol.mass[0:self.pol.nat_qm] * np.sum(self.pol.pnac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(np.sum(self.pol.pnac[self.rstate_old, self.rstate] * self.pol.vel[0:self.pol.nat_qm], axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c

                # Default: hopping is allowed
                self.l_reject = False

                # Velocities cannot be adjusted when zero kinetic energy is given
                if (self.hop_rescale == "energy" and self.pol.ekin_qm < eps):
                    self.l_reject = True
                # Clasically forbidden hop due to lack of kinetic energy
                if (self.pol.ekin_qm < pot_diff):
                    self.l_reject = True
                # Kinetic energy is enough, but there is no solution for scaling factor
                if (det < 0.):
                    self.l_reject = True
                # When kinetic energy is enough, velocities are always rescaled in 'augment' case
                if (self.hop_rescale == "augment" and self.pol.ekin_qm > pot_diff):
                    self.l_reject = False

                if (self.l_reject):
                    # Record event for frustrated hop
                    if (self.pol.ekin_qm < pot_diff):
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
                    pol_list[0] = self.rstate
                else:
                    if (self.hop_rescale == "energy" or (det < 0. and self.hop_rescale == "augment")):
                        if (det < 0.):
                            self.event["HOP"].append("Accept hopping: no solution to find rescale factor, but velocity is simply rescaled")
                        x = np.sqrt(1. - pot_diff / self.pol.ekin_qm)
                    else:
                        if (b < 0.):
                            x = 0.5 * (- b - np.sqrt(det)) / a
                        else:
                            x = 0.5 * (- b + np.sqrt(det)) / a

                # Rescale velocities for QM atoms
                if (not (self.hop_reject == "keep" and self.l_reject)):
                    if (self.hop_rescale == "energy"):
                        self.pol.vel[0:self.pol.nat_qm] *= x

                    elif (self.hop_rescale == "velocity"):
                        self.pol.vel[0:self.pol.nat_qm] += x * self.pol.pnac[self.rstate_old, self.rstate]

                    elif (self.hop_rescale == "momentum"):
                        self.pol.vel[0:self.pol.nat_qm] += x * self.pol.pnac[self.rstate_old, self.rstate] / \
                            self.pol.mass[0:self.pol.nat_qm].reshape((-1, 1))

                    elif (self.hop_rescale == "augment"):
                        if (det > 0. or self.pol.ekin_qm < pot_diff):
                            self.pol.vel[0:self.pol.nat_qm] += x * self.pol.pnac[self.rstate_old, self.rstate] / \
                                self.pol.mass[0:self.pol.nat_qm].reshape((-1, 1))
                        else:
                            self.pol.vel[0:self.pol.nat_qm] *= x

                # Update kinetic energy
                self.pol.update_kinetic()

        # Record hopping event
        if (self.rstate != self.rstate_old):
            if (not qed.l_trivial):
                if (self.force_hop):
                    self.event["HOP"].append(f"Accept hopping: force hop {self.rstate_old} -> {self.rstate}")
                else:
                    self.event["HOP"].append(f"Accept hopping: hop {self.rstate_old} -> {self.rstate}")
            else:
                self.event["HOP"].append(f"Trivial crossing hopping: hop {self.rstate_old} -> {self.rstate}")


    def calculate_force(self):
        """ Routine to calculate the forces
        """
        self.rforce = np.copy(self.pol.pol_states[self.rstate].force)

    def update_energy(self):
        """ Routine to update the energy of molecules in surface hopping dynamics
        """
        # Update kinetic energy
        self.pol.update_kinetic()
        self.pol.epot = self.pol.pol_states[self.rstate].energy
        self.pol.etot = self.pol.epot + self.pol.ekin

    def check_decoherence(self):
        """ Routine to check if the polaritonic state is decohered
        """
        if (self.l_hop):
            if (True in self.l_coh):
                self.event["DECO"].append(f"Destroy auxiliary trajectories: hopping occurs")
            self.l_coh = [False] * self.pol.pst
            self.l_first = [False] * self.pol.pst
            self.l_fix = [False] * self.pol.pst
        else:
            for ist in range(self.pol.pst):
                if (self.l_coh[ist]):
                    rho = self.pol.rho_a.real[ist, ist]
                    if (rho > self.upper_th):
                        self.set_decoherence(ist)
                        return

    def check_coherence(self):
        """ Routine to check coherence among polaritonic states
        """
        count = 0
        tmp_st = ""
        for ist in range(self.pol.pst):
            rho = self.pol.rho_a.real[ist, ist]
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
            self.l_coh = [False] * self.pol.pst
            self.l_first = [False] * self.pol.pst
            tmp_st = ""

        if (len(tmp_st) >= 1):
            tmp_st = tmp_st.rstrip(', ')
            self.event["DECO"].append(f"Generate auxiliary trajectory on {tmp_st} state")

    def set_decoherence(self, one_st):
        """ Routine to reset coefficient/density if the state is decohered

            :param integer one_st: State index that its population is one
        """
        self.phase = np.zeros((self.pol.pst, self.aux.nat, self.aux.ndim))
        self.pol.rho_a = np.zeros((self.pol.pst, self.pol.pst), dtype=np.complex128)
        self.pol.rho_a[one_st, one_st] = 1. + 0.j

        self.l_coh = [False] * self.pol.pst
        self.l_first = [False] * self.pol.pst
        self.l_fix = [False] * self.pol.pst

        self.event["DECO"].append(f"Destroy auxiliary trajectories: decohered to {one_st} state")

        if (self.elec_object == "coefficient"):
            for ist in range(self.pol.pst):
                if (ist == one_st):
                    self.pol.pol_states[ist].coef_a /= np.absolute(self.pol.pol_states[ist].coef_a).real
                else:
                    self.pol.pol_states[ist].coef_a = 0. + 0.j

    def aux_propagator(self):
        """ Routine to propagate auxiliary molecule
        """
        # Get auxiliary position
        for ist in range(self.pol.pst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    self.aux.pos[ist] = self.pol.pos[0:self.aux.nat]
                else:
                    if (ist == self.rstate):
                        self.aux.pos[ist] = self.pol.pos[0:self.aux.nat]
                    else:
                        self.aux.pos[ist] += self.aux.vel[ist] * self.dt

        self.pos_0 = np.copy(self.aux.pos[self.rstate])

        # Get auxiliary velocity
        self.l_collapse = False
        self.aux.vel_old = np.copy(self.aux.vel)
        for ist in range(self.pol.pst):
            # Calculate propagation factor alpha
            if (self.l_coh[ist]):
                if (self.l_fix[ist]):
                    alpha = 0.
                else:
                    if (ist == self.rstate):
                        alpha = self.pol.ekin_qm
                    else:
                        if (self.l_first[ist]):
                            alpha = self.pol.ekin_qm
                            if (self.l_econs_state):
                                alpha += self.pol.pol_states[self.rstate].energy - self.pol.pol_states[ist].energy
                        else:
                            ekin_old = np.sum(0.5 * self.aux.mass * np.sum(self.aux.vel_old[ist] ** 2, axis=1))
                            alpha = ekin_old + self.pol.pol_states[ist].energy_old - self.pol.pol_states[ist].energy
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
                alpha /= self.pol.ekin_qm
                self.aux.vel[ist] = self.pol.vel[0:self.aux.nat] * np.sqrt(alpha)

    def collapse(self, cstate):
        """ Routine to collapse coefficient/density of a state to zero
        """
        fac = 1. - self.pol.rho_a.real[cstate, cstate]

        if (self.elec_object == "coefficient"):
            for ist in range(self.pol.pst):
                if (ist == cstate):
                    self.pol.pol_states[ist].coef_a = 0. + 0.j
                else:
                    self.pol.pol_states[ist].coef_a /= np.sqrt(fac)

        self.pol.rho_a[cstate,:] = 0. + 0.j
        self.pol.rho_a[:,cstate] = 0. + 0.j
        self.pol.rho_a /= fac
         
    def get_phase(self):
        """ Routine to calculate phase term
        """
        for ist in range(self.pol.pst):
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
        tmp = f'{istep + 1:9d}' + "".join([f'{self.prob[ist]:15.8f}' for ist in range(self.pol.pst)])
        typewriter(tmp, unixmd_dir, "SHPROB", "a")

    def write_dotpop(self, unixmd_dir, istep):
        """ Write time-derivative BO population

            :param string unixmd_dir: PyUNIxMD directory
            :param integer istep: Current MD step
        """
        if (self.verbosity >= 1):
            # Write NAC term in DOTPOPNACD
            tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopnac_d])
            typewriter(tmp, unixmd_dir, "DOTPOPNACD", "a")

            # Write decoherence term in DOTPOPDECD
            tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopdec_d])
            typewriter(tmp, unixmd_dir, "DOTPOPDECD", "a")

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
            for ist in range(self.pol.pst):
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

    def print_init(self, qed, qm, mm, restart):
        """ Routine to print the initial information of dynamics

            :param object qed: QED object containing cavity-molecule interaction
            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string restart: Option for controlling dynamics restarting
        """
        # Print initial information about polariton, qed, qm, mm and thermostat
        super().print_init(qed, qm, mm, restart)

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
        ctemp = self.pol.ekin * 2. / float(self.pol.ndof) * au_to_K
        norm = 0.
        for ist in range(self.pol.pst):
            norm += self.pol.rho_a.real[ist, ist]

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d}{self.rstate:>5d}"
        INFO += f"{self.pol.ekin:16.8f}{self.pol.epot:15.8f}{self.pol.etot:15.8f}"
        INFO += f"{ctemp:13.6f}"
        INFO += f"{norm:11.5f}"
        print (INFO, flush=True)

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            DEBUG1 += f"{self.rand:11.5f}"
            for ist in range(self.pol.pst):
                DEBUG1 += f"{self.acc_prob[ist]:12.5f} ({self.rstate}->{ist})"
            print (DEBUG1, flush=True)

        # Print event in SHXF
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event["HOP"] = []
        self.event["DECO"] = []


