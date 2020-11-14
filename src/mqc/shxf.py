from __future__ import division
from build.el_propagator_xf import el_run
from mqc.mqc import MQC
from misc import eps, au_to_K, au_to_A, call_name, typewriter
import random, os, shutil, textwrap
import numpy as np

class Auxiliary_Molecule(object):
    """ Class for auxiliary molecule that is used for the calculation of decoherence term

        :param object molecule: molecule object
    """
    def __init__(self, molecule, one_dim):
        # Initialize auxiliary molecule
        if (one_dim):

            self.nat = 1
            self.nsp = 1
            self.symbols = ['XX']

            self.mass = np.zeros((self.nat))
            self.mass[0] = 1. / np.sum(1. / molecule.mass[0:molecule.nat_qm])

        else:

            self.nat = molecule.nat_qm
            self.nsp = molecule.nsp
            self.symbols = molecule.symbols

            self.mass = np.copy(molecule.mass)
        
        self.pos = np.zeros((molecule.nst, self.nat, self.nsp))
        self.vel = np.zeros((molecule.nst, self.nat, self.nsp))
        self.vel_old = np.copy(self.vel)


class SHXF(MQC):
    """ Class for DISH-XF dynamics

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
        :param string vel_rescale: velocity rescaling method after hop
        :param double threshold: electronic density threshold for decoherence term calculation
        :param wsigma: width of nuclear wave packet of auxiliary trajectory
        :type wsigma: double or double,list
        :param coefficient: initial BO coefficient
        :type coefficient: double, list or complex, list
        :param boolean l_state_wise: logical to use state-wise total energies for auxiliary trajectories
        :param string unit_dt: unit of time step (fs = femtosecond, au = atomic unit)
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", solver="rk4", l_pop_print=False, l_adjnac=True, vel_rescale="momentum", \
        threshold=0.01, wsigma=None, one_dim=False, coefficient=None, l_state_wise=False, unit_dt="fs"):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, nesteps, \
            propagation, solver, l_pop_print, l_adjnac, coefficient, unit_dt)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.rand = 0.
        self.prob = np.zeros(self.mol.nst)
        self.acc_prob = np.zeros(self.mol.nst + 1)

        self.l_hop = False
        self.force_hop = False

        self.vel_rescale = vel_rescale
        self.l_state_wise = l_state_wise

        if (self.vel_rescale == "energy"):
            pass
        elif (self.vel_rescale == "velocity"):
            if (self.mol.l_nacme): 
                raise ValueError (f"( {self.md_type}.{call_name()} ) Nonadiabatic coupling vectors are not available! l_nacme: {self.mol.l_nacme}")
        elif (self.vel_rescale == "momentum"):
            if (self.mol.l_nacme): 
                raise ValueError (f"( {self.md_type}.{call_name()} ) Nonadiabatic coupling vectors are not available! l_nacme: {self.mol.l_nacme}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'vel_rescale'! {self.vel_rescale}")

        # Initialize XF related variables
        self.one_dim = one_dim
        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst
        self.threshold = threshold

        self.wsigma = wsigma
        if (self.wsigma == None):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Sigma values should be provided in input arguments! {self.wsigma}")

        if (isinstance(self.wsigma, float)):
            # uniform value for wsigma
            pass
        elif (isinstance(self.wsigma, list)):
            # atom-resolved values for wsigma
            if (len(self.wsigma) != self.mol.nat_qm):
                raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong number of elements of sigma given! {self.wsigma}")
            if (self.one_dim):
                raise ValueError (f"( {self.md_type}.{call_name()} ) SHXF1D requires only 1 float number for sigma! {self.wsigma}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong type for sigma given! {self.wsigma}")

        self.upper_th = 1. - self.threshold
        self.lower_th = self.threshold

        # Initialize auxiliary molecule object
        self.aux = Auxiliary_Molecule(self.mol, self.one_dim)
        self.pos_0 = np.zeros((self.aux.nat, self.aux.nsp))
        self.phase = np.array(np.zeros((self.mol.nst, self.aux.nat, self.aux.nsp)))

        # Debug variables
        self.dotpopd = np.zeros(self.mol.nst)

        # Initialize event to print
        self.event = {"HOP": [], "DECO": []}

    def run(self, qm, mm=None, input_dir="./", \
        save_QMlog=False, save_MMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to decoherence-induced surface hopping dynamics

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
            shutil.move(unixmd_dir, unixmd_dir + "_old_" + str(os.getpid()))
        os.makedirs(unixmd_dir)

        QMlog_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(QMlog_dir)):
            shutil.move(QMlog_dir, QMlog_dir + "_old_" + str(os.getpid()))
        if (save_QMlog):
            os.makedirs(QMlog_dir)

        if (self.mol.qmmm and mm != None):
            MMlog_dir = os.path.join(base_dir, "MMlog")
            if (os.path.exists(MMlog_dir)):
                shutil.move(MMlog_dir, MMlog_dir + "_old_" + str(os.getpid()))
            if (save_MMlog):
                os.makedirs(MMlog_dir)

        if ((self.mol.qmmm and mm == None) or (not self.mol.qmmm and mm != None)):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Both self.mol.qmmm and mm object is necessary! {self.mol.qmmm} and {mm}")

        # Check compatibility for QM and MM objects
        if (self.mol.qmmm and mm != None):
            self.check_qmmm(qm, mm)

        # Initialize UNI-xMD
        os.chdir(base_dir)
        bo_list = [self.rstate]
        qm.calc_coupling = True

        self.touch_file(unixmd_dir)
        self.print_init(qm, mm, debug)

        # Initialize decoherence variables
        self.append_wsigma()

        # Calculate initial input geometry at t = 0.0 s
        self.mol.reset_bo(qm.calc_coupling)
        qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=-1, calc_force_only=False)
        if (self.mol.qmmm and mm != None):
            mm.get_data(self.mol, base_dir, bo_list, istep=-1, calc_force_only=False)
        if (not self.mol.l_nacme):
            self.mol.get_nacme()

        self.hop_prob(istep=-1)
        self.hop_check(bo_list)
        self.evaluate_hop(bo_list, istep=-1)
        if (qm.re_calc and self.l_hop):
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=-1, calc_force_only=True)
            if (self.mol.qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep=-1, calc_force_only=True)

        self.update_energy()

        self.check_decoherence()
        self.check_coherence()
        self.aux_propagator()
        self.get_phase()

        self.write_md_output(unixmd_dir, istep=-1)
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

            self.hop_prob(istep=istep)
            self.hop_check(bo_list)
            self.evaluate_hop(bo_list, istep=istep)
            if (qm.re_calc and self.l_hop):
                qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=istep, calc_force_only=True)
                if (self.mol.qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, istep=istep, calc_force_only=True)

            if (self.thermo != None):
                self.thermo.run(self)

            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()

            self.write_md_output(unixmd_dir, istep=istep)
            self.print_step(debug, istep=istep)
            if (istep == self.nsteps - 1):
                self.write_final_xyz(unixmd_dir, istep=istep)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            if (self.mol.qmmm and mm != None):
                tmp_dir = os.path.join(unixmd_dir, "scr_mm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def hop_prob(self, istep):
        """ Routine to calculate hopping probabilities

            :param integer istep: current MD step
        """
        # Reset surface hopping variables
        self.rstate_old = self.rstate

        self.prob = np.zeros(self.mol.nst)
        self.acc_prob = np.zeros(self.mol.nst + 1)

        self.l_hop = False
        self.force_hop = False

        accum = 0.

        if (self.mol.rho.real[self.rstate, self.rstate] < self.threshold):
            self.force_hop = True

        for ist in range(self.mol.nst):
            if (ist != self.rstate):
                if (self.force_hop):
                    self.prob[ist] = self.mol.rho.real[ist, ist] / (1. - self.threshold)
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

            :param integer,list bo_list: list of BO states for BO calculation
        """
        self.rand = random.random()
        for ist in range(self.mol.nst):
            if (ist == self.rstate):
                continue
            if (self.rand > self.acc_prob[ist] and self.rand <= self.acc_prob[ist + 1]):
                self.l_hop = True
                self.rstate = ist
                bo_list[0] = self.rstate

    def evaluate_hop(self, bo_list, istep):
        """ Routine to evaluate hopping and velocity rescaling

            :param integer,list bo_list: list of BO states for BO calculation
            :param integer istep: current MD step
        """
        if (self.l_hop):
            pot_diff = self.mol.states[self.rstate].energy - self.mol.states[self.rstate_old].energy
            if (self.mol.ekin_qm < pot_diff):
                self.l_hop = False
                self.force_hop = False
                self.event["HOP"].append(f"Reject hopping: smaller kinetic energy than potential energy difference between {self.rstate} and {self.rstate_old}")
                self.rstate = self.rstate_old
                bo_list[0] = self.rstate
            else:
                if (self.mol.ekin_qm < eps):
                    raise ValueError (f"( {self.md_type}.{call_name()} ) Too small kinetic energy! {self.mol.ekin_qm}")

                if (self.vel_rescale == "energy"):
                    fac = 1. - pot_diff / self.mol.ekin_qm
                    # Rescale velocities for QM atoms
                    self.mol.vel[0:self.mol.nat_qm] *= np.sqrt(fac)

                elif (self.vel_rescale == "velocity"):
                    a = np.sum(self.mol.mass * np.sum(self.mol.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(self.mol.mass * np.sum(self.mol.nac[self.rstate_old, self.rstate] * self.mol.vel, axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c

                    if (det < 0.):
                        self.l_hop = False
                        self.force_hop = False
                        self.rstate = self.rstate_old
                        bo_list[0] = self.rstate
                        self.event["HOP"].append("Reject hopping: no solution to find rescale factor")
                    else:
                        if (b < 0.):
                            x = 0.5 * (- b - np.sqrt(det)) / a
                        else:
                            x = 0.5 * (- b + np.sqrt(det)) / a

                        # Rescale velocities for QM atoms
                        self.mol.vel[0:self.mol.nat_qm] += x * self.mol.nac[self.rstate_old, self.rstate, 0:self.mol.nat_qm]

                elif (self.vel_rescale == "momentum"):
                    a = np.sum(1. / self.mol.mass * np.sum(self.mol.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(np.sum(self.mol.nac[self.rstate_old, self.rstate] * self.mol.vel, axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c

                    if (det < 0.):
                        self.l_hop = False
                        self.force_hop = False
                        self.rstate = self.rstate_old
                        bo_list[0] = self.rstate
                        self.event["HOP"].append("Reject hopping: no solution to find rescale factor")
                    else:
                        if (b < 0.):
                            x = 0.5 * (- b - np.sqrt(det)) / a
                        else:
                            x = 0.5 * (- b + np.sqrt(det)) / a

                        # Rescale velocities for QM atoms
                        self.mol.vel[0:self.mol.nat_qm] += x * self.mol.nac[self.rstate_old, self.rstate, 0:self.mol.nat_qm] /\
                            self.mol.mass[0:self.mol.nat_qm].reshape((-1,1))

                # Update kinetic energy
                self.mol.update_kinetic()

        # Record event
        if (self.rstate != self.rstate_old):
            if (self.force_hop):
                self.event["HOP"].append(f"Force hop {self.rstate_old} -> {self.rstate}")
            else:
                self.event["HOP"].append(f"Hopping {self.rstate_old} -> {self.rstate}")

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

            :param integer one_st: state index that its population is one
        """
        self.phase = np.zeros((self.mol.nst, self.aux.nat, self.aux.nsp))
        self.mol.rho = np.zeros((self.mol.nst, self.mol.nst), dtype=np.complex_)
        self.mol.rho[one_st, one_st] = 1. + 0.j

        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst

        if (self.propagation == "coefficient"):
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
                    if (self.one_dim):
                        self.aux.pos[ist] = np.zeros((self.aux.nat, self.aux.nsp))
                    else:
                        self.aux.pos[ist] = self.mol.pos[0:self.aux.nat]
                else:
                    if (self.one_dim):
                        self.aux.pos[ist] += self.aux.vel[ist] * self.dt
                    else:
                        if (ist == self.rstate):
                            self.aux.pos[ist] = self.mol.pos[0:self.aux.nat]
                        else:
                            self.aux.pos[ist] += self.aux.vel[ist] * self.dt

        self.pos_0 = np.copy(self.aux.pos[self.rstate])

        # Get auxiliary velocity
        self.aux.vel_old = np.copy(self.aux.vel)
        for ist in range(self.mol.nst):
            # Calculate propagation factor alpha
            if (self.l_coh[ist]):
                if (ist == self.rstate):
                    alpha = self.mol.ekin_qm
                else:
                    if (self.l_first[ist]):
                        alpha = self.mol.ekin_qm
                        if (not self.l_state_wise):
                            alpha += self.mol.states[self.rstate].energy - self.mol.states[ist].energy
                    else:
                        ekin_old = np.sum(0.5 * self.aux.mass * np.sum(self.aux.vel_old[ist] ** 2, axis=1))
                        alpha = ekin_old + self.mol.states[ist].energy_old - self.mol.states[ist].energy
                if (alpha < 0.):
                    alpha = 0.

                # Calculate auxiliary velocity from alpha
                if (self.one_dim):
                    alpha /= 0.5 * self.aux.mass[0]
                    self.aux.vel[ist] = np.sqrt(alpha)
                else:
                    alpha /= self.mol.ekin_qm
                    self.aux.vel[ist] = self.mol.vel[0:self.aux.nat] * np.sqrt(alpha) 

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

    def append_wsigma(self):
        """ Routine to append sigma values when single float number is provided
        """
        # Create a list from single float number
        if (isinstance(self.wsigma, float)):
            sigma = self.wsigma
            self.wsigma = self.aux.nat * [sigma]

    def write_md_output(self, unixmd_dir, istep): 
        """ Write output files

            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, istep)

        # Write hopping-related quantities
        self.write_sh(unixmd_dir, istep)

        # Write decoherence information
        self.write_deco(unixmd_dir, istep)

    def write_sh(self, unixmd_dir, istep): 
        """ Write hopping-related quantities into files

            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        # Write SHSTATE file
        tmp = f'{istep + 1:9d}{"":14s}{self.rstate}'
        typewriter(tmp, unixmd_dir, "SHSTATE")

        # Write SHPROB file
        tmp = f'{istep + 1:9d}' + "".join([f'{self.prob[ist]:15.8f}' for ist in range(self.mol.nst)])
        typewriter(tmp, unixmd_dir, "SHPROB")

    def write_deco(self, unixmd_dir, istep):
        """ Write XF-based decoherence information

            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        # Write time-derivative density matrix elements in DOTPOTD
        tmp = f'{istep + 1:9d}' + "".join([f'{self.dotpopd[ist]:15.8f}' for ist in range(self.mol.nst)])
        typewriter(tmp, unixmd_dir, "DOTPOPD")        

        # Write auxiliary trajectories
        for ist in range(self.mol.nst):
            if (self.l_coh[ist]):
                self.write_aux_movie(unixmd_dir, ist, istep=istep)

    def write_aux_movie(self, unixmd_dir, ist, istep):
        """ Write auxiliary trajecoty movie file    

            :param string unixmd_dir: unixmd directory
            :param integer ist: current adiabatic state
            :param integer istep: current MD step
        """
        # Write auxiliary trajectory movie files
        tmp = f'{self.aux.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)'
        typewriter(tmp, unixmd_dir, f"AUX_MOVIE_{ist}.xyz")
        for iat in range(self.aux.nat):
            tmp = f'{self.aux.symbols[iat]:5s}' + \
                "".join([f'{self.aux.pos[ist, iat, isp] * au_to_A:15.8f}' for isp in range(self.aux.nsp)]) \
                + "".join([f"{self.aux.vel[ist, iat, isp]:15.8f}" for isp in range(self.aux.nsp)])
            typewriter(tmp, unixmd_dir, f"AUX_MOVIE_{ist}.xyz")

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
        INIT = f" #INFO{'STEP':>8s}{'State':>7s}{'Max. Prob.':>14s}{'Rand.':>12s}{'Kinetic(H)':>15s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'Norm.':>8s}"
        dynamics_step_info += INIT

        # Print DEBUG1 for each step
        if (debug >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        # Print DEBUG2 for each step
        if (debug >= 2):
            DEBUG2 = f" #DEBUG2{'STEP':>6s}{'Acc. Hopping Prob.':>22s}"
            dynamics_step_info += "\n" + DEBUG2

        print (dynamics_step_info, flush=True)

    def print_step(self, debug, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer debug: verbosity level for standard output
            :param integer istep: current MD step
        """
        if (istep == -1):
            max_prob = 0.
            hstate = self.rstate
        else:
            max_prob = max(self.prob)
            hstate = np.where(self.prob == max_prob)[0][0]

        ctemp = self.mol.ekin * 2. / float(self.mol.dof) * au_to_K
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
        if (debug >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{self.mol.states[ist].energy:17.8f} "
            print (DEBUG1, flush=True)

        # Print DEBUG2 for each step
        if (debug >= 2):
            DEBUG2 = f" DEBUG2{istep + 1:>7d}"
            for ist in range(self.mol.nst):
                DEBUG2 += f"{self.acc_prob[ist]:12.5f}({self.rstate}->{ist})"
            print (DEBUG2, flush=True)

        # Print event in SHXF
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event["HOP"] = []
        self.event["DECO"] = []
