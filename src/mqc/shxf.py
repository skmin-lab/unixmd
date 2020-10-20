from __future__ import division
from build.el_propagator import *
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz, typewriter
from misc import eps, au_to_K, call_name
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

            self.mass = np.zeros((self.nat))
            self.mass[0] = 1. / np.sum(1. / molecule.mass[0:molecule.nat_qm])

        else:

            self.nat = molecule.nat_qm
            self.nsp = molecule.nsp

            self.mass = np.copy(molecule.mass)
        
        self.pos = np.zeros((molecule.nst, self.nat, self.nsp))
        self.vel = np.zeros((molecule.nst, self.nat, self.nsp))
        self.vel_old = np.copy(self.vel)


class SHXF(MQC):
    """ Class for DISH-XF dynamics

        :param object molecule: molecule object
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string propagation: propagation scheme
        :param boolean l_pop_print: logical to print BO population and coherence
        :param boolean l_adjnac: logical to adjust nonadiabatic coupling
        :param string vel_rescale: velocity rescaling method after hop
        :param double threshold: electronic density threshold for decoherence term calculation
        :param wsigma: width of nuclear wave packet of auxiliary trajectory
        :type wsigma: double or double,list
        :param coefficient: initial BO coefficient
        :type coefficient: double, list or complex, list
        :param boolean l_state_wise: logical to use state-wise total energies for auxiliary trajectories
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_pop_print=False, l_adjnac=True, vel_rescale="momentum", \
        threshold=0.01, wsigma=None, one_dim=False, coefficient=None, l_state_wise=False):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, l_pop_print, l_adjnac, coefficient)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.rand = 0.
        self.prob = np.zeros(molecule.nst)
        self.acc_prob = np.zeros(molecule.nst + 1)

        self.l_hop = False
        self.force_hop = False

        self.vel_rescale = vel_rescale
        self.l_state_wise = l_state_wise
        
        if (self.vel_rescale == "energy"):
            pass
        elif (self.vel_rescale == "velocity"):
            if (molecule.l_nacme): 
                raise ValueError (f"( {self.md_type}.{call_name()} ) Nonadiabatic coupling vectors are not available! l_nacme: {molecule.l_nacme}")
        elif (self.vel_rescale == "momentum"):
            if (molecule.l_nacme): 
                raise ValueError (f"( {self.md_type}.{call_name()} ) Nonadiabatic coupling vectors are not available! l_nacme: {molecule.l_nacme}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'vel_rescale'! {self.vel_rescale}")

        # Initialize XF related variables
        self.one_dim = one_dim
        self.l_coh = [False] * molecule.nst
        self.l_first = [False] * molecule.nst
#        self.tot_E = np.array(np.zeros((molecule.nst)))
        self.threshold = threshold

        self.wsigma = wsigma
        if (self.wsigma == None):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Sigma values should be provided in input arguments! {self.wsigma}")

        if (isinstance(self.wsigma, float)):
            # uniform value for wsigma
            pass
        elif (isinstance(self.wsigma, list)):
            # atom-resolved values for wsigma
            if (len(self.wsigma) != molecule.nat_qm):
                raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong number of elements of sigma given! {self.wsigma}")
            if (self.one_dim):
                raise ValueError (f"( {self.md_type}.{call_name()} ) SHXF1D requires only 1 float number for sigma! {self.wsigma}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong type for sigma given! {self.wsigma}")

        self.upper_th = 1. - self.threshold
        self.lower_th = self.threshold

        # Initialize auxiliary molecule object
        self.aux = Auxiliary_Molecule(molecule, self.one_dim)
        self.pos_0 = np.zeros((self.aux.nat, self.aux.nsp))
        self.phase = np.array(np.zeros((molecule.nst, self.aux.nat, self.aux.nsp)))

    def run(self, molecule, qm, mm=None, thermostat=None, input_dir="./", \
        save_QMlog=False, save_MMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to decoherence-induced surface hopping dynamics

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

        # Initialize UNI-xMD
        os.chdir(base_dir)
        bo_list = [self.rstate]
        qm.calc_coupling = True

        touch_file(molecule, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, SH_chk=True)
        self.print_init(molecule, qm, mm, thermostat, debug)

        # Initialize decoherence variables
        self.append_wsigma()

        # Calculate initial input geometry at t = 0.0 s
        molecule.reset_bo(qm.calc_coupling)
        qm.get_data(molecule, base_dir, bo_list, self.dt, istep=-1, calc_force_only=False)
        if (molecule.qmmm and mm != None):
            mm.get_data(molecule, base_dir, bo_list, istep=-1, calc_force_only=False)
        if (not molecule.l_nacme):
            molecule.get_nacme()

        self.hop_prob(molecule, unixmd_dir, istep=-1)
        self.hop_check(molecule, bo_list)
        self.evaluate_hop(molecule, bo_list, unixmd_dir, istep=-1)
        if (qm.re_calc and self.l_hop):
            qm.get_data(molecule, base_dir, bo_list, self.dt, istep=-1, calc_force_only=True)
            if (molecule.qmmm and mm != None):
                mm.get_data(molecule, base_dir, bo_list, istep=-1, calc_force_only=True)

        self.update_energy(molecule)

        self.check_decoherence(molecule)
        self.check_coherence(molecule)
        self.aux_propagator(molecule)
        self.get_phase(molecule)

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

            self.el_propagator(molecule)

            self.hop_prob(molecule, unixmd_dir, istep=istep)
            self.hop_check(molecule, bo_list)
            self.evaluate_hop(molecule, bo_list, unixmd_dir, istep=istep)
            if (qm.re_calc and self.l_hop):
                qm.get_data(molecule, base_dir, bo_list, self.dt, istep=istep, calc_force_only=True)
                if (molecule.qmmm and mm != None):
                    mm.get_data(molecule, base_dir, bo_list, istep=istep, calc_force_only=True)

            if (thermostat != None):
                thermostat.run(molecule, self)

            self.update_energy(molecule)

            self.check_decoherence(molecule)
            self.check_coherence(molecule)
            self.aux_propagator(molecule)
            self.get_phase(molecule)

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

    def hop_prob(self, molecule, unixmd_dir, istep):
        """ Routine to calculate hopping probabilities

            :param object molecule: molecule object
            :param string unixmd_dir: md directory
            :param integer istep: current MD step
        """
        # Reset surface hopping variables
        self.rstate_old = self.rstate

        self.prob = np.zeros(molecule.nst)
        self.acc_prob = np.zeros(molecule.nst + 1)

        self.l_hop = False
        self.force_hop = False

        accum = 0.

        if (molecule.rho.real[self.rstate, self.rstate] < eps):
            self.force_hop = True

        for ist in range(molecule.nst):
            if (ist != self.rstate):
                if (self.force_hop):
                    self.prob[ist] = molecule.rho.real[ist, ist] / (1. - eps)
                else:
                    self.prob[ist] = - 2. * molecule.rho.real[ist, self.rstate] * \
                        molecule.nacme[ist, self.rstate] * self.dt / molecule.rho.real[self.rstate, self.rstate]

                if (self.prob[ist] < 0.):
                    self.prob[ist] = 0.
                accum += self.prob[ist]
            self.acc_prob[ist + 1] = accum
        psum = self.acc_prob[molecule.nst]
 
        if (psum > 1.):
            self.prob /= psum
            self.acc_prob /= psum

        # Write SHPROB file
        tmp = f'{istep + 1:9d}' + "".join([f'{self.prob[ist]:15.8f}' for ist in range(molecule.nst)])
        typewriter(tmp, unixmd_dir, "SHPROB")

    def hop_check(self, molecule, bo_list):
        """ Routine to check hopping occurs with random number

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
        """
        self.rand = random.random()
        for ist in range(molecule.nst):
            if (ist == self.rstate):
                continue
            if (self.rand > self.acc_prob[ist] and self.rand <= self.acc_prob[ist + 1]):
                self.l_hop = True
                self.rstate = ist
                bo_list[0] = self.rstate

    def evaluate_hop(self, molecule, bo_list, unixmd_dir, istep):
        """ Routine to evaluate hopping and velocity rescaling

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        if (self.l_hop):
            pot_diff = molecule.states[self.rstate].energy - molecule.states[self.rstate_old].energy
            if (molecule.ekin_qm < pot_diff):
                if (not self.force_hop):
                    self.l_hop = False
                    self.rstate = self.rstate_old
                    bo_list[0] = self.rstate
            else:
                if (molecule.ekin_qm < eps):
                    raise ValueError (f"( {self.md_type}.{call_name()} ) Too small kinetic energy! {molecule.ekin_qm}")

                if (self.vel_rescale == "energy"):
                    fac = 1. - pot_diff / molecule.ekin_qm
                    # Rescale velocities for QM atoms
                    molecule.vel[0:molecule.nat_qm] *= np.sqrt(fac)

                elif (self.vel_rescale == "velocity"):
                    a = np.sum(molecule.mass * np.sum(molecule.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(molecule.mass * np.sum(molecule.nac[self.rstate_old, self.rstate] * molecule.vel, axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c

                    if (det < 0.):
                        self.l_hop = False
                        self.force_hop = False
                        self.rstate = self.rstate_old
                        bo_list[0] = self.rstate
                    else:
                        if (b < 0.):
                            x = 0.5 * (- b - np.sqrt(det)) / a
                        else:
                            x = 0.5 * (- b + np.sqrt(det)) / a

                        # Rescale velocities for QM atoms
                        molecule.vel[0:molecule.nat_qm] += x * molecule.nac[self.rstate_old, self.rstate, 0:molecule.nat_qm]

                elif (self.vel_rescale == "momentum"):
                    a = np.sum(1. / molecule.mass * np.sum(molecule.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(np.sum(molecule.nac[self.rstate_old, self.rstate] * molecule.vel, axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c

                    if (det < 0.):
                        self.l_hop = False
                        self.force_hop = False
                        self.rstate = self.rstate_old
                        bo_list[0] = self.rstate
                    else:
                        if (b < 0.):
                            x = 0.5 * (- b - np.sqrt(det)) / a
                        else:
                            x = 0.5 * (- b + np.sqrt(det)) / a

                        # Rescale velocities for QM atoms
                        molecule.vel[0:molecule.nat_qm] += x * molecule.nac[self.rstate_old, self.rstate, 0:molecule.nat_qm] /\
                            molecule.mass[0:molecule.nat_qm].reshape((-1,1))

                # Update kinetic energy
                molecule.update_kinetic()

        # Write SHSTATE file
        tmp = f'{istep + 1:9d}{"":14s}{self.rstate}'
        typewriter(tmp, unixmd_dir, "SHSTATE")

    def calculate_force(self, molecule):
        """ Routine to calculate the forces

            :param object molecule: molecule object
        """
        self.rforce = np.copy(molecule.states[self.rstate].force)

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in surface hopping dynamics

            :param object molecule: molecule object
        """
        # Update kinetic energy
        molecule.update_kinetic()
        molecule.epot = molecule.states[self.rstate].energy
        molecule.etot = molecule.epot + molecule.ekin

    def check_decoherence(self, molecule):
        """ Routine to check if the electronic state is decohered

            :param object molecule: molecule object
        """
        if (self.l_hop):
            self.l_coh = [False] * molecule.nst
            self.l_first = [False] * molecule.nst
        else:
            for ist in range(molecule.nst):
                if (self.l_coh[ist]):
                    rho = molecule.rho.real[ist, ist]
                    if (rho > self.upper_th):
                        self.set_decoherence(molecule, ist)
                        return

    def check_coherence(self, molecule):
        """ Routine to check coherence among BO states

            :param object molecule: molecule object
        """
        count = 0
        for ist in range(molecule.nst):
            rho = molecule.rho.real[ist, ist]
            if (rho > self.upper_th or rho < self.lower_th):
                self.l_coh[ist] = False
            else:
                if (self.l_coh[ist]):
                    self.l_first[ist] = False
                else:
                    self.l_first[ist] = True
                self.l_coh[ist] = True
                count += 1

        if (count < 2):
            self.l_coh = [False] * molecule.nst
            self.l_first = [False] * molecule.nst

    def set_decoherence(self, molecule, one_st):
        """ Routine to reset coefficient/density if the state is decohered

            :param object molecule: molecule object
            :param integer one_st: state index that its population is one
        """
        self.phase = np.zeros((molecule.nst, self.aux.nat, self.aux.nsp))
        molecule.rho = np.zeros((molecule.nst, molecule.nst), dtype=np.complex_)
        molecule.rho[one_st, one_st] = 1. + 0.j

        self.l_coh = [False] * molecule.nst
        self.l_first = [False] * molecule.nst

        if (self.propagation == "coefficient"):
            for ist in range(molecule.nst):
                if (ist == one_st):
                    molecule.states[ist].coef /= np.absolute(molecule.states[ist].coef).real
                else:
                    molecule.states[ist].coef = 0. + 0.j
 
    def aux_propagator(self, molecule):
        """ Routine to propagate auxiliary molecule

            :param object molecule: molecule object
        """
        # Get auxiliary position
        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    if (self.one_dim):
                        self.aux.pos[ist] = np.zeros((self.aux.nat, self.aux.nsp))
                    else:
                        self.aux.pos[ist] = molecule.pos[0:self.aux.nat]
                else:
                    if (self.one_dim):
                        self.aux.pos[ist] += self.aux.vel[ist] * self.dt
                    else:
                        if (ist == self.rstate):
                            self.aux.pos[ist] = molecule.pos[0:self.aux.nat]
                        else:
                            self.aux.pos[ist] += self.aux.vel[ist] * self.dt

        self.pos_0 = np.copy(molecule.pos[0:self.aux.nat])

        # Get auxiliary velocity
        self.aux.vel_old = np.copy(self.aux.vel)
        for ist in range(molecule.nst):
            # Calculate propagation factor alpha
            if (self.l_coh[ist]):
                if (ist == self.rstate):
                    alpha = molecule.ekin_qm
                else:
                    if (self.l_first[ist]):
                        alpha = molecule.ekin_qm
                        if (not self.l_state_wise):
                            alpha += molecule.states[self.rstate].energy - molecule.states[ist].energy
                    else:
                        ekin_old = np.sum(0.5 * self.aux.mass * np.sum(self.aux.vel_old[ist] ** 2, axis=1))
                        alpha = ekin_old + molecule.states[ist].energy_old - molecule.states[ist].energy
                if (alpha < 0.):
                    alpha = 0.
                
                # Calculate auxiliary velocity from alpha
                if (self.one_dim):
                    alpha /= 0.5 * self.aux.mass[0]
                    self.aux.vel[ist] = np.sqrt(alpha)
                else:
                    alpha /= molecule.ekin_qm
                    self.aux.vel[ist] = molecule.vel[0:self.aux.nat] * np.sqrt(alpha) 

    def get_phase(self, molecule):
        """ Routine to calculate phase term

            :param object molecule: molecule object
        """
        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    self.phase[ist] = 0.
                else:
                    for iat in range(self.aux.nat):
                        self.phase[ist, iat] += self.aux.mass[iat] * \
                            (self.aux.vel[ist, iat] - self.aux.vel_old[ist, iat])

    def el_propagator(self, molecule):
        """ Routine to propagate BO coefficients or density matrix

            :param object molecule: molecule object
        """
        if (self.propagation == "coefficient"):
            el_coef_xf(self, molecule)
        elif (self.propagation == "density"):
            el_rho_xf(self, molecule)
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Other propagator not implemented! {self.propagation}")

    def append_wsigma(self):
        """ Routine to append sigma values when single float number is provided
        """
        # Create a list from single float number
        if (isinstance(self.wsigma, float)):
            sigma = self.wsigma
            self.wsigma = self.aux.nat * [sigma]

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
        INIT = f" #INFO{'STEP':>8s}{'State':>7s}{'Max. Prob.':>14s}{'Rand.':>12s}{'Kinetic(H)':>15s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'Norm.':>8s}"
        dynamics_step_info += INIT

        # Print DEBUG1 for each step
        if (debug >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}"
            for ist in range(molecule.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        # Print DEBUG2 for each step
        if (debug >= 2):
            DEBUG2 = f" #DEBUG2{'STEP':>6s}{'Acc. Hopping Prob.':>22s}"
            dynamics_step_info += "\n" + DEBUG2

        print (dynamics_step_info, flush=True)

    def print_step(self, molecule, debug, istep):
        """ Routine to print each steps infomation about dynamics

            :param object molecule: molecule object
            :param integer debug: verbosity level for standard output
            :param integer istep: current MD step
        """
        if (istep == -1):
            max_prob = 0.
            hstate = self.rstate
        else:
            max_prob = max(self.prob)
            hstate = np.where(self.prob == max_prob)[0][0]

        ctemp = molecule.ekin * 2. / float(molecule.dof) * au_to_K
        norm = 0.
        for ist in range(molecule.nst):
            norm += molecule.rho.real[ist, ist]

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d}{self.rstate:>5d}{max_prob:11.5f} ({self.rstate}->{hstate}){self.rand:11.5f}"
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

        # Print DEBUG2 for each step
        if (debug >= 2):
            DEBUG2 = f" DEBUG2{istep + 1:>7d}"
            for ist in range(molecule.nst):
                DEBUG2 += f"{self.acc_prob[ist]:12.5f}({self.rstate}->{ist})"
            print (DEBUG2, flush=True)

        # Print event in surface hopping
        if (self.rstate != self.rstate_old):
            print (f" Hopping {self.rstate_old} -> {self.rstate}", flush=True)

        if (self.force_hop):
            print (f" Force hop {self.rstate_old} -> {self.rstate}", flush=True)


