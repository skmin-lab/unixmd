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
        
        self.one_dim = one_dim

        self.pos = []
        self.vel = []
        self.vel_old = []
        
        if (self.one_dim):
            
            self.nat = 1
            self.nsp = 1

            self.mass = np.zeros((self.nat))
            self.pos = np.zeros((molecule.nst, self.nat, self.nsp))
            self.vel = np.zeros((molecule.nst, self.nat, self.nsp))
            self.vel_old = np.copy(self.vel)
            self.mass[0] = 1. / np.sum(1. / molecule.mass)
        
        else:
            
            self.nat = molecule.nat
            self.nsp = molecule.nsp

            for ist in range(molecule.nst):
                self.pos.append(molecule.pos)
                self.vel.append(molecule.vel)
                self.vel_old.append(molecule.vel)

            self.pos = np.array(self.pos)
            self.vel = np.array(self.vel)
            self.vel_old = np.array(self.vel_old)

            self.mass = np.copy(molecule.mass)


class SHXF(MQC):
    """ Class for DISH-XF dynamics

        :param object molecule: molecule object
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string propagation: propagation scheme
        :param boolean l_adjnac: logical to adjust nonadiabatic coupling
        :param string vel_rescale: velocity rescaling method after hop
        :param double threshold: electronic density threshold for decoherence term calculation
        :param double wsigma: width of nuclear wave packet of auxiliary trajectory
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True, vel_rescale="simple", threshold=0.01, wsigma=0.1, one_dim=False):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, l_adjnac)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.rand = 0.
        self.prob = np.zeros(molecule.nst)
        self.acc_prob = np.zeros(molecule.nst + 1)

        self.l_hop = False
        self.force_hop = False
        
        if (vel_rescale == "simple"):
            self.vel_rescale = vel_rescale
        elif (vel_rescale == "nac"):
            if (molecule.l_nacme): 
                raise ValueError (f"( {self.md_type}.{call_name()} ) Nonadiabatic coupling vectors are not available! l_nacme: {molecule.l_nacme}")
            else:
                self.vel_rescale = vel_rescale
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'vel_rescale'! {self.vel_rescale}")

        # Initialize XF related variables
        self.one_dim = one_dim
        self.l_coh = []
        self.l_first = []
        for ist in range(molecule.nst):
            self.l_coh.append(False)
            self.l_first.append(False)
        self.phase = np.array(np.zeros((molecule.nst, molecule.nat, molecule.nsp)))
        self.tot_E = np.array(np.zeros((molecule.nst)))
        self.threshold = threshold
        self.wsigma = wsigma

        self.upper_th = 1. - self.threshold
        self.lower_th = self.threshold

        # Initialize auxiliary molecule object
        self.aux = Auxiliary_Molecule(molecule, one_dim)
        # TODO: pos_i\0 to aux obj??? 1D, it seems right, but full dimension, unnecessary
        self.pos_0 = np.zeros((self.aux.nat, self.aux.nsp))

    def run(self, molecule, theory, thermostat=None, input_dir="./", \
        save_QMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to decoherence-induced surface hopping dynamics

            :param object molecule: molecule object
            :param object theory: theory object containing on-the-fly calculation infomation
            :param object thermostat: thermostat type
            :param string input_dir: location of input directory
            :param boolean save_QMlog: logical for saving QM calculation log
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

        # Initialize UNI-xMD
        os.chdir(base_dir)
        bo_list = [self.rstate]
        theory.calc_coupling = True

        touch_file(molecule, theory.calc_coupling, self.propagation, \
            unixmd_dir, SH_chk=True)
        self.print_init(molecule, theory, thermostat, debug)

        # Calculate initial input geometry at t = 0.0 s
        theory.get_bo(molecule, base_dir, -1, bo_list, self.dt, calc_force_only=False)
        if (not molecule.l_nacme):
            molecule.get_nacme()

        self.hop_prob(molecule, -1, unixmd_dir)
        self.hop_check(molecule, bo_list)
        self.evaluate_hop(molecule, bo_list, -1, unixmd_dir)
        if (theory.re_calc and self.l_hop):
            theory.get_bo(molecule, base_dir, -1, bo_list, self.dt, calc_force_only=True)

        self.update_energy(molecule)

        self.check_decoherence(molecule)
        self.check_coherence(molecule)
        self.aux_propagator(molecule)
        self.get_phase(molecule)

        write_md_output(molecule, theory.calc_coupling, -1, \
            self.propagation, unixmd_dir)
        self.print_step(molecule, -1, debug)

        # Main MD loop
        for istep in range(self.nsteps):

            self.cl_update_position(molecule)

            molecule.backup_bo()
            theory.get_bo(molecule, base_dir, istep, bo_list, self.dt, calc_force_only=False)

            if (not molecule.l_nacme):
                molecule.adjust_nac()

            self.cl_update_velocity(molecule)

            if (not molecule.l_nacme):
                molecule.get_nacme()

            self.el_propagator(molecule)

            self.hop_prob(molecule, istep, unixmd_dir)
            self.hop_check(molecule, bo_list)
            self.evaluate_hop(molecule, bo_list, istep, unixmd_dir)
            if (theory.re_calc and self.l_hop):
                theory.get_bo(molecule, base_dir, istep, bo_list, self.dt, calc_force_only=True)

            if (thermostat != None):
                thermostat.run(molecule, self)

            self.update_energy(molecule)

            self.check_decoherence(molecule)
            self.check_coherence(molecule)
            self.aux_propagator(molecule)
            self.get_phase(molecule)

            write_md_output(molecule, theory.calc_coupling, istep, \
                self.propagation, unixmd_dir)
            self.print_step(molecule, istep, debug)
            if (istep == self.nsteps - 1):
                write_final_xyz(molecule, istep, unixmd_dir)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

    def hop_prob(self, molecule, istep, unixmd_dir):
        """ Routine to calculate hopping probabilities

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param string unixmd_dir: md directory
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

    def evaluate_hop(self, molecule, bo_list, istep, unixmd_dir):
        """ Routine to evaluate hopping and velocity rescaling

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param integer istep: current MD step
            :param string unixmd_dir: unixmd directory
        """
        if (self.l_hop):        
            pot_diff = molecule.states[self.rstate].energy - molecule.states[self.rstate_old].energy
            if (molecule.ekin < pot_diff):
                if (not self.force_hop):
                    self.l_hop = False
                    self.rstate = self.rstate_old
                    bo_list[0] = self.rstate
            else:
                if (molecule.ekin < eps):
                    raise ValueError (f"( {self.md_type}.{call_name()} ) Too small kinetic energy! {molecule.ekin}")
                 
                if (self.vel_rescale == "simple"):
                    fac = 1. - pot_diff / molecule.ekin
                    molecule.vel *= np.sqrt(fac)
                 
                elif (self.vel_rescale == "nac"):
                    
                    a = np.sum(molecule.mass * np.sum(molecule.nac[self.rstate_old, self.rstate] ** 2., axis=1))
                    b = 2. * np.sum(molecule.mass * np.sum(molecule.nac[self.rstate_old, self.rstate] * molecule.vel, axis=1))
                    c = 2. * pot_diff
                    det = b ** 2. - 4. * a * c
                    
                    if (det < 0.):
                        self.l_hop = False
                        self.rstate = self.rstate_old
                        bo_list[0] = self.rstate
                    else:
                        if(b < 0.):
                            x = 0.5 * (- b - np.sqrt(det)) / a 
                    
                        else:
                            x = 0.5 * (- b + np.sqrt(det)) / a 
                    
                        molecule.vel += x * molecule.nac[self.rstate_old, self.rstate]
                
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
            for ist in range(molecule.nst):
                self.l_coh[ist] = False
                self.l_first[ist] = False


    def check_decoherence(self, molecule):
        """ Routine to check if the electronic state is decohered

            :param object molecule: molecule object
        """
        if (self.l_hop):
            for ist in range(molecule.nst):
                self.l_coh[ist] = False
                self.l_first[ist] = False
        else:
            for ist in range(molecule.nst):
                if (self.l_coh[ist]):
                    rho = molecule.rho.real[ist, ist]
                    if (rho > self.upper_th):
                        self.set_decoherence(molecule, ist)
                        return

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
                        self.aux.pos[ist] = np.copy(molecule.pos)
                else:
                    if (self.one_dim):
                        self.aux.pos[ist] += self.aux.vel[ist] * self.dt
                    else:
                        if (ist == self.rstate):
                            self.aux.pos[ist] = np.copy(molecule.pos)
                        else:
                            self.aux.pos[ist] += self.aux.vel[ist] * self.dt

        if (self.one_dim):
            self.pos_0 = np.copy(self.aux.pos[self.rstate])
        else:
            self.pos_0 = np.copy(molecule.pos)

        # Get auxiliary velocity
        
        self.aux.vel_old = np.copy(self.aux.vel)
        
        if (self.one_dim):
            self.aux.vel[self.rstate] = np.sqrt(2. * molecule.ekin / self.aux.mass[0])
        else:
            self.aux.vel[self.rstate] = np.copy(molecule.vel)

        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    self.tot_E[ist] = molecule.ekin + molecule.states[ist].energy
                    self.aux.vel[ist] = np.copy(self.aux.vel[self.rstate])
                else:
                    alpha = self.tot_E[ist] - molecule.states[ist].energy
                    if (alpha < eps):
                        self.aux.vel[ist] *= 0.
                    else:
                        if (self.one_dim):
                            alpha /= 0.5 * self.aux.mass[0]
                            self.aux.vel[ist] = np.sqrt(alpha)
                        else:
                            alpha /= molecule.ekin
                            alpha = np.sqrt(alpha)
                            self.aux.vel[ist] = molecule.vel * alpha
            else:
                if (self.one_dim):
                    self.aux.vel[ist] = np.copy(self.aux.vel[self.rstate])
                else:
                    self.aux.vel[ist] = np.copy(molecule.vel)

    def set_decoherence(self, molecule, one_st):
        """ Routine to reset coefficient/density if the state is decohered

            :param object molecule: molecule object
            :param integer one_st: state index that its population is one
        """
        self.phase = np.zeros((molecule.nst, self.aux.nat, self.aux.nsp))
        molecule.rho = np.zeros((molecule.nst, molecule.nst), dtype=np.complex_)
        molecule.rho[one_st, one_st] = 1. + 0.j

        for ist in range(molecule.nst):
            self.l_coh[ist] = False
            self.l_first[ist] = False
        
        if (self.propagation == "coefficient"):
            for ist in range(molecule.nst):
                if (ist == one_st):
                    molecule.states[ist].coef /= np.absolute(molecule.states[ist].coef).real
                else:
                    molecule.states[ist].coef = 0. + 0.j
 
    def get_phase(self, molecule):
        """ Routine to calculate phase term

            :param object molecule: molecule object
        """
        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    self.phase[ist] = 0.
                else:
                    self.phase[ist] += molecule.mass * \
                        (self.aux.vel[ist] - self.aux.vel_old[ist])

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

    def print_init(self, molecule, theory, thermostat, debug):
        """ Routine to print the initial information of dynamics

            :param object molecule: molecule object
            :param object theory: theory object containing on-the-fly calculation infomation
            :param object thermostat: thermostat type
            :param integer debug: verbosity level for standard output
        """
        # Print initial information about molecule, theory and thermostat
        super().print_init(molecule, theory, thermostat, debug)

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

    def print_step(self, molecule, istep, debug):
        """ Routine to print each steps infomation about dynamics

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param integer debug: verbosity level for standard output
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


