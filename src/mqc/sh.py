from __future__ import division
from build.el_propagator import *
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz, typewriter
from misc import eps, au_to_K, call_name
import random, os, shutil, textwrap
import numpy as np
 
class SH(MQC):
    """ Class for surface hopping dynamics

        :param object molecule: molecule object
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string propagation: propagation scheme
        :param boolean l_adjnac: logical to adjust nonadiabatic coupling
        :param string vel_rescale: velocity rescaling method after hop
        :param double threshold: electronic density threshold for unphysical density
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", l_adjnac=True, vel_rescale="simple", threshold=0.01):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, l_adjnac)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.threshold = threshold

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

    def run(self, molecule, qm, mm=None, thermostat=None, input_dir="./", \
        save_QMlog=False, save_MMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to surface hopping dynamics

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

        touch_file(molecule, qm.calc_coupling, self.propagation, unixmd_dir, SH_chk=True)
        self.print_init(molecule, qm, mm, thermostat, debug)

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

        write_md_output(molecule, qm.calc_coupling, self.propagation, unixmd_dir, istep=-1)
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

            write_md_output(molecule, qm.calc_coupling, self.propagation, unixmd_dir, istep=istep)
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

        if (molecule.rho.real[self.rstate, self.rstate] < self.threshold):
            self.force_hop = True

        for ist in range(molecule.nst):
            if (ist != self.rstate):
                if (self.force_hop):
                    self.prob[ist] = molecule.rho.real[ist, ist] / (1. - self.threshold)
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

                if (self.vel_rescale == "simple"):
                    fac = 1. - pot_diff / molecule.ekin_qm
                    # Rescale velocities for QM atoms
                    molecule.vel[0:molecule.nat_qm] *= np.sqrt(fac)

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

                        # Rescale velocities for QM atoms
                        molecule.vel[0:molecule.nat_qm] += x * molecule.nac[self.rstate_old, self.rstate, 0:molecule.nat_qm]

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

    def el_propagator(self, molecule):
        """ Routine to propagate BO coefficients or density matrix

            :param object molecule: molecule object
        """
        if (self.propagation == "coefficient"):
            el_coef(self.nesteps, self.dt, molecule)
        elif (self.propagation == "density"):
            el_rho(self.nesteps, self.dt, molecule)
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Other propagator not implemented! {self.propagation}")

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


