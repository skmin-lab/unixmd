from __future__ import division
from build.el_propagator import el_run
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz, typewriter
from misc import eps, au_to_K, call_name
import random, os, shutil, textwrap
import numpy as np
 
class SH(MQC):
    """ Class for surface hopping dynamics

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
        :param coefficient: initial BO coefficient
        :param string decoherence_method: simple decoherence correction schemes 
        :param double c_edc: energy constant for rescaling coefficients
        :type coefficient: double, list or complex, list
        :param string unit_dt: unit of time step (fs = femtosecond, au = atomic unit)
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", solver="rk4", l_pop_print=False, l_adjnac=True, \
        vel_rescale="momentum", coefficient=None, decoherence_method="IDC", c_edc=0.1 , unit_dt="fs"):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, nesteps, \
            propagation, solver, l_pop_print, l_adjnac, coefficient, decoherence_method, c_edc, unit_dt)

        # Initialize SH variables
        self.rstate = istate
        self.rstate_old = self.rstate

        self.rand = 0.
        self.prob = np.zeros(self.mol.nst)
        self.acc_prob = np.zeros(self.mol.nst + 1)
        self.decoherence_method = decoherence_method
        self.c_edc = c_edc

        self.l_hop = False

        self.vel_rescale = vel_rescale
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

        # Initialize event to print
        self.event = {"HOP": []}

    def run(self, qm, mm=None, input_dir="./", \
        save_QMlog=False, save_MMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to surface hopping dynamics

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
            shutil.rmtree(unixmd_dir)
        os.makedirs(unixmd_dir)

        QMlog_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(QMlog_dir)):
            shutil.rmtree(QMlog_dir)
        if (save_QMlog):
            os.makedirs(QMlog_dir)

        if (self.mol.qmmm and mm != None):
            MMlog_dir = os.path.join(base_dir, "MMlog")
            if (os.path.exists(MMlog_dir)):
                shutil.rmtree(MMlog_dir)
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

        touch_file(self.mol, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, SH_chk=True, XF_chk=False)
        self.print_init(qm, mm, debug)

        # Calculate initial input geometry at t = 0.0 s
        self.mol.reset_bo(qm.calc_coupling)
        qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=-1, calc_force_only=False)
        if (self.mol.qmmm and mm != None):
            mm.get_data(self.mol, base_dir, bo_list, istep=-1, calc_force_only=False)
        if (not self.mol.l_nacme):
            self.mol.get_nacme()

        self.hop_prob(unixmd_dir, istep=-1)
        self.hop_check(bo_list)
        self.evaluate_hop(bo_list, unixmd_dir, istep=-1)
        if (qm.re_calc and self.l_hop):
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=-1, calc_force_only=True)
            if (self.mol.qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep=-1, calc_force_only=True)

        self.update_energy()

        write_md_output(self.mol, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=-1)
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

            self.hop_prob(unixmd_dir, istep=istep)
            self.hop_check(bo_list)
            self.evaluate_hop(bo_list, unixmd_dir, istep=istep)
            if (qm.re_calc and self.l_hop):
                qm.get_data(self.mol, base_dir, bo_list, self.dt, istep=istep, calc_force_only=True)
                if (self.mol.qmmm and mm != None):
                    mm.get_data(self.mol, base_dir, bo_list, istep=istep, calc_force_only=True)

            self.decoherence_scheme()

            if (self.thermo != None):
                self.thermo.run(self)

            self.update_energy()

            write_md_output(self.mol, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=istep)
            self.print_step(debug, istep=istep)
            if (istep == self.nsteps - 1):
                write_final_xyz(self.mol, unixmd_dir, istep=istep)

        # Delete scratch directory
        if (not save_scr):
            tmp_dir = os.path.join(unixmd_dir, "scr_qm")
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)

            if (self.mol.qmmm and mm != None):
                tmp_dir = os.path.join(unixmd_dir, "scr_mm")
                if (os.path.exists(tmp_dir)):
                    shutil.rmtree(tmp_dir)

    def hop_prob(self, unixmd_dir, istep):
        """ Routine to calculate hopping probabilities

            :param string unixmd_dir: md directory
            :param integer istep: current MD step
        """
        # Reset surface hopping variables
        self.rstate_old = self.rstate

        self.prob = np.zeros(self.mol.nst)
        self.acc_prob = np.zeros(self.mol.nst + 1)

        self.l_hop = False

        accum = 0.

        for ist in range(self.mol.nst):
            if (ist != self.rstate):
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

        # Write SHPROB file
        tmp = f'{istep + 1:9d}' + "".join([f'{self.prob[ist]:15.8f}' for ist in range(self.mol.nst)])
        typewriter(tmp, unixmd_dir, "SHPROB")

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

    def evaluate_hop(self, bo_list, unixmd_dir, istep):
        """ Routine to evaluate hopping and velocity rescaling

            :param integer,list bo_list: list of BO states for BO calculation
            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        if (self.l_hop):
            pot_diff = self.mol.states[self.rstate].energy - self.mol.states[self.rstate_old].energy
            if (self.mol.ekin_qm < pot_diff):
                self.l_hop = False
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
            self.event["HOP"].append(f"Hopping {self.rstate_old} -> {self.rstate}")

        # Write SHSTATE file
        tmp = f'{istep + 1:9d}{"":14s}{self.rstate}'
        typewriter(tmp, unixmd_dir, "SHSTATE")

    def decoherence_scheme(self):
        """ Routine to simple decoherence correction
        """
        if (self.decoherence_method == "IDC" and self.l_hop == True):
            if (self.propagation == "density"):
#                self.mol.rho = np.zeros((self.mol.nst, self.mol.nst), dtype=np.complex_)
                self.mol.rho[:] = 0.
                self.mol.rho.real[self.rstate, self.rstate] = 1.0 
            elif (self.propagation == "coefficient"):
                for ist in range(self.mol.nst):
                    self.mol.states[ist].coef = 0. 
                self.mol.states[self.rstate].coef = 1. + 0.j 
        elif (self.decoherence_method == "EDC"):
            raise ValueError (f"EDC is not yet implemented")

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

        # Print event in surface hopping
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event["HOP"] = []

