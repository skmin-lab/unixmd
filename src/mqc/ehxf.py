from __future__ import division
from build.el_propagator_xf import el_run
from mqc.mqc import MQC
from misc import eps, au_to_K, au_to_A, call_name, typewriter
import os, shutil, textwrap
import numpy as np
import pickle

class Auxiliary_Molecule(object):
    """ Class for auxiliary molecule that is used for the calculation of decoherence term

        :param object molecule: molecule object
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
    """ Class for Ehrenfest-XF dynamics

        :param object molecule: molecule object
        :param object thermostat: thermostat type
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string obj: Representation for electronic state
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: logical to print BO population and coherence
        :param boolean l_adj_nac: logical to adjust nonadiabatic coupling
        :param double threshold: electronic density threshold for decoherence term calculation
        :param sigma: width of nuclear wave packet of auxiliary trajectory
        :type sigma: double or double,list
        :param init_coef: initial BO coefficient
        :type init_coef: double, list or complex, list
        :param boolean l_econs_state: logical to use identical total energies for auxiliary trajectories
        :param string unit_dt: unit of time step (fs = femtosecond, au = atomic unit)
        :param integer out_freq: frequency of printing output
        :param integer verbosity: verbosity of output
    """
    def __init__(self, molecule, thermostat=None, istate=0, dt=0.5, nsteps=1000, nesteps=20, \
        obj="density", propagator="rk4", l_print_dm=True, l_adj_nac=True, \
        threshold=0.01, sigma=None, l_deco_force=False, init_coef=None, \
        l_econs_state=True, unit_dt="fs", out_freq=1, verbosity=0):
        # Initialize input values
        super().__init__(molecule, thermostat, istate, dt, nsteps, nesteps, \
            obj, propagator, l_print_dm, l_adj_nac, init_coef, unit_dt, out_freq, verbosity)

        # Initialize XF related variables
        self.l_coh = []
        self.l_first = []
        for ist in range(self.mol.nst):
            self.l_coh.append(False)
            self.l_first.append(False)
        self.threshold = threshold
        self.sigma = sigma

        if (isinstance(self.sigma, float)):
            # uniform value for sigma
            pass
        elif (isinstance(self.sigma, list)):
            # atom-resolved values for sigma
            if (len(self.sigma) != self.mol.nat_qm):
                raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong number of elements of sigma given! {self.sigma}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong type for sigma given! {self.sigma}")

        self.upper_th = 1. - self.threshold
        self.lower_th = self.threshold

        self.l_deco_force = l_deco_force

        self.l_econs_state = l_econs_state

        # Initialize auxiliary molecule object
        self.aux = Auxiliary_Molecule(self.mol)
        self.pos_0 = np.zeros((self.aux.nat, self.aux.ndim))
        self.phase = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))

        # Debug variables
        self.dotpopd = np.zeros(self.mol.nst)
        self.qmom = np.zeros((self.aux.nat, self.aux.ndim))

        # Initialize event to print
        self.event = {"DECO": []}

    def run(self, qm, mm=None, output_dir="./", l_save_qm_log=False, l_save_mm_log=False, l_save_scr=True, restart=None):
        """ Run MQC dynamics according to Ehrenfest-XF dynamics

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param string output_dir: Name of directory where outputs to be saved.
            :param boolean l_save_qm_log: logical for saving QM calculation log
            :param boolean l_save_mm_log: logical for saving MM calculation log
            :param boolean l_save_scr: logical for saving scratch directory
            :param string restart: option for controlling dynamics restarting
        """
        # Initialize UNI-xMD
        base_dir, unixmd_dir, qm_log_dir, mm_log_dir =\
             self.run_init(qm, mm, output_dir, l_save_qm_log, l_save_mm_log, l_save_scr, restart)
        bo_list = [ist for ist in range(self.mol.nst)]
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

            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()

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

            self.cl_update_position()

            self.mol.backup_bo()
            self.mol.reset_bo(qm.calc_coupling)
            qm.get_data(self.mol, base_dir, bo_list, self.dt, istep, calc_force_only=False)
            if (self.mol.l_qmmm and mm != None):
                mm.get_data(self.mol, base_dir, bo_list, istep, calc_force_only=False)

            self.mol.adjust_nac()

            self.cl_update_velocity()

            self.mol.get_nacme()

            el_run(self)

            if (self.thermo != None):
                self.thermo.run(self)

            self.update_energy()

            self.check_decoherence()
            self.check_coherence()
            self.aux_propagator()
            self.get_phase()

            if ((istep + 1) % self.out_freq == 0):
                self.write_md_output(unixmd_dir, istep)
            if ((istep + 1) % self.out_freq == 0 or len(self.event["DECO"]) > 0):
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

    def calculate_force(self):
        """ Calculate the Ehrenfest-XF force
        """
        self.rforce = np.zeros((self.mol.nat, self.mol.ndim))

        for ist, istate in enumerate(self.mol.states):
            self.rforce += istate.force * self.mol.rho.real[ist, ist]

        for ist in range(self.mol.nst):
            for jst in range(ist + 1, self.mol.nst):
                self.rforce += 2. * self.mol.nac[ist, jst] * self.mol.rho.real[ist, jst] \
                    * (self.mol.states[ist].energy - self.mol.states[jst].energy)

        if (self.l_deco_force):
            # Calculate quantum momentum
            qmom = np.zeros((self.aux.nat, self.mol.ndim))
            for ist in range(self.mol.nst):
                for iat in range(self.aux.nat):
                    qmom[iat] += 0.5 * self.mol.rho.real[ist, ist] * (self.pos_0[iat] - self.aux.pos[ist, iat]) \
                        / self.sigma[iat] ** 2 / self.mol.mass[iat]

            # Calculate XF force
            for ist in range(self.mol.nst):
                for jst in range(self.mol.nst):
                    self.rforce -= 2. * self.mol.rho.real[ist, ist] * self.mol.rho.real[jst, jst] \
                        * np.sum(qmom * (self.phase[ist] - self.phase[jst])) * self.phase[jst]

    def update_energy(self):
        """ Routine to update the energy of molecules in Ehrenfest-XF dynamics
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
        self.phase = np.zeros((self.mol.nst, self.aux.nat, self.aux.ndim))
        self.mol.rho = np.zeros((self.mol.nst, self.mol.nst), dtype=np.complex_)
        self.mol.rho[one_st, one_st] = 1. + 0.j

        self.l_coh = [False] * self.mol.nst
        self.l_first = [False] * self.mol.nst

        if (self.obj == "coefficient"):
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
                    self.aux.pos[ist] += self.aux.vel[ist] * self.dt
            else:
                self.aux.pos[ist] = self.mol.pos[0:self.aux.nat]

        self.pos_0 = np.copy(self.mol.pos[0:self.aux.nat])

        # Get auxiliary velocity
        self.aux.vel_old = np.copy(self.aux.vel)
        for ist in range(self.mol.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    alpha = self.mol.ekin_qm
                    if (not self.l_econs_state):
                        alpha += self.mol.epot - self.mol.states[ist].energy
                else:
                    ekin_old = np.sum(0.5 * self.aux.mass * np.sum(self.aux.vel_old[ist] ** 2, axis=1))
                    alpha = ekin_old + self.mol.states[ist].energy_old - self.mol.states[ist].energy
                if (alpha < 0.):
                    alpha = 0.

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

    def append_sigma(self):
        """ Routine to append sigma values when single float number is provided
        """
        # Create a list from single float number
        if (isinstance(self.sigma, float)):
            sigma = self.sigma
            self.sigma = self.aux.nat * [sigma]

    def write_md_output(self, unixmd_dir, istep):
        """ Write output files

            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        # Write the common part
        super().write_md_output(unixmd_dir, istep)

        # Write decoherence information
        self.write_deco(unixmd_dir, istep)

    def write_deco(self, unixmd_dir, istep):
        """ Write XF-based decoherence information

            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        # Write time-derivative density matrix elements in DOTPOTD
        tmp = f'{istep + 1:9d}' + "".join([f'{pop:15.8f}' for pop in self.dotpopd])
        typewriter(tmp, unixmd_dir, "DOTPOPD", "a")

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
                    tmp = f'{self.aux.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Momentum (au)' + \
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

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
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

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        print (dynamics_step_info, flush=True)

    def print_step(self, istep):
        """ Routine to print each steps infomation about dynamics

            :param integer istep: current MD step
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

        # Print DEBUG1 for each step
        if (self.verbosity >= 1):
            DEBUG1 = f" DEBUG1{istep + 1:>7d}"
            for ist in range(self.mol.nst):
                DEBUG1 += f"{self.mol.states[ist].energy:17.8f} "
            print (DEBUG1, flush=True)

        # Print event in EhXF
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event["DECO"] = []
