from __future__ import division
from build.el_propagator_xf import el_run
from mqc.mqc import MQC
from fileio import touch_file, write_md_output, write_final_xyz, write_aux_movie, typewriter
from misc import eps, au_to_K, call_name
import os, shutil, textwrap
import numpy as np

class Auxiliary_Molecule(object):
    """ Class for auxiliary molecule that is used for the calculation of decoherence term

        :param object molecule: molecule object
    """
    def __init__(self, molecule):
        # Initialize auxiliary molecule
        self.nat = molecule.nat_qm
        self.nsp = molecule.nsp

        self.mass = np.copy(molecule.mass)
        
        self.pos = np.zeros((molecule.nst, self.nat, self.nsp))
        self.vel = np.zeros((molecule.nst, self.nat, self.nsp))
        self.vel_old = np.copy(self.vel)


class EhXF(MQC):
    """ Class for Ehrenfest-XF dynamics

        :param object molecule: molecule object
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string propagation: propagation scheme
        :param string solver: propagation solver
        :param boolean l_pop_print: logical to print BO population and coherence
        :param boolean l_adjnac: logical to adjust nonadiabatic coupling
        :param double threshold: electronic density threshold for decoherence term calculation
        :param wsigma: width of nuclear wave packet of auxiliary trajectory
        :type wsigma: double or double,list
        :param coefficient: initial BO coefficient
        :type coefficient: double, list or complex, list
        :param boolean l_state_wise: logical to use state-wise total energies for auxiliary trajectories
        :param string unit_dt: unit of time step (fs = femtosecond, au = atomic unit)
    """
    def __init__(self, molecule, istate=0, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", solver="rk4", l_pop_print=False, l_adjnac=True, \
        threshold=0.01, wsigma=None, l_qmom_force=False, coefficient=None, \
        l_state_wise=False, unit_dt="fs"):
        # Initialize input values
        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, solver, l_pop_print, l_adjnac, coefficient, unit_dt)

        # Initialize XF related variables
        self.l_coh = []
        self.l_first = []
        for ist in range(molecule.nst):
            self.l_coh.append(False)
            self.l_first.append(False)
        self.threshold = threshold
        self.wsigma = wsigma

        if (isinstance(self.wsigma, float)):
            # uniform value for wsigma
            pass
        elif (isinstance(self.wsigma, list)):
            # atom-resolved values for wsigma
            if (len(self.wsigma) != molecule.nat_qm):
                raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong number of elements of sigma given! {self.wsigma}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Wrong type for sigma given! {self.wsigma}")

        self.upper_th = 1. - self.threshold
        self.lower_th = self.threshold

        self.l_qmom_force = l_qmom_force

        self.l_state_wise = l_state_wise

        # Initialize auxiliary molecule object
        self.aux = Auxiliary_Molecule(molecule)
        self.pos_0 = np.zeros((self.aux.nat, self.aux.nsp))
        self.phase = np.zeros((molecule.nst, self.aux.nat, self.aux.nsp))

        # Debug variables
        self.dotpopd = np.zeros(molecule.nst)
        
        # Initialize event to print
        self.event = {"DECO": []}

    def run(self, molecule, qm, mm=None, thermostat=None, input_dir="./", \
        save_QMlog=False, save_MMlog=False, save_scr=True, debug=0):
        """ Run MQC dynamics according to Ehrenfest-XF dynamics

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
        bo_list = [ist for ist in range(molecule.nst)]
        qm.calc_coupling = True

        touch_file(molecule, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, SH_chk=False, XF_chk=True)
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

        self.update_energy(molecule)

        self.check_decoherence(molecule)
        self.check_coherence(molecule)
        self.aux_propagator(molecule)
        self.get_phase(molecule)
        self.print_deco(molecule, unixmd_dir, istep=-1)

        write_md_output(molecule, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=-1)
        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                write_aux_movie(self.aux, unixmd_dir, ist, istep=-1) 
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

            el_run(self, molecule)

            if (thermostat != None):
                thermostat.run(molecule, self)

            self.update_energy(molecule)

            self.check_decoherence(molecule)
            self.check_coherence(molecule)
            self.aux_propagator(molecule)
            self.get_phase(molecule)
            self.print_deco(molecule, unixmd_dir, istep=istep)

            write_md_output(molecule, qm.calc_coupling, self.propagation, self.l_pop_print, unixmd_dir, istep=istep)
            for ist in range(molecule.nst):
                if (self.l_coh[ist]):
                    write_aux_movie(self.aux, unixmd_dir, ist, istep=istep) 
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

    def calculate_force(self, molecule):
        """ Calculate the Ehrenfest-XF force

            :param object molecule: molecule object
        """
        self.rforce = np.zeros((molecule.nat, molecule.nsp))

        for ist, istate in enumerate(molecule.states):
            self.rforce += istate.force * molecule.rho.real[ist, ist]

        for ist in range(molecule.nst):
            for jst in range(ist + 1, molecule.nst):
                self.rforce += 2. * molecule.nac[ist, jst] * molecule.rho.real[ist, jst] \
                    * (molecule.states[ist].energy - molecule.states[jst].energy)

        if (self.l_qmom_force):
            # Calculate quantum momentum
            qmom = np.zeros((self.aux.nat, molecule.nsp))
            for ist in range(molecule.nst):
                for iat in range(self.aux.nat):
                    qmom[iat] += 0.5 * molecule.rho.real[ist, ist] * (self.pos_0[iat] - self.aux.pos[ist, iat]) \
                        / self.wsigma[iat] ** 2 / molecule.mass[iat]

            # Calculate XF force
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    self.rforce -= 2. * molecule.rho.real[ist, ist] * molecule.rho.real[jst, jst] \
                        * np.sum(qmom * (self.phase[ist] - self.phase[jst])) * self.phase[jst]

    def update_energy(self, molecule):
        """ Routine to update the energy of molecules in Ehrenfest-XF dynamics

            :param object molecule: molecule object
        """
        # Update kinetic energy
        molecule.update_kinetic()
        molecule.epot = 0.
        for ist, istate in enumerate(molecule.states):
            molecule.epot += molecule.rho.real[ist, ist] * molecule.states[ist].energy
        molecule.etot = molecule.epot + molecule.ekin

    def check_decoherence(self, molecule):
        """ Routine to check if the electronic state is decohered

            :param object molecule: molecule object
        """
        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                rho = molecule.rho.real[ist, ist]
                if (rho > self.upper_th):
                    self.set_decoherence(molecule, ist)
                    self.event["DECO"].append(f"Destroy auxiliary trajectories: decohered to {ist} state")
                    return

    def check_coherence(self, molecule):
        """ Routine to check coherence among BO states

            :param object molecule: molecule object
        """
        count = 0
        tmp_st = ""
        for ist in range(molecule.nst):
            rho = molecule.rho.real[ist, ist]
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
            self.l_coh = [False] * molecule.nst
            self.l_first = [False] * molecule.nst

        if (len(tmp_st) >= 1):
            tmp_st = tmp_st.rstrip(', ')
            self.event["DECO"].append(f"Generate auxiliary trajectory on {tmp_st} state")

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
                    self.aux.pos[ist] = molecule.pos[0:self.aux.nat]
                else:
                    self.aux.pos[ist] += self.aux.vel[ist] * self.dt
            else:
                self.aux.pos[ist] = molecule.pos[0:self.aux.nat]

        self.pos_0 = np.copy(molecule.pos[0:self.aux.nat])

        # Get auxiliary velocity
        self.aux.vel_old = np.copy(self.aux.vel)
        for ist in range(molecule.nst):
            if (self.l_coh[ist]):
                if (self.l_first[ist]):
                    alpha = molecule.ekin_qm
                    if (not self.l_state_wise):
                        alpha += molecule.epot - molecule.states[ist].energy
                else:
                    ekin_old = np.sum(0.5 * self.aux.mass * np.sum(self.aux.vel_old[ist] ** 2, axis=1))
                    alpha = ekin_old + molecule.states[ist].energy_old - molecule.states[ist].energy
                if (alpha < 0.):
                    alpha = 0.
                
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

    def append_wsigma(self):
        """ Routine to append sigma values when single float number is provided
        """
        # Create a list from single float number
        if (isinstance(self.wsigma, float)):
            sigma = self.wsigma
            self.wsigma = self.aux.nat * [sigma]

    def print_deco(self, molecule, unixmd_dir, istep):
        """ Routine to print decoherence information

            :param object molecule: molecule object
            :param string unixmd_dir: unixmd directory
            :param integer istep: current MD step
        """
        tmp = f'{istep + 1:9d}' + "".join([f'{self.dotpopd[ist]:15.8f}' for ist in range(molecule.nst)])
        typewriter(tmp, unixmd_dir, "DOTPOPD")

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
        INIT = f" #INFO{'STEP':>8s}{'Kinetic(H)':>15s}{'Potential(H)':>15s}{'Total(H)':>13s}{'Temperature(K)':>17s}{'norm':>8s}"
        dynamics_step_info += INIT

        # Print DEBUG1 for each step
        if (debug >= 1):
            DEBUG1 = f" #DEBUG1{'STEP':>6s}"
            for ist in range(molecule.nst):
                DEBUG1 += f"{'Potential_':>14s}{ist}(H)"
            dynamics_step_info += "\n" + DEBUG1

        print (dynamics_step_info, flush=True)

    def print_step(self, molecule, debug, istep):
        """ Routine to print each steps infomation about dynamics

            :param object molecule: molecule object
            :param integer debug: verbosity level for standard output
            :param integer istep: current MD step
        """
        ctemp = molecule.ekin * 2. / float(molecule.dof) * au_to_K
        norm = 0.
        for ist in range(molecule.nst):
            norm += molecule.rho.real[ist, ist]

        # Print INFO for each step
        INFO = f" INFO{istep + 1:>9d} "
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

        # Print event in EhXF
        for category, events in self.event.items():
            if (len(events) != 0):
                for ievent in events:
                    print (f" {category}{istep + 1:>9d}  {ievent}", flush=True)
        self.event = {"DECO": []}

