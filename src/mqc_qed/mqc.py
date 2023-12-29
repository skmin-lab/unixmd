from __future__ import division
from misc import fs_to_au, au_to_A, call_name, typewriter
import textwrap, datetime
import numpy as np
import os, shutil


class MQC_QED(object):
    """ Class for nuclear/electronic propagator used in MQC dynamics coupled to confined cavity mode

        :param object polariton: Polariton object
        :param object thermostat: Thermostat type
        :param integer istate: Initial adiabatic state
        :param double dt: Time interval
        :param integer nsteps: Nuclear step
        :param integer nesteps: Electronic step
        :param string elec_object: Electronic equation of motions
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: Logical to print BO population and coherence
        :param boolean l_adj_nac: Logical to adjust nonadiabatic coupling
        :param boolean l_adj_tdp: Adjust transition dipole moments to align the phases
        :param init_coef: Initial BO coefficient
        :type init_coef: Double, list or complex, list
        :param string unit_dt: Unit of time step (fs = femtosecond, au = atomic unit)
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, polariton, thermostat, istate, dt, nsteps, nesteps, elec_object, \
        propagator, l_print_dm, l_adj_nac, l_adj_tdp, init_coef, unit_dt, out_freq, verbosity):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize Polariton object
        self.pol = polariton

        # Initialize Thermostat object
        self.thermo = thermostat

        # Initialize input values
        self.istate = istate
        self.nsteps = nsteps
        self.nesteps = nesteps

        # Initialize time step
        self.istep = -1
        self.fstep = -1

        # Decide unit of time step
        self.unit_dt = unit_dt.lower()
        if (self.unit_dt == 'au'):
            self.dt = dt
        elif (self.unit_dt == 'fs'):
            self.dt = dt * fs_to_au
        else:
            error_message = "Invalid unit for time step!"
            error_vars = f"unit_dt = {unit_dt}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # None for BOMD case
        self.elec_object = elec_object
        if (self.elec_object != None):
            self.elec_object = self.elec_object.lower()

        if not (self.elec_object in [None, "coefficient", "density"]):
            error_message = "Invalid electronic object!"
            error_vars = f"elec_object = {self.elec_object}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        if not (self.elec_object in [None, "coefficient"]):
            error_message = "Only coefficient object is available for QED calculation!"
            error_vars = f"elec_object = {self.elec_object}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.propagator = propagator
        if (self.propagator != None):
            self.propagator = self.propagator.lower()

        if not (self.propagator in [None, "rk4"]):
            error_message = "Invalid electronic propagator!"
            error_vars = f"propagator = {self.propagator}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.l_print_dm = l_print_dm

        self.l_adj_nac = l_adj_nac
        self.l_adj_tdp = l_adj_tdp

        self.rforce = np.zeros((self.pol.nat, self.pol.ndim))

        self.out_freq = out_freq
        self.verbosity = verbosity

        # Initialize coefficients and densities
        self.pol.get_coefficient(init_coef, self.istate)

    def run_init(self, qed, qm, mm, output_dir, l_save_qed_log, l_save_qm_log, l_save_mm_log, \
        l_save_scr, restart):
        """ Initialize MQC dynamics

            :param object qed: QED object containing cavity-molecule interaction
            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string output_dir: Location of input directory
            :param boolean l_save_qed_log: Logical for saving QED calculation log
            :param boolean l_save_qm_log: Logical for saving QM calculation log
            :param boolean l_save_mm_log: Logical for saving MM calculation log
            :param boolean l_save_scr: Logical for saving scratch directory
            :param string restart: Option for controlling dynamics restarting
        """
        # Check whether the restart option is right
        if (restart != None):
            restart = restart.lower()

        if not (restart in [None, "write", "append"]):
            error_message = "Invalid restart option!"
            error_vars = f"restart = {restart}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Check if NACVs are calculated for Ehrenfest dynamics
        if (self.md_type == "Eh" and self.pol.l_nacme):
            error_message = "Ehrenfest dynamics needs evaluation of NACVs, check your QM object!"
            error_vars = f"(QM) qm_prog.qm_method = {qm.qm_prog}.{qm.qm_method}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Exception for QED with QM/MM dynamics
        # TODO : not yet tested
        if (mm != None):
            error_message = "QED using polariton object is not compatible with QM/MM now!"
            error_vars = f"mm = {mm}"
            raise NotImplementedError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Exception for QED with thermostat
        # TODO : thermostat is not yet tested with QED
        if (self.thermo != None):
            error_message = "QED using polariton object is not compatible with thermostat now!"
            error_vars = f"thermostat = {self.thermo}"
            raise NotImplementedError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Check compatibility of variables for QM and MM calculation
        if ((self.pol.l_qmmm and mm == None) or (not self.pol.l_qmmm and mm != None)):
            error_message = "Both logical for QM/MM and MM object is necessary!"
            error_vars = f"Polariton.l_qmmm = {self.pol.l_qmmm}, mm = {mm}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
        if (self.pol.l_qmmm and mm != None):
            self.check_qmmm(qm, mm)

        # Exception for CTMQC/Ehrenfest with QM/MM
#        if ((self.md_type in ["CT", "Eh"]) and (mm != None)):
#            error_message = "QM/MM calculation is not compatible with CTMQC or Ehrenfest now!"
#            error_vars = f"mm = {mm}"
#            raise NotImplementedError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Exception for CTMQC/Ehrenfest
        # TODO : not yet tested
        if (self.md_type in ["CT", "Eh"]):
            error_message = "Ehrenfest and CTMQC are not available for QED!"
            error_vars = f"mqc_qed = {self.md_type}"
            raise NotImplementedError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Warning for SH/SHXF
        if (self.md_type in ["SH", "SHXF"]):
            if (not self.pol.l_pnacme and self.hop_rescale == "energy"):
                print ("\n\n WARNING: Isotropic rescaling for accepted hop does not require the computation of pNACVs! \n\n", flush=True)
            if (not self.pol.l_pnacme and self.hop_reject == "keep"):
                print ("\n\n WARNING: Keep rescaling for frustrated hop does not require the computation of pNACVs! \n\n", flush=True)

        # Check compatibility between QED and QM objects
        self.check_qed(qed, qm)

        # Set directory information
        output_dir = os.path.expanduser(output_dir)
        base_dir = []
        unixmd_dir = []
        qed_log_dir = []
        qm_log_dir = []
        mm_log_dir = [None]

        if (self.pol.l_qmmm and mm != None):
            mm_log_dir = []

        dir_tmp = os.path.join(os.getcwd(), output_dir)
        if (self.md_type != "CT"):
            base_dir.append(dir_tmp)
        else:
            for itraj in range(self.ntrajs):
                itraj_dir = os.path.join(dir_tmp, f"TRAJ_{itraj + 1:0{self.digit}d}")
                base_dir.append(itraj_dir)

        for idir in base_dir:
            unixmd_dir.append(os.path.join(idir, "md"))
            qed_log_dir.append(os.path.join(idir, "qed_log"))
            qm_log_dir.append(os.path.join(idir, "qm_log"))
            if (self.pol.l_qmmm and mm != None):
                mm_log_dir.append(os.path.join(idir, "mm_log"))

        # Check and make directories
        if (restart == "append"):
            # For MD output directory
            for md_idir in unixmd_dir:
                if (not os.path.exists(md_idir)):
                    error_message = f"Directory {md_idir} to be appended for restart not found!"
                    error_vars = f"restart = {restart}, output_dir = {output_dir}"
                    raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

            # For QED output directory
            if (l_save_qed_log):
                for qed_idir in qed_log_dir:
                    if (not os.path.exists(qed_idir)):
                        os.makedirs(qed_idir)

            # For QM output directory
            if (l_save_qm_log):
                for qm_idir in qm_log_dir:
                    if (not os.path.exists(qm_idir)):
                        os.makedirs(qm_idir)

            # For MM output directory
            if (self.pol.l_qmmm and mm != None):
                if (l_save_mm_log):
                    for mm_idir in mm_log_dir:
                        if (not os.path.exists(mm_idir)):
                            os.makedirs(mm_idir)
        else:
            # For MD output directory
            for md_idir in unixmd_dir:
                if (os.path.exists(md_idir)):
                    shutil.move(md_idir, md_idir + "_old_" + str(os.getpid()))
                os.makedirs(md_idir)

                self.touch_file(md_idir)

            # For QED output directory
            for qed_idir in qed_log_dir:
                if (os.path.exists(qed_idir)):
                    shutil.move(qed_idir, qed_idir + "_old_" + str(os.getpid()))
                if (l_save_qed_log):
                    os.makedirs(qed_idir)

            # For QM output directory
            for qm_idir in qm_log_dir:
                if (os.path.exists(qm_idir)):
                    shutil.move(qm_idir, qm_idir + "_old_" + str(os.getpid()))
                if (l_save_qm_log):
                    os.makedirs(qm_idir)

            # For MM output directory
            for mm_idir in mm_log_dir:
                if (self.pol.l_qmmm and mm != None):
                    if (os.path.exists(mm_idir)):
                        shutil.move(mm_idir, mm_idir + "_old_" + str(os.getpid()))
                    if (l_save_mm_log):
                        os.makedirs(mm_idir)

        os.chdir(base_dir[0])

        if (self.md_type != "CT"):
            return base_dir[0], unixmd_dir[0], qed_log_dir[0], qm_log_dir[0], mm_log_dir[0]
        else:
            return base_dir, unixmd_dir, qed_log_dir, qm_log_dir, mm_log_dir

    def cl_update_position(self):
        """ Routine to update nuclear positions
        """
        self.pol.vel += 0.5 * self.dt * self.rforce / np.column_stack([self.pol.mass] * self.pol.ndim)
        self.pol.pos += self.dt * self.pol.vel

    def cl_update_velocity(self):
        """ Routine to update nuclear velocities
        """
        self.pol.vel += 0.5 * self.dt * self.rforce / np.column_stack([self.pol.mass] * self.pol.ndim)
        self.pol.update_kinetic()

#    def calculate_temperature(self):
#        """ Routine to calculate current temperature
#        """
#        pass
#        #self.temperature = self.pol.ekin * 2 / float(self.pol.ndof) * au_to_K

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self, qed, qm, mm, restart):
        """ Routine to print the initial information of dynamics

            :param object qed: QED object containing cavity-molecule interaction
            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string restart: Option for controlling dynamics restarting
        """
        # Print PyUNIxMD version
        cur_time = datetime.datetime.now()
        cur_time = cur_time.strftime("%Y-%m-%d %H:%M:%S")
        prog_info = textwrap.dedent(f"""\
        {"-" * 68}

        {"PyUNIxMD version 20.1":>43s}

        {"< Developers >":>40s}
        {" " * 4}Seung Kyu Min,  In Seong Lee,  Jong-Kwon Ha,  Daeho Han,
        {" " * 4}Kicheol Kim,  Tae In Kim,  Sung Wook Moon

        {"-" * 68}

        {" " * 4}Please cite PyUNIxMD as follows:
        {" " * 4}I. S. Lee, J.-K. Ha, D. Han, T. I. Kim, S. W. Moon, & S. K. Min.
        {" " * 4}PyUNIxMD: A Python-based excited state molecular dynamics package.
        {" " * 4}Journal of Computational Chemistry, 42:1755-1766. 2021

        {" " * 4}PyUNIxMD begins on {cur_time}
        """)
        print (prog_info, flush=True)

        # Print restart info
        if (restart != None):
            restart_info = textwrap.indent(textwrap.dedent(f"""\
            Dynamics is restarted from the last step of a previous dynamics.
            Restart Mode: {restart}
            """), "    ")
            print (restart_info, flush=True)

        # Print self.pol information: coordinate, velocity
        if (self.md_type != "CT"):
            self.pol.print_init(mm)
        else:
            for itraj, pol in enumerate(self.pols):
                pol.print_init(mm)

        # Print QED information
        qed.print_init()

        # Print dynamics information
        dynamics_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Dynamics Information":>43s}
        {"-" * 68}
        """)

        # Print QM information
        dynamics_info += textwrap.indent(textwrap.dedent(f"""\
          QM Program               = {qm.qm_prog:>16s}
          QM Method                = {qm.qm_method:>16s}
        """), "  ")
        if (self.pol.l_qmmm and mm != None):
            dynamics_info += textwrap.indent(textwrap.dedent(f"""\
              MM Program               = {mm.mm_prog:>16s}
              QMMM Scheme              = {mm.scheme:>16s}
            """), "  ")
            # Print charge embedding in MM program
            if (mm.embedding != None):
                dynamics_info += f"  Charge Embedding         = {mm.embedding:>16s}\n"
            else:
                dynamics_info += f"  Charge Embedding         = {'No':>16s}\n"
            # Print vdw interaction in MM program
            if (mm.vdw != None):
                dynamics_info += f"  VDW Interaction          = {mm.vdw:>16s}\n"
            else:
                dynamics_info += f"  VDW Interaction          = {'No':>16s}\n"

        dynamics_info += textwrap.indent(textwrap.dedent(f"""\

          MQC Method               = {self.md_type:>16s}
          Time Interval (fs)       = {self.dt / fs_to_au:16.6f}
          Initial State (0:GS)     = {self.istate:>16d}
          Nuclear Step             = {self.nsteps:>16d}
        """), "  ")

        if (self.md_type != "BOMD"):
            dynamics_info += f"  Electronic Step          = {self.nesteps:>16d}\n"
            dynamics_info += f"  Propagation Scheme       = {self.elec_object:>16s}\n"

        # Print surface hopping variables
        if (self.md_type in ["SH", "SHXF"]):
            dynamics_info += f"\n  Rescaling after Hop      = {self.hop_rescale:>16s}\n"
            dynamics_info += f"  Rescaling after Reject   = {self.hop_reject:>16s}\n"

        # Print ad-hoc decoherence variables
        if (self.md_type in ["SH"]):
            if (self.dec_correction != None):
                dynamics_info += f"\n  Decoherence Scheme       = {self.dec_correction:>16s}\n"
                if (self.dec_correction == "edc"):
                    dynamics_info += f"  Energy Constant          = {self.edc_parameter:>16.6f}\n"

        # Print XF variables
        if (self.md_type == "SHXF"):
            # Print density threshold used in decoherence term
            dynamics_info += f"\n  Density Threshold        = {self.rho_threshold:>16.6f}"
            # Print sigma values
            if (isinstance(self.sigma, float)):
                dynamics_info += f"\n  Sigma                    = {self.sigma:16.3f}\n"
            elif (isinstance(self.sigma, list)):
                dynamics_info += f"\n  Sigma (1:N)              =\n"
                nlines = int(self.aux.nat / 6)
                if (self.aux.nat % 6 != 0):
                    nlines += 1
                sigma_info = ""
                for iline in range(nlines):
                    iline1 = iline * 6
                    iline2 = (iline + 1) * 6
                    if (iline2 > self.aux.nat):
                        iline2 = self.aux.nat
                    sigma_info += f"  {iline1 + 1:>3d}:{iline2:<3d};"
                    sigma_info += "".join([f'{sigma:7.3f}' for sigma in self.sigma[iline1:iline2]])
                    sigma_info += "\n"
                dynamics_info += sigma_info
            # Print auxiliary trajectory setting
            if (self.l_econs_state):
                dynamics_info += f"\n  Aux. Total E             = {'Real Total E':>16s}\n"
            else:
                dynamics_info += f"\n  Aux. Total E             = {'Real Total E':>16s}" \
                    + f"\n                             {'+ Pot. E Diff.':>16s}\n"

        # Print system information
        dynamics_info += f"\n  Output Frequency         = {self.out_freq:>16d}\n"
        dynamics_info += f"  Verbosity Level          = {self.verbosity:>16d}\n"

        print (dynamics_info, flush=True)

        # Print thermostat information
        if (self.thermo != None):
            self.thermo.print_init()
        else:
            thermostat_info = "  No Thermostat: Total energy is conserved!\n"
            print (thermostat_info, flush=True)

    def touch_file(self, unixmd_dir):
        """ Routine to write PyUNIxMD output files

            :param string unixmd_dir: Directory where MD output files are written
        """
        # Energy information file header
        tmp = f'{"#":5s}{"Step":9s}{"Kinetic(H)":15s}{"Potential(H)":15s}{"Total(H)":15s}' + \
            "".join([f'E({ist})(H){"":8s}' for ist in range(self.pol.pst)])
        typewriter(tmp, unixmd_dir, "MDENERGY", "w")

        if (self.md_type != "BOMD"):
            # polaritonic state coefficents, densities file header
            if (self.elec_object == "coefficient"):
                tmp = f'{"#":5s} Polaritonic State Coefficients: state Re-Im; see the manual for detail orders'
                typewriter(tmp, unixmd_dir, "QEDCOEFA", "w")
                typewriter(tmp, unixmd_dir, "QEDCOEFD", "w")
                if (self.l_print_dm):
                    tmp = f'{"#":5s} Density Matrix: population Re; see the manual for detail orders'
                    typewriter(tmp, unixmd_dir, "QEDPOPA", "w")
                    tmp = f'{"#":5s} Density Matrix: coherence Re-Im; see the manual for detail orders'
                    typewriter(tmp, unixmd_dir, "QEDCOHA", "w")
                    tmp = f'{"#":5s} Density Matrix: population Re; see the manual for detail orders'
                    typewriter(tmp, unixmd_dir, "QEDPOPD", "w")
                    tmp = f'{"#":5s} Density Matrix: coherence Re-Im; see the manual for detail orders'
                    typewriter(tmp, unixmd_dir, "QEDCOHD", "w")

            # NACME file header
            tmp = f'{"#":5s}Non-Adiabatic Coupling Matrix Elements: off-diagonal'
            typewriter(tmp, unixmd_dir, "PNACME", "w")

            # DOTPOPNACD file header
            if (self.verbosity >= 1):
                tmp = f'{"#":5s} Time-derivative Density Matrix by NAC: population; see the manual for detail orders'
                typewriter(tmp, unixmd_dir, "DOTPOPNACD", "w")

        # file header for SH-based methods
        if (self.md_type in ["SH", "SHXF"]):
            tmp = f'{"#":5s}{"Step":8s}{"Running State":10s}'
            typewriter(tmp, unixmd_dir, "SHSTATE", "w")

            tmp = f'{"#":5s}{"Step":12s}' + "".join([f'Prob({ist}){"":8s}' for ist in range(self.pol.pst)])
            typewriter(tmp, unixmd_dir, "SHPROB", "w")

        # file header for XF-based methods
        if (self.md_type in ["SHXF"]):
            if (self.verbosity >= 1):
                tmp = f'{"#":5s} Time-derivative Density Matrix by decoherence: population; see the manual for detail orders'
                typewriter(tmp, unixmd_dir, "DOTPOPDECD", "w")

    def write_md_output(self, unixmd_dir, istep):
        """ Write output files

            :param string unixmd_dir: Directory where MD output files are written
            :param integer istep: Current MD step
        """
        # Write MOVIE.xyz file including positions and velocities
        tmp = f'{self.pol.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)' + \
            "".join(["\n" + f'{self.pol.symbols[iat]:5s}' + \
            "".join([f'{self.pol.pos[iat, isp] * au_to_A:15.8f}' for isp in range(self.pol.ndim)]) + \
            "".join([f"{self.pol.vel[iat, isp]:15.8f}" for isp in range(self.pol.ndim)]) for iat in range(self.pol.nat)])
        typewriter(tmp, unixmd_dir, "MOVIE.xyz", "a")

        # Write MDENERGY file including several energy information
        tmp = f'{istep + 1:9d}{self.pol.ekin:15.8f}{self.pol.epot:15.8f}{self.pol.etot:15.8f}' \
            + "".join([f'{states.energy:15.8f}' for states in self.pol.pol_states])
        typewriter(tmp, unixmd_dir, "MDENERGY", "a")

        if (self.md_type != "BOMD"):
            # Write QEDCOEF, QEDPOP, QEDCOH files
            if (self.elec_object == "coefficient"):
                # Adiabatic quantities
                tmp = f'{istep + 1:9d}' + "".join([f'{states.coef_a.real:15.8f}{states.coef_a.imag:15.8f}' \
                    for states in self.pol.pol_states])
                typewriter(tmp, unixmd_dir, "QEDCOEFA", "a")
                if (self.l_print_dm):
                    tmp = f'{istep + 1:9d}' + "".join([f'{self.pol.rho_a.real[ist, ist]:15.8f}' for ist in range(self.pol.pst)])
                    typewriter(tmp, unixmd_dir, "QEDPOPA", "a")
                    tmp = f'{istep + 1:9d}' + "".join([f"{self.pol.rho_a.real[ist, jst]:15.8f}{self.pol.rho_a.imag[ist, jst]:15.8f}" \
                        for ist in range(self.pol.pst) for jst in range(ist + 1, self.pol.pst)])
                    typewriter(tmp, unixmd_dir, "QEDCOHA", "a")
                # Diabatic quantities
                tmp = f'{istep + 1:9d}' + "".join([f'{states.coef_d.real:15.8f}{states.coef_d.imag:15.8f}' \
                    for states in self.pol.pol_states])
                typewriter(tmp, unixmd_dir, "QEDCOEFD", "a")
                if (self.l_print_dm):
                    tmp = f'{istep + 1:9d}' + "".join([f'{self.pol.rho_d.real[ist, ist]:15.8f}' for ist in range(self.pol.pst)])
                    typewriter(tmp, unixmd_dir, "QEDPOPD", "a")
                    tmp = f'{istep + 1:9d}' + "".join([f"{self.pol.rho_d.real[ist, jst]:15.8f}{self.pol.rho_d.imag[ist, jst]:15.8f}" \
                        for ist in range(self.pol.pst) for jst in range(ist + 1, self.pol.pst)])
                    typewriter(tmp, unixmd_dir, "QEDCOHD", "a")

            # Write NACME file
            tmp = f'{istep + 1:10d}' + "".join([f'{self.pol.pnacme[ist, jst]:15.8f}' \
                for ist in range(self.pol.pst) for jst in range(ist + 1, self.pol.pst)])
            typewriter(tmp, unixmd_dir, "PNACME", "a")

            # Write NACV file
            if (not self.pol.l_pnacme and self.verbosity >= 2):
                for ist in range(self.pol.pst):
                    for jst in range(ist + 1, self.pol.pst):
                        tmp = f'{self.pol.nat_qm:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}NACV' + \
                            "".join(["\n" + f'{self.pol.symbols[iat]:5s}' + \
                            "".join([f'{self.pol.pnac[ist, jst, iat, isp]:15.8f}' for isp in range(self.pol.ndim)]) for iat in range(self.pol.nat_qm)])
                        typewriter(tmp, unixmd_dir, f"PNACV_{ist}_{jst}", "a")

    def write_final_xyz(self, unixmd_dir, istep):
        """ Write final positions and velocities

            :param string unixmd_dir: Directory where MD output files are written
            :param integer istep: Current MD step
        """
        # Write FINAL.xyz file including positions and velocities
        tmp = f'{self.pol.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)'
        for iat in range(self.pol.nat):
            tmp += "\n" + f'{self.pol.symbols[iat]:5s}' + \
                "".join([f'{self.pol.pos[iat, isp] * au_to_A:15.8f}' for isp in range(self.pol.ndim)]) \
                + "".join([f"{self.pol.vel[iat, isp]:15.8f}" for isp in range(self.pol.ndim)])

        typewriter(tmp, unixmd_dir, "FINAL.xyz", "w")

    def check_qmmm(self, qm, mm):
        """ Routine to check compatibility between QM and MM objects

            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
        """
        # Now check MM object
        if (mm.mm_prog == "Tinker"):
            # Now check QM object
            if (qm.qm_prog == "dftbplus"):
                if (qm.qm_method == "SSR"):
                    do_qmmm = True
                else:
                    do_qmmm = False
            else:
                do_qmmm = False
        else:
            do_qmmm = False

        if (do_qmmm):
            if (qm.embedding != mm.embedding):
                error_message = "Inconsistent charge embedding between QM and MM objects!"
                error_vars = f"(QM) embedding = {qm.embedding}, (MM) embedding = {mm.embedding}"
                raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Incompatible QM and MM objects for QM/MM calculation!"
            error_vars = f"(QM) qm_prog.qm_method = {qm.qm_prog}.{qm.qm_method}, (MM) mm_prog = {mm.mm_prog}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

    def check_qed(self, qed, qm):
        """ Routine to check compatibility between QED and QM objects

            :param object qed: QED object containing cavity-molecule interaction
            :param object qm: QM object containing on-the-fly calculation infomation
        """
        # Now check QED object
        if (qed.qed_method == "Jaynes_Cummings"):
            # Now check QM object
            if (qm.qm_prog == "dftbplus"):
                if (qm.qm_method == "SSR"):
                    do_qed = True
                else:
                    do_qed = False
            else:
                do_qed = False
        else:
            do_qed = False

        if (not do_qed):
            error_message = "Incompatible QED and QM objects!"
            error_vars = f"(QED) qed_method = {qed.qed_method}, (QM) qm_prog.qm_method = {qm.qm_prog}.{qm.qm_method}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")


