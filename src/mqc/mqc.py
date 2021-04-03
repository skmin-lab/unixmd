from __future__ import division
from misc import fs_to_au, au_to_A, call_name, typewriter
import textwrap, datetime
import numpy as np
import os, shutil


class MQC(object):
    """ Class for nuclear/electronic propagator used in MQC dynamics

        :param object molecule: Molecule object
        :param object thermostat: Thermostat type
        :param integer istate: Initial adiabatic state
        :param double dt: Time interval
        :param integer nsteps: Nuclear step
        :param integer nesteps: Electronic step
        :param string obj: Representation for electronic state
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: Logical to print BO population and coherence
        :param boolean l_adj_nac: Logical to adjust nonadiabatic coupling
        :param initial_coef: Initial BO coefficient
        :type initial_coef: Double, list or complex, list
        :param string unit_dt: Unit of time step (fs = femtosecond, au = atomic unit)
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, thermostat, istate, dt, nsteps, nesteps, \
        obj, propagator, l_print_dm, l_adj_nac, init_coef, unit_dt, out_freq, verbosity):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize Molecule object
        self.mol = molecule

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
        if (unit_dt == 'au'):
            self.dt = dt
        elif (unit_dt == 'fs'):
            self.dt = dt * fs_to_au
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid unit for time step! {unit_dt}")

        # Check number of state and initial state
        if (self.istate >= self.mol.nst):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Index for initial state must be smaller than number of states! {self.istate}")

        # None for BOMD case
        self.obj = obj
        if not (self.obj in [None, "coefficient", "density"]):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'obj'! {self.obj}")
        self.propagator = propagator
        if not (self.propagator in [None, "rk4"]):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'propagator'! {self.propagator}")

        self.l_print_dm = l_print_dm

        self.l_adj_nac = l_adj_nac

        self.rforce = np.zeros((self.mol.nat, self.mol.ndim))

        self.out_freq = out_freq
        self.verbosity = verbosity

        # Initialize coefficients and densities
        self.mol.get_coefficient(init_coef, self.istate)

    def run_init(self, qm, mm, output_dir, save_qm_log, save_mm_log, save_scr, restart):
        """ Initialize MQC dynamics

            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
            :param string output_dir: Location of input directory
            :param boolean save_qm_log: Logical for saving QM calculation log
            :param boolean save_mm_log: Logical for saving MM calculation log
            :param boolean save_scr: Logical for saving scratch directory
            :param string restart: Option for controlling dynamics restarting
        """
        # Check whether the restart option is right
        if not (restart in [None, "write", "append"]):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'restart'! {restart}")

        # Check if NACVs are calculated for Ehrenfest dynamics
        if (self.md_type in ["Eh", "EhXF"] and self.mol.l_nacme):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Ehrenfest dynamics needs NACV! {self.mol.l_nacme}")

        # Check compatibility of variables for QM and MM calculation
        if ((self.mol.qmmm and mm == None) or (not self.mol.qmmm and mm != None)):
            raise ValueError (f"( {self.md_type}.{call_name()} ) Both self.mol.qmmm and mm object is necessary! {self.mol.qmmm} and {mm}")
        if (self.mol.qmmm and mm != None):
            self.check_qmmm(qm, mm)

        # Set directory information
        output_dir = os.path.expanduser(output_dir)
        base_dir = os.path.join(os.getcwd(), output_dir)
        unixmd_dir = os.path.join(base_dir, "md")
        qm_log_dir = os.path.join(base_dir, "qm_log")
        mm_log_dir = None
        if (self.mol.qmmm and mm != None):
            mm_log_dir = os.path.join(base_dir, "mm_log")

        # Check and make directories
        if (restart == "append"):
            if (not os.path.exists(unixmd_dir)):
                raise ValueError (f"( {self.md_type}.{call_name()} ) Directory to be appended for restart not found! {restart} and {unixmd_dir}")
            if (not os.path.exists(unixmd_dir) and save_qm_log):
                os.makedirs(qm_log_dir)
            if (self.mol.qmmm and mm != None):
                if (not os.path.exists(mm_log_dir) and save_mm_log):
                    os.makedirs(mm_log_dir)
        else:
            if (os.path.exists(unixmd_dir)):
                shutil.move(unixmd_dir, unixmd_dir + "_old_" + str(os.getpid()))
            os.makedirs(unixmd_dir)

            if (os.path.exists(qm_log_dir)):
                shutil.move(qm_log_dir, qm_log_dir + "_old_" + str(os.getpid()))
            if (save_qm_log):
                os.makedirs(qm_log_dir)

            if (self.mol.qmmm and mm != None):
                if (os.path.exists(mm_log_dir)):
                    shutil.move(mm_log_dir, mm_log_dir + "_old_" + str(os.getpid()))
                if (save_mm_log):
                    os.makedirs(mm_log_dir)

            self.touch_file(unixmd_dir)

        os.chdir(base_dir)

        return base_dir, unixmd_dir, qm_log_dir, mm_log_dir

    def cl_update_position(self):
        """ Routine to update nuclear positions

        """
        self.calculate_force()

        self.mol.vel += 0.5 * self.dt * self.rforce / np.column_stack([self.mol.mass] * self.mol.ndim)
        self.mol.pos += self.dt * self.mol.vel

    def cl_update_velocity(self):
        """ Routine to update nuclear velocities

        """
        self.calculate_force()

        self.mol.vel += 0.5 * self.dt * self.rforce / np.column_stack([self.mol.mass] * self.mol.ndim)
        self.mol.update_kinetic()

#    def calculate_temperature(self):
#        """ Routine to calculate current temperature
#        """
#        pass
#        #self.temperature = self.mol.ekin * 2 / float(self.mol.ndof) * au_to_K

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self, qm, mm, restart):
        """ Routine to print the initial information of dynamics

            :param object qm: QM object containing on-the-fly calculation infomation
            :param object mm: MM object containing MM calculation infomation
        """
        # Print UNI-xMD version
        cur_time = datetime.datetime.now()
        cur_time = cur_time.strftime("%Y-%m-%d %H:%M:%S")
        prog_info = textwrap.dedent(f"""\
        {"-" * 68}

        {"UNI-xMD version 20.1":>43s}

        {"< Developers >":>40s}
        {" " * 4}Seung Kyu Min,  In Seong Lee,  Jong-Kwon Ha,  Daeho Han,
        {" " * 4}Kicheol Kim,  Tae In Kim,  Sung Wook Moon

        {"-" * 68}

        {" " * 4}Please cite UNI-xMD as follows:
        {" " * 4}This is article

        {" " * 4}UNI-xMD begins on {cur_time}
        """)
        print (prog_info, flush=True)

        # Print restart info
        if (restart != None):
            restart_info = textwrap.indent(textwrap.dedent(f"""\
            Dynamics is restarted from the last step of a previous dynamics.
            Restart Mode: {restart}
            """), "    ")
            print (restart_info, flush=True)

        # Print molecule information: coordinate, velocity
        self.mol.print_init(mm)

        # Print dynamics information
        dynamics_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Dynamics Information":>43s}
        {"-" * 68}
          QM Program               = {qm.qm_prog:>16s}
          QM Method                = {qm.qm_method:>16s}
        """)
        if (self.mol.qmmm and mm != None):
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
            dynamics_info += f"  Propagation Scheme       = {self.obj:>16s}\n"

        # Print surface hopping variables
        if (self.md_type == "SH" or self.md_type == "SHXF"):
            dynamics_info += f"\n  Velocity Rescale in Hop  = {self.hop_rescale:>16s}\n"
            dynamics_info += f"  Rescale when Hop Reject  = {self.hop_reject:>16s}\n"

        # Print XF variables
        if (self.md_type == "SHXF" or self.md_type == "EhXF"):
            # Print density threshold used in decoherence term
            dynamics_info += f"\n  Density Threshold        = {self.threshold:>16.6f}"
            if (self.md_type == "SHXF" and self.l_xf1d):
                # Print reduced mass
                dynamics_info += f"\n  Reduced Mass             = {self.aux.mass[0]:16.6f}"
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
            "".join([f'E({ist})(H){"":8s}' for ist in range(self.mol.nst)])
        typewriter(tmp, unixmd_dir, "MDENERGY", "w")

        if (self.md_type != "BOMD"):
            # BO coefficents, densities file header
            if (self.obj == "density"):
                tmp = f'{"#":5s} Density Matrix: population Re; see the manual for detail orders'
                typewriter(tmp, unixmd_dir, "BOPOP", "w")
                tmp = f'{"#":5s} Density Matrix: coherence Re-Im; see the manual for detail orders'
                typewriter(tmp, unixmd_dir, "BOCOH", "w")
            elif (self.obj == "coefficient"):
                tmp = f'{"#":5s} BO State Coefficients: state Re-Im; see the manual for detail orders'
                typewriter(tmp, unixmd_dir, "BOCOEF", "w")
                if (self.l_print_dm):
                    tmp = f'{"#":5s} Density Matrix: population Re; see the manual for detail orders'
                    typewriter(tmp, unixmd_dir, "BOPOP", "w")
                    tmp = f'{"#":5s} Density Matrix: coherence Re-Im; see the manual for detail orders'
                    typewriter(tmp, unixmd_dir, "BOCOH", "w")

            # NACME file header
            tmp = f'{"#":5s}Non-Adiabatic Coupling Matrix Elements: off-diagonal'
            typewriter(tmp, unixmd_dir, "NACME", "w")

        # file header for SH-based methods
        if (self.md_type == "SH" or self.md_type == "SHXF"):
            tmp = f'{"#":5s}{"Step":8s}{"Running State":10s}'
            typewriter(tmp, unixmd_dir, "SHSTATE", "w")

            tmp = f'{"#":5s}{"Step":12s}' + "".join([f'Prob({ist}){"":8s}' for ist in range(self.mol.nst)])
            typewriter(tmp, unixmd_dir, "SHPROB", "w")

        # file header for XF-based methods
        if (self.md_type == "SHXF" or self.md_type == "EhXF"):
            tmp = f'{"#":5s} Time-derivative Density Matrix by decoherence: population; see the manual for detail orders'
            typewriter(tmp, unixmd_dir, "DOTPOPD", "w")

    def write_md_output(self, unixmd_dir, istep):
        """ Write output files

            :param string unixmd_dir: Directory where MD output files are written
            :param integer istep: Current MD step
        """
        # Write MOVIE.xyz file including positions and velocities
        tmp = f'{self.mol.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)' + \
            "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
            "".join([f'{self.mol.pos[iat, isp] * au_to_A:15.8f}' for isp in range(self.mol.ndim)]) + \
            "".join([f"{self.mol.vel[iat, isp]:15.8f}" for isp in range(self.mol.ndim)]) for iat in range(self.mol.nat)])
        typewriter(tmp, unixmd_dir, "MOVIE.xyz", "a")

        # Write MDENERGY file including several energy information
        tmp = f'{istep + 1:9d}{self.mol.ekin:15.8f}{self.mol.epot:15.8f}{self.mol.etot:15.8f}' \
            + "".join([f'{states.energy:15.8f}' for states in self.mol.states])
        typewriter(tmp, unixmd_dir, "MDENERGY", "a")

        if (self.md_type != "BOMD"):
            # Write BOCOEF, BOPOP, BOCOH files
            if (self.obj == "density"):
                tmp = f'{istep + 1:9d}' + "".join([f'{self.mol.rho.real[ist, ist]:15.8f}' for ist in range(self.mol.nst)])
                typewriter(tmp, unixmd_dir, "BOPOP", "a")
                tmp = f'{istep + 1:9d}' + "".join([f"{self.mol.rho.real[ist, jst]:15.8f}{self.mol.rho.imag[ist, jst]:15.8f}" \
                    for ist in range(self.mol.nst) for jst in range(ist + 1, self.mol.nst)])
                typewriter(tmp, unixmd_dir, "BOCOH", "a")
            elif (self.obj == "coefficient"):
                tmp = f'{istep + 1:9d}' + "".join([f'{states.coef.real:15.8f}{states.coef.imag:15.8f}' \
                    for states in self.mol.states])
                typewriter(tmp, unixmd_dir, "BOCOEF", "a")
                if (self.l_print_dm):
                    tmp = f'{istep + 1:9d}' + "".join([f'{self.mol.rho.real[ist, ist]:15.8f}' for ist in range(self.mol.nst)])
                    typewriter(tmp, unixmd_dir, "BOPOP", "a")
                    tmp = f'{istep + 1:9d}' + "".join([f"{self.mol.rho.real[ist, jst]:15.8f}{self.mol.rho.imag[ist, jst]:15.8f}" \
                        for ist in range(self.mol.nst) for jst in range(ist + 1, self.mol.nst)])
                    typewriter(tmp, unixmd_dir, "BOCOH", "a")

            # Write NACME file
            tmp = f'{istep + 1:10d}' + "".join([f'{self.mol.nacme[ist, jst]:15.8f}' \
                for ist in range(self.mol.nst) for jst in range(ist + 1, self.mol.nst)])
            typewriter(tmp, unixmd_dir, "NACME", "a")

            # Write NACV file
            if (not self.mol.l_nacme and self.verbosity >= 2):
                for ist in range(self.mol.nst):
                    for jst in range(ist + 1, self.mol.nst):
                        tmp = f'{self.mol.nat_qm:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}NACV' + \
                            "".join(["\n" + f'{self.mol.symbols[iat]:5s}' + \
                            "".join([f'{self.mol.nac[ist, jst, iat, isp]:15.8f}' for isp in range(self.mol.ndim)]) for iat in range(self.mol.nat_qm)])
                        typewriter(tmp, unixmd_dir, f"NACV_{ist}_{jst}", "a")

    def write_final_xyz(self, unixmd_dir, istep):
        """ Write final positions and velocities

            :param string unixmd_dir: Directory where MD output files are written
            :param integer istep: Current MD step
        """
        # Write FINAL.xyz file including positions and velocities
        tmp = f'{self.mol.nat:6d}\n{"":2s}Step:{istep + 1:6d}{"":12s}Position(A){"":34s}Velocity(au)'
        for iat in range(self.mol.nat):
            tmp += "\n" + f'{self.mol.symbols[iat]:5s}' + \
                "".join([f'{self.mol.pos[iat, isp] * au_to_A:15.8f}' for isp in range(self.mol.ndim)]) \
                + "".join([f"{self.mol.vel[iat, isp]:15.8f}" for isp in range(self.mol.ndim)])

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
                raise ValueError (f"( {self.md_type}.{call_name()} ) Inconsistent charge embedding options between QM and MM objects! {qm.embedding} and {mm.embedding}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Not compatible objects in QMMM! {qm.qm_prog}.{qm.qm_method} and {mm.mm_prog}")


