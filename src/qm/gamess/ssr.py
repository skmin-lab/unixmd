from __future__ import division
from qm.gamess.gamess import GAMESS
from misc import data, call_name, au_to_A, eps
import os, shutil, re, textwrap
import numpy as np

class SSR(GAMESS):
    """ Class for SSR method of GAMESS

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string functional: Exchange-correlation functional information
        :param integer active_space: Active space for SSR calculation
        :param string guess: Initial guess for REKS SCF iterations
        :param integer reks_rho_tol: Density convergence for REKS SCF iterations
        :param integer reks_max_iter: Maximum number of REKS SCF iterations
        :param boolean l_reks_diis: Use DIIS algorithm for REKS SCF iterations
        :param double shift: Level shifting value in REKS SCF iterations
        :param boolean l_state_interactions: Include state-interaction terms to SA-REKS
        :param double cpreks_grad_tol: Gradient tolerance for CP-REKS equations
        :param integer cpreks_max_iter: Maximum number of CP-REKS iterations
        :param string qm_path: Path for QM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of GAMESS, check $VERNO
    """
    def __init__(self, molecule, basis_set="6-31g*", memory="50", functional="bhhlyp", \
        active_space=2, guess="dft", reks_rho_tol=5, reks_max_iter=30, l_reks_diis=True, \
        shift=0.3, l_state_interactions=False, cpreks_grad_tol=1E-6, cpreks_max_iter=100, \
        qm_path="./", nthreads=1, version="00"):
        # Initialize GAMESS common variables
        super(SSR, self).__init__(basis_set, memory, qm_path, nthreads, version)

        # Initialize GAMESS SSR variables
        self.functional = functional

        self.active_space = active_space
        if not (self.active_space in [2]):
            error_message = "Invalid active space for SSR!"
            error_vars = f"active_space = {self.active_space}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (self.active_space == 2):
            if (molecule.nst > 2):
                error_message = "SSR(3,2) calculation is not implemented in current version of pyUNIxMD!"
                error_vars = f"Molecule.nstates = {molecule.nst}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
            elif (molecule.nst == 1):
                error_message = "GAMESS does not support single-state REKS calculation!"
                error_vars = f"Molecule.nstates = {molecule.nst}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set initial guess for REKS SCF iterations
        self.guess = guess.lower()
        if not (self.guess in ["dft", "read"]):
            error_message = "Invalid initial guess for SSR!"
            error_vars = f"guess = {self.guess}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.reks_rho_tol = reks_rho_tol
        self.reks_max_iter = reks_max_iter
        self.l_reks_diis = l_reks_diis
        self.shift = shift
        self.l_state_interactions = l_state_interactions

        self.cpreks_grad_tol = cpreks_grad_tol
        self.cpreks_max_iter = cpreks_max_iter

        # Check the closed shell for systems
        if (not int(molecule.nelec) % 2 == 0):
            error_message = "Only closed shell configuration supported, check charge!"
            error_vars = f"Molecule.nelec = {int(molecule.nelec)}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set 'l_nacme' with respect to the computational method
        # SSR can produce NACs, so we do not need to get NACME from CIoverlap
        # SSR can compute the gradient of several states simultaneously.
        molecule.l_nacme = False
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from SSR method

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer,list bo_list: List of BO states for BO calculation
            :param double dt: Time interval
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        self.copy_files(istep)
        super().get_data(base_dir, calc_force_only)
        self.get_input(molecule, istep, bo_list)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_QM(molecule, bo_list)
        self.move_dir(base_dir)

    def copy_files(self, istep):
        """ Copy necessary scratch files in previous step

            :param integer istep: Current MD step
        """
        # Copy required files to read initial guess
        if (self.guess == "read" and istep >= 0):
            # After T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "./userscr/gamess.dat"), \
                os.path.join(self.scr_qm_dir, "../gamess.dat"))

    def get_input(self, molecule, istep, bo_list):
        """ Generate GAMESS input files: gamess.inp.xx

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Geometry Block
        input_geom = ""
        for iat in range(molecule.nat_qm):
            if (iat > 0):
                input_geom += " " * 8
            input_geom += f"{molecule.symbols[iat]:4s}"
            sym_index = list(data.keys()).index(molecule.symbols[iat])
            input_geom += f"{sym_index:4.1f}"
            input_geom += "".join([f"{i:14.8f}" for i in molecule.pos[iat] * au_to_A])
            if (iat < molecule.nat_qm - 1):
                input_geom += "\n"

        input_data = textwrap.indent(textwrap.dedent(f"""\
        $data
        gamess
        C1
        {input_geom}
        $end
        """), " ")

        # Basis Block
        # TODO : Currently, only 6-31g* is possible
        if not (self.basis_set in ["6-31g*"]):
            error_message = "Invalid basis sets for GAMESS!"
            error_vars = f"basis_set = {self.basis_set}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        input_basis = textwrap.indent(textwrap.dedent(f"""\
        $basis
        gbasis=n31 ngauss=6 ndfunc=1
        $end
        """), " ")

        # System Block
        input_system = textwrap.indent(textwrap.dedent(f"""\
        $system
        mwords={self.memory}
        $end
        """), " ")

        # Guess options
        if (self.guess == "read"):
            if (istep == -1):
                self.do_dft = True
            elif (istep >= 0):
                # Move previous file to currect directory
                os.rename("../gamess.dat", os.path.join(self.scr_qm_dir, "gamess.dat"))
                self.do_dft = False
        elif (self.guess == "dft"):
            self.do_dft = True

        if (self.do_dft):
            # Control Block: Run DFT method to obtain initial guess of REKS calculation
            input_control_dft = textwrap.indent(textwrap.dedent(f"""\
            $contrl
            scftyp=rhf runtyp=energy dfttyp={self.functional} pltorb=.true.
            units=angs maxit={self.reks_max_iter} icharg={molecule.charge} mult=1
            $end
            """), " ")

        # SCF Block
        input_scf = textwrap.indent(textwrap.dedent(f"""\
        $scf
        nconv={self.reks_rho_tol}
        $end
        """), " ")

        if (self.do_dft):
            # Make 'gamess.inp.0' file
            input_gamess = ""

            input_gamess += input_control_dft
            input_gamess += input_system
            input_gamess += input_basis
            input_gamess += input_scf
            input_gamess += input_data

            # Write 'gamess.inp.0' file
            file_name = "gamess.inp.0"
            with open(file_name, "w") as f:
                f.write(input_gamess)

        # Control Block: calculate gradient option in REKS method
        input_control_reks = textwrap.indent(textwrap.dedent(f"""\
        $contrl
        scftyp=REKS runtyp=gradient dfttyp={self.functional} pltorb=.true.
        units=angs maxit={self.reks_max_iter} icharg={molecule.charge} mult=1
        $end
        """), " ")

        # Guess Block
        input_guess = textwrap.indent(textwrap.dedent(f"""\
        $guess
        guess=moread norb={int(molecule.nelec / 2 + 1)}
        $end
        """), " ")

        # Energy functional and NAC options
        if (self.l_state_interactions):
            # SSR state
            reks_type = 1
            if (self.calc_coupling):
                # SHXF, SH, Eh
                self.nac = "Yes"
            else:
                # BOMD
                self.nac = "No"
        else:
            # SA-REKS state
            reks_type = 0
            if (self.calc_coupling):
                # SHXF, SH, Eh
                error_message = "SA-REKS calculation does not provide NACVs!"
                error_vars = f"l_state_interactions = {self.l_state_interactions}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                # BOMD
                self.nac = "No"

        # Set DIIS option in REKS
        if (self.l_reks_diis):
            reks_diis = "yes"
        else:
            reks_diis = "no"

        # First calculation
        # REKS Block: calculate gradient for target SA-REKS or SSR state
        input_reks = textwrap.indent(textwrap.dedent(f"""\
        $reks
        rexType={reks_type} rexTarget={bo_list[0] + 1} rexShift={self.shift}
        rexDIIS={reks_diis} rxCGit={self.cpreks_max_iter} rxCGth={self.cpreks_grad_tol}
        $end
        """), " ")

        # Make 'gamess.inp.1' file
        input_gamess = ""

        input_gamess += input_control_reks
        input_gamess += input_system
        input_gamess += input_basis
        input_gamess += input_scf
        input_gamess += input_reks
        input_gamess += input_guess
        input_gamess += input_data

        # Write 'gamess.inp.1' file
        file_name = "gamess.inp.1"
        with open(file_name, "w") as f:
            f.write(input_gamess)

        # For evaluation of NACVs, two additional calculations are required
        if (self.nac == "Yes"):
            # Second calculation
            # REKS Block: calculate gradient for target SA-REKS state
            input_reks = textwrap.indent(textwrap.dedent(f"""\
            $reks
            rexType=0 rexTarget={bo_list[0] + 1} rexShift={self.shift}
            rexDIIS={reks_diis} rxCGit={self.cpreks_max_iter} rxCGth={self.cpreks_grad_tol}
            $end
            """), " ")

            # Make 'gamess.inp.2' file
            input_gamess = ""

            input_gamess += input_control_reks
            input_gamess += input_system
            input_gamess += input_basis
            input_gamess += input_scf
            input_gamess += input_reks
            input_gamess += input_guess
            input_gamess += input_data

            # Write 'gamess.inp.2' file
            file_name = "gamess.inp.2"
            with open(file_name, "w") as f:
                f.write(input_gamess)

            # Third calculation
            # REKS Block: calculate gradient for averaged functional
            input_reks = textwrap.indent(textwrap.dedent(f"""\
            $reks
            rexType=0 rexTarget=0 rexShift={self.shift}
            rexDIIS={reks_diis} rxCGit={self.cpreks_max_iter} rxCGth={self.cpreks_grad_tol}
            $end
            """), " ")

            # Make 'gamess.inp.3' file
            input_gamess = ""

            input_gamess += input_control_reks
            input_gamess += input_system
            input_gamess += input_basis
            input_gamess += input_scf
            input_gamess += input_reks
            input_gamess += input_guess
            input_gamess += input_data

            # Write 'gamess.inp.3' file
            file_name = "gamess.inp.3"
            with open(file_name, "w") as f:
                f.write(input_gamess)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run SSR calculation and save the output files to qm_log directory

            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Environmental variable setting
        user_scr_dir = os.path.join(self.scr_qm_dir, "userscr")
        os.makedirs(user_scr_dir)
        os.environ["USERSCR"] = user_scr_dir
        tmp_dir = os.path.join(self.scr_qm_dir, "tmp")
        os.makedirs(tmp_dir)
        os.environ["TMPDIR"] = tmp_dir
        os.environ["GMSPATH"] = self.qm_path

        # Set run command
        qm_command = os.path.join(self.qm_path, "rungms")
        command = f"{qm_command} gamess {self.version} {self.nthreads} >& gamess.log"
        # OpenMP setting
        os.environ["OMP_NUM_THREADS"] = "1"

        # Run GAMESS method
        # Preparation: Run DFT to obtain initial guess
        if (self.do_dft):
            shutil.copy("gamess.inp.0", "gamess.inp")
            os.system(command)
            shutil.copy("gamess.log", "gamess.log.0")

            # Save the molecular orbitals
            guess_file = os.path.join(user_scr_dir, "gamess.dat")
            shutil.copy(guess_file, "gamess.dat")

            if (os.path.exists(user_scr_dir)):
                shutil.rmtree(user_scr_dir)
            os.makedirs(user_scr_dir)
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)
            os.makedirs(tmp_dir)

        if (istep == -1):
            # Read the number of basis functions
            # Read 'gamess.log.0' file
            file_name = "gamess.log.0"
            with open(file_name, "r") as f:
                log_out = f.read()

            nbasis = re.findall('BASIS FUNCTIONS =\s+(\S+)', log_out)
            nbasis = np.array(nbasis[0], dtype=np.int32)

            # Read 'gamess.dat' file
            if (nbasis % 5 == 0):
                ncol = nbasis // 5
            else:
                ncol = nbasis // 5 + 1
            # Calculate total number of lines for the molecular orbitals
            self.total_nrow = nbasis * ncol + 1

        # Add the initial guess to input file
        mo_command = f"grep 'VEC' -A {self.total_nrow} gamess.dat > gamess.vec"
        os.system(mo_command)

        cat_command = f"cat gamess.inp.1 gamess.vec > tmp.inp"
        os.system(cat_command)
        os.rename("tmp.inp", "gamess.inp.1")

        # First calculation: Run REKS using initial guess obtained from DFT or previous result
        shutil.copy("gamess.inp.1", "gamess.inp")
        os.system(command)
        shutil.copy("gamess.log", "gamess.log.1")

        if (self.nac == "Yes"):
            # Save the molecular orbitals
            guess_file = os.path.join(user_scr_dir, "gamess.dat")
            shutil.copy(guess_file, "gamess.dat")

            # Add the initial guess to input file
            mo_command = f"grep 'VEC' -A {self.total_nrow} gamess.dat > gamess.vec"
            os.system(mo_command)

            cat_command = f"cat gamess.inp.2 gamess.vec > tmp.inp"
            os.system(cat_command)
            os.rename("tmp.inp", "gamess.inp.2")
            cat_command = f"cat gamess.inp.3 gamess.vec > tmp.inp"
            os.system(cat_command)
            os.rename("tmp.inp", "gamess.inp.3")

            if (os.path.exists(user_scr_dir)):
                shutil.rmtree(user_scr_dir)
            os.makedirs(user_scr_dir)
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)
            os.makedirs(tmp_dir)
            # Second calculation: Run REKS using initial guess obtained from REKS
            shutil.copy("gamess.inp.2", "gamess.inp")
            os.system(command)
            shutil.copy("gamess.log", "gamess.log.2")

            if (os.path.exists(user_scr_dir)):
                shutil.rmtree(user_scr_dir)
            os.makedirs(user_scr_dir)
            if (os.path.exists(tmp_dir)):
                shutil.rmtree(tmp_dir)
            os.makedirs(tmp_dir)
            # Third calculation: Run REKS using initial guess obtained from REKS
            shutil.copy("gamess.inp.3", "gamess.inp")
            os.system(command)
            shutil.copy("gamess.log", "gamess.log.3")

        # Copy the output file to 'qm_log' directory
        tmp_dir = os.path.join(base_dir, "qm_log")
        if (os.path.exists(tmp_dir)):
            if (self.do_dft):
                command = f"cat gamess.log.0 gamess.log.1 gamess.log.2 gamess.log.3 > gamess.log"
            else:
                command = f"cat gamess.log.1 gamess.log.2 gamess.log.3 > gamess.log"
            os.system(command)
            log_step = f"gamess.log.{istep + 1}.{bo_list[0]}"
            shutil.copy("gamess.log", os.path.join(tmp_dir, log_step))

    def extract_QM(self, molecule, bo_list):
        """ Read the output files to get BO information

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Read 'gamess.log.1' file
        file_name = "gamess.log.1"
        with open(file_name, "r") as f:
            log_out = f.read()

        # Energy
        if (self.l_state_interactions):
            # SSR state
            energy = re.findall('SSR state\s\d\s+([-]\S+)', log_out)
            energy = np.array(energy, dtype=np.float64)
            for ist in range(molecule.nst):
                molecule.states[ist].energy = energy[ist]
        else:
            # SA-REKS state
            tmp_e = '\(PPS\):\S+\s+([-]\S+)'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy, dtype=np.float64)
            molecule.states[0].energy = energy[0]
            tmp_e = '\(OSS\):\S+\s+([-]\S+)'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy, dtype=np.float64)
            molecule.states[1].energy = energy[0]

        # Force
        tmp_g = 'GRADIENT OF THE ENERGY\n\s+' + '-' * 22 + '\n\n UNITS ARE HARTREE/BOHR' + \
            "\s+E'X\s+E'Y\s+E'Z \n" + "\s+\d+ \S+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n" * molecule.nat_qm
        grad = re.findall(tmp_g, log_out)
        grad = np.array(grad[0], dtype=np.float64)
        grad = grad.reshape(molecule.nat_qm, 3, order='C')

        if (self.nac == "No"):
            # SSR or SA-REKS force
            molecule.states[bo_list[0]].force = - np.copy(grad)
        else:
            # Read coefficients
            ssr_coef = np.zeros((molecule.nst, molecule.nst))
            coef = re.findall('SSR state\s\d\s+[-]\S+\s+([-]*\S+)\s+([-]*\S+)', log_out)
            coef = np.array(coef, dtype=np.float64)
            for ist in range(molecule.nst):
                ssr_coef[ist] = coef[ist]

            # Save SSR gradient
            ssr_grad = np.zeros((molecule.nst, molecule.nat, molecule.ndim))
            ssr_grad[bo_list[0]] = np.copy(grad)

            # Read 'gamess.log.2' file
            file_name = "gamess.log.2"
            with open(file_name, "r") as f:
                log_out = f.read()

            # Read SA-REKS gradient
            sa_grad = np.zeros((molecule.nst, molecule.nat, molecule.ndim))

            grad = re.findall(tmp_g, log_out)
            grad = np.array(grad[0], dtype=np.float64)
            grad = grad.reshape(molecule.nat_qm, 3, order='C')
            sa_grad[bo_list[0]] = np.copy(grad)

            # Read 'gamess.log.3' file
            file_name = "gamess.log.3"
            with open(file_name, "r") as f:
                log_out = f.read()

            # Read averaged functional gradient
            avg_grad = np.zeros((molecule.nat, molecule.ndim))

            grad = re.findall(tmp_g, log_out)
            grad = np.array(grad[0], dtype=np.float64)
            grad = grad.reshape(molecule.nat_qm, 3, order='C')
            avg_grad = np.copy(grad)

            ssr_grad[1 - bo_list[0]] = 2. * avg_grad - ssr_grad[bo_list[0]]
            sa_grad[1 - bo_list[0]] = 2. * avg_grad - sa_grad[bo_list[0]]

            for ist in range(molecule.nst):
                molecule.states[ist].force = - ssr_grad[ist]

            # Calculate G and g vectors from SSR and SA-REKS gradients
            G_vec = np.zeros((molecule.nat, molecule.ndim))
            g_vec = np.zeros((molecule.nat, molecule.ndim))

            G_vec = 0.5 * (ssr_grad[0] - ssr_grad[1])
            g_vec = 0.5 * (sa_grad[0] - sa_grad[1])

            # Calculate NACVs from G, g and coefficients
            cos_theta = ssr_coef[0, 0] ** 2. - ssr_coef[0, 1] ** 2.
            sin_theta = 2. * ssr_coef[0, 0] * ssr_coef[0, 1]
            if (abs(sin_theta) <= 0.0001):
                if (sin_theta <= eps):
                    sin_theta = - 0.0001
                else:
                    sin_theta = 0.0001

            H_vec = np.zeros((molecule.nat, molecule.ndim))
            H_vec = 0.5 * (cos_theta / sin_theta * (G_vec - cos_theta * g_vec) + sin_theta * g_vec) \
                / (molecule.states[1].energy - molecule.states[0].energy)

            molecule.nac[0, 1] = np.copy(H_vec)
            molecule.nac[1, 0] = - np.copy(H_vec)


