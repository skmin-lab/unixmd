from __future__ import division
from bo.terachem.terachem import TeraChem
from misc import call_name
import os, shutil, re, textwrap
import numpy as np

class SSR(TeraChem):
    """ Class for SSR method of TeraChem program

        :param object molecule: molecule object
        :param string basis_set: basis set information
        :param string functional: level of DFT theory
        :param string precision: precision in the calculations
        :param double scf_rho_tol: wavefunction convergence for SCF iterations
        :param integer scf_max_iter: maximum number of SCF iterations
        :param boolean ssr22: use REKS(2,2) calculation?
        :param string guess: initial guess for REKS SCF iterations
        :param string guess_file: initial guess file
        :param double reks_rho_tol: DIIS error for REKS SCF iterations
        :param integer reks_max_iter: maximum number of REKS SCF iterations
        :param double shift: level shifting value in REKS SCF iterations
        :param boolean use_ssr_state: calculate SSR state, if not, treat SA-REKS
        :param double cpreks_grad_tol: gradient tolerance for CP-REKS equations
        :param integer cpreks_max_iter: maximum number of CP-REKS iterations
        :param string qm_path: path for QM binary
        :param integer ngpus: number of GPUs
        :param string gpu_id: ID of used GPUs
        :param double version: version of TeraChem program
    """
    def __init__(self, molecule, ngpus=1, gpu_id="1", precision="dynamic", \
        version=1.92, functional="hf", basis_set="sto-3g", scf_rho_tol=1E-2, \
        scf_max_iter=300, ssr22=True, guess="dft", guess_file="./c0", \
        reks_rho_tol=1E-6, reks_max_iter=1000, shift=0.3, use_ssr_state=True, \
        cpreks_grad_tol=1E-6, cpreks_max_iter=1000, qm_path="./"):
        # Initialize TeraChem common variables
        super(SSR, self).__init__(functional, basis_set, qm_path, ngpus, \
            gpu_id, precision, version)

        # Initialize TeraChem SSR variables
        self.scf_rho_tol = scf_rho_tol
        self.scf_max_iter = scf_max_iter

        self.ssr22 = ssr22
        if (self.ssr22):
            if (molecule.nst > 2):
                raise ValueError (f"( {self.qm_method}.{call_name()} ) 3-state REKS not implemented! {molecule.nst}")
            self.reks_rho_tol = reks_rho_tol
            self.reks_max_iter = reks_max_iter
            self.shift = shift
            self.use_ssr_state = use_ssr_state

            # Set initial guess for REKS SCF iterations
            self.guess = guess
            self.guess_file = guess_file
            if (not (self.guess == "dft" or self.guess == "read")):
                raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong input for initial guess option! {self.guess}")

            if (molecule.nst > 1):
                self.cpreks_grad_tol = cpreks_grad_tol
                self.cpreks_max_iter = cpreks_max_iter
        else:
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Other active space not implemented! {self.ssr22}")

        # Set 'l_nacme' with respect to the computational method
        # SSR can produce NACs, so we do not need to get NACME from CIoverlap
        # SSR can compute the gradient of several states simultaneously.
        molecule.l_nacme = False
        self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from SSR method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        self.copy_files(istep)
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list)
        self.move_dir(base_dir)

    def copy_files(self, istep):
        """ Copy necessary scratch files in previous step

            :param integer istep: current MD step
        """
        # Copy required files to read initial guess
        if (self.guess == "read" and istep >= 0):
            # After T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "./scr/c0"), \
                os.path.join(self.scr_qm_dir, "../c0"))

    def get_input(self, molecule, bo_list):
        """ Generate TeraChem input files: input.tcin

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Make 'input.tcin' file
        input_terachem = ""

        # Guess Block
        if (self.guess == "read"):
            c0_dir = os.path.join(self.scr_qm_dir, "scr")
            os.makedirs(c0_dir)
            if (istep == -1):
                if (os.path.isfile(self.guess_file)):
                    # Copy guess file to currect directory
                    shutil.copy(self.guess_file, os.path.join(c0_dir, "c0"))
                    restart = 1
                    dft = False
                else:
                    restart = 0
                    dft = True
            elif (istep >= 0):
                # Move previous file to currect directory
                os.rename("../c0", os.path.join(c0_dir, "c0"))
                restart = 1
                dft = False
        elif (self.guess == "dft"):
            restart = 0
            dft = True

        # Control Block
        input_control = \
        input_control = textwrap.dedent(f"""\

        run gradient

        coordinates geometry.xyz

        """)
        input_terachem += input_control

        # System Block
        input_system = textwrap.dedent(f"""\
        precision {self.precision}
        gpus {self.ngpus}   {self.gpu_id}

        """)
        input_terachem += input_system

        # Setting Block
        input_setting = textwrap.dedent(f"""\
        method {self.functional}
        basis {self.basis_set}
        charge {molecule.charge}

        """)
        input_terachem += input_setting

        # DFT Block
        if (dft):
            input_dft = textwrap.dedent(f"""\
            convthre {self.scf_rho_tol}
            maxit {self.scf_max_iter}

            """)
            input_terachem += input_dft

        # REKS Block
        if (self.ssr22):

            # Energy functional options
            if (molecule.nst == 1):
                sa_reks = 0
            elif (molecule.nst == 2):
                if (self.use_ssr_state):
                    sa_reks = 2
                else:
                    sa_reks = 1

            # NAC calculation options
            if (self.calc_coupling and sa_reks == 2):
                # SHXF, SH, Eh : SSR state
                reks_target = 12
                self.nac = "Yes"
            else:
                # BOMD : SSR state, SA-REKS state or single-state REKS
                reks_target = bo_list[0] + 1
                self.nac = "No"

            # TODO: pointcharges? in qmmm?

            input_reks_basic = textwrap.dedent(f"""\
            reks22 yes
            reks_convthre {self.reks_rho_tol}
            reks_maxit {self.reks_max_iter}
            reks_diis yes
            reks_shift {self.shift}
            sa_reks {sa_reks}
            reks_target {reks_target}
            """)
            input_terachem += input_reks_basic

            if (restart == 1):
                input_reks_guess = textwrap.dedent(f"""\
                reks_guess {restart}
                guess ./scr/c0
                """)
                input_terachem += input_reks_guess

            if (molecule.nst > 1):
                input_cpreks = textwrap.dedent(f"""\
                cpreks_thresh {self.cpreks_grad_tol}
                cpreks_maxit {self.cpreks_max_iter}
                """)
                input_terachem += input_cpreks

        # Write 'input.tcin' file
        file_name = "input.tcin"
        with open(file_name, "w") as f:
            f.write(input_terachem)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run SSR calculation and save the output files to QMlog directory

            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Run TeraChem method
        qm_command = os.path.join(self.qm_path, "terachem")
        # OpenMP setting
        os.environ["OMP_NUM_THREADS"] = "1"
        command = f"{qm_command} input.tcin > log"
        os.system(command)
        # Copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Read 'log' file
        file_name = "log"
        with open(file_name, "r") as f:
            log_out = f.read()

        # Energy
        for states in molecule.states:
            states.energy = 0.

        if (molecule.nst == 1):
            # Single-state REKS
            tmp_e = 'REKS energy:\s+([-]\S+)'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy)
            energy = energy.astype(float)
            molecule.states[0].energy = energy[0]
        else:
            if (self.use_ssr_state):
                # SSR state
                energy = re.findall('SSR state\s\d\s+([-]\S+)', log_out)
                energy = np.array(energy)
                energy = energy.astype(float)
                for ist in range(molecule.nst):
                    molecule.states[ist].energy = energy[ist]
            else:
                # SA-REKS state
                tmp_e = '\(GSS\):\s+([-]\S+)'
                energy = re.findall(tmp_e, log_out)
                energy = np.array(energy)
                energy = energy.astype(float)
                molecule.states[0].energy = energy[0]
                tmp_e = '\(OSS\):\s+([-]\S+)'
                energy = re.findall(tmp_e, log_out)
                energy = np.array(energy)
                energy = energy.astype(float)
                molecule.states[1].energy = energy[0]

        # Force
        for states in molecule.states:
            states.force = np.zeros((molecule.nat, molecule.nsp))

        if (self.nac == "Yes"):
            # SHXF, SH, Eh : SSR state
            for ist in range(molecule.nst):
                tmp_f = f'Eigen state {ist + 1} gradient\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
                    '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                force = re.findall(tmp_f, log_out)
                force = np.array(force[0])
                force = force.astype(float)
                force = force.reshape(molecule.nat, 3, order='C')
                molecule.states[ist].force = - np.copy(force)
        else:
            # BOMD : SSR state, SA-REKS state or single-state REKS
            tmp_f = 'Gradient units are Hartree/Bohr\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	              '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = - np.copy(force)

        # NAC
        if (self.nac == "Yes"):

            # 1.92 version do not show H vector
            if (self.version == 1.92):
                # Zeroing for G, h and H vectors
                Gvec = np.zeros((molecule.nat, molecule.nsp))
                hvec = np.zeros((molecule.nat, molecule.nsp))
                ssr_coef = np.zeros((molecule.nst, molecule.nst))
                Hvec = np.zeros((molecule.nat, molecule.nsp))
                # Calculate G vector, G vector is difference gradient so minus sign is needed
                Gvec = - 0.5 * (molecule.states[0].force - molecule.states[1].force)
                # Read h vector
                tmp_c = 'Coupling gradient\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	                  '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                hvec = re.findall(tmp_c, log_out)
                hvec = np.array(hvec[0])
                hvec = hvec.astype(float)
                hvec = hvec.reshape(molecule.nat, 3, order='C')
                # Read coefficients of SSR state
                ssr_coef = re.findall('SSR state\s\d\s+[-]\S+\s+([-]*\S+)\s+([-]*\S+)', log_out)
                ssr_coef = np.array(ssr_coef)
                ssr_coef = ssr_coef.astype(float)
                # Calculate H vector from G, h vector
                Hvec = 1. / (molecule.states[1].energy - molecule.states[0].energy) * \
                    ((2. * ssr_coef[0, 0] * ssr_coef[1, 0] * Gvec + hvec) / \
                    (ssr_coef[0, 0] * ssr_coef[1, 1] + ssr_coef[1, 0] * ssr_coef[0, 1]))

            kst = 0
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.nac[ist, jst] = np.zeros((molecule.nat, molecule.nsp))
                    elif (ist < jst):
                        if (self.version == 1.92):
                            molecule.nac[ist, jst] = Hvec
                        elif (self.version == 1.93):
                            tmp_c = '>\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	                              '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                            nac = re.findall(tmp_c, log_out)
                            nac = np.array(nac[kst])
                            nac = nac.astype(float)
                            nac = nac.reshape(molecule.nat, 3, order='C')
                            molecule.nac[ist, jst] = np.copy(nac)
                        kst += 1
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]



