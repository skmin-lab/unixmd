from __future__ import division
from qm.terachem.terachem import TeraChem
from misc import call_name
import os, shutil, re, textwrap
import numpy as np

class SSR(TeraChem):
    """ Class for SSR method of TeraChem

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param string functional: Exchange-correlation functional information
        :param string precision: Precision in the calculations
        :param double scf_wf_tol: Wavefunction convergence for SCF iterations
        :param integer scf_max_iter: Maximum number of SCF iterations
        :param integer active_space: Active space for SSR calculation
        :param string guess: Initial guess for REKS SCF iterations
        :param string guess_file: Initial guess file
        :param double reks_diis_tol: DIIS error for REKS SCF iterations
        :param integer reks_max_iter: Maximum number of REKS SCF iterations
        :param double shift: Level shifting value in REKS SCF iterations
        :param boolean l_state_interactions: Include state-interaction terms to SA-REKS
        :param double cpreks_grad_tol: Gradient tolerance for CP-REKS equations
        :param integer cpreks_max_iter: Maximum number of CP-REKS iterations
        :param string root_path: Path for TeraChem root directory
        :param integer ngpus: Number of GPUs
        :param integer,list gpu_id: ID of used GPUs
        :param string version: Version of TeraChem
    """
    def __init__(self, molecule, ngpus=1, gpu_id=None, precision="dynamic", \
        version="1.93", functional="hf", basis_set="sto-3g", scf_wf_tol=1E-2, \
        scf_max_iter=300, active_space=2, guess="dft", guess_file="./c0", \
        reks_diis_tol=1E-6, reks_max_iter=1000, shift=0.3, l_state_interactions=False, \
        cpreks_grad_tol=1E-6, cpreks_max_iter=1000, root_path="./"):
        # Initialize TeraChem common variables
        super(SSR, self).__init__(functional, basis_set, root_path, ngpus, \
            gpu_id, precision, version)

        # Initialize TeraChem SSR variables
        self.scf_wf_tol = scf_wf_tol
        self.scf_max_iter = scf_max_iter

        self.active_space = active_space
        if not (self.active_space in [2]):
            error_message = "Invalid active space for SSR!"
            error_vars = f"active_space = {self.active_space}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (self.active_space == 2):
            if (molecule.nst > 2):
                error_message = "SSR(2,2) can calculate up to 2 electronic states!"
                error_vars = f"Molecule.nstates = {molecule.nst}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set initial guess for REKS SCF iterations
        self.guess = guess.lower()
        self.guess_file = guess_file
        if not (self.guess in ["dft", "read"]):
            error_message = "Invalid initial guess for SSR!"
            error_vars = f"guess = {self.guess}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.reks_diis_tol = reks_diis_tol
        self.reks_max_iter = reks_max_iter
        self.shift = shift
        self.l_state_interactions = l_state_interactions

        self.cpreks_grad_tol = cpreks_grad_tol
        self.cpreks_max_iter = cpreks_max_iter

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
        self.write_xyz(molecule)
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
            shutil.copy(os.path.join(self.scr_qm_dir, "./scr/c0"), \
                os.path.join(self.scr_qm_dir, "../c0"))

    def get_input(self, molecule, istep, bo_list):
        """ Generate TeraChem input files: input.tcin

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
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
        gpu_ids = ""
        for gpu_cur_id in self.gpu_id:
            gpu_ids += f" {gpu_cur_id} "

        input_system = textwrap.dedent(f"""\
        precision {self.precision}
        gpus {self.ngpus}   {gpu_ids}

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
            convthre {self.scf_wf_tol}
            maxit {self.scf_max_iter}

            """)
            input_terachem += input_dft

        # REKS Block
        if (self.active_space == 2):

            # Energy functional options
            space = "reks22 yes"
            if (molecule.nst == 1):
                sa_reks = 0
            elif (molecule.nst == 2):
                if (self.l_state_interactions):
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
        {space}
        reks_convthre {self.reks_diis_tol}
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
        """ Run SSR calculation and save the output files to qm_log directory

            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Run TeraChem method
        qm_command = os.path.join(self.qm_path, "terachem")
        # OpenMP setting
        os.environ["OMP_NUM_THREADS"] = "1"
        command = f"{qm_command} input.tcin > log"
        os.system(command)
        # Copy the output file to 'qm_log' directory
        tmp_dir = os.path.join(base_dir, "qm_log")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_QM(self, molecule, bo_list):
        """ Read the output files to get BO information

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Read 'log' file
        file_name = "log"
        with open(file_name, "r") as f:
            log_out = f.read()

        # Energy
        if (molecule.nst == 1):
            # Single-state REKS
            tmp_e = 'REKS energy:\s+([-]\S+)'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy, dtype=np.float64)
            molecule.states[0].energy = energy[0]
        else:
            if (self.l_state_interactions):
                # SSR state
                energy = re.findall('SSR state\s\d\s+([-]\S+)', log_out)
                energy = np.array(energy, dtype=np.float64)
                for ist in range(molecule.nst):
                    molecule.states[ist].energy = energy[ist]
            else:
                # SA-REKS state
                tmp_e = '\(GSS\):\s+([-]\S+)'
                energy = re.findall(tmp_e, log_out)
                energy = np.array(energy, dtype=np.float64)
                molecule.states[0].energy = energy[0]
                tmp_e = '\(OSS\):\s+([-]\S+)'
                energy = re.findall(tmp_e, log_out)
                energy = np.array(energy, dtype=np.float64)
                molecule.states[1].energy = energy[0]

        # Force
        if (self.nac == "Yes"):
            # SHXF, SH, Eh : SSR state
            for ist in range(molecule.nst):
                tmp_f = f'Eigen state {ist + 1} gradient\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
                    '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
                force = re.findall(tmp_f, log_out)
                force = np.array(force[0], dtype=np.float64)
                force = force.reshape(molecule.nat_qm, 3, order='C')
                molecule.states[ist].force = - np.copy(force)
        else:
            # BOMD : SSR state, SA-REKS state or single-state REKS
            tmp_f = 'Gradient units are Hartree/Bohr\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	              '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0], dtype=np.float64)
            force = force.reshape(molecule.nat_qm, 3, order='C')
            molecule.states[bo_list[0]].force = - np.copy(force)

        # NAC
        if (self.nac == "Yes"):

            ## 1.99 version do not show H vector
            #if (self.version == "1.99"):
            #    # Zeroing for G, h and H vectors

            #    Gvec = np.zeros((molecule.nat_qm, molecule.ndim))
            #    hvec = np.zeros((molecule.nat_qm, molecule.ndim))
            #    ssr_coef = np.zeros((molecule.nst, molecule.nst))
            #    Hvec = np.zeros((molecule.nat_qm, molecule.ndim))

            #    # Calculate G vector, G vector is difference gradient so minus sign is needed
            #    Gvec = - 0.5 * (molecule.states[0].force - molecule.states[1].force)
            #    # Read h vector
            #    tmp_c = 'Coupling gradient\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	          #        '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
            #    hvec = re.findall(tmp_c, log_out)
            #    hvec = np.array(hvec[0], dtype=np.float64)
            #    hvec = hvec.reshape(molecule.nat_qm, 3, order='C')
            #    # Read coefficients of SSR state
            #    ssr_coef = re.findall('SSR state\s\d\s+[-]\S+\s+([-]*\S+)\s+([-]*\S+)', log_out)
            #    ssr_coef = np.array(ssr_coef, dtype=np.float64)
            #    # Calculate H vector from G, h vector
            #    Hvec = 1. / (molecule.states[1].energy - molecule.states[0].energy) * \
            #        ((2. * ssr_coef[0, 0] * ssr_coef[1, 0] * Gvec + hvec) / \
            #        (ssr_coef[0, 0] * ssr_coef[1, 1] + ssr_coef[1, 0] * ssr_coef[0, 1]))

            #    molecule.nac[0, 1] = Hvec
            #    molecule.nac[1, 0] = - Hvec

            #elif (self.version == "1.93"):
            kst = 0
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    tmp_c = '>\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
                        '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
                    nac = re.findall(tmp_c, log_out)
                    nac = np.array(nac[kst], dtype=np.float64)
                    nac = nac.reshape(molecule.nat_qm, 3, order='C')
                    molecule.nac[ist, jst] = nac
                    molecule.nac[jst, ist] = - nac
                    kst += 1


