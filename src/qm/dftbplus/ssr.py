from __future__ import division
from qm.dftbplus.dftbplus import DFTBplus
from qm.dftbplus.dftbpar import spin_w, spin_w_lc, max_l
from misc import au_to_A, call_name
import os, shutil, re, textwrap
import numpy as np

class SSR(DFTBplus):
    """ Class for SSR method of DFTB+

        :param object molecule: Molecule object
        :param boolean scc: Include self-consistent charge (SCC) scheme
        :param double scc_tol: Stopping criteria for the SCC iterations
        :param integer scc_max_iter: Maximum number of SCC iterations
        :param boolean ocdftb: Include onsite correction to SCC term
        :param boolean lcdftb: Include long-range corrected functional
        :param string lc_method: Algorithms for LC-DFTB
        :param boolean ssr22: Use SSR(2,2) calculation?
        :param boolean ssr44: Use SSR(4,4) calculation?
        :param string guess: Initial guess method for SCC scheme
        :param string guess_file: Initial guess file for eigenvetors
        :param boolean state_interactions: Include state-interaction terms to SA-REKS
        :param double shift: Level shifting value in SCC iterations
        :param double,list tuning: Scaling factor for atomic spin constants
        :param string cpreks_grad_alg: Algorithms used in CP-REKS equations
        :param double cpreks_grad_tol: Tolerance used in the conjugate-gradient based algorithms
        :param boolean save_memory: Save memory in cache used in CP-REKS equations
        :param string embedding: Charge-charge embedding options
        :param boolean periodic: Use periodicity in the calculations
        :param double,list cell_length: The lattice vectors of periodic unit cell
        :param string sk_path: Path for slater-koster files
        :param string install_path: Path for DFTB+ install directory
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of DFTB+
    """
    def __init__(self, molecule, scc=True, scc_tol=1E-6, scc_max_iter=1000, ocdftb=False, \
        lcdftb=False, lc_method="MatrixBased", ssr22=False, ssr44=False, guess="h0", \
        guess_file="./eigenvec.bin", state_interactions=False, shift=0.3, tuning=None, \
        cpreks_grad_alg="pcg", cpreks_grad_tol=1E-8, save_memory=False, embedding=None, \
        periodic=False, cell_length=[0., 0., 0., 0., 0., 0., 0., 0., 0.], sk_path="./", \
        install_path="./", nthreads=1, version="20.1"):
        # Initialize DFTB+ common variables
        super(SSR, self).__init__(molecule, sk_path, install_path, nthreads, version)

        # Initialize DFTB+ SSR variables
        self.scc = scc
        self.scc_tol = scc_tol
        self.scc_max_iter = scc_max_iter

        self.ocdftb = ocdftb
        if (self.ocdftb):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Onsite-correction not implemented! {self.ocdftb}")

        self.lcdftb = lcdftb
        self.lc_method = lc_method

        self.ssr22 = ssr22
        self.ssr44 = ssr44
        if (not (self.ssr22 or self.ssr44)):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Other active space not implemented! {self.ssr22 or self.ssr44}")
        if (self.ssr44):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) REKS(4,4) not implemented! {self.ssr44}")

        # Set initial guess for eigenvectors
        self.guess = guess
        self.guess_file = guess_file
        if (not (self.guess == "h0" or self.guess == "read")):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong input for initial guess option! {self.guess}")

        self.state_interactions = state_interactions
        self.shift = shift

        # Set scaling factor for atomic spin constants
        self.tuning = tuning
        if (self.tuning != None):
            if (len(self.tuning) != len(self.atom_type)):
                raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong number of elements for scaling factors! {self.tuning}")

        self.cpreks_grad_alg = cpreks_grad_alg
        self.cpreks_grad_tol = cpreks_grad_tol
        self.save_memory = save_memory

        self.embedding = embedding
        if (self.embedding != None):
            if (not (self.embedding == "mechanical" or self.embedding == "electrostatic")):
                raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong charge embedding given! {self.embedding}")

        self.periodic = periodic
        self.a_axis = np.copy(cell_length[0:3])
        self.b_axis = np.copy(cell_length[3:6])
        self.c_axis = np.copy(cell_length[6:9])

        # Set 'l_nacme' and 're_calc' with respect to the computational method
        # DFTB/SSR can produce NACs, so we do not need to get NACME from CIoverlap
        # DFTB/SSR can compute the gradient of several states simultaneously.
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
            shutil.copy(os.path.join(self.scr_qm_dir, "eigenvec.bin"), \
                os.path.join(self.scr_qm_dir, "../eigenvec.bin.pre"))

    def get_input(self, molecule, istep, bo_list):
        """ Generate DFTB+ input files: geometry.gen, dftb_in.hsd

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Make 'geometry.gen' file
        os.system("xyz2gen geometry.xyz")
        if (self.periodic):
            # Substitute C to S in first line
            file_be = open('geometry.gen', 'r')
            file_af = open('tmp.gen', 'w')
            first_row = True
            for row in file_be:
                if (first_row):
                    row = f'{molecule.nat_qm} S\n'
                    first_row = False
                file_af.write(row)
            # Add gamma-point and cell lattice information
            geom_periodic = textwrap.dedent(f"""\
            {0.0:15.8f} {0.0:15.8f} {0.0:15.8f}
            {self.a_axis[0]:15.8f} {self.a_axis[1]:15.8f} {self.a_axis[2]:15.8f}
            {self.b_axis[0]:15.8f} {self.b_axis[1]:15.8f} {self.b_axis[2]:15.8f}
            {self.c_axis[0]:15.8f} {self.c_axis[1]:15.8f} {self.c_axis[2]:15.8f}
            """)
            file_af.write(geom_periodic)
            file_be.close()
            file_af.close()
            os.rename('tmp.gen', 'geometry.gen')

        # Make 'point_charges.xyz' file used in electrostatic charge embedding of QM/MM
        if (self.embedding == "electrostatic"):
            # Make 'point_charges.xyz' file
            input_geom_pc = ""
            for iat in range(molecule.nat_qm, molecule.nat):
                input_geom_pc += "".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A])
                input_geom_pc += f"  {molecule.mm_charge[iat - molecule.nat_qm]:8.4f}\n"

            # Write 'point_charges.xyz' file
            file_name = "point_charges.xyz"
            with open(file_name, "w") as f:
                f.write(input_geom_pc)

        # Make 'dftb_in.hsd' file
        input_dftb = ""

        # Geometry Block
        input_geom = textwrap.dedent(f"""\
        Geometry = GenFormat{{
          <<< 'geometry.gen'
        }}
        """)
        input_dftb += input_geom

        # Hamiltonian Block
        input_ham_init = textwrap.dedent(f"""\
        Hamiltonian = DFTB{{
        """)
        input_dftb += input_ham_init

        # SCC-DFTB option
        if (self.scc):
            input_ham_scc = textwrap.indent(textwrap.dedent(f"""\
              SCC = Yes
              SCCTolerance = {self.scc_tol}
              MaxSCCIterations = {self.scc_max_iter}
            """), "  ")
            input_dftb += input_ham_scc

            # Read atomic spin constants used in DFTB/SSR
            if (self.lcdftb):
                spin_constant = ("\n" + " " * 14).join([f"  {itype} = {{ {spin_w_lc[f'{itype}']} }}" for itype in self.atom_type])
            else:
                spin_constant = ("\n" + " " * 14).join([f"  {itype} = {{ {spin_w[f'{itype}']} }}" for itype in self.atom_type])
            input_ham_spin = textwrap.indent(textwrap.dedent(f"""\
              SpinConstants = {{
                ShellResolvedSpin = Yes
              {spin_constant}
              }}
            """), "  ")
            input_dftb += input_ham_spin

            # Long-range corrected DFTB (LC-DFTB) option
            if (self.lcdftb):
                input_ham_lc = textwrap.indent(textwrap.dedent(f"""\
                  RangeSeparated = LC{{
                    Screening = {self.lc_method}{{}}
                  }}
                """), "  ")
                input_dftb += input_ham_lc

            # Add point charges to Hamiltonian used in electrostatic charge embedding of QM/MM
            if (self.embedding == "electrostatic"):
                input_ham_pc = textwrap.indent(textwrap.dedent(f"""\
                  ElectricField = PointCharges{{
                    CoordsAndCharges [Angstrom] = DirectRead{{
                      Records = {molecule.nat_mm}
                      File = "point_charges.xyz"
                    }}
                  }}
                """), "  ")
                input_dftb += input_ham_pc

        if (self.periodic):
            input_ham_periodic = textwrap.indent(textwrap.dedent(f"""\
              KPointsAndWeights = {{
                0.0 0.0 0.0 1.0
              }}
            """), "  ")
            input_dftb += input_ham_periodic

        angular_momentum = ("\n" + " " * 10).join([f"  {itype} = '{max_l[f'{itype}']}'" for itype in self.atom_type])
        input_ham_basic = textwrap.dedent(f"""\
          Charge = {molecule.charge}
          MaxAngularMomentum = {{
          {angular_momentum}
          }}
          SlaterKosterFiles = Type2FileNames{{
            Prefix = '{self.sk_path}'
            Separator = '-'
            Suffix = '.skf'
            LowerCaseTypeName = No
          }}
        }}
        """)
        input_dftb += input_ham_basic

        # Analysis Block
        input_analysis = textwrap.dedent(f"""\
        Analysis = {{
          CalculateForces = Yes
          WriteBandOut = Yes
          WriteEigenvectors = Yes
          MullikenAnalysis = Yes
        }}
        """)
        input_dftb += input_analysis

        # Options Block
        input_options = textwrap.dedent(f"""\
        Options = {{
          WriteDetailedXml = Yes
          WriteDetailedOut = Yes
          TimingVerbosity = -1
        }}
        """)
        input_dftb += input_options

        # REKS Block

        # Energy functional options
        if (self.ssr22):
            # Set active space and energy functional used in SA-REKS
            space = "SSR22"
            if (molecule.nst == 1):
                energy_functional = "{ 'PPS' }"
                all_states = "No"
            elif (molecule.nst == 2):
                energy_functional = "{ 'PPS' 'OSS' }"
                all_states = "No"
            elif (molecule.nst == 3):
                energy_functional = "{ 'PPS' 'OSS' }"
                all_states = "Yes"
            else:
                raise ValueError (f"( {self.qm_method}.{call_name()} ) Too many electrnoic states! {molecule.nst}")

            # Include state-interaction terms to SA-REKS; SI-SA-REKS
            if (self.state_interactions):
                do_ssr = "Yes"
            else:
                do_ssr = "No"

            # Read 'eigenvec.bin' from previous step
            if (self.guess == "read"):
                if (istep == -1):
                    if (os.path.isfile(self.guess_file)):
                        # Copy guess file to currect directory
                        shutil.copy(self.guess_file, os.path.join(self.scr_qm_dir, "eigenvec.bin"))
                        restart = "Yes"
                    else:
                        restart = "No"
                elif (istep >= 0):
                    # Move previous file to currect directory
                    os.rename("../eigenvec.bin.pre", "./eigenvec.bin")
                    restart = "Yes"
            elif (self.guess == "h0"):
                restart = "No"

            # Scale the atomic spin constants
            if (self.tuning != None):
                spin_tuning = ""
                spin_tuning += "{"
                for scale_W in self.tuning:
                    spin_tuning += f" {scale_W} "
                spin_tuning += "}"
            else:
                spin_tuning = "{}"

        # CP-REKS algorithm options
        if (self.cpreks_grad_alg == "pcg"):
            cpreks_alg = "ConjugateGradient"
            preconditioner = "Yes"
        elif (self.cpreks_grad_alg == "cg"):
            cpreks_alg = "ConjugateGradient"
            preconditioner = "No"
        elif (self.cpreks_grad_alg == "direct"):
            cpreks_alg = "Direct"
            preconditioner = "No"
        else:
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong CP-REKS algorithm given! {self.cpreks_grad_alg}")

        # Save memory in cache to reduce computational cost in CP-REKS
        if (self.save_memory):
            memory = "Yes"
        else:
            memory = "No"
 
        # NAC calculation options
        if (molecule.nst == 1 or not self.state_interactions):
            # Single-state REKS or SA-REKS state
            self.nac = "No"
        else:
            # SSR state
            if (self.calc_coupling):
                # SHXF, SH, Eh need NAC calculations
                self.nac = "Yes"
            else:
                # BOMD do not need NAC calculations
                self.nac = "No"

        # Relaxed density options; It is determined automatically
        if (molecule.qmmm and self.embedding == "electrostatic"):
            relaxed_density = "Yes"
        else:
            relaxed_density = "No"

        input_reks = textwrap.dedent(f"""\
        REKS = {space}{{
          Energy = {{
            Functional = {energy_functional}
            StateInteractions = {do_ssr}
            IncludeAllStates = {all_states}
          }}
          TargetState = {bo_list[0] + 1}
          ReadEigenvectors = {restart}
          FONmaxIter = 50
          shift = {self.shift}
          SpinTuning = {spin_tuning}
          Gradient = {cpreks_alg}{{
            CGmaxIter = 100
            Tolerance = {self.cpreks_grad_tol}
            Preconditioner = {preconditioner}
            SaveMemory = {memory}
          }}
          RelaxedDensity = {relaxed_density}
          NonAdiabaticCoupling = {self.nac}
          VerbosityLevel = 1
        }}
        """)
        input_dftb += input_reks

        # ParserOptions Block
        if (self.version == "19.1"):
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) SSR not implemented in this version! {self.version}")
        elif (self.version == "20.1"):
            parser_version = 8

        input_parseroptions = textwrap.dedent(f"""\
        ParserOptions = {{
          ParserVersion = {parser_version}
        }}
        """)
        input_dftb += input_parseroptions

        # Write 'dftb_in.hsd' file
        file_name = "dftb_in.hsd"
        with open(file_name, "w") as f:
            f.write(input_dftb)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run SSR calculation and save the output files to QMlog directory

            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Set run command
        qm_command = os.path.join(self.qm_path, "dftb+")
        # OpenMP setting
        os.environ["OMP_NUM_THREADS"] = f"{self.nthreads}"
        command = f"{qm_command} > log"
        # Run DFTB+ method for molecular dynamics
        os.system(command)

        # Copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
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

        # Read 'detailed.out' file
        file_name = "detailed.out"
        with open(file_name, "r") as f:
            detailed_out = f.read()

        # Energy
        if (molecule.nst == 1):
            # Single-state REKS
            tmp_e = 'Spin' + '\n\s+\w+\s+([-]\S+)(?:\s+\S+){3}' * molecule.nst
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy)
        else:
            if (self.state_interactions):
                # SSR state
                energy = re.findall('SSR state\s+\S+\s+([-]\S+)', log_out)
                energy = np.array(energy)
            else:
                # SA-REKS state
                tmp_e = 'Spin' + '\n\s+\w+\s+([-]\S+)(?:\s+\S+){3}' * molecule.nst
                energy = re.findall(tmp_e, log_out)
                energy = np.array(energy[0])

        energy = energy.astype(float)
        for ist in range(molecule.nst):
            molecule.states[ist].energy = energy[ist]

        # Force
        if (self.nac == "Yes"):
            # SHXF, SH, Eh : SSR state
            if (molecule.qmmm):
                for ist in range(molecule.nst):
                    tmp_g = f' {ist + 1} st state \(SSR\)' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
                    grad = re.findall(tmp_g, log_out)
                    grad = np.array(grad[0])
                    grad = grad.astype(float)
                    grad = grad.reshape(molecule.nat_qm, 3, order='C')
                    molecule.states[ist].force[0:molecule.nat_qm] = - grad
                    if (self.embedding == "electrostatic"):
                        tmp_f = 'Forces on external charges' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_mm
                        force = re.findall(tmp_f, detailed_out)
                        force = np.array(force[0])
                        force = force.astype(float)
                        force = force.reshape(molecule.nat_mm, 3, order='C')
                        molecule.states[ist].force[molecule.nat_qm:molecule.nat] = force
            else:
                for ist in range(molecule.nst):
                    tmp_g = f' {ist + 1} st state \(SSR\)' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
                    grad = re.findall(tmp_g, log_out)
                    grad = np.array(grad[0])
                    grad = grad.astype(float)
                    grad = grad.reshape(molecule.nat_qm, 3, order='C')
                    molecule.states[ist].force = - np.copy(grad)
        else:
            # BOMD : SSR state, SA-REKS state or single-state REKS
            if (molecule.qmmm):
                tmp_g = f' {bo_list[0] + 1} state \(\w+[-]*\w+\)' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
                grad = re.findall(tmp_g, log_out)
                grad = np.array(grad[0])
                grad = grad.astype(float)
                grad = grad.reshape(molecule.nat_qm, 3, order='C')
                molecule.states[bo_list[0]].force[0:molecule.nat_qm] = - grad
                if (self.embedding == "electrostatic"):
                    tmp_f = 'Forces on external charges' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_mm
                    force = re.findall(tmp_f, detailed_out)
                    force = np.array(force[0])
                    force = force.astype(float)
                    force = force.reshape(molecule.nat_mm, 3, order='C')
                    molecule.states[bo_list[0]].force[molecule.nat_qm:molecule.nat] = force
            else:
                tmp_g = f' {bo_list[0] + 1} state \(\w+[-]*\w+\)' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
                grad = re.findall(tmp_g, log_out)
                grad = np.array(grad[0])
                grad = grad.astype(float)
                grad = grad.reshape(molecule.nat_qm, 3, order='C')
                molecule.states[bo_list[0]].force = - np.copy(grad)

        # NAC
        if (self.nac == "Yes"):
            kst = 0
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    tmp_c = 'non-adiabatic coupling' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
                    nac = re.findall(tmp_c, log_out)
                    nac = np.array(nac[kst])
                    nac = nac.astype(float)
                    nac = nac.reshape(molecule.nat_qm, 3, order='C')
                    molecule.nac[ist, jst] = nac
                    molecule.nac[jst, ist] = - nac
                    kst += 1


