from __future__ import division
from build.cioverlap import *
from bo.dftbplus.dftbplus import DFTBplus
from bo.dftbplus.dftbpar import spin_w, max_l
from misc import eV_to_au, call_name
import os, shutil, re, textwrap
import numpy as np

class DFTB(DFTBplus):
    """ Class for (TD)DFTB method of DFTB+ program

        :param object molecule: molecule object
        :param boolean scc: include SCC scheme
        :param double scc_tol: energy convergence for SCC iterations
        :param integer scc_max_iter: maximum number of SCC iterations
        :param string guess: initial guess for SCC iterations
        :param string guess_file: initial guess file
        :param boolean sdftb: include spin-polarisation scheme
        :param double unpaired_e: number of unpaired electrons
        :param double e_temp: electronic temperature for Fermi-Dirac scheme
        :param string mixer: charge mixing method used in DFTB
        :param string ex_symmetry: symmetry (singlet) in TDDFTB
        :param string sk_path: path for slater-koster files
        :param boolean periodic: use periodicity in the calculations
        :param double,list cell_length: the length of cell lattice
        :param string qm_path: path for QM binary
        :param string script_path: path for DFTB+ python script (dptools)
        :param integer nthreads: number of threads in the calculations
        :param boolean mpi: use MPI parallelization
        :param string mpi_path: path for MPI binary
        :param double version: version of DFTB+ program
    """
    def __init__(self, molecule, scc=True, scc_tol=1E-6, scc_max_iter=100, \
        guess="h0", guess_file="./charges.bin", sdftb=False, unpaired_e=0., e_temp=0., mixer="Broyden", \
        ex_symmetry="singlet", sk_path="./", periodic=False, cell_length=[0., 0., 0., 0., 0., 0., 0., 0., 0.,], \
        qm_path="./", script_path="./", nthreads=1, mpi=False, mpi_path="./", version=19.1):
        # Initialize DFTB+ common variables
        super(DFTB, self).__init__(molecule, sk_path, qm_path, script_path, nthreads, version)

        # Initialize DFTB+ DFTB variables
        self.scc = scc
        self.scc_tol = scc_tol
        self.scc_max_iter = scc_max_iter

        # Set initial guess for SCC term
        self.guess = guess
        self.guess_file = guess_file
        if (not (self.guess == "h0" or self.guess == "read")):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong input for initial guess option! {self.guess}")

        self.sdftb = sdftb
        self.unpaired_e = unpaired_e

        self.e_temp = e_temp
        self.mixer = mixer

        self.ex_symmetry = ex_symmetry

        self.mpi = mpi
        self.mpi_path = mpi_path

        self.periodic = periodic
        self.a_axis = np.zeros(3)
        self.b_axis = np.zeros(3)
        self.c_axis = np.zeros(3)
        self.a_axis = cell_length[0:3]
        self.b_axis = cell_length[3:6]
        self.c_axis = cell_length[6:9]

        # Check excitation symmetry in TDDFTB
        # TODO : Currently, disable triplet excited states with TDDFTB
#        if (not (self.ex_symmetry == "singlet" or self.ex_symmetry == "triplet")):
        if (not (self.ex_symmetry == "singlet")):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong symmetry for excited state! {self.ex_symmetry}")

        # Set 'l_nacme' and 're_calc' with respect to the computational method
        # TDDFTB do not produce NACs, so we should get NACME from CIoverlap
        # TDDFTB cannot compute the gradient of several states simultaneously.
        molecule.l_nacme = True
        self.re_calc = True

        # Calculate number of basis for current system
        # Set new variable to decide the position of basis functions in terms of atoms
        # DFTB method considers only valence electrons, so core electrons should be removed
        core_elec = 0.
        self.nbasis = 0
        self.check_atom = [0]
        for iat in range(molecule.nat):
            max_ang = max_l[molecule.symbols[iat]]
            if (max_ang == 's'):
                self.nbasis += 1
            elif (max_ang == 'p'):
                self.nbasis += 4
                core_elec += 2.
            self.check_atom.append(self.nbasis)

        # Set new variable to decide the position of atoms in terms of basis functions
        self.check_basis = []
        for ibasis in range(self.nbasis):
            for iat in range(molecule.nat):
                ind_a = self.check_atom[iat] + 1
                ind_b = self.check_atom[iat + 1]
                if (ibasis + 1 >= ind_a and ibasis + 1 <= ind_b):
                    self.check_basis.append(iat + 1)

        # Initialize NACME variables
        # There is no core orbitals in TDDFTB (fixed occupations)
        # nocc is number of occupied orbitals and nvirt is number of virtual orbitals
        self.norb = self.nbasis
        self.nocc = int(int(molecule.nelec - core_elec) / 2)
        self.nvirt = self.norb - self.nocc

        self.ao_overlap = np.zeros((self.nbasis, self.nbasis))
        self.mo_coef_old = np.zeros((self.norb, self.nbasis))
        self.mo_coef_new = np.zeros((self.norb, self.nbasis))
        self.ci_coef_old = np.zeros((molecule.nst, self.nocc, self.nvirt))
        self.ci_coef_new = np.zeros((molecule.nst, self.nocc, self.nvirt))

    def get_data(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFTB method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        self.copy_files(molecule, istep, calc_force_only)
        super().get_data(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, istep, bo_list, calc_force_only)
        self.run_QM(molecule, base_dir, istep, bo_list, calc_force_only)
        self.extract_BO(molecule, base_dir, istep, bo_list, dt, calc_force_only)
        self.move_dir(base_dir)

    def copy_files(self, molecule, istep, calc_force_only):
        """ Copy necessary scratch files in previous step

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Copy required files for NACME
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):
            # After T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "geometry.xyz"), \
                os.path.join(self.scr_qm_dir, "../geometry.xyz.pre"))
            if (istep == 0):
                shutil.copy(os.path.join(self.scr_qm_dir, "eigenvec.bin"), \
                    os.path.join(self.scr_qm_dir, "../eigenvec.bin.pre"))
                shutil.copy(os.path.join(self.scr_qm_dir, "SPX.DAT"), \
                    os.path.join(self.scr_qm_dir, "../SPX.DAT.pre"))
                shutil.copy(os.path.join(self.scr_qm_dir, "XplusY.DAT"), \
                    os.path.join(self.scr_qm_dir, "../XplusY.DAT.pre"))

        # Copy required files to read initial guess
        if (self.guess == "read" and istep >= 0):
            # After T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "charges.bin"), \
                os.path.join(self.scr_qm_dir, "../charges.bin.pre"))

    def get_input(self, molecule, istep, bo_list, calc_force_only):
        """ Generate DFTB+ input files: geometry.gen, dftb_in.hsd

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
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
                    row = f'{molecule.nat} S\n'
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

        # Make 'double.gen' file for CIoverlap in TDDFTB
        # In this case, we do not need to consider periodicity
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):
            # Move previous files to currect directory
            os.rename('../geometry.xyz.pre', './geometry.xyz.pre')
            if (istep == 0):
                os.rename('../eigenvec.bin.pre', './eigenvec.bin.pre')
                os.rename('../SPX.DAT.pre', './SPX.DAT.pre')
                os.rename('../XplusY.DAT.pre', './XplusY.DAT.pre')
            # Open 'geometry.xyz.pre'
            file_af = open('double.xyz', 'w')
            file_be = open('geometry.xyz.pre', 'r')
            first_row = True
            for row in file_be:
                if (first_row):
                    row = f'{molecule.nat * 2}\n'
                    first_row = False
                file_af.write(row)
            file_be.close()
            # Open 'geometry.xyz'
            file_be = open('geometry.xyz', 'r')
            iline = 1
            for row in file_be:
                if (iline > 2):
                    file_af.write(row)
                iline += 1
            file_be.close()
            file_af.close()
            os.system("xyz2gen double.xyz")

        # Generate input files for multiple evaluation of force at same geometry
        irun = 1
        for ist in bo_list:
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

            if (self.scc):
                input_ham_scc = textwrap.indent(textwrap.dedent(f"""\
                  SCC = Yes
                  SCCTolerance = {self.scc_tol}
                  MaxSCCIterations = {self.scc_max_iter}
                  Mixer = {self.mixer}{{}}
                """), "  ")
                input_dftb += input_ham_scc

                if (self.sdftb and molecule.nst == 1):
                    input_ham_spin = textwrap.dedent(f"""\
                    SpinPolarisation = Colinear{{
                      UnpairedElectrons = {self.unpaired_e}
                    }}
                    """)
                    input_dftb += input_ham_spin

                # TODO : Currently, disable triplet excited states with TDDFTB
#                if (self.sdftb or self.ex_symmetry == "triplet"):
                if (self.sdftb and molecule.nst == 1):
                    spin_constant = ("\n" + " " * 18).join([f"  {itype} = {{ {spin_w[f'{itype}']} }}" for itype in self.atom_type])
                    input_ham_spin_w = textwrap.indent(textwrap.dedent(f"""\
                      SpinConstants = {{
                        ShellResolvedSpin = Yes
                      {spin_constant}
                      }}
                    """), "  ")
                    input_dftb += input_ham_spin_w

                # Read 'charges.bin' from previous step
                if (irun == 1):
                    if (self.guess == "read"):
                        if (istep == -1):
                            if (os.path.isfile(self.guess_file)):
                                # Copy guess file to currect directory
                                shutil.copy(self.guess_file, os.path.join(self.scr_qm_dir, "charges.bin"))
                                restart = "Yes"
                            else:
                                restart = "No"
                        elif (istep >= 0):
                            # Move previous file to currect directory
                            os.rename("../charges.bin.pre", "./charges.bin")
                            restart = "Yes"
                    elif (self.guess == "h0"):
                        restart = "No"

                # Read 'charges.bin' for Ehrenfest or surface hopping when hop occurs
                if (calc_force_only or irun > 1):
                    restart = "Yes"

                input_ham_restart = textwrap.indent(textwrap.dedent(f"""\
                  ReadInitialCharges = {restart}
                """), "  ")
                input_dftb += input_ham_restart

            # TODO: for QM/MM, point_charge??

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
              Filling = Fermi{{
                Temperature[K] = {self.e_temp}
              }}
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

            # ExcitedState Block
            if (len(bo_list) == 1 or irun > 1):
                do_excited = True
            else:
                do_excited = False

            if (molecule.nst > 1 and do_excited):
                # Calculate excited state force for target state
                if (ist > 0):
                    ex_force = "Yes"
                    rst = ist
                else:
                    ex_force = "No"
                    rst = ist + 1

                # Set number of excitations in TDDFTB
                # This part can be modified by users
                if (molecule.nat <= 5):
                    num_ex = molecule.nst + 2
                elif (molecule.nat > 5 and molecule.nat <= 15):
                    num_ex = 2 * molecule.nst + 2
                else:
                    num_ex = 3 * molecule.nst + 2

                # Write XplusY data?
                if (self.calc_coupling):
                    xpy = "Yes"
                else:
                    xpy = "No"

                input_excited = textwrap.dedent(f"""\
                ExcitedState = Casida{{
                  NrOfExcitations = {num_ex}
                  StateOfInterest = {rst}
                  Symmetry = {self.ex_symmetry}
                  WriteTransitions = Yes
                  WriteSPTransitions = {xpy}
                  WriteMulliken = Yes
                  WriteXplusY = {xpy}
                  ExcitedStateForces = {ex_force}
                }}
                """)
                input_dftb += input_excited

            # ParserOptions Block
            if (self.version == 19.1):
                parser_version = 7
            elif (self.version == 20.1):
                parser_version = 8

            input_parseroptions = textwrap.dedent(f"""\
            ParserOptions = {{
              ParserVersion = {parser_version}
            }}
            """)
            input_dftb += input_parseroptions

            # Parallel Block
            if (self.mpi):
                if (self.sdftb and self.nthreads > 1):
                    groups = 2
                else:
                    groups = 1
                input_parallel = textwrap.dedent(f"""\
                Parallel = {{
                  Groups = {groups}
                  UseOmpThreads = No
                  Blacs = BlockSize {{ 32 }}
                }}
                """)
                input_dftb += input_parallel

            # Write 'dftb_in.hsd.geom' file
            file_name = f"dftb_in.hsd.geom.{ist}"
            with open(file_name, "w") as f:
                f.write(input_dftb)

            irun += 1

        # Write 'dftb_in.hsd.double' file
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):
            # New input for dftb
            input_dftb = ""

            # Geometry Block
            input_geom = textwrap.dedent(f"""\
            Geometry = GenFormat{{
              <<< 'double.gen'
            }}
            """)
            input_dftb += input_geom
            input_dftb += input_ham_init
            input_dftb += input_ham_basic

            # Options Block
            input_options = textwrap.dedent(f"""\
            Options = {{
              WriteDetailedXml = Yes
              WriteDetailedOut = Yes
              WriteHS = Yes
              TimingVerbosity = -1
            }}
            """)
            input_dftb += input_options

            file_name = "dftb_in.hsd.double"
            with open(file_name, "w") as f:
                f.write(input_dftb)

    def run_QM(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Run (TD)DFTB calculation and save the output files to QMlog directory

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Set run command
        qm_command = os.path.join(self.qm_path, "dftb+")
        if (self.mpi):
            # MPI setting
            os.environ["OMP_NUM_THREADS"] = "1"
            mpi_command = os.path.join(self.mpi_path, "mpirun")
            command = f"{mpi_command} -np {self.nthreads} {qm_command} > log"
        else:
            # OpenMP setting
            os.environ["OMP_NUM_THREADS"] = f"{self.nthreads}"
            command = f"{qm_command} > log"

        # Run DFTB+ for calculation of overlap matrix
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):
            shutil.copy("dftb_in.hsd.double", "dftb_in.hsd")
            os.system(command)

        # Run DFTB+ method for molecular dynamics
        tmp_dir = os.path.join(base_dir, "QMlog")
        for ist in bo_list:
            # Copy dftb_in.hsd for target state
            file_name = f"dftb_in.hsd.geom.{ist}"
            shutil.copy(file_name, "dftb_in.hsd")
            os.system(command)

            # Copy detailed.out for target state
            file_name = f"detailed.out.{ist}"
            shutil.copy("detailed.out", file_name)

            # Copy the output file to 'QMlog' directory
            if (os.path.exists(tmp_dir)):
                detailed_out_step = f"detailed.out.{istep + 1}.{ist}"
                shutil.copy("detailed.out", os.path.join(tmp_dir, detailed_out_step))
                log_step = f"log.{istep + 1}.{ist}"
                shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Read 'detailed.out' file
        # TODO: the qmmm information is written in this file
        file_name = f"detailed.out.{bo_list[0]}"
        with open(file_name, "r") as f:
            detailed_out = f.read()

        # Read 'EXC.DAT' file
        if (molecule.nst > 1):
            file_name = "EXC.DAT"
            with open(file_name, "r") as f:
                exc_out = f.read()

        # Energy
        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            energy = re.findall('Total energy:\s+([-]\S+) H', detailed_out)
            energy = np.array(energy[0])
            energy = energy.astype(float)
            molecule.states[0].energy = energy

            if (molecule.nst > 1):
                tmp_e = f'[=]+\n' + ('\s+([-]*\S+)\s+\S+\s+\d+\s+->\s+\d+\s+\S+\s+\S+\s+[ST]') * molecule.nst
                energy = re.findall(tmp_e, exc_out)
                energy = np.array(energy[0])
                energy = energy.astype(float)
                energy *= eV_to_au
                for ist in range(1, molecule.nst):
                    molecule.states[ist].energy = molecule.states[0].energy + energy[ist - 1]

        # Force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        for ist in bo_list:
            file_name = f"detailed.out.{ist}"
            with open(file_name, "r") as f:
                detailed_out = f.read()

            tmp_f = 'Total Forces' + '\n\s+\d*\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
            force = re.findall(tmp_f, detailed_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[ist].force = np.copy(force)

        # NACME
        if (self.calc_coupling and not calc_force_only):
            molecule.nacme = np.zeros((molecule.nst, molecule.nst))
            if (istep >= 0):
                self.CI_overlap(molecule, istep, dt)

    def CI_overlap(self, molecule, istep, dt):
        """ Read the necessary files and calculate NACME from tdnac.c routine,
            note that only reading of several files is required in this method

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param double dt: time interval
        """
        # Read upper right block of 'oversqr.dat' file (< t | t+dt >)
        file_name_in = "oversqr.dat"

        self.ao_overlap = np.zeros((self.nbasis, self.nbasis))
        with open(file_name_in, "r") as f_in:
            lines = f_in.readlines()
            row = 0
            iline = 0
            for line in lines:
                # Skip first five lines and read upper block
                if (iline in range(5, 5 + self.nbasis)):
                    col = 0
                    count = False
                    field = line.split()
                    for element in field:
                        # Read right block
                        if (count):
                            ind_a = self.check_basis[row]
                            ind_b = self.check_basis[col]
                            if (ind_a == ind_b):
                                # Choose onsite (same-atom) block
                                # Sometimes NaN or too large values appear in the onsite block due to the slater-koster file
                                # The values set to 1 or 0 regardless of original elements
                                if (row == col):
                                    # Diagonal element in onsite block
                                    new_val = 1.
                                else:
                                    # Off-diagonal element in onsite block
                                    new_val = 0.
                            else:
                                # Choose offsite (different-atom) block
                                new_val = float(element)
                            # Set overlap matrix element
                            self.ao_overlap[row, col] = new_val
                        col += 1
                        # Read right block
                        if (col > self.nbasis - 1):
                            col -= self.nbasis
                            count = True
                    row += 1
                iline += 1
#        np.savetxt("test-over", self.ao_overlap, fmt=f"%6.3f")

        # Read 'eigenvec.bin.pre' file at time t
        if (istep == 0):
            file_name_in = "eigenvec.bin.pre"

            self.mo_coef_old = np.zeros((self.norb, self.nbasis))
            with open(file_name_in, "rb") as f_in:
                dummy = np.fromfile(f_in, dtype=np.integer, count=1)
                for iorb in range(self.norb):
                    dummy = np.fromfile(f_in, dtype=np.integer, count=1)
                    data = np.fromfile(f_in, dtype=np.float64, count=self.nbasis)
                    self.mo_coef_old[iorb] = data
#            np.savetxt("test-mo1", self.mo_coef_old, fmt=f"%12.6f")

        # Read 'eigenvec.bin' file at time t + dt
        file_name_in = "eigenvec.bin"

        self.mo_coef_new = np.zeros((self.norb, self.nbasis))
        with open(file_name_in, "rb") as f_in:
            dummy = np.fromfile(f_in, dtype=np.integer, count=1)
            for iorb in range(self.norb):
                dummy = np.fromfile(f_in, dtype=np.integer, count=1)
                data = np.fromfile(f_in, dtype=np.float64, count=self.nbasis)
                self.mo_coef_new[iorb] = data
#        np.savetxt("test-mo2", self.mo_coef_new, fmt=f"%12.6f")

        # Dimension for CI coefficients
        nmat = self.nocc * self.nvirt

        # The CI coefficients are arranged in order of single-particle excitations
        # Read 'SPX.DAT.pre' file at time t
        if (istep == 0):
            file_name_in = "SPX.DAT.pre"

            get_wij_ind_old = np.zeros((nmat, 2), dtype=np.int_)
            with open(file_name_in, "r") as f_in:
                lines = f_in.readlines()
                iline = 0
                for line in lines:
                    # Skip first five lines
                    if (iline in range(5, 5 + nmat)):
                        # Column information: 1st = index, 4th = occ(i), 6th = virt(a)
                        field = line.split()
                        get_wij_ind_old[int(field[0]) - 1] = [int(field[3]), int(field[5])]
                    iline += 1

        # Read 'SPX.DAT' file at time t + dt
        file_name_in = "SPX.DAT"

        get_wij_ind_new = np.zeros((nmat, 2), dtype=np.int_)
        with open(file_name_in, "r") as f_in:
            lines = f_in.readlines()
            iline = 0
            for line in lines:
                # Skip first five lines
                if (iline in range(5, 5 + nmat)):
                    # Column information: 1st = index, 4th = occ(i), 6th = virt(a)
                    field = line.split()
                    get_wij_ind_new[int(field[0]) - 1] = [int(field[3]), int(field[5])]
                iline += 1

        # Read 'XplusY.DAT.pre' file at time t
        if (istep == 0):
            file_name_in = "XplusY.DAT.pre"

            self.ci_coef_old = np.zeros((molecule.nst, self.nocc, self.nvirt))
            with open(file_name_in, "r") as f_in:
                lines = f_in.readlines()
                iline = 0
                for line in lines:
                    if (iline == 0):
                        field = line.split()
                        assert (int(field[0]) == nmat)
                        assert (int(field[1]) >= molecule.nst - 1)
                        # nxply is number of lines for each excited state in 'XplusY.dat'
                        nxply = int(nmat / 6) + 1
                        if (nmat % 6 != 0):
                            nxply += 1
                    else:
                        field = line.split()
                        if (iline % nxply == 1):
                            ind = 0
                            ist = int(field[0]) - 1
                            # In general, TDDFTB calculate the excited states more than molecule.nst,
                            # so we do not need to read all data for 'XplusY.DAT'
                            if (ist == molecule.nst - 1):
                                break
                        else:
                            # Currently, elements for CI coefficients for S0 state are zero (not used values)
                            for element in field:
                                ind_occ = get_wij_ind_old[ind, 0] - 1
                                ind_virt = get_wij_ind_old[ind, 1] - self.nocc - 1
                                self.ci_coef_old[ist + 1, ind_occ, ind_virt] = float(element)
                                ind += 1
                    iline += 1
#            np.savetxt("test-ci1", self.ci_coef_old[1], fmt=f"%12.6f")

        # Read 'XplusY.DAT' file at time t + dt
        file_name_in = "XplusY.DAT"

        self.ci_coef_new = np.zeros((molecule.nst, self.nocc, self.nvirt))
        with open(file_name_in, "r") as f_in:
            lines = f_in.readlines()
            iline = 0
            for line in lines:
                if (iline == 0):
                    field = line.split()
                    assert (int(field[0]) == nmat)
                    assert (int(field[1]) >= molecule.nst - 1)
                    # nxply is number of lines for each excited state in 'XplusY.dat'
                    nxply = int(nmat / 6) + 1
                    if (nmat % 6 != 0):
                        nxply += 1
                else:
                    field = line.split()
                    if (iline % nxply == 1):
                        ind = 0
                        ist = int(field[0]) - 1
                        # In general, TDDFTB calculate the excited states more than molecule.nst,
                        # so we do not need to read all data for 'XplusY.DAT'
                        if (ist == molecule.nst - 1):
                            break
                    else:
                        # Currently, elements for CI coefficients for S0 state are zero (not used values)
                        for element in field:
                            ind_occ = get_wij_ind_new[ind, 0] - 1
                            ind_virt = get_wij_ind_new[ind, 1] - self.nocc - 1
                            self.ci_coef_new[ist + 1, ind_occ, ind_virt] = float(element)
                            ind += 1
                iline += 1
#        np.savetxt("test-ci2", self.ci_coef_new[1], fmt=f"%12.6f")

        # Calculate wavefunction overlap with orbital scheme
        # Reference: J. Phys. Chem. Lett. 2015, 6, 4200-4203
        wf_overlap(self, molecule, istep, dt)


