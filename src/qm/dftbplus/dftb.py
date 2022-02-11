from __future__ import division
from build.cioverlap import *
from qm.dftbplus.dftbplus import DFTBplus
from qm.dftbplus.dftbpar import spin_w, spin_w_lc, onsite_uu, onsite_ud, max_l
from misc import data, eps, eV_to_au, call_name
import os, shutil, re, textwrap
import numpy as np

class DFTB(DFTBplus):
    """ Class for (TD)DFTB method of DFTB+

        :param object molecule: Molecule object
        :param boolean l_scc: Include self-consistent charge (SCC) scheme
        :param double scc_tol: Stopping criteria for the SCC iterations
        :param integer scc_max_iter: Maximum number of SCC iterations
        :param boolean l_onsite: Include onsite correction to SCC term
        :param boolean l_range_sep: Include long-range corrected functional
        :param string lc_method: Algorithms for LC-DFTB
        :param boolean l_spin_pol: Include spin-polarisation scheme
        :param double unpaired_elec: Number of unpaired electrons
        :param string guess: Initial guess method for SCC scheme
        :param string guess_file: Initial guess file for charges
        :param double elec_temp: Electronic temperature in Fermi-Dirac scheme
        :param string mixer: Charge mixing method used in DFTB
        :param string ex_symmetry: Symmetry of excited state in TDDFTB
        :param double e_window: Energy window for TDDFTB. Increases efficiency of NACME calculation.
        :param integer,list k_point: Number of k-point samplings
        :param boolean l_periodic: Use periodicity in the calculations
        :param double,list cell_length: The lattice vectors of periodic unit cell
        :param string sk_path: Path for Slater-Koster files
        :param string install_path: Path for DFTB+ install directory
        :param boolean mpi: Use MPI parallelization
        :param string mpi_path: Path for MPI binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of DFTB+
    """
    def __init__(self, molecule, l_scc=True, scc_tol=1E-6, scc_max_iter=100, l_onsite=False, \
        l_range_sep=False, lc_method="MatrixBased", l_spin_pol=False, unpaired_elec=0., guess="h0", \
        guess_file="./charges.bin", elec_temp=0., mixer="Broyden", ex_symmetry="singlet", e_window=0., \
        k_point=[1, 1, 1], l_periodic=False, cell_length=[0., 0., 0., 0., 0., 0., 0., 0., 0.,], \
        sk_path="./", install_path="./", mpi=False, mpi_path="./", nthreads=1, version="20.1"):
        # Initialize DFTB+ common variables
        super(DFTB, self).__init__(molecule, sk_path, install_path, nthreads, version)

        # Initialize DFTB+ DFTB variables
        self.l_scc = l_scc
        self.scc_tol = scc_tol
        self.scc_max_iter = scc_max_iter

        self.l_onsite = l_onsite

        self.l_range_sep = l_range_sep
        self.lc_method = lc_method.lower()

        self.l_spin_pol = l_spin_pol
        self.unpaired_elec = unpaired_elec

        # Set initial guess for SCC term
        self.guess = guess.lower()
        self.guess_file = guess_file
        if not (self.guess in ["h0", "read"]):
            error_message = "Invalid initial guess for DFTB!"
            error_vars = f"guess = {self.guess}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.elec_temp = elec_temp
        self.mixer = mixer.lower()

        self.ex_symmetry = ex_symmetry.lower()
        self.e_window = e_window

        self.k_point = k_point
        self.l_periodic = l_periodic
        self.a_axis = np.copy(cell_length[0:3])
        self.b_axis = np.copy(cell_length[3:6])
        self.c_axis = np.copy(cell_length[6:9])

        # Check excitation symmetry in TDDFTB
        # TODO : Currently, allows only singlet excited states with TDDFTB
#        if not (self.ex_symmetry in ["singlet", "triplet"]):
        if (not self.ex_symmetry == "singlet"):
            error_message = "Invalid symmetry of excited states for TDDFTB given!"
            error_vars = f"ex_symmetry = {self.ex_symmetry}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.mpi = mpi
        self.mpi_path = mpi_path

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
        for iat in range(molecule.nat_qm):
            # Check number of basis functions with respect to maximum angular momentum
            max_ang = max_l[molecule.symbols[iat]]
            if (max_ang == 's'):
                self.nbasis += 1
            elif (max_ang == 'p'):
                self.nbasis += 4
            elif (max_ang == 'd'):
                self.nbasis += 9
            else:
                error_message = "Number of basis for f orbital not implemented, see '$PYUNIXMDHOME/src/qm/dftbplus/dftb.py'!"
                error_vars = f"current atom = {molecule.symbols[iat]}, max_ang = {max_ang}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
            self.check_atom.append(self.nbasis)
            # Check number of core electrons with respect to atomic number
            sym_index = list(data.keys()).index(molecule.symbols[iat])
            if (sym_index > 0 and sym_index <= 2):
                core_elec += 0.
            elif (sym_index > 2 and sym_index <= 10):
                core_elec += 2.
            elif (sym_index > 10 and sym_index <= 18):
                core_elec += 10.
            elif (sym_index > 18 and sym_index <= 36):
                core_elec += 18.
            elif (sym_index > 36 and sym_index <= 54):
                core_elec += 36.
            else:
                error_message = "Core electrons for current element not implemented, see '$PYUNIXMDHOME/src/qm/dftbplus/dftb.py'!"
                error_vars = f"current atom = {molecule.symbols[iat]}, sym_index = {sym_index}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set new variable to decide the position of atoms in terms of basis functions
        self.check_basis = []
        for ibasis in range(self.nbasis):
            for iat in range(molecule.nat_qm):
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
    
        # Replace norb by arrays containing the limits of the for loops.
        # For energy window calculations loops will not go from (0 to nocc/nvirt) or (0 to norb)
        # but from (nocc_min to nocc/0 to nvirt_max) or (nocc_min to norb).
        self.orb_ini = np.zeros(1, dtype=np.int32)
        self.orb_final = np.zeros(1, dtype=np.int32)
        self.orb_final[0] = self.norb

        if (self.e_window > eps):
            # Swap minimal/maximal values to replace them in reading of SPX.DAT by the minimal/maximal values.
            self.orb_ini[0] = self.norb
            self.orb_final[0] = 0

        self.ao_overlap = np.zeros((self.nbasis, self.nbasis))
        self.mo_coef_old = np.zeros((self.norb, self.nbasis))
        self.mo_coef_new = np.zeros((self.norb, self.nbasis))
        self.ci_coef_old = np.zeros((molecule.nst, self.nocc, self.nvirt))
        self.ci_coef_new = np.zeros((molecule.nst, self.nocc, self.nvirt))

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFTB method

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer,list bo_list: List of BO states for BO calculation
            :param double dt: Time interval
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        self.copy_files(molecule, istep, calc_force_only)
        super().get_data(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, istep, bo_list, calc_force_only)
        self.run_QM(molecule, base_dir, istep, bo_list, calc_force_only)
        self.extract_QM(molecule, base_dir, istep, bo_list, dt, calc_force_only)
        self.move_dir(base_dir)

    def copy_files(self, molecule, istep, calc_force_only):
        """ Copy necessary scratch files in previous step

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
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

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Make 'geometry.gen' file
        os.system("xyz2gen geometry.xyz")
        if (self.l_periodic):
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
                    row = f'{molecule.nat_qm * 2}\n'
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
        if (self.l_scc):
            input_ham_scc = textwrap.indent(textwrap.dedent(f"""\
              SCC = Yes
              SCCTolerance = {self.scc_tol}
              MaxSCCIterations = {self.scc_max_iter}
              Mixer = {self.mixer}{{}}
            """), "  ")
            input_dftb += input_ham_scc

            # Onsite-corrected DFTB (OC-DFTB) option
            if (self.l_onsite):
                onsite_const_uu = ("\n" + " " * 18).join([f"  {itype}uu = {{ {onsite_uu[f'{itype}']} }}" for itype in self.atom_type])
                onsite_const_ud = ("\n" + " " * 18).join([f"  {itype}ud = {{ {onsite_ud[f'{itype}']} }}" for itype in self.atom_type])
                input_ham_oc = textwrap.indent(textwrap.dedent(f"""\
                  OnsiteCorrection = {{
                  {onsite_const_uu}
                  {onsite_const_ud}
                  }}
                """), "  ")
                input_dftb += input_ham_oc

            # Long-range corrected DFTB (LC-DFTB) option
            if (self.l_range_sep):
                input_ham_lc = textwrap.indent(textwrap.dedent(f"""\
                  RangeSeparated = LC{{
                    Screening = {self.lc_method}{{}}
                  }}
                """), "  ")
                input_dftb += input_ham_lc

            # Spin-polarized DFTB option
            if (self.l_spin_pol and molecule.nst == 1):
                input_ham_spin = textwrap.dedent(f"""\
                SpinPolarisation = Colinear{{
                  UnpairedElectrons = {self.unpaired_elec}
                }}
                """)
                input_dftb += input_ham_spin

            # Read atomic spin constants used in spin-polarized DFTB or TDDFTB
            # TODO : Currently, allows only singlet excited states with TDDFTB
#            if (self.l_spin_pol or self.ex_symmetry == "triplet"):
            if (self.l_spin_pol and molecule.nst == 1):
                if (self.l_range_sep):
                    spin_constant = ("\n" + " " * 18).join([f"  {itype} = {{ {spin_w_lc[f'{itype}']} }}" for itype in self.atom_type])
                else:
                    spin_constant = ("\n" + " " * 18).join([f"  {itype} = {{ {spin_w[f'{itype}']} }}" for itype in self.atom_type])
                input_ham_spin_w = textwrap.indent(textwrap.dedent(f"""\
                  SpinConstants = {{
                    ShellResolvedSpin = Yes
                  {spin_constant}
                  }}
                """), "  ")
                input_dftb += input_ham_spin_w

            # Read 'charges.bin' from previous step
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

            # Read 'charges.bin' for surface hopping dynamics when hop occurs
            if (calc_force_only):
                restart = "Yes"

            input_ham_restart = textwrap.indent(textwrap.dedent(f"""\
              ReadInitialCharges = {restart}
            """), "  ")
            input_dftb += input_ham_restart

        # TODO: for QM/MM, point_charge??

        if (self.l_periodic):
            num_k_point = np.sum(self.k_point)
            if (num_k_point == 3):
                # gamma-point sampling
                input_ham_periodic = textwrap.indent(textwrap.dedent(f"""\
                  KPointsAndWeights = {{
                    0.0 0.0 0.0 1.0
                  }}
                """), "  ")
            else:
                # K-point sampling
                shift_vector = [0.5 if (ik % 2 == 0) else 0 for ik in self.k_point]
                input_ham_periodic = textwrap.indent(textwrap.dedent(f"""\
                  KPointsAndWeights = SupercellFolding{{
                    {self.k_point[0]} 0.0 0.0
                    0.0 {self.k_point[1]} 0.0
                    0.0 0.0 {self.k_point[2]}
                    {shift_vector[0]} {shift_vector[1]} {shift_vector[2]}
                  }}
                """), "  ")
            input_dftb += input_ham_periodic

        angular_momentum = ("\n" + " " * 10).join([f"  {itype} = '{max_l[f'{itype}']}'" for itype in self.atom_type])
        input_ham_basic = textwrap.dedent(f"""\
          Charge = {molecule.charge}
          Filling = Fermi{{
            Temperature[K] = {self.elec_temp}
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
        if (molecule.nst > 1):
            # Calculate excited state force for target state
            if (bo_list[0] > 0):
                ex_force = "Yes"
                rst = bo_list[0]
            else:
                ex_force = "No"
                rst = bo_list[0] + 1

            # Set number of excitations in TDDFTB
            # This part can be modified by users
            if (molecule.nat_qm <= 5):
                num_ex = molecule.nst + 2
            elif (molecule.nat_qm > 5 and molecule.nat_qm <= 15):
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
              EnergyWindow [eV] = {self.e_window}
              ExcitedStateForces = {ex_force}
            }}
            """)

            input_dftb += input_excited

        # ParserOptions Block
        if (self.version == "19.1"):
            parser_version = 7
        elif (self.version == "20.1"):
            parser_version = 8
        elif (self.version == "21.1"):
            parser_version = 9

        input_parseroptions = textwrap.dedent(f"""\
        ParserOptions = {{
          ParserVersion = {parser_version}
        }}
        """)
        input_dftb += input_parseroptions

        # Parallel Block
        if (self.mpi):
            if (self.l_spin_pol and self.nthreads > 1):
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
        file_name = f"dftb_in.hsd.geom.{bo_list[0]}"
        with open(file_name, "w") as f:
            f.write(input_dftb)

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
        """ Run (TD)DFTB calculation and save the output files to qm_log directory

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
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

        # Copy dftb_in.hsd for target state
        file_name = f"dftb_in.hsd.geom.{bo_list[0]}"
        shutil.copy(file_name, "dftb_in.hsd")

        # Run DFTB+ method for molecular dynamics
        os.system(command)

        # Copy detailed.out for target state
        file_name = f"detailed.out.{bo_list[0]}"
        shutil.copy("detailed.out", file_name)

        # Copy the output file to 'qm_log' directory
        tmp_dir = os.path.join(base_dir, "qm_log")
        if (os.path.exists(tmp_dir)):
            detailed_out_step = f"detailed.out.{istep + 1}.{bo_list[0]}"
            shutil.copy("detailed.out", os.path.join(tmp_dir, detailed_out_step))
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_QM(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
            :param double dt: Time interval
            :param boolean calc_force_only: Logical to decide whether calculate force only
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
            energy = re.findall('Total energy:\s+([-]\S+) H', detailed_out)
            energy = np.array(energy[0], dtype=np.float64)
            molecule.states[0].energy = energy

            if (molecule.nst > 1):
                tmp_e = f'[=]+\n' + ('\s+([-]*\S+)\s+\S+\s+\d+\s+->\s+\d+\s+\S+\s+\S+\s+[ST]') * molecule.nst
                energy = re.findall(tmp_e, exc_out)
                energy = np.array(energy[0], dtype=np.float64)
                energy *= eV_to_au
                for ist in range(1, molecule.nst):
                    molecule.states[ist].energy = molecule.states[0].energy + energy[ist - 1]

        # Force
        tmp_f = 'Total Forces' + '\n\s+\d*\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
        force = re.findall(tmp_f, detailed_out)
        force = np.array(force[0], dtype=np.float64)
        force = force.reshape(molecule.nat_qm, 3, order='C')
        molecule.states[bo_list[0]].force = np.copy(force)

        # NACME
        if (self.calc_coupling and not calc_force_only):
            if (istep >= 0):
                self.CI_overlap(molecule, istep, dt)

    def CI_overlap(self, molecule, istep, dt):
        """ Read the necessary files and calculate NACME from tdnac.c routine,
            note that only reading of several files is required in this method

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param double dt: Time interval
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
        dum_byte = 4
        if (istep == 0):
            file_name_in = "eigenvec.bin.pre"

            self.mo_coef_old = np.zeros((self.norb, self.nbasis))
            with open(file_name_in, "rb") as f_in:
                f_in.seek(dum_byte, 1)
                dummy = np.fromfile(f_in, dtype=np.int32, count=1)
                f_in.seek(dum_byte, 1)
                for iorb in range(self.norb):
                    f_in.seek(dum_byte, 1)
                    data = np.fromfile(f_in, dtype=np.float64, count=self.nbasis)
                    f_in.seek(dum_byte, 1)
                    self.mo_coef_old[iorb] = data
#            np.savetxt("test-mo1", self.mo_coef_old, fmt=f"%12.6f")

        # Read 'eigenvec.bin' file at time t + dt
        file_name_in = "eigenvec.bin"

        self.mo_coef_new = np.zeros((self.norb, self.nbasis))
        with open(file_name_in, "rb") as f_in:
            f_in.seek(dum_byte, 1)
            dummy = np.fromfile(f_in, dtype=np.int32, count=1)
            f_in.seek(dum_byte, 1)
            for iorb in range(self.norb):
                f_in.seek(dum_byte, 1)
                data = np.fromfile(f_in, dtype=np.float64, count=self.nbasis)
                f_in.seek(dum_byte, 1)
                self.mo_coef_new[iorb] = data
#        np.savetxt("test-mo2", self.mo_coef_new, fmt=f"%12.6f")

        # The CI coefficients are arranged in order of single-particle excitations
        # Read 'SPX.DAT.pre' file at time t
        if (istep == 0):
            file_name_in = "SPX.DAT.pre"

            with open(file_name_in, "r") as f_in:
                lines = f_in.readlines()
                # Dimension for CI coefficients (number of excitations)
                ndim_old = int(lines[-2].strip().split()[0])
                get_wij_ind_old = np.zeros((ndim_old, 2), dtype=np.int32)
                iline = 0
                for line in lines:
                    # Skip first five lines
                    if (iline in range(5, 5 + ndim_old)):
                        # Column information: 1st = index, 4th = occ(i), 6th = virt(a)
                        field = line.split()
                        # Determine new limits for for-loops
                        if (int(field[5]) > self.orb_final[0]):
                            self.orb_final[0] = int(field[5])
                        if (int(field[3]) < self.orb_ini[0] + 1):
                            self.orb_ini[0] = int(field[3]) - 1
                        get_wij_ind_old[int(field[0]) - 1] = [int(field[3]), int(field[5])]
                    iline += 1

        # Read 'SPX.DAT' file at time t + dt
        file_name_in = "SPX.DAT"

        with open(file_name_in, "r") as f_in:
            lines = f_in.readlines()
            # Dimension for CI coefficients (number of excitations)
            ndim = int(lines[-2].strip().split()[0])
            get_wij_ind_new = np.zeros((ndim, 2), dtype=np.int32)
            iline = 0
            for line in lines:
                # Skip first five lines
                if (iline in range(5, 5 + ndim)):
                    # Column information: 1st = index, 4th = occ(i), 6th = virt(a)
                    field = line.split()
                    if (int(field[5]) > self.orb_final[0]):
                        self.orb_final[0] = int(field[5])
                    if (int(field[3]) < self.orb_ini[0] + 1):
                        self.orb_ini[0] = int(field[3]) - 1
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
                        assert (int(field[0]) == ndim_old)
                        assert (int(field[1]) >= molecule.nst - 1)
                        # nxply is number of lines for each excited state in 'XplusY.dat'
                        nxply = int(ndim_old / 6) + 1
                        if (ndim_old % 6 != 0):
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
                    assert (int(field[0]) == ndim)
                    assert (int(field[1]) >= molecule.nst - 1)
                    # nxply is number of lines for each excited state in 'XplusY.dat'
                    nxply = int(ndim / 6) + 1
                    if (ndim % 6 != 0):
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


