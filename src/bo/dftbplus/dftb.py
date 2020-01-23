from __future__ import division
import os, shutil, re, struct
from bo.dftbplus.dftbplus import DFTBplus
import numpy as np
import textwrap
from misc import eV_to_au

spin_w = {"H":"-0.072", "C":"-0.031 -0.025 -0.025 -0.023", "N":"-0.033 -0.027 -0.027 -0.026", \
    "O":"-0.035 -0.030 -0.030 -0.028"}

max_l = {"H":"s", "C":"p", "N":"p", "O":"p"}

class DFTB(DFTBplus):
    """ Class for (TD)DFTB method of DFTB+ program

        :param object molecule: molecule object
        :param boolean scc: include SCC scheme
        :param double scc_tol: energy convergence for SCC iterations
        :param integer max_scc_iter: maximum number of SCC iterations
        :param boolean sdftb: include spin-polarisation scheme
        :param double unpaired_e: number of unpaired electrons
        :param double e_temp: electronic temperature for Fermi-Dirac scheme
        :param string mixer: charge mixing method used in SCC-DFTB
        :param string ex_symmetry: symmetry (singlet or triplet) in TD-DFTB
        :param string sk_path: path for slater-koster files
        :param boolean periodic: use periodicity in the calculations
        :param double a(b, c)_axis: the length of cell lattice
        :param string qm_path: path for QM binary
        :param integer nthreads: number of threads in the calculations
        :param boolean mpi: use MPI parallelization
        :param string mpi_path: path for MPI binary
        :param double version: version of DFTB+ program
    """
    def __init__(self, molecule, scc=True, scc_tol=1E-6, max_scc_iter=100, \
        sdftb=False, unpaired_e=0., e_temp=0., mixer="Broyden", \
        ex_symmetry="S", sk_path="./", periodic=False, a_axis=0., b_axis=0., c_axis=0., \
        qm_path="./", nthreads=1, mpi=False, mpi_path="./", version=19.1):
        # Initialize DFTBplus common variables
        super().__init__(molecule, sk_path, qm_path, nthreads, version)

        # Initialize DFTBplus DFTB variables
        self.scc = scc
        self.scc_tol = scc_tol
        self.max_scc_iter = max_scc_iter

        self.sdftb = sdftb
        self.unpaired_e = unpaired_e

        self.e_temp = e_temp
        self.mixer = mixer

        self.ex_symmetry = ex_symmetry

        self.mpi = mpi
        self.mpi_path = mpi_path

        self.periodic = periodic
        self.a_axis = a_axis
        self.b_axis = b_axis
        self.c_axis = c_axis

        # set 'l_nacme' with respect to the computational method
        # TD-DFTB do not produce NACs, so we should get NACME from CIoverlap
        molecule.l_nacme = True

        # re-calculation of excited state forces is not needed for ground state dynamics
        if (molecule.nst > 1):
            self.re_calc = True
        else:
            self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFTB method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        self.copy_files(molecule, istep, calc_force_only)
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, istep, bo_list, calc_force_only)
        self.run_QM(molecule, base_dir, istep, bo_list, calc_force_only)
        self.extract_BO(molecule, base_dir, istep, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def copy_files(self, molecule, istep, calc_force_only):
        """ Copy necessary scratch files in previous step

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):
            # after T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "geometry.xyz"), \
                os.path.join(self.scr_qm_dir, "../geometry.xyz.pre"))
            shutil.copy(os.path.join(self.scr_qm_dir, "eigenvec.bin"), \
                os.path.join(self.scr_qm_dir, "../eigenvec.bin.pre"))
            shutil.copy(os.path.join(self.scr_qm_dir, "XplusY.DAT"), \
                os.path.join(self.scr_qm_dir, "../XplusY.DAT.pre"))

    def get_input(self, molecule, istep, bo_list, calc_force_only):
        """ Generate DFTB+ input files: geometry.gen, dftb_in.hsd

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # TODO : currently, CIoverlap is not correct -> only BOMD possible with TD-DFTB
        if (self.calc_coupling):
            raise ValueError ("only BOMD possible with TD-DFTB")

        # make 'geometry.gen' file
        os.system("xyz2gen geometry.xyz")
        if (self.periodic):
            # substitute C to S in first line
            file_be = open('geometry.gen', 'r')
            file_af = open('tmp.gen', 'w')
            first_row = True
            for row in file_be:
                if (first_row):
                    row = f'{molecule.nat} S\n'
                    first_row = False
                file_af.write(row)
            # add gamma-point and cell lattice information
            geom_periodic = textwrap.dedent(f"""\
            {0.0:15.8f} {0.0:15.8f} {0.0:15.8f}
            {self.a_axis:15.8f} {0.0:15.8f} {0.0:15.8f}
            {0.0:15.8f} {self.b_axis:15.8f} {0.0:15.8f}
            {0.0:15.8f} {0.0:15.8f} {self.c_axis:15.8f}
            """)
            file_af.write(geom_periodic)
            file_be.close()
            file_af.close()
            os.rename('tmp.gen', 'geometry.gen')

        # make 'double.gen' file for overlap in TD-DFTB
        # In this case, we do not need to consider periodicity
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):
            # move previous files to currect directory
            os.rename('../geometry.xyz.pre', './geometry.xyz.pre')
            os.rename('../eigenvec.bin.pre', './eigenvec.bin.pre')
            os.rename('../XplusY.DAT.pre', './XplusY.DAT.pre')
            # open geometry.xyz.pre
            file_af = open('double.xyz', 'w')
            file_be = open('geometry.xyz.pre', 'r')
            first_row = True
            for row in file_be:
                if (first_row):
                    row = f'{molecule.nat * 2}\n'
                    first_row = False
                file_af.write(row)
            file_be.close()
            # open geometry.xyz
            file_be = open('geometry.xyz', 'r')
            iline = 1
            for row in file_be:
                if (iline > 2):
                    file_af.write(row)
                iline += 1
            file_be.close()
            file_af.close()
            os.system("xyz2gen double.xyz")

        # make 'dftb_in.hsd' file
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
              MaxSCCIterations = {self.max_scc_iter}
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
            if (self.sdftb or self.ex_symmetry == "T"):
                spin_constant = ("\n" + " " * 18).join([f"  {itype} = {{ {spin_w[f'{itype}']} }}" for itype in self.atom_type])
                input_ham_spin_w = textwrap.indent(textwrap.dedent(f"""\
                  SpinConstants = {{
                    ShellResolvedSpin = Yes
                  {spin_constant}
                  }}
                """), "  ")
                input_dftb += input_ham_spin_w
            # TODO : read information from previous step
#            if (calc_force_only):
#                input_ham_restart = textwrap.indent(textwrap.dedent(f"""\
#                  ReadInitialCharges = Yes
#                """), "  ")
#                input_dftb += input_ham_restart
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
        if (molecule.nst > 1):

            # calculate excited state force and root state?
            if (bo_list[0] > 0):
                ex_force = "Yes"
                rst = bo_list[0]
            else:
                ex_force = "No"
                rst = bo_list[0] + 1

            # set symmetry in TD-DFTB
            if (self.ex_symmetry == "S"):
                symmetry = "singlet"
            elif (self.ex_symmetry == "T"):
                symmetry = "triplet"
            else:
                raise ValueError (f"wrong input given in {self.ex_symmetry}")

            # set number of excitations in TD-DFTB
            if (molecule.nat <= 5):
                num_ex = molecule.nst + 2
            elif (molecule.nat > 5 and molecule.nat <= 15):
                num_ex = 2 * molecule.nst + 2
            else:
                num_ex = 3 * molecule.nst + 2

            # write XplusY data?
            if (self.calc_coupling):
                xpy = "Yes"
            else:
                xpy = "No"

            input_excited = textwrap.dedent(f"""\
            ExcitedState = Casida{{
              NrOfExcitations = {num_ex}
              StateOfInterest = {rst}
              Symmetry = {symmetry}
              WriteTransitions = Yes
              WriteMulliken = Yes
              WriteXplusY = {xpy}
              ExcitedStateForces = {ex_force}
            }}
            """)
            input_dftb += input_excited

        # ParserOptions Block
        if (self.version == 19.1):
            parser_version = 7
        else:
            raise ValueError ("Other Versions Not Implemented")
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

        # write 'dftb_in.hsd.geom' file
        file_name = "dftb_in.hsd.geom"
        with open(file_name, "w") as f:
            f.write(input_dftb)

        # write 'dftb_in.hsd.double' file
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):

            # new input for dftb
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
        # set run command
        qm_command = os.path.join(self.qm_path, "dftb+")
        if (self.mpi):
            # mpi setting
            os.environ["OMP_NUM_THREADS"] = "1"
            mpi_command = os.path.join(self.mpi_path, "mpirun")
            command = f"{mpi_command} -np {self.nthreads} {qm_command} > log"
        else:
            # openmp setting
            os.environ["OMP_NUM_THREADS"] = f"{self.nthreads}"
            command = f"{qm_command} > log"

        # run DFTBplus for calculation of overlap matrix
        if (self.calc_coupling and not calc_force_only and istep >= 0 and molecule.nst > 1):
            shutil.copy("dftb_in.hsd.double", "dftb_in.hsd")
            os.system(command)

        # run DFTBplus method for molecular dynamics
        shutil.copy("dftb_in.hsd.geom", "dftb_in.hsd")
        os.system(command)

        # copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # read 'detailed.out' file
        # TODO: the qmmm information is written in this file
        file_name = "detailed.out"
        with open(file_name, "r") as f:
            detailed_out = f.read()
        # read 'EXC.DAT' file
        if (molecule.nst > 1):
            file_name = "EXC.DAT"
            with open(file_name, "r") as f:
                exc_out = f.read()

        # energy
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

        # force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        tmp_f = 'Total Forces' + '\n\s+\d*\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
        force = re.findall(tmp_f, detailed_out)
        force = np.array(force[0])
        force = force.astype(float)
        force = force.reshape(molecule.nat, 3, order='C')
        molecule.states[bo_list[0]].force = np.copy(force)

        # nacme
        if (self.calc_coupling and not calc_force_only):
            molecule.nacme = np.zeros((molecule.nst, molecule.nst))
            if (istep >= 0):
                self.CIoverlap(molecule, base_dir)
                # read 'NACME.DAT'
                ist = 0
                jst = 0
                nline = 1
                file_name_in = "NACME.DAT"
                with open(file_name_in, "r") as f_in:
                    lines = f_in.readlines()
                    for line in lines:
                        field = line.split()
                        if (nline % 2 == 0):
                            # TODO : current TDNAC gives too large values
                            #molecule.nacme[ist, jst] = field[0]
                            jst += 1
                            if (jst == molecule.nst):
                                ist += 1
                                jst = 0
                        nline += 1

    def CIoverlap(self, molecule, base_dir):
        """ Read the necessary files and generate NACME file and
            this is an experimental feature and not used

            :param object molecule: molecule object
            :param string base_dir: base directory
        """
        # set new variable to decide the number of basis functions for atoms
        check_atom = [0]
        num_basis = 0
        core_elec = 0.
        for iat in range(molecule.nat):
            max_ang = max_l[molecule.symbols[iat]]
            if (max_ang == 'p'):
                num_basis += 4
                core_elec += 2.
                check_atom.append(num_basis)
            elif (max_ang == 's'):
                num_basis += 1
                check_atom.append(num_basis)

        # set new variable to decide the position of atoms in basis functions
        check_basis = []
        for ibasis in range(num_basis):
            for iat in range(molecule.nat):
                ind_a = check_atom[iat] + 1
                ind_b = check_atom[iat + 1]
                if (ibasis + 1 >= ind_a and ibasis + 1 <= ind_b):
                    check_basis.append(iat + 1)

        # write 'INPUT' file
        ncore = 0
        nocc = int(int(molecule.nelec - core_elec) / 2) - ncore
        nvirt = num_basis - nocc - ncore

        file_name_out = "INPUT"
        f_out = open(file_name_out, "w")

        f_print = f"{num_basis:5d} {ncore:4d} {nocc:4d} {nvirt:4d} {molecule.nst:3d} 5.16767" + "\n"
        f_out.write(f_print)

        f_out.close()

        # write 'AOVERLAP' file
        file_name_out = "AOVERLAP"
        f_out = open(file_name_out, "w")

        file_name_in = "oversqr.dat"
#        over_mat = np.loadtxt(file_name_in, skiprows=5, dtype=np.float)
#        nan_ind = np.argwhere(np.isnan(over_mat))
#        for row, col in nan_ind:
#            ind_a = check_basis[row%num_basis]
#            ind_b = check_basis[col%num_basis]
#            if (ind_a == ind_b):
#                if (row%num_basis == col%num_basis):
#                    over_mat[row, col] = 1.
#                else:
#                    over_mat[row, col] = 0.
#        #np.savetxt("test", over_mat, fmt=f"%24.15e")
#        np.savetxt("test", over_mat, fmt=f"%23.15e")
        with open(file_name_in, "r") as f_in:
            lines = f_in.readlines()
            nline = 1
            nrow = 1
            for line in lines:
                if (nline >= 6):
                    field = line.split()
                    if ("NaN" in field):
                        ncolumn = 1
                        for element in field:
                            if (element == 'NaN'):
                                ind_a = check_basis[nrow - 1]
                                ind_b = check_basis[ncolumn - 1]
                                if (ind_a == ind_b):
                                    if (nrow == ncolumn):
                                        f_print = f"{1.0:24.15e}"
                                    else:
                                        f_print = f"{0.0:24.15e}"
                                    f_out.write(f_print)
                            else:
                                f_print = f"{float(element):24.15e}"
                                f_out.write(f_print)
                            if (field.index(element) == 2 * num_basis - 1):
                                f_out.write("\n")
                            ncolumn += 1
                            if (ncolumn > num_basis):
                                ncolumn -= num_basis
                    else:
                        f_print = f"{line}"
                        f_out.write(f_print)
                    nrow += 1
                    if (nrow > num_basis):
                        nrow -= num_basis
                nline += 1

        f_out.close()

        # write 'MOCOEF' file
        file_name_out = "MOCOEF"
        f_out = open(file_name_out, "w")

        file_name_in = "eigenvec.bin"
        mocoef = []
        with open(file_name_in, "rb") as f_in:
            dummy = np.fromfile(f_in, dtype=np.integer, count = 1)
            for ibasis in range(num_basis):
                dummy = np.fromfile(f_in, dtype=np.integer, count = 1)
                data = np.fromfile(f_in, dtype=np.float64, count = num_basis)
                mocoef.append(data)
            mocoef = np.array(mocoef)
            mocoef = np.transpose(mocoef)
            mocoef = [val for sublist in mocoef for val in sublist]
            for ibasis in range(num_basis):
                ind_a = num_basis * ibasis
                ind_b = num_basis * (ibasis + 1)
                f_print = " ".join([f"{mocoef[ind]:13.8f}" for ind in range(ind_a, ind_b)]) + "\n"
                f_out.write(f_print)

        f_out.close()

        # write 'MOCOEFOLD' file
        file_name_out = "MOCOEFOLD"
        f_out = open(file_name_out, "w")

        file_name_in = "eigenvec.bin.pre"
        mocoef = []
        with open(file_name_in, "rb") as f_in:
            dummy = np.fromfile(f_in, dtype=np.integer, count = 1)
            for ibasis in range(num_basis):
                dummy = np.fromfile(f_in, dtype=np.integer, count = 1)
                data = np.fromfile(f_in, dtype=np.float64, count = num_basis)
                mocoef.append(data)
            mocoef = np.array(mocoef)
            mocoef = np.transpose(mocoef)
            mocoef = [val for sublist in mocoef for val in sublist]
            for ibasis in range(num_basis):
                ind_a = num_basis * ibasis
                ind_b = num_basis * (ibasis + 1)
                f_print = " ".join([f"{mocoef[ind]:13.8f}" for ind in range(ind_a, ind_b)]) + "\n"
                f_out.write(f_print)

        f_out.close()

        # write 'CICOEF' file
        file_name_out = "CICOEF"
        f_out = open(file_name_out, "w")

        file_name_in = "XplusY.DAT"
        with open(file_name_in, "r") as f_in:
            lines = f_in.readlines()
            nline = 1
            ind_b = -1
            for line in lines:
                if (nline == 1):
                    field = line.split()
                    nmat = int(field[0])
                    nexc = int(field[1])
                    nstd = int(nmat / 6) + 1
                    if (nmat % 6 != 0):
                        nstd += 1
                    xply = np.zeros((nocc, nvirt, nexc))
                else:
                    if ((nline - 1) % nstd == 1):
                        ind_occ = nocc - 1
                        ind_virt = 0
                        ind_b += 1
                    else:
                        field = line.split()
                        for element in field:
                            xply[ind_occ, ind_virt, ind_b] = float(element)
                            if (ind_occ == 0):
                                ind_occ = nocc - 1
                                ind_virt += 1
                            else:
                                ind_occ -= 1
                nline += 1

        for iexc in range(nexc):

            # normalize the CI coefficients
            norm_val = 0.
            for iocc in range(nocc):
                for ivirt in range(nvirt):
                    norm_val += xply[iocc, ivirt, iexc] ** 2
            norm_val = np.sqrt(norm_val)
            for iocc in range(nocc):
                for ivirt in range(nvirt):
                    xply[iocc, ivirt, iexc] /= norm_val

            f_print = f"{iexc + 1:4d}" + "\n"
            f_out.write(f_print)

            for iocc in range(nocc):
                for ivirt in range(nvirt):
                    f_print = f"{xply[iocc, ivirt, iexc]:13.8f}"
                    if (ivirt == nvirt - 1):
                        f_print += "\n"
                    f_out.write(f_print)

        f_out.close()

        # write 'CICOEFOLD' file
        file_name_out = "CICOEFOLD"
        f_out = open(file_name_out, "w")

        file_name_in = "XplusY.DAT.pre"
        with open(file_name_in, "r") as f_in:
            lines = f_in.readlines()
            nline = 1
            ind_b = -1
            for line in lines:
                if (nline == 1):
                    field = line.split()
                    nmat = int(field[0])
                    nexc = int(field[1])
                    nstd = int(nmat / 6) + 1
                    if (nmat % 6 != 0):
                        nstd += 1
                    xply = np.zeros((nocc, nvirt, nexc))
                else:
                    if ((nline - 1) % nstd == 1):
                        ind_occ = nocc - 1
                        ind_virt = 0
                        ind_b += 1
                    else:
                        field = line.split()
                        for element in field:
                            xply[ind_occ, ind_virt, ind_b] = float(element)
                            if (ind_occ == 0):
                                ind_occ = nocc - 1
                                ind_virt += 1
                            else:
                                ind_occ -= 1
                nline += 1

        for iexc in range(nexc):

            # normalize the CI coefficients
            norm_val = 0.
            for iocc in range(nocc):
                for ivirt in range(nvirt):
                    norm_val += xply[iocc, ivirt, iexc] ** 2
            norm_val = np.sqrt(norm_val)
            for iocc in range(nocc):
                for ivirt in range(nvirt):
                    xply[iocc, ivirt, iexc] /= norm_val

            f_print = f"{iexc + 1:4d}" + "\n"
            f_out.write(f_print)

            for iocc in range(nocc):
                for ivirt in range(nvirt):
                    f_print = f"{xply[iocc, ivirt, iexc]:13.8f}"
                    if (ivirt == nvirt - 1):
                        f_print += "\n"
                    f_out.write(f_print)

        f_out.close()

        # TODO: this is temporary path, the directory for tdnac.x can be changed
#        tdnac_command = os.path.join(base_dir, "../tdnac.x")
#        command = f"{tdnac_command}"
#        os.system(command)


