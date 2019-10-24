from __future__ import division
import os, shutil, re
from bo.dftbplus.dftbplus import DFTBplus
import numpy as np
import textwrap

spin_w = {"H":"-0.072", "C":"-0.031 -0.025 -0.025 -0.023", "N":"-0.033 -0.027 -0.027 -0.026", \
    "O":"-0.035 -0.030 -0.030 -0.028"}

onsite_uu = {"H":"0.00000", "C":"0.00000 0.04973 0.04973 -0.01203", "N":"0.00000 0.06816 0.06816 -0.00879", \
    "O":"0.00000 0.08672 0.08672 -0.00523"}
onsite_ud = {"H":"0.00000", "C":"0.00000 0.10512 0.10512 0.02643", "N":"0.00000 0.12770 0.12770 0.03246", \
    "O":"0.00000 0.14969 0.14969 0.03834"}

max_l = {"H":"s", "C":"p", "N":"p", "O":"p"}

class SSR(DFTBplus):
    """ Class for density-functional tight-binding (DFTB) / SI-SA-REKS(SSR) method
                            bomd | sh | eh | nac | re_calc
        single-state REKS :  o     x    x     F      F
        SA-REKS           :  o     o    x     F      T
        SI-SA-REKS(SSR)   :  o     o    o     T      F
    """
    def __init__(self, molecule, scc=True, scc_tol=1E-6, max_scc_iter=1000, \
        sdftb=True, lcdftb=True, lc_method="NB", ocdftb=False, \
        ssr22=True, use_ssr_state=1, state_l=0, guess=1, shift=0.3, tuning=1., \
        grad_level=1, grad_tol=1E-8, mem_level=2, sk_path="./", periodic=False, \
        a_axis=0., b_axis=0., c_axis=0., qm_path="./", nthreads=1, version=19.1):
        # Initialize DFTBplus common variables
        super().__init__(molecule, sk_path, qm_path, nthreads, version)

        # Initialize DFTBplus SSR variables
        self.scc = scc
        self.scc_tol = scc_tol
        self.max_scc_iter = max_scc_iter

        self.sdftb = sdftb

        self.lcdftb = lcdftb
        self.lc_method = lc_method

        self.ocdftb = ocdftb

        self.ssr22 = ssr22
        self.use_ssr_state = use_ssr_state
        self.state_l = state_l
        self.guess = guess
        self.shift = shift
        self.tuning = tuning
        self.grad_level = grad_level
        self.grad_tol = grad_tol
        self.mem_level = mem_level

        self.periodic = periodic
        self.a_axis = a_axis
        self.b_axis = b_axis
        self.c_axis = c_axis

        # set 'l_nacme' with respect to the computational method
        # DFTB/SSR can produce NACs, so we do not need to get NACME from CIoverlap
        # when we calculate SA-REKS state, NACME can be obtained directly from diabetic Hamiltonian
        if (self.use_ssr_state == 1):
            molecule.l_nacme = False
        else:
            molecule.l_nacme = True

        if (molecule.nst > 1 and self.use_ssr_state == 0):
            # SA-REKS state with sh
            self.re_calc = True
        else:
            # SSR state or single-state REKS
            self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Get/Extract BO information from DFTBplus
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list):
        """ Generate DFTBplus input files: geometry.gen, dftb_in.hsd
        """
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
            """), "  ")
            input_dftb += input_ham_scc
            if (self.sdftb):
                spin_constant = ("\n" + " " * 18).join([f"  {itype} = {{ {spin_w[f'{itype}']} }}" for itype in self.atom_type])
                input_ham_spin = textwrap.indent(textwrap.dedent(f"""\
                  SpinConstants = {{
                    ShellResolvedSpin = Yes
                  {spin_constant}
                  }}
                """), "  ")
                input_dftb += input_ham_spin
            if (self.lcdftb):
                if (self.lc_method == "MM"):
                    lc_method = "MatrixBased"
                elif (self.lc_method == "NB"):
                    lc_method = "NeighbourBased"
                else:
                    raise ValueError("Other LC screening Not Compatible with SSR")
                input_ham_lc = textwrap.indent(textwrap.dedent(f"""\
                  RangeSeparated = LC{{
                    Screening = {lc_method}{{}}
                  }}
                """), "  ")
                input_dftb += input_ham_lc
            if (self.ocdftb):
                onsite_const_uu = ("\n" + " " * 18).join([f"  {itype}uu = {{ {onsite_uu[f'{itype}']} }}" for itype in self.atom_type])
                onsite_const_ud = ("\n" + " " * 18).join([f"  {itype}ud = {{ {onsite_ud[f'{itype}']} }}" for itype in self.atom_type])
                input_ham_oc = textwrap.indent(textwrap.dedent(f"""\
                  OnsiteCorrection = {{
                  {onsite_const_uu}
                  {onsite_const_ud}
                  }}
                """), "  ")
                input_dftb += input_ham_oc
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
        if (not self.ssr22):
            raise ValueError ("Other active spaces Not Implemented")

        # energy functional options
        if (molecule.nst == 1):
            energy_functional = 1
            energy_level = 1
        elif (molecule.nst == 2):
            energy_functional = 2
            energy_level = 1
        else:
            energy_functional = 2
            energy_level = 2

        # NAC calculation options
        if (molecule.nst == 1 or self.use_ssr_state == 0):
            # single-state REKS or SA-REKS state
            self.nac = "No"
        else:
            # SSR state
            if (self.calc_coupling):
                # sh, eh need NAC calculations
                self.nac = "Yes"
            else:
                # bomd do not need NAC calculations
                self.nac = "No"

        # TODO: rd will be determined automatically
        # qm => do not use rd, qmmm => use rd + external pc
        rd = "No"

        # options for SCF optimization
        if (self.guess == 2):
            raise ValueError("read external guess not implemented")

        # TODO : read information from previous step
#        if (calc_force_only):
#            guess = 2
#        else:
#            guess = self.guess
        guess = self.guess

        input_reks = textwrap.dedent(f"""\
        REKS = SSR22{{
          EnergyFunctional = {energy_functional}
          EnergyLevel = {energy_level}
          useSSRstate = {self.use_ssr_state}
          TargetState = {bo_list[0] + 1}
          TargetStateL = {self.state_l}
          InitialGuess = {guess}
          FONmaxIter = 50
          shift = {self.shift}
          GradientLevel = {self.grad_level}
          CGmaxIter = 100
          GradientTolerance = {self.grad_tol}
          RelaxedDensity = {rd}
          NonAdiabaticCoupling = {self.nac}
          PrintLevel = 1
          MemoryLevel = {self.mem_level}
        }}
        """)
        input_dftb += input_reks

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

        # write 'dftb_in.hsd' file
        file_name = "dftb_in.hsd"
        with open(file_name, "w") as f:
            f.write(input_dftb)

    def run_QM(self, base_dir, istep, bo_list):
        """ run DFTBplus calculation and save the output files
        """
        # run DFTBplus method
        qm_command = os.path.join(self.qm_path, "dftb+")
        # openmp setting
        os.environ["OMP_NUM_THREADS"] = f"{self.nthreads}"
        command = f"{qm_command} > log"
        os.system(command)
        # copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ read the output files to get BO data
        """
        # read 'log' file
        file_name = "log"
        with open(file_name, "r") as f:
            log_out = f.read()
        # read 'detailed.out' file
        # TODO: the qmmm information is written in this file
#        file_name = "detailed.out"
#        with open(file_name, "r") as f:
#            detailed_out = f.read()

        # energy
        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            if (molecule.nst == 1):
                # single-state REKS
                tmp_e = 'Spin' + '\n\s+\w+\s+([-]\S+)(?:\s+\S+){3}' * molecule.nst
                energy = re.findall(tmp_e, log_out)
                energy = np.array(energy)
            else:
                if (self.use_ssr_state == 1):
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

        # force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        if (self.nac == "Yes"):
            # SSR state with sh, eh
            for ist in range(molecule.nst):
                tmp_f = f' {ist + 1} st state \(SSR\)' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                force = re.findall(tmp_f, log_out)
                force = np.array(force[0])
                force = force.astype(float)
                force = force.reshape(molecule.nat, 3, order='C')
                molecule.states[ist].force = - np.copy(force)
        else:
            # sh : SA-REKS state
            # bomd : SSR state, SA-REKS state or single-state REKS
            tmp_f = f' {bo_list[0] + 1} state \(\w+[-]*\w+\)' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = - np.copy(force)

        # NAC
        if (not calc_force_only and self.nac == "Yes"):
            kst = 0
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.nac[ist, jst, :, :] = 0.
                    elif (ist < jst):
                        tmp_c = 'non-adiabatic coupling' + '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                        nac = re.findall(tmp_c, log_out)
                        nac = np.array(nac[kst])
                        nac = nac.astype(float)
                        nac = nac.reshape(molecule.nat, 3, order='C')
                        molecule.nac[ist, jst] = np.copy(nac)
                        kst += 1
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]




