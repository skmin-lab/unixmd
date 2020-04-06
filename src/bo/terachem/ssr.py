from __future__ import division
from bo.terachem.terachem import TeraChem
import os, shutil, re, textwrap
import numpy as np

class SSR(TeraChem):
    """ Class for SSR method of TeraChem program

        :param object molecule: molecule object
        :param string basis_set: basis set information
        :param string functional: level of DFT theory
        :param string precision: precision in the calculations
        :param double scf_tol: energy convergence for SCF iterations
        :param integer max_scf_iter: maximum number of SCF iterations
        :param string reks22: use REKS(2,2) calculation?
        :param double reks_scf_tol: energy convergence for REKS SCF iterations
        :param integer reks_max_scf_iter: maximum number of REKS SCF iterations
        :param string reks_diis: DIIS acceleration in REKS SCF iterations
        :param double shift: level shifting value in REKS SCF iterations
        :param integer use_ssr_state: calculate SSR state, if not, treat SA-REKS
        :param double cpreks_max_tol: gradient tolerance for CP-REKS equations
        :param integer cpreks_max_iter: maximum number of CP-REKS iterations
        :param string qm_path: path for QM binary
        :param integer ngpus: number of GPUs
        :param string gpu_id: ID of used GPUs
        :param double version: version of TeraChem program
    """
    def __init__(self, molecule, ngpus=1, gpu_id="1", precision="dynamic", \
        version=1.92, functional="hf", basis_set="sto-3g", scf_tol=1E-2, \
        max_scf_iter=300, reks22="yes", reks_scf_tol=1E-6, \
        reks_max_scf_iter=1000, reks_diis="yes", shift=0.3, use_ssr_state=1, \
        cpreks_max_tol=1E-6, cpreks_max_iter=1000, qm_path="./"):
        # Initialize TeraChem common variables
        super(SSR, self).__init__(functional, basis_set, qm_path, ngpus, \
            gpu_id, precision, version)

        # Initialize TeraChem SSR variables
        self.scf_tol = scf_tol
        self.max_scf_iter = max_scf_iter

        self.reks22 = reks22
        if (self.reks22 == "yes"):
            if (molecule.nst > 2):
                raise ValueError ("3state REKS with gradient Not Implemented")
            self.reks_scf_tol = reks_scf_tol
            self.reks_max_scf_iter = reks_max_scf_iter
            self.reks_diis = reks_diis
            self.shift = shift
            self.use_ssr_state = use_ssr_state
            if (molecule.nst > 1):
                self.cpreks_max_tol = cpreks_max_tol
                self.cpreks_max_iter = cpreks_max_iter
        else:
            raise ValueError("reks22 should be switched on in our interface")

        # Set 'l_nacme' with respect to the computational method
        # SSR can produce NACs, so we do not need to get NACME from CIoverlap
        # When we calculate SA-REKS state, NACME can be obtained directly from diabetic Hamiltonian
        if (self.use_ssr_state == 1):
            molecule.l_nacme = False
        else:
            # TODO : diabatic SH?
            molecule.l_nacme = True

        if (molecule.nst > 1 and self.use_ssr_state == 0):
            # SA-REKS state with SH
            # TODO : diabatic SH?
            self.re_calc = True
        else:
            # SSR state or single-state REKS
            self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from SSR method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, calc_force_only):
        """ Generate TeraChem input files: input.tcin

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Make 'input.tcin' file
        input_terachem = ""

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

        # DFT Block
        input_dft = textwrap.dedent(f"""\
        method {self.functional}
        basis {self.basis_set}
        convthre {self.scf_tol}
        maxit {self.max_scf_iter}
        charge {molecule.charge}

        """)
        input_terachem += input_dft

        # Options for SCF initial guess
        # TODO : read information from previous step
#        if (calc_force_only):
#            guess = 1
#            input_dft_guess = textwrap.dedent(f"""\
#            guess scr/c0
#            """)
#            input_terachem += input_dft_guess

        # REKS Block
        if (self.reks22 == "yes"):

            # Energy functional options
            if (molecule.nst == 1):
                sa_reks = 0
            elif (molecule.nst == 2):
                if (self.use_ssr_state == 1):
                    sa_reks = 2
                else:
                    sa_reks = 1

            # NAC calculation options
            if (self.calc_coupling and sa_reks == 2):
                # SSR state with SH, Eh
                reks_target = 12
                self.nac = "Yes"
            else:
                # Any method with BOMD / SA-REKS state with SH
                # TODO : diabatic SH?
                reks_target = bo_list[0] + 1
                self.nac = "No"

            # TODO: pointcharges? in qmmm?

            # Options for REKS SCF initial guess
            # TODO : read information from previous step
#            if (restart):
#                reks_guess = 1
#            else:
#                reks_guess = 0
            reks_guess = 0

            input_reks_basic = textwrap.dedent(f"""\
            reks22 {self.reks22}
            reks_convthre {self.reks_scf_tol}
            reks_maxit {self.reks_max_scf_iter}
            reks_diis {self.reks_diis}
            reks_shift {self.shift}
            sa_reks {sa_reks}
            reks_target {reks_target}
            reks_guess {reks_guess}
            """)
            input_terachem += input_reks_basic
            if (molecule.nst > 1):
                input_cpreks = textwrap.dedent(f"""\
                cpreks_thresh {self.cpreks_max_tol}
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

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Read 'log' file
        file_name = "log"
        with open(file_name, "r") as f:
            log_out = f.read()

        # Energy
        if (not calc_force_only):
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
                if (self.use_ssr_state == 1):
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
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        if (self.nac == "Yes"):
            # SSR state with SH, Eh
            for ist in range(molecule.nst):
                tmp_f = f'Eigen state {ist + 1} gradient\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
                    '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                force = re.findall(tmp_f, log_out)
                force = np.array(force[0])
                force = force.astype(float)
                force = force.reshape(molecule.nat, 3, order='C')
                molecule.states[ist].force = - np.copy(force)
        else:
            # SH : SA-REKS state
            # TODO : diabatic SH?
            # BOMD : SSR state, SA-REKS state or single-state REKS
            tmp_f = 'Gradient units are Hartree/Bohr\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	              '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = - np.copy(force)

        # NAC
        # TODO : NACME - diabatic SH?
        if (not calc_force_only and self.nac == "Yes"):

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
                        molecule.nac[ist, jst, :, :] = 0.
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
                            molecule.nac[ist, jst] = nac
                        kst += 1
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]



