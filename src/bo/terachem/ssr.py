from __future__ import division
import os, shutil, re
from bo.terachem.terachem import TeraChem
import numpy as np
import textwrap

class SSR(TeraChem):
    """ Class for TeraChem program
                            bomd | sh | eh | nac | re_calc
        single-state REKS :  o     x    x     F      F
        SA-REKS           :  o     o    x     F      T
        SI-SA-REKS(SSR)   :  o     o    o     T      F
    """
    def __init__(self, molecule, ngpus=1, gpu_id="1", precision="dynamic", \
        version=1.92, functional="hf", basis_set="sto-3g", scf_tol=1E-2, \
	      max_scf_iter=300, reks22="yes", reks_scf_tol=1E-6, \
	      reks_max_scf_iter=1000, reks_diis="yes", shift=0.3, use_ssr_state=1, \
        cpreks_max_tol=1E-6, cpreks_max_iter=1000, qm_path="./"):
        # Initialize TeraChem common variables
        super().__init__(functional, basis_set, qm_path, ngpus, \
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

        # set 'l_nacme' with respect to the computational method
        # SSR can produce NACs, so we do not need to get NACME from CIoverlap
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
        """ Get/Extract BO information from DFTB
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, calc_force_only):
        """ Generate TeraChem input files: input.tcin
        """
        # make 'input.tcin' file
        input_terachem = ""

        # control Block
        input_control = \
        input_control = textwrap.dedent(f"""\

        run gradient

        coordinates geometry.xyz

        """)
        input_terachem += input_control

        # system Block
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

        # options for SCF initial guess
        # TODO : read information from previous step
#        if (calc_force_only):
#            guess = 1
#            input_dft_guess = textwrap.dedent(f"""\
#            guess scr/c0
#            """)
#            input_terachem += input_dft_guess

        # REKS Block
        if (self.reks22 == "yes"):

            # energy functional options
            if (molecule.nst == 1):
                sa_reks = 0
            elif (molecule.nst == 2):
                if (self.use_ssr_state == 1):
                    sa_reks = 2
                else:
                    sa_reks = 1

            # NAC calculation options
            if (self.calc_coupling and sa_reks == 2):
                # SSR state with sh, eh
                reks_target = 12
                self.nac = "Yes"
            else:
                # any method with bomd / SA-REKS state with sh
                reks_target = bo_list[0] + 1
                self.nac = "No"

            # TODO: pointcharges? in qmmm?

            # options for REKS SCF initial guess
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

        # write 'input.tcin' file
        file_name = "input.tcin"
        with open(file_name, "w") as f:
            f.write(input_terachem)

    def run_QM(self, base_dir, istep, bo_list):
        """ run TeraChem calculation and save the output files
        """
        # run TeraChem method
        qm_command = os.path.join(self.qm_path, "terachem")
        # openmp setting
        os.environ["OMP_NUM_THREADS"] = "1"
        command = f"{qm_command} input.tcin > log"
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

        # energy
        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            if (molecule.nst == 1):
                # single-state REKS
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

        # force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        if (self.nac == "Yes"):
            # SSR state with sh, eh
            for ist in range(molecule.nst):
                tmp_f = f'Eigen state {ist + 1} gradient\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
                    '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                force = re.findall(tmp_f, log_out)
                force = np.array(force[0])
                force = force.astype(float)
                force = force.reshape(molecule.nat, 3, order='C')
                molecule.states[ist].force = - np.copy(force)
        else:
            # sh : SA-REKS state
            # bomd : SSR state, SA-REKS state or single-state REKS
            tmp_f = 'Gradient units are Hartree/Bohr\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	              '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = - np.copy(force)

        # NAC
        if (not calc_force_only and self.nac == "Yes"):

            # 1.92 version do not show H vector
            if (self.version == 1.92):
                # zeroing for G, h and H vectors
                Gvec = np.zeros((molecule.nat, molecule.nsp))
                hvec = np.zeros((molecule.nat, molecule.nsp))
                ssr_coef = np.zeros((molecule.nst, molecule.nst))
                Hvec = np.zeros((molecule.nat, molecule.nsp))
                # calculate G vector, G vector is difference gradient so minus sign is needed
                Gvec = - 0.5 * (molecule.states[0].force - molecule.states[1].force)
                # read h vector
                tmp_c = 'Coupling gradient\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	                  '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat
                hvec = re.findall(tmp_c, log_out)
                hvec = np.array(hvec[0])
                hvec = hvec.astype(float)
                hvec = hvec.reshape(molecule.nat, 3, order='C')
                # read coefficients of SSR state
                ssr_coef = re.findall('SSR state\s\d\s+[-]\S+\s+([-]*\S+)\s+([-]*\S+)', log_out)
                ssr_coef = np.array(ssr_coef)
                ssr_coef = ssr_coef.astype(float)
                # calculate H vector from G, h vector
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



