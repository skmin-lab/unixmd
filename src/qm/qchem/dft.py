from __future__ import division
from qm.qchem.qchem import QChem
from misc import au_to_A, eV_to_au
import os, shutil, re, textwrap, subprocess
import numpy as np

class DFT(QChem):
    """ Class for DFT method of QChem5.2 program

        :param object molecule: molecule object
        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param integer nthreads: number of threads in the calculation
        :param string functional: xc functional 
        :param integer scf_max_iter: maximum number of SCF iterations
        :param integer scf_rho_tol: density convergence for SCF iterations
        :param integer cis_max_iter: maximum number of CIS iterations
        :param integer cis_en_tol: energy convergence for CIS iterations
        :param integer cpscf_max_iter: maximum number of CP iterations
        :param integer cpscf_grad_tol: gradient convergence of CP iterations
        :param string qm_path: path for QChem
        :param double version: QChem version
    """
    def __init__(self, molecule, basis_set="sto-3g", memory="500m", nthreads=1, \
        functional="blyp", scf_max_iter=50, scf_rho_tol=5, cis_max_iter=30, cis_en_tol=6, \
        cpscf_max_iter=30, cpscf_grad_tol=6, qm_path="./", version=5.2):
        # Initialize QChem common variables
        super(DFT, self).__init__(basis_set, memory, qm_path, nthreads, version)

        self.functional = functional
        self.scf_max_iter = scf_max_iter
        self.nthreads = nthreads
        self.scf_rho_tol = scf_rho_tol
        self.cis_max_iter = cis_max_iter
        self.cis_en_tol = cis_en_tol
        self.cpscf_max_iter = cpscf_max_iter
        self.cpscf_grad_tol = cpscf_grad_tol

        # QChem can provide NACs
        molecule.l_nacme = False
        self.re_calc = True

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFT method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        super().get_data(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, calc_force_only):
        """ Generate QChem input files: qchem.in

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Make QChem input file
        input_qc = ""

        # Molecular information such as charge, geometry
        input_molecule = textwrap.dedent(f"""\
        $molecule
        {int(molecule.charge)}  1
        """)

        for iat in range(molecule.nat):
            input_molecule += f"{molecule.symbols[iat]}"\
                + "".join([f"{i:15.8f}" for i in molecule.pos[iat]]) + "\n"
        input_molecule += "$end\n\n"
        input_qc += input_molecule

        # Job control to calculate NAC
        if (not calc_force_only and self.calc_coupling):
            # Arguments about SCF, xc functional and basis set
            input_nac = textwrap.dedent(f"""\
            $rem
            JOBTYPE SP
            INPUT_BOHR TRUE
            METHOD {self.functional}
            BASIS {self.basis_set}
            SCF_CONVERGENCE {self.scf_rho_tol}
            SYMMETRY FALSE
            SYM_IGNORE TRUE
            """)

            # Arguments about TDDFT and NAC
            input_nac += textwrap.dedent(f"""\
            CIS_N_ROOTS {molecule.nst-1}
            CIS_TRIPLETS FALSE
            CIS_CONVERGENCE {self.cis_en_tol}
            MAX_CIS_CYCLES {self.cis_max_iter}
            CALC_NAC TRUE
            CIS_DER_NUMSTATE {molecule.nst}
            SET_ITER {self.cpscf_max_iter}
            SET_CONV {self.cpscf_grad_tol}
            $end

            $derivative_coupling
            This is comment line
            """)

            for ist in range(molecule.nst):
                input_nac += f"{ist}  "
            input_nac += "\n$end\n\n"
            input_qc += input_nac

        # Job control to calculate force
        input_force = ""

        # BOMD: calc_force_only = F, self_calc_coupling = F
        # In BOMD, read in molecule section, scf_guess and skip_scf are not valid
        guess = "SAD"; skip = "FALSE"
        for ist in bo_list:
            if (not calc_force_only and self.calc_coupling):
                guess = "READ"; skip = "TRUE"
                input_force = textwrap.dedent(f"""\
                @@@

                $molecule
                read
                $end

                """)

            input_force += textwrap.dedent(f"""\
            $rem
            JOBTYPE force
            INPUT_BOHR TRUE
            METHOD {self.functional}
            BASIS {self.basis_set}
            SCF_GUESS {guess}
            SYMMETRY FALSE
            SYM_IGNORE TRUE
            """)

            # When ground state force is calculated, QChem doesn't need CIS option.
            if (ist != 0):
                input_force += textwrap.dedent(f"""\
                CIS_N_ROOTS {molecule.nst-1}
                CIS_STATE_DERIV {ist}
                CIS_TRIPLETS FALSE
                CIS_CONVERGENCE {self.cis_en_tol}
                MAX_CIS_CYCLES {self.cis_max_iter}
                SET_ITER {self.cpscf_max_iter}
                SET_CONV {self.cpscf_grad_tol}
                """)

                # CIS solution isn't saved in scratch.
                if (not calc_force_only and self.calc_coupling):
                    input_force += textwrap.dedent(f"""\
                    CIS_GUESS_DISK TRUE
                    CIS_GUESS_DISK_TYPE 2
                    SKIP_CIS_RPA TRUE
                    """)
            input_force += "$end\n\n"

            input_qc += input_force

        file_name = "qchem.in"
        with open(file_name, "w") as f:
            f.write(input_qc)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run (TD)DFT calculation and save the output files to QMlog directory

            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Set environment variable 
        os.environ["QC"] = self.qm_path
        path_qcenv = os.path.join(self.qm_path, "qcenv.sh")
        command = f'env -i sh -c "source {path_qcenv} && env"'
        for line in subprocess.getoutput(command).split("\n"):
            key, value = line.split("=")
            os.environ[key] = value
        os.environ["QCSCRATCH"] = self.scr_qm_dir
        os.environ["QCLOCALSCR"] = self.scr_qm_dir

        #TODO: MPI binary
        qm_exec_command = f"$QC/bin/qchem -nt {self.nthreads} qchem.in log save > qcprog.info "

        # Run QChem
        os.system(qm_exec_command)

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
        file_name = "log"
        with open(file_name, "r") as f:
            log = f.read()

        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            # Ground state energy
            energy = re.findall('Total energy in the final basis set =\s*([-]*\S*)', log)
            energy = np.array(energy)
            energy = energy.astype(float)
            molecule.states[0].energy = energy[0]

            # Excited state energy
            if (molecule.nst > 1):
                energy = re.findall('Total energy for state\s*\d*.\s*([-]*\S*)', log)
                energy = np.array(energy)
                energy = energy.astype(float)

                for ist, en in enumerate(energy):
                    if ist < molecule.nst - 1:
                        molecule.states[ist + 1].energy = en

        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        # Adiabatic force 
        tmp_f = "Gradient of\D*\s*" 
        num_line = int(molecule.nat / 6)
        if (num_line >= 1):
            tmp_f += ("\s*\d*\s*\d*\s*\d*\s*\d*\s*\d*\s*\d*"
                 + ("\s*\d?\s*" + "([-]*\S*)\s*" * 6) * 3) * num_line

        dnum = molecule.nat % 6
        tmp_f += "\s*\d*" * dnum
        tmp_f += ("\s*\d?\s*" + "([-]*\S*)\s*" * dnum) * 3

        force = re.findall(tmp_f, log)
        force = np.array(force)
        force = force.astype(float)

        # QChem provides energy gradient not force
        force = -force

        for index, ist in enumerate(bo_list):
            iline = 0; iiter = 0
            for iiter in range(num_line):
                tmp_force = np.transpose(force[index][18 * iiter:18 * (iiter + 1)].reshape(3, 6, order="C"))
                for iat in range(6):
                    molecule.states[ist].force[6 * iline + iat] = np.copy(tmp_force[iat])
                iline += 1

            if (dnum != 0):
                if (num_line != 0):
                    tmp_force = np.transpose(force[index][18 * (iiter + 1):].reshape(3, dnum, order="C"))
                else:
                    tmp_force = np.transpose(force[index][0:].reshape(3, dnum, order="C"))

                for iat in range(dnum):
                    molecule.states[ist].force[6 * iline + iat] = np.copy(tmp_force[iat])

        # NACs
        if (not calc_force_only and self.calc_coupling):
            tmp_nac = "with ETF[:]*\s*Atom\s*X\s*Y\s*Z\s*[-]*" + ("\s*\d*\s*" + "([-]*\S*)\s*"*3) * molecule.nat
            nac = re.findall(tmp_nac, log)
            nac = np.array(nac)
            nac = nac.astype(float)

            num = 0
            molecule.nac = np.zeros((molecule.nst, molecule.nst, molecule.nat, molecule.nsp))
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    molecule.nac[ist, jst] = np.copy(nac[num].reshape(molecule.nat, 3, order='C'))
                    molecule.nac[jst, ist] = - molecule.nac[ist, jst]
                    num += 1
