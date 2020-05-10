from __future__ import division
from bo.qchem.qchem import QChem
from misc import au_to_A, eV_to_au
import os, shutil, re, textwrap, subprocess
import numpy as np

class DFT(QChem):
    """ Class for DFT method of QChem5.2 program
    """
    def __init__(self, molecule, basis_set="sto-3g", memory="500m", \
        functional="b3lyp", scf_max_iter=50, \
        qm_path="/opt/qchem", nthreads=1, version=5.2):
        # Initialize QChem common variables
        super(DFT, self).__init__(basis_set, memory, qm_path, nthreads, version)

        self.functional = functional
        self.scf_max_iter = scf_max_iter
        self.nthreads = nthreads

    def get_bo(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFT method
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, calc_force_only):
        """ Generate QChem input files: qchem.in
        """

        # Make QChem5.2 input file
        input_qc = ""
          
        # Molecular information such as charge, geometry
        input_molecule = textwrap.dedent(f"""\
        $molecule
        {int(molecule.charge)}  1
        """)
        
        for iat in range(molecule.nat):
            list_pos = list(molecule.pos[iat])
            input_molecule += \
                f"{molecule.symbols[iat]}{list_pos[0]:15.8f}{list_pos[1]:15.8f}{list_pos[2]:15.8f}\n"
        input_molecule += "$end\n\n$rem\n"
        input_qc += input_molecule

        if (not calc_force_only):
            # Job control to calculate NAC
            input_nac = textwrap.dedent(f"""\
            JOBTYPE  SP
            METHOD {self.functional}
            BASIS {self.basis_set}
            SCF_MAX_CYCLES {self.scf_max_iter}
            CIS_N_ROOTS  {molecule.nst}
            CIS_TRIPLETS FALSE
            CALC_NAC TRUE
            CIS_DER_NUMSTATE {molecule.nst}
            $end

            $derivative_coupling
            running state is {bo_list[0]}
            """)

            for ist in range(molecule.nst):
                input_nac += f"{ist}  "
            input_nac += "\n$end\n\n"

            input_nac += textwrap.dedent(f"""\
            @@@

            $molecule
            read
            $end

            $rem
            SCF_GUESS read
            """)

            input_qc += input_nac
      
        # Job control to calculate force
        input_force = textwrap.dedent(f"""\
        JOBTYPE force
        METHOD {self.functional}
        BASIS {self.basis_set}
        CIS_N_ROOTS  {molecule.nst}
        CIS_TRIPLETS FALSE
        CIS_STATE_DERIV {bo_list[0]}
        $end
        """)

        input_qc += input_force

        file_name = "qchem.in"
        with open(file_name, "w") as f:
            f.write(input_qc)
        
    def run_QM(self, base_dir, istep, bo_list):
        """ Run (TD)DFT calculation and save the output files to QMlog directory
        """
        # Set environment variable 
        os.environ["QC"] = self.qm_path
        path_qcenv = os.path.join(self.qm_path, "qcenv.sh")
        command = f'env -i sh -c "source {path_qcenv} && env"'
        for line in subprocess.getoutput(command).split("\n"):
            key, value = line.split("=")
            os.environ[key] = value
        os.environ["QCLOCALSCR"] = self.scr_qm_dir

        qm_exec_command = f"$QC/bin/qchem -mpi -np {self.nthreads} qchem.in > log"
        
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

        # Run QChem
        os.system(qm_exec_command)

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information
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

        # QChem provides gradient, not force
        force = -force

        nline = 0
        for iiter in range(num_line):
            tmp_force = np.transpose(force[0][18 * iiter:18 * (iiter + 1)].reshape(3, 6, order="C"))
            for iat in range(6):
                states.force[6 * nline + iat] = tmp_force[iat]
            nline += 1
            if (iiter == num_line -1):
                tmp_force = np.transpose(force[0][18 * (iiter + 1):].reshape(3, dnum, order="C"))
                for iat in range(dnum):
                    states.force[6 * nline + iat] = tmp_force[iat]

        if (not calc_force_only):
            # Non-adiabatic coupling vector
            tmp_nac = "with ETF[:]*\s*Atom\s*X\s*Y\s*Z\s*[-]*" + ("\s*\d*\s*" + "([-]*\S*)\s*"*3) * molecule.nat
            nac = re.findall(tmp_nac, log)
            nac = np.array(nac)
            nac = nac.astype(float)

            num = 0
            molecule.nac[:, :, :, :] = 0.
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    molecule.nac[ist, jst] = nac[num].reshape(molecule.nat, 3, order='C')
                    molecule.nac[jst, ist] = -molecule.nac[ist, jst]
                    num += 1
