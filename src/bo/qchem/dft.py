from __future__ import division
from bo.qchem.qchem import QChem
from misc import au_to_A, eV_to_au
import os, shutil, re, textwrap, subprocess
import numpy as np

class DFT(QChem):
    """ Class for DFT method of QChem5.2 program
    """
    def __init__(self, molecule, basis_set="sto-3g", memory="500m", \
        functional="b3lyp", scf_max_iter=20, \
        qm_path="/opt/qchem", nthreads=1, version=5.2):
        # Initialize Molpro common variables
        super(DFT, self).__init__(basis_set, memory, qm_path, nthreads, version)

        self.functional = functional
        self.scf_max_iter = scf_max_iter
        self.nthreads = nthreads

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFT method
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
       # self.extract_BO(molecule, bo_list, calc_force_only)
       # self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, calc_force_only):
        """ Generate QChem input files: qchem.in
        """

        # Make QC5.2 input file
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
        input_molecule += "$end\n\n"
        input_qc += input_molecule

        # Job control
        input_rem = textwrap.dedent(f"""\
        $rem
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

        for ist in range(molecule.nst + 1):
            input_rem += f"{ist}  "
        input_rem += "\n$end\n\n"

        input_qc += input_rem
  
        file_name = "qchem.in"
        with open(file_name, "w") as f:
            f.write(input_qc)
        
    def run_QM(self, base_dir, istep, bo_list):
        """ Run (TD)DFT calculation and save the output files to QMlog directory
        """
        # Set environment variable 
        os.environ["QC"] = self.qm_path
        os.system(f"echo $QC")
        os.system(f"source $QC/qcenv.sh")
        os.system(f"echo $QCSCRATCH")
        os.environ["QCLOCALSCR"] = self.scr_qm_dir

        qm_exec_command = f"$QC/bin/qchem -mpi -np {self.nthreads} qchem.in > log"
        
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

        # run QChem
        os.system(qm_exec_command)

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information
        """
