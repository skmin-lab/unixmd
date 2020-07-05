from __future__ import division
from bo.turbomole.turbomole import Turbomole
from misc import call_name
import os, shutil, re, textwrap
import numpy as np

class DFT(Turbomole):
    """ Class for (TD)DFT method of Turbomole program

        :param object molecule: molecule object
        :param string functional: level of DFT theory
        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param integer scf_max_iter: maximum number of SCF iterations
        :param double scf_en_tol: energy convergence for SCF iterations
        :param string qm_path: path for QM turbomole
        :param integer nthreads: number of threads in the calculations
        :param double version: version of Turbomole program
    """
    def __init__(self, molecule, functional="b-lyp", basis_set="SV(P)", memory="", \
        scf_max_iter=50, scf_en_tol=1E-8, qm_path="./", nthreads=1, version=6.4):
        # Initialize Turbomole common variables
        super(DFT, self).__init__(functional, basis_set, memory, qm_path, nthreads, version)

        self.scf_max_iter = scf_max_iter
        self.scf_en_tol = scf_en_tol

        # Set 'l_nacme' with respect to the computational method
        # TDDFT cannot produce NAC between ground and first excited state,
        # so we need to get NACME from CIoverlap
        molecule.l_nacme = False

        # Re-calculation of excited state forces is not needed for ground state dynamics
        self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFT method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)

        self.get_input(molecule, bo_list)
        self.run_QM(molecule, base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list):
        """ Generate Turbomole input files: define.in, control, etc

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
        """
        if (self.calc_coupling):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) BOMD is only valid! {self.calc_coupling}")

        x2t_command = os.path.join(self.qm_scripts_path, "x2t")
        command = f"{x2t_command} geometry.xyz > coord"
        os.system(command)

        input_define = ""

        # Job title
        input_title = textwrap.dedent(f"""\


        """)
        input_define += input_title

        # Molcule geometry
        input_geom = textwrap.dedent(f"""\
        a coord
        *
        no
        """)
        input_define += input_geom

        # Basis set
        input_bss = textwrap.dedent(f"""\
        b all {self.basis_set}
        *
        """)
        input_define += input_bss

        # Occupation number and MO
        input_initMO = textwrap.dedent(f"""\
        eht

        {int(molecule.charge)}

        """)
        input_define += input_initMO

        # Functional
        input_functional = textwrap.dedent(f"""\
        dft
        func
        {self.functional}
        on
        *
        """)
        input_define += input_functional

        # TODO: TDA or RPA (rpas -> ciss)
        # Excited state calculation
        input_ES = textwrap.dedent(f"""\
        ex
        rpas
        q
        a {molecule.nst - 1}
        q
        q

        *
        """)
        input_define += input_ES

        file_name = "define.in"
        with open(file_name, "w") as f:
            f.write(input_define)

        define_command = os.path.join(self.qm_bin_path, "define")
        command = f"{define_command} < {file_name} >& define_log"
        os.system(command)

        file_name = "control"
        with open(file_name, "r") as f:
            control_prev = f.readlines()

        control = ""
        control += f"$exopt {bo_list[0]}\n"

        iline = 0
        while "$scfiterlimit" not in control_prev[iline]:
            control += control_prev[iline]
            iline += 1
        scfiter = f"$scfiterlimit   {self.scf_max_iter}\n"
        control += scfiter
        iline += 1

        while "$dft" not in control_prev[iline]:
            control += control_prev[iline]
            iline += 1
        control += control_prev[iline]
        control += "  weight derivatives\n"
        iline += 1

        while iline != len(control_prev):
            control += control_prev[iline]
            iline += 1

        with open(file_name, "w") as f:
            f.write(control)

    def run_QM(self, molecule, base_dir, istep, bo_list):
        """ Run (TD)DFT calculation and save the output files to QMlog directory

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Run dscf
        scf_command = os.path.join(self.qm_bin_path, "dscf")
        # OpenMP setting
        # TODO : parallel implementation
        #os.environ["OMP_NUM_THREADS"] = f"{self.nthreads}"
        command = f"{scf_command} >& dscf.out"
        os.system(command)
        
        if (bo_list[0] == 0):
            grad_command = os.path.join(self.qm_bin_path, "grad")
            command = f"{grad_command} >& grad.out"
            os.system(command)
            if (molecule.nst > 1):
                grad_command = os.path.join(self.qm_bin_path, "escf")
                command = f"{grad_command} >& escf.out"
                os.system(command)
        else:
            egrad_command = os.path.join(self.qm_bin_path, "egrad")
            command = f"{egrad_command} >& egrad.out"
            os.system(command)

        # Copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            shutil.copy("gradient", os.path.join(tmp_dir, f"grad.{istep + 1}"))

    def extract_BO(self, molecule, bo_list):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
        """
        file_name = "gradient"
        with open(file_name, "r") as f:
            bo_out = f.read()
        bo_out = bo_out.replace('D', 'E')

        # Energy of running state
        for states in molecule.states:
            states.energy = 0.

        find_e = "energy =\s+([-]\d+[.]\d+)"
        energy = re.findall(find_e, bo_out)
        energy = np.array(energy)
        energy = energy.astype(float)
        molecule.states[bo_list[0]].energy = energy[0]

        # Force of running state
        for states in molecule.states:
            states.force = np.zeros((molecule.nat, molecule.nsp))

        find_grad = "\s+([-]*\d*[.]\d+[E][-|+]\d+)"
        grad = re.findall(find_grad, bo_out)
        grad = np.array(grad)
        grad = grad.astype(float)
        grad = grad.reshape(molecule.nat, 3, order='C')
        molecule.states[bo_list[0]].force = - np.copy(grad)

        # Energy of other states (except running state)
        if (molecule.nst > 1):
            if (bo_list[0] != 0):
                file_name = "egrad.out"
            else:
                file_name = "escf.out"

            with open(file_name, "r") as f:
                bo_out = f.read()
            find_e = 'Total energy:\s+([-]\d+[.]\d+)'
            energy = re.findall(find_e, bo_out)
            energy = np.array(energy)
            energy = energy.astype(float)
            for ist in range(molecule.nst):
                if (ist != bo_list[0]):
                    molecule.states[ist].energy = energy[ist]

        # NACME
        # Turbomole cannot provides NACVs between excited states
