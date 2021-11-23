from __future__ import division
from qm.turbomole.turbomole import Turbomole
from misc import call_name
import os, shutil, re, textwrap
import numpy as np

class DFT(Turbomole):
    """ Class for (TD)DFT method of Turbomole

        :param object molecule: Molecule object
        :param string functional: Exchange-correlation functional information
        :param string basis_set: Basis set information
        :param integer memory: Allocatable memory in the calculations
        :param integer scf_max_iter: Maximum number of SCF iterations
        :param integer scf_en_tol: Energy convergence for SCF iterations
        :param integer cis_max_iter: Maximum number of CIS iterations
        :param integer cis_en_tol: Energy convergence for CIS iterations
        :param string root_path: Path for Turbomole root directory
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of Turbomole
    """
    def __init__(self, molecule, functional="b-lyp", basis_set="SV(P)", memory=50, \
        scf_max_iter=50, scf_en_tol=6, cis_max_iter=25, cis_en_tol=6, \
        root_path="./", nthreads=1, version="6.4"):
        # Initialize Turbomole common variables
        super(DFT, self).__init__(functional, basis_set, memory, root_path, nthreads, version)

        self.scf_max_iter = scf_max_iter
        self.scf_en_tol = scf_en_tol
        self.cis_max_iter = cis_max_iter
        self.cis_en_tol = cis_en_tol

        # Set 'l_nacme' with respect to the computational method
        # TDDFT cannot produce NAC between excited states,
        # so we need to get NACME from CIoverlap but Turbomole does not provide AO overlap.
        # Hence, CIoverlap is not valid yet.
        molecule.l_nacme = False

        # Re-calculation of excited state forces is not needed for ground state dynamics
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from (TD)DFT method

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer,list bo_list: List of BO states for BO calculation
            :param double dt: Time interval
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        super().get_data(base_dir, calc_force_only)
        self.write_xyz(molecule)

        self.get_input(molecule, bo_list)
        self.run_QM(molecule, base_dir, istep, bo_list)
        self.extract_QM(molecule, bo_list)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list):
        """ Generate Turbomole input files: define.in, control, etc

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
        """
        if (self.calc_coupling):
            error_message = "Turbomole supports only BOMD!"
            error_vars = f"qm_prog.qm_method = {qm.qm_prog}.{qm.qm_method}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        x2t_command = os.path.join(self.scripts_path, "x2t")
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

        define_command = os.path.join(self.qm_path, "define")
        command = f"{define_command} < {file_name} >& define_log"
        os.system(command)

        file_name = "control"
        with open(file_name, "r") as f:
            control_prev = f.readlines()

        control = ""

        # Root state to calculate gradient
        control += f"$exopt {bo_list[0]}\n"

        # Memory to use
        control += f"$maxcor {self.memory}\n"

        iline = 0
        # SCF options such as iteration and energy tolerance
        while "$scfiterlimit" not in control_prev[iline]:
            control += control_prev[iline]
            iline += 1

        control += f"$scfiterlimit   {self.scf_max_iter}\n"
        iline += 1
        control += f"$scfconv {self.scf_en_tol}\n"

        if (molecule.nst > 1):
            control += f"$rpacor {self.memory}\n"
            control += f"$rpaconv {self.cis_en_tol}\n"
            control += f"$escfiterlimit {self.cis_max_iter}\n"
            
        # Calculate energy gradient
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
        """ Run (TD)DFT calculation and save the output files to qm_log directory

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Run dscf
        scf_command = os.path.join(self.qm_path, "dscf")
        command = f"{scf_command} >& dscf.out"
        os.system(command)

        if (bo_list[0] == 0):
            grad_command = os.path.join(self.qm_path, "grad")
            command = f"{grad_command} >& grad.out"
            os.system(command)
            if (molecule.nst > 1):
                grad_command = os.path.join(self.qm_path, "escf")
                command = f"{grad_command} >& escf.out"
                os.system(command)
        else:
            egrad_command = os.path.join(self.qm_path, "egrad")
            command = f"{egrad_command} >& egrad.out"
            os.system(command)

        # Copy the output file to 'qm_log' directory
        tmp_dir = os.path.join(base_dir, "qm_log")
        if (os.path.exists(tmp_dir)):
            shutil.copy("dscf.out", os.path.join(tmp_dir, f"dscf.out.{istep + 1}"))
            if (bo_list[0] == 0):
                shutil.copy("grad.out", os.path.join(tmp_dir, f"grad.out.{istep + 1}"))
                if (molecule.nst > 1):
                    shutil.copy("escf.out", os.path.join(tmp_dir, f"escf.out.{istep + 1}"))
            else:
                shutil.copy("egrad.out", os.path.join(tmp_dir, f"egrad.out.{istep + 1}"))

    def extract_QM(self, molecule, bo_list):
        """ Read the output files to get BO information

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
        """
        file_name = "gradient"
        with open(file_name, "r") as f:
            bo_out = f.read()

        bo_out = bo_out.replace('D', 'E')

        # Energy of running state
        find_e = "energy =\s+([-]\d+[.]\d+)"
        energy = re.findall(find_e, bo_out)
        energy = np.array(energy, dtype=np.float64)
        molecule.states[bo_list[0]].energy = energy[0]

        # Force of running state
        find_grad = "\s+([-]*\d*[.]\d+[E][-|+]\d+)"
        grad = re.findall(find_grad, bo_out)
        grad = np.array(grad, dtype=np.float64)
        grad = grad.reshape(molecule.nat_qm, 3, order='C')
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
            energy = np.array(energy, dtype=np.float64)
            for ist in range(molecule.nst):
                if (ist != bo_list[0]):
                    molecule.states[ist].energy = energy[ist]

        # NACME
        # Turbomole cannot provides NACVs between excited states
        pass


