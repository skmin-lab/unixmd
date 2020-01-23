from __future__ import division
from bo.turbomole.turbomole import Turbomole
import os, shutil, re
import numpy as np
import textwrap

class DFT(Turbomole):
    """ Class for TURBOMOLE TD-DFT method

        :param object molecule: molecule object
        :param string functional: xc functional
        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param integer max_iter: maximum number of SCF iterations
        :param double scf_en_tol: energy convergence for SCF iterations
        :param string qm_path: path for QM turbomole
        :param string qm_bin_path: path for QM binary
        :param string qm_scripts_path: path for QM scripts
        :param integer nthreads: number of threads in the calculations
        :param double version: version of Turbomole program
    """
    def __init__(self, molecule, functional="b-lyp", basis_set="STO-3g", memory="", \
        max_iter=50, scf_en_tol=1E-8, qm_path="./", qm_bin_path="./", qm_scripts_path="./", \
        nthreads=1, version=6.4):
        #Initialize
        super().__init__(functional, basis_set, memory, qm_path, nthreads, version)

        self.max_iter = max_iter
        self.scf_en_tol = scf_tol
        
        #self.qm_bin_path = os.path.join(self.qm_path, "bin/em64t-unknown-linux-gnu_smp/")
        self.qm_bin_path = os.path.join(self.qm_path, "bin/em64t-unknown-linux-gnu/")
        self.qm_scripts_path = os.path.join(self.qm_path, "scripts/")
        
        #os.environ["PARA_ARCH"] = "SMP"
        #os.environ["PARNODES"] = f"{self.nthreads}"
        os.environ["TURBODIR"] = qm_path
        
        # TURBOMOLE can provide NAC except NAC between ground state and first excited state.
        #
        molecule.l_nacme = False
        self.re_calc = True
    
    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Get/Extract BO information from TURBOMOLE

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        if (self.calc_coupling):
            raise ValueError ("Turbomole is available only for BOMD.")
        
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        x2t_command = os.path.join(self.qm_scripts_path, "x2t")
        command = f"{x2t_command} geometry.xyz > coord" 
        os.system(command)      
        #os.system("cp coord coord.ref")
        
        
        self.get_input(molecule, bo_list)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)
    
    def get_input(self, molecule, bo_list):
        """ Generate TURBOMOLE input files

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
        """
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
        
        # Occupation # and MO
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
        
        # Excited state calculation
        input_ES = textwrap.dedent(f"""\
        ex
        rpas
        q
        a {molecule.nst-1}
        q
        q

        *
        """)
        input_define += input_ES

        file_name = "define.in"
        with open(file_name, "w") as f:
            f.write(input_define)
        
        define_command = os.path.join(self.qm_bin_path, "define")
        command = f"{define_command} < {file_name} > define_log"
        os.system(command)
        
        file_name = "control"
        with open(file_name, "r") as f:
            control_prev = f.readlines()

        control = ""
        iline = 0
        while "$scfiterlimit" not in control_prev[iline]:
            control += control_prev[iline]
            iline += 1
        scfiter = f"$scfiterlimit   {self.max_iter}\n"
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

    def run_QM(self, base_dir, istep, bo_list):
        """ run TURBOMOLE calculation and save the output files

            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # run dscf
        scf_command = os.path.join(self.qm_bin_path, "dscf")
        # openmp setting
        # TODO
        #os.environ["OMP_NUM_THREADS"] = f"{self.nthreads}"
        command = f"{scf_command} > dscf.out"
        os.system(command)

        if (bo_list[0] == 0):
            grad_command = os.path.join(self.qm_bin_path, "grad")
            command = f"{grad_command} > grad.out"
        else:
            egrad_command = os.path.join(self.qm_bin_path, "egrad")
            command = f"{egrad_command} > egrad.out"
        os.system(command)

        #  copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            shutil.copy("gradient", os.path.join(tmp_dir, f"grad.{istep + 1}"))

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ read output file and extract data

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """

        file_name = "gradient"
        with open(file_name, "r") as f:
            bo_out = f.read()
        bo_out = bo_out.replace('D', 'E')

        # force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        find_grad = "\s+([-]*\d*[.]\d+[E][-|+]\d+)"
        grad = re.findall(find_grad, bo_out)
        grad = np.array(grad)
        grad = grad.astype(float)
        grad = grad.reshape(molecule.nat, 3, order='C')

        molecule.states[bo_list[0]].force = - np.copy(grad)

        # energy
        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            if (bo_list[0] == 0):
                find_e = "energy =\s+([-]\d+[.]\d+)"

                energy = re.findall(find_e, bo_out)
                energy = np.array(energy)
                energy = energy.astype(float)

                molecule.states[0].energy = energy[0]
            else:
                file_name = "egrad.out"
                with open(file_name, "r") as f:
                    bo_out = f.read()
                find_e ='Total energy:\s+([-]\d+[.]\d+)'

                energy = re.findall(find_e, bo_out)
                energy = np.array(energy)
                energy = energy.astype(float)
                for ist in range(molecule.nst):
                    molecule.states[ist].energy = energy[ist]
