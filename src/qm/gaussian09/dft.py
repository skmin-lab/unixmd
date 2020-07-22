from __future__ import division
from qm.gaussian09.gaussian09 import Gaussian09
from misc import au_to_A, eV_to_au, call_name
import os, shutil, re, textwrap, subprocess
import numpy as np

class DFT(Gaussian09):
    """ Class for the (TD)DFT method of Gaussian09 program

        :param object molecule: molecule object
        :param integer nthreads: number of threads in the calculations
        :param string memory: allocatable memory in the calculations
        :param string functional: level of DFT theory
        :param string basis_set: basis set information
        :param string guess: initial guess type
        :param string guess_file: initial guess file
        :param string g09_root_path: path for Gaussian09 root
        :param string version: version of Gaussian09 program
    """
    def __init__(self, molecule, nthreads=1, memory="1gb", \
        functional="BLYP", basis_set="STO-3G", \
        guess="Harris", guess_file="./g09.chk", \
        g09_root_path="/opt/gaussian/", version="Revision A.02"):
        # Initialize Gaussian09 common variables
        super(DFT, self).__init__(basis_set, memory, nthreads, g09_root_path, version)

        # Initialize Gaussian09 DFT variables
        self.functional = functional

        # Set initial guess for DFT calculation
        self.guess = guess
        self.guess_file = os.path.abspath(guess_file)
        if (not (self.guess == "Harris" or self.guess == "read")):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) Wrong input for initial guess option! {self.guess}")

        # Set 'l_nacme' with respect to the computational method
        molecule.l_nacme = True

        # Re-calculation of excited state forces is not needed for ground state dynamics
        if (molecule.nst > 1):
            self.re_calc = True
        else:
            self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient from (TD)DFT method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        self.copy_files(istep)
        super().get_data(base_dir, calc_force_only)
        self.get_input(molecule, istep, bo_list)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_QM(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def copy_files(self, istep):
        """ Copy necessary scratch files in previous step

            :param integer istep: current MD step
        """
        # Copy required files to read initial guess
        if (self.guess == "read" and istep >= 0):
            # After T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "g09.chk"), \
                os.path.join(self.scr_qm_dir, "../g09.chk.pre"))

    def get_input(self, molecule, istep, bo_list):
        """ Generate Gaussian09 input files: g09.inp

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # TODO : currently, CIoverlap is not correct -> only BOMD possible with TDDFT
        if (self.calc_coupling):
            raise ValueError ("only BOMD possible with TDDFT")

        # Read check-point file from previous step
        if (self.guess == "read"):
            if (istep == -1):
                if (os.path.isfile(self.guess_file)):
                    # Copy guess file to currect directory
                    shutil.copy(self.guess_file, os.path.join(self.scr_qm_dir, "g09.chk"))
                    restart = True
                else:
                    # TODO : Printout about reading a checkpoint file for the initial guess
                    # print(f"( {self.qm_method}.{call_name()} ) Make the initial guess of density only for the 1st step.\n", flush=True)
                    restart = False
            elif (istep >= 0):
                # Move previous file to currect directory
                os.rename("../g09.chk.pre", "./g09.chk")
                restart = True
        elif (self.guess == "Harris"):
            restart = False

        # Make 'g09.inp' file
        input_g09 = ""

        # Route section block
        input_route = textwrap.dedent(f"""\
        %nproc={self.nthreads}
        %mem={self.memory}
        %chk=g09.chk
        # {self.functional}/{self.basis_set} force""")

        if (restart):
            input_route += f" Guess=Read"

        if (bo_list[0] != 0):
            input_route += f" td(Root={bo_list[0]}, Nstates={molecule.nst - 1})\n\n"
        else:
            input_route += f"\n\n"
        input_g09 += input_route

        # Title section block
        input_title = f"g09 input\n\n"
        input_g09 += input_title

        # Molecule specification block
        input_molecule = textwrap.dedent(f"""\
        {int(molecule.charge)} 1
        """)
        for iat in range(molecule.nat):
            list_pos = list(molecule.pos[iat] * au_to_A)
            input_molecule += \
                f"{molecule.symbols[iat]}{list_pos[0]:15.8f}{list_pos[1]:15.8f}{list_pos[2]:15.8f}\n"
        input_molecule += "\n"
        input_g09 += input_molecule

        # Write 'g09.inp' file
        file_name = "g09.inp"
        with open(file_name, "w") as f:
            f.write(input_g09)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run (TD)DFT calculation and save the output files to QMlog directory

            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Set environment variables
        # TODO : move to gaussian09.py
        os.environ["g09root"] = self.g09_root_path
        os.environ["GAUSS_SCDIR"] = self.scr_qm_dir
        path_profile = os.path.join(self.g09_root_path, "g09/bsd/g09.profile")
        command = f'env -i sh -c "source {path_profile} && env"'
        for line in subprocess.getoutput(command).split("\n"):
            key, value = line.split("=")
            os.environ[key] = value

        # Set run command
        qm_command = os.path.join(self.g09_root_path, "g09/g09")
        command = f"{qm_command} < g09.inp > log"

        # Run Gaussian09
        os.system(command)

        # Copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_QM(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        file_name = "log"
        with open(file_name, "r") as f:
            log = f.read()

        # Energy
        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            # Read ground energy
            energy = re.findall('SCF Done:\s+E\(\w+\)\s+=\s+([-]\S+)\s+A.U.', log)
            energy = np.array(energy[0])
            energy = energy.astype(float)
            molecule.states[0].energy = energy

            if (molecule.nst > 1 and bo_list[0] != 0):
                energy = re.findall('Excited\sState\s+\w+:\s+\w+-\w+\s+(\S+)\s+eV', log)
                energy = np.array(energy)
                energy = energy.astype(float)
                energy *= eV_to_au
                for ist in range(1, molecule.nst):
                    molecule.states[ist].energy = molecule.states[0].energy + energy[ist - 1]

        # Force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        tmp_f = "Forces\s+\(Hartrees\/Bohr\)\n.+\n.+" \
            + "\n\s+\d*\s+\d*\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)" * molecule.nat
        force = re.findall(tmp_f, log)
        force = np.array(force[0])
        force = force.astype(float)
        force = force.reshape(molecule.nat, 3, order='C')
        molecule.states[bo_list[0]].force = np.copy(force)

        # NACME
        pass


