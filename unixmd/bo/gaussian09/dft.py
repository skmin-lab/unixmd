from __future__ import division
import os, shutil, re, subprocess
import numpy as np
import textwrap
from bo.gaussian09.gaussian09 import Gaussian09
from misc import au_to_A, eV_to_au

class DFT(Gaussian09):
    """ Class for Gaussian09 (TD-)DFT method
                  bomd | sh | eh | nac | re_calc
        DFT    :   o     x    x     F      F
        TD-DFT :   o     x    x     F      T
    """
    def __init__(self, molecule, nthreads=8, memory="4gb",\
        functional="BLYP", basis_set="STO-3G", \
        g09_root_path="/opt/gaussian/", version="Revision A.02"):
        # Initialize Gaussian09 common variables
        super().__init__(basis_set, memory, nthreads, g09_root_path, version)

        # Initialize Gaussian09 DFT variables
        self.functional = functional

        # Set 'l_nacme' with respect to the computational method
        molecule.l_nacme = True

        # Re-calculation of excited state forces is not needed for ground state dynamics
        if (molecule.nst > 1):
            self.re_calc = True
        else:
            self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Get/Extract BO information from Gaussian09
        """
        super().get_bo(base_dir, calc_force_only)
        self.get_input(molecule, istep, bo_list)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, istep, bo_list):
        """ Generate Gaussian09 input files: g09.inp
        """
        if (len(bo_list) > 1):
            raise ValueError(f"Ehrenfest dynamics with g09 is not implemented yet!")

        # make 'g09.inp' file
        input_g09 = ""

        # Route section block
        input_route = textwrap.dedent(f"""\
        %nproc={self.nthreads}
        %mem={self.memory}
        %chk=g{istep}.chk
        # {self.functional}/{self.basis_set} force""")

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

        # write 'g09.inp' file
        file_name = "g09.inp"
        with open(file_name, "w") as f:
            f.write(input_g09)

    def run_QM(self, base_dir, istep, bo_list):
        """ run g09 calculation and save the output files
        """
        # set environment variables
        os.environ["g09root"] = self.g09_root_path
        os.environ["GAUSS_SCDIR"] = self.scr_qm_dir
        path_profile = os.path.join(self.g09_root_path, "g09/bsd/g09.profile")
        command = f'env -i sh -c "source {path_profile} && env"'
        for line in subprocess.getoutput(command).split("\n"):
            key, value = line.split("=")
            os.environ[key] = value

        # set run command
        qm_command = os.path.join(self.g09_root_path, "g09/g09")
        command = f"{qm_command} < g09.inp > log"

        # Run Gaussian09
        os.system(command)

        # copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ read the output files to get BO data
        """
        file_name = "log"
        with open(file_name, "r") as f:
            log = f.read()

        # energy
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

        # force
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

        # nacme
        pass


