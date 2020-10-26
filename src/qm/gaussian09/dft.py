from __future__ import division
from build.cioverlap import wf_overlap 
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

        # MO dimension, initialized later by reading Gaussian log
        self.nbasis = 0
        self.norb = 0
        self.nfc = 0
        self.nocc = 0
        self.nvirt = 0
        
        # Temporaries for NACME calculation, also initialized later if not allocated
        # ao_overlap - the nr. of AOs, the nr. of AOs
        # mo_coef - the nr. of MOs, the nr. of AOs
        # ci_coef - the nr. BO states, the nr. of occ, the nr. of virt
        self.pos_old = []
        self.ao_overlap = []
        self.mo_coef_old = []
        self.mo_coef_new = []
        self.ci_coef_old = []
        self.ci_coef_new = []

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient from (TD)DFT method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        self.copy_files(molecule, istep, calc_force_only)
        super().get_data(base_dir, calc_force_only)
        self.get_input(molecule, istep, bo_list, calc_force_only)
        self.run_QM(molecule, base_dir, istep, bo_list, calc_force_only)
        self.extract_QM(molecule, istep, bo_list, dt, calc_force_only)
        self.move_dir(base_dir)

    def copy_files(self, molecule, istep, calc_force_only):
        """ Copy necessary scratch files in previous step

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Copy required files for NACME
        if (molecule.nst > 1 and not calc_force_only and istep >= 0):
            if (istep == 0):
                shutil.copy(os.path.join(self.scr_qm_dir, "g09.rwf"), \
                    os.path.join(self.scr_qm_dir, "../g09.rwf.pre"))

        # Copy required files to read initial guess
        if (self.guess == "read" and istep >= 0):
            # After T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "g09.chk"), \
                os.path.join(self.scr_qm_dir, "../g09.chk.pre"))

    def get_input(self, molecule, istep, bo_list, calc_force_only):
        """ Generate Gaussian09 input files: g09.inp

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """

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
        %chk=g09.chk\n""")
        
        if (molecule.nst > 1 and not calc_force_only):
            input_route = input_route + textwrap.dedent(f"""\
            %rwf=g09.rwf\n""")
        
        input_route = input_route + textwrap.dedent(f"""\
        # {self.functional}/{self.basis_set}""")

        if (restart):
            input_route += f" Guess=Read"

        if (molecule.nst >1):
            if (self.calc_coupling):
                if (calc_force_only and bo_list[0] == 0):
                    input_route += f" force nosymm"
                else:
                    input_route += f" td(Root={bo_list[0]}, Nstates={molecule.nst - 1}) nosymm"
            # BOMD
            else:
                input_route += f" td(Root={bo_list[0]}, Nstates={molecule.nst - 1}) nosymm"

            if (bo_list[0] != 0):
                input_route += f" force\n\n"
            # If the running state is the ground state, an additional input is created later.
            else:
                input_route += f"\n\n"
        else:
            if (bo_list[0] != 0):
                input_route += f" td(Root={bo_list[0]}, Nstates={molecule.nst - 1}) force nosymm\n\n"
            else:
                input_route += f" force nosymm\n\n"
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

        # In the nonadiabatic case where the running state is the ground state, another static 
        # calculation without the td option must be done to provide the BO(ground-state) force.
        if (molecule.nst > 1 and bo_list[0] == 0 and not calc_force_only):
            input_g09 = ""

            input_route = textwrap.dedent(f"""\
            %nproc={self.nthreads}
            %mem={self.memory}
            %chk=g09_g.chk
            # {self.functional}/{self.basis_set}""")

            if (restart):
                input_route += f" Guess=Read"

            input_route += f" force nosymm\n\n"
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
    
            # Write 'g09_g.inp' file
            file_name = "g09_g.inp"
            with open(file_name, "w") as f:
                f.write(input_g09)
        
        # Write "doubled molecule" input
        if (self.calc_coupling and molecule.nst > 1 and not calc_force_only and istep >= 0):
            if (istep == 0):
                os.rename('../g09.rwf.pre', './g09.rwf.pre')
            input_g09 = ""
            
            # Stop the run after L302 calculating overlap
            # Keep running the job regardless of interatomic distances; IOp(2/12=3)
            input_route = textwrap.dedent(f"""\
            %nproc={self.nthreads}
            %mem={self.memory}
            %kjob l302
            %rwf=g09_double.rwf
            # {self.functional}/{self.basis_set} IOp(2/12=3) nosymm\n\n""")
            
            input_g09 += input_route
            
            # Title section block
            input_title = f"g09 double input\n\n"
            input_g09 += input_title
            
            # Molecule specification block
            input_molecule = textwrap.dedent(f"""\
            {2 * int(molecule.charge)} 1
            """)
            for iat in range(molecule.nat):
                list_pos = list(self.pos_old[iat] * au_to_A)
                input_molecule += \
                    f"{molecule.symbols[iat]}{list_pos[0]:15.8f}{list_pos[1]:15.8f}{list_pos[2]:15.8f}\n"
            
            for iat in range(molecule.nat):
                list_pos = list(molecule.pos[iat] * au_to_A)
                input_molecule += \
                    f"{molecule.symbols[iat]}{list_pos[0]:15.8f}{list_pos[1]:15.8f}{list_pos[2]:15.8f}\n"
            input_molecule += "\n"
            input_g09 += input_molecule
            
            file_name = "g09_double.inp"
            with open(file_name, "w") as f:
                f.write(input_g09)

    def run_QM(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Run (TD)DFT calculation and save the output files to QMlog directory

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
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

        # Run Gaussian for the ground-state force
        if (molecule.nst > 1 and bo_list[0] == 0 and not calc_force_only):
           qm_command = os.path.join(self.g09_root_path, "g09/g09")
           command = f"{qm_command} < g09_g.inp > log_g"
           os.system(command)

        # Run Gaussian09 for the overlap matrix
        if (self.calc_coupling and molecule.nst > 1 and not calc_force_only and istep >= 0):
           qm_command = os.path.join(self.g09_root_path, "g09/g09")
           command = f"{qm_command} < g09_double.inp > log_double"
           os.system(command)

    def extract_QM(self, molecule, istep, bo_list, dt, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        file_name = "log"
        with open(file_name, "r") as f:
            log = f.read()

        # Energy
        if (not calc_force_only):
            # Read ground energy
            energy = re.findall('SCF Done:\s+E\(\S+\)\s+=\s+([-]\S+)\s+A.U.', log)
            energy = np.array(energy[0], dtype=np.float)
            molecule.states[0].energy = energy

            if (molecule.nst > 1):
                energy = re.findall('Excited\sState\s+\w+:\s+\w+-\S+\s+(\S+)\s+eV', log)
                energy = np.array(energy, dtype=np.float)
                energy *= eV_to_au
                for ist in range(1, molecule.nst):
                    molecule.states[ist].energy = molecule.states[0].energy + energy[ist - 1]

        if (molecule.nst == 1 or bo_list[0] != 0 or calc_force_only):
            # Force
            tmp_f = "Forces\s+\(Hartrees\/Bohr\)\n.+\n.+" \
                + "\n\s+\d*\s+\d*\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)" * molecule.nat
            force = re.findall(tmp_f, log)
            force = np.array(force[0], dtype=np.float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = np.copy(force)
        else:
            file_name = "log_g"
            with open(file_name, "r") as f:
                log = f.read()

            # Force
            tmp_f = "Forces\s+\(Hartrees\/Bohr\)\n.+\n.+" \
                + "\n\s+\d*\s+\d*\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)" * molecule.nat
            force = re.findall(tmp_f, log)
            force = np.array(force[0], dtype=np.float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = np.copy(force)

        # NACME
        if (self.calc_coupling and molecule.nst > 1 and not calc_force_only):
            molecule.nacme = np.zeros((molecule.nst, molecule.nst))
            if (istep == -1):
                self.pos_old = np.zeros((molecule.nat, molecule.nsp))
                self.init_buffer(molecule)
            else:
                self.CI_overlap(molecule, istep, dt)
        
            # Save geometry in the buffer
            self.pos_old = np.copy(molecule.pos)

    def init_buffer(self, molecule):
        """ Initialize buffer variables to get NACME
    
            :param object molecule: molecule object
        """
        file_name = "log"
        with open(file_name, "r") as f:
            log = f.read()
    
        self.nbasis = re.findall('NBasis=\s+(\d+)\s+', log)
        self.nbasis = int(self.nbasis[0])
        self.nfc = re.findall('NFC=\s+(\d+)\s+', log)
        self.nfc = int(self.nfc[0])
        self.nocc = re.findall('NOA=\s+(\d+)\s+', log)
        self.nocc = int(self.nocc[0])
        self.nvirt = re.findall('NVA=\s+(\d+)\s+', log)
        self.nvirt = int(self.nvirt[0])
        self.norb = self.nocc + self.nvirt
    
        self.ao_overlap = np.zeros((self.nbasis, self.nbasis))
        self.mo_coef_old = np.zeros((self.norb, self.nbasis))
        self.mo_coef_new = np.zeros((self.norb, self.nbasis))
        self.ci_coef_old = np.zeros((molecule.nst, self.nocc, self.nvirt))
        self.ci_coef_new = np.zeros((molecule.nst, self.nocc, self.nvirt))

    def CI_overlap(self, molecule, istep, dt):
        """ Read the necessary files and calculate NACME from tdnac.c routine
            note that only reading of several files is required in this method
    
            :param object molecule: molecule object
            :param integer istep: current MD step
            :param double dt: time interval
        """
        # Read overlap
        path_rwfdump = os.path.join(self.g09_root_path, "g09/rwfdump")
        os.system(path_rwfdump+' g09_double.rwf ao_overlap.dat 514R')
    
        with open('ao_overlap.dat', "r") as f:
            log = f.read()
    
        tmp = re.findall('[-]?\d+\.\d+D[+-]\d\d', log)
        tmp = [float(x.replace('D', 'e')) for x in tmp]
    
        tmp_ovr = np.zeros((self.nbasis * 2, self.nbasis * 2))
    
        cnt = 0
        for ibasis in range(self.nbasis * 2):
            for jbasis in range(ibasis + 1):
                tmp_ovr[ibasis, jbasis] = tmp[cnt]
                cnt += 1
    
        tmp_ovr = tmp_ovr + np.transpose(tmp_ovr) - np.diag(np.diag(tmp_ovr))
    
        # Slicing the components between t and t+dt
        self.ao_overlap = np.copy(tmp_ovr[:self.nbasis, self.nbasis:])
    
        # Read mo coefficients
        if (istep == 0):
            os.system(path_rwfdump+' g09.rwf.pre mo_coef.dat 524R')
        
            with open('mo_coef.dat', "r") as f:
                log = f.read()
        
            tmp = re.findall('[-]?\d+\.\d+D[+-]\d\d', log)
            tmp = np.array([x.replace('D','e') for x in tmp], dtype=np.float)
        
            tmp_mo = tmp.reshape(self.nbasis, self.nbasis)
        
            tmp_mo = tmp_mo[self.nfc:self.nbasis]
            self.mo_coef_old = np.copy(tmp_mo)
    
        os.system(path_rwfdump+' g09.rwf mo_coef.dat 524R')
        
        with open('mo_coef.dat', "r") as f:
            log = f.read()
        
        tmp = re.findall('[-]?\d+\.\d+D[+-]\d\d', log)
        tmp = np.array([x.replace('D','e') for x in tmp], dtype=np.float)
        
        tmp_mo = tmp.reshape(self.nbasis, self.nbasis)
        
        tmp_mo = tmp_mo[self.nfc:self.nbasis]
        self.mo_coef_new = np.copy(tmp_mo)
    
        # Read CI coefficients
        if (istep == 0):
            os.system(path_rwfdump+' g09.rwf.pre xy_coef.dat 635R')
        
            with open('xy_coef.dat', "r") as f:
                log = f.read()
        
            tmp = re.findall('[-]?\d+\.\S+[+-]\d+', log)
        
            # Drop the first 12 dummy elements
            tmp = tmp[12:]
        
            # Gaussian deals with 4 times as much roots as the input NStates value. 
            # the nr. of excitation function => nocc \times nvirt
            # spin degrees of freedom => 2
            # X+Y, X-Y => 2
            roots = (molecule.nst - 1) * 4
            num_coef = 4 * (self.nocc * self.nvirt) * roots
        
            tmp = tmp[:num_coef]
            tmp = [t.replace('D','e') for t in tmp]
            tmp = np.array(tmp, dtype=np.float)
            xpy, xmy = tmp.reshape(2, roots, 2, -1)
            x = 0.5 * (xpy + xmy)
        
            # Drop beta part and unrequested excited states
            x = x[:(molecule.nst - 1), 0, :]
        
            self.ci_coef_old[1:] = x.reshape(-1, self.nocc, self.nvirt)
        
        os.system(path_rwfdump+' g09.rwf xy_coef.dat 635R')
        
        with open('xy_coef.dat', "r") as f:
            log = f.read()
        
        tmp = re.findall('[-]?\d+\.\S+[+-]\d+', log)
        
        # Drop the first 12 dummy elements
        tmp = tmp[12:]
    
        # Gaussian deals with 4 times as much roots as the input NStates value.
        # the nr. of excitation function => nocc \times nvirt
        # spin degrees of freedom => 2
        # X+Y, X-Y => 2
        roots = (molecule.nst - 1) * 4
        num_coef = 4 * (self.nocc * self.nvirt) * roots
        
        tmp = tmp[:num_coef]
        tmp = np.array([x.replace('D','e') for x in tmp], dtype=np.float)
        xpy, xmy = tmp.reshape(2, roots, 2, -1)
        x = 0.5 * (xpy + xmy)
        
        # Drop beta part and unrequested excited states
        x = x[:(molecule.nst - 1), 0, :]
        
        self.ci_coef_new[1:] = x.reshape(-1, self.nocc, self.nvirt)
        
        # Calculate wavefunction overlap with orbital scheme
        # Reference: J. Phys. Chem. Lett. 2015, 6, 4200-4203
        wf_overlap(self, molecule, istep, dt)

