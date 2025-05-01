from __future__ import division
from qm.terachem.terachem import TeraChem
from misc import call_name
import os, shutil, re, textwrap
import numpy as np

class DFT(TeraChem):
    """ Class for (TD)DFT method of TeraChem

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param string functional: Exchange-correlation functional information
        :param string precision: Precision in the calculations
        :param double scf_wf_tol: Wavefunction convergence for SCF iterations
        :param integer scf_max_iter: Maximum number of SCF iterations
        :param boolean l_spin_pol: Include spin-polarisation scheme
        :param integer spin_mult: Spin multiplicity of wavefunction
        :param string guess: Initial guess for SCF iterations
        :param string guess_file: Initial guess file
        :param string root_path: Path for TeraChem root directory
        :param integer ngpus: Number of GPUs
        :param integer,list gpu_id: ID of used GPUs
        :param string version: Version of TeraChem
    """
    def __init__(self, molecule, ngpus=1, gpu_id=None, precision="dynamic", \
        version="1.93", functional="hf", basis_set="sto-3g", scf_wf_tol=3E-5, \
        scf_max_iter=100, l_spin_pol=False, spin_mult=1, guess="dft", \
        guess_file="./c0", root_path="./"):
        # Initialize TeraChem common variables
        super(DFT, self).__init__(functional, basis_set, root_path, ngpus, \
            gpu_id, precision, version)

        # Initialize TeraChem DFT variables
        self.scf_wf_tol = scf_wf_tol
        self.scf_max_iter = scf_max_iter

        self.l_spin_pol = l_spin_pol
        self.spin_mult = spin_mult
        if (self.l_spin_pol):
            if ((int(molecule.nelec) % 2 == 0) and (self.spin_mult % 2 == 0)):
                error_message = "Spin multiplicity for DFT must be singlet, triplet, and so on!"
                error_vars = f"spin_mult = {self.spin_mult}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
            elif ((int(molecule.nelec) % 2 == 1) and (self.spin_mult % 2 == 1)):
                error_message = "Spin multiplicity for DFT must be doublet, quartet, and so on!"
                error_vars = f"spin_mult = {self.spin_mult}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            if (self.spin_mult != 1):
                error_message = "Invalid spin multiplicity for DFT!"
                error_vars = f"spin_mult = {self.spin_mult}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set initial guess for SCF iterations
        self.guess = guess.lower()
        self.guess_file = guess_file
        if not (self.guess in ["read"]):
            error_message = "Invalid initial guess for DFT!"
            error_vars = f"guess = {self.guess}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Error for excited state dynamics
        if (molecule.nst > 1):
            error_message = "Excited state dynamics using TDDFT not implemented!"
            error_vars = f"Molecule.nstates = {molecule.nst}"
            raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set 'l_nacme' with respect to the computational method
        molecule.l_nacme = False

        # Re-calculation of excited state forces is not needed for ground state dynamics
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only, traj=None):
        """ Extract energy, gradient and nonadiabatic couplings from DFT method

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer,list bo_list: List of BO states for BO calculation
            :param double dt: Time interval
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
            :param object traj: Trajectory object containing the calculator and trajectory
        """
        self.copy_files(istep)
        super().get_data(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, istep)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_QM(molecule)
        self.move_dir(base_dir)

    def copy_files(self, istep):
        """ Copy necessary scratch files in previous step

            :param integer istep: Current MD step
        """
        # Copy required files to read initial guess
        if (self.guess == "read" and istep >= 0):
            # After T = 0.0 s
            if (self.l_spin_pol):
                shutil.copy(os.path.join(self.scr_qm_dir, "./scr/ca0"), \
                    os.path.join(self.scr_qm_dir, "../ca0"))
                shutil.copy(os.path.join(self.scr_qm_dir, "./scr/cb0"), \
                    os.path.join(self.scr_qm_dir, "../cb0"))
            else:
                shutil.copy(os.path.join(self.scr_qm_dir, "./scr/c0"), \
                    os.path.join(self.scr_qm_dir, "../c0"))

    def get_input(self, molecule, istep):
        """ Generate TeraChem input files: input.tcin

            :param object molecule: Molecule object
            :param integer istep: Current MD step
        """
        # Make 'input.tcin' file
        input_terachem = ""

        # Guess Block
        if (self.guess == "read"):
            if (istep == -1):
                if (os.path.isfile(self.guess_file)):
                    # Copy guess file to currect directory
                    shutil.copy(self.guess_file, os.path.join(self.scr_qm_dir, "c0"))
                    restart = True
                else:
                    restart = False
            elif (istep >= 0):
                # Move previous file to currect directory
                if (self.l_spin_pol):
                    os.rename("../ca0", os.path.join(self.scr_qm_dir, "ca0"))
                    os.rename("../cb0", os.path.join(self.scr_qm_dir, "cb0"))
                else:
                    os.rename("../c0", os.path.join(self.scr_qm_dir, "c0"))
                restart = True
        elif (self.guess == "dft"):
            restart = False

        # Control Block
        input_control = textwrap.dedent(f"""\

        run gradient

        coordinates geometry.xyz

        """)
        input_terachem += input_control

        # System Block
        gpu_ids = ""
        for gpu_cur_id in self.gpu_id:
            gpu_ids += f" {gpu_cur_id} "

        input_system = textwrap.dedent(f"""\
        scrdir ./scr/
        precision {self.precision}
        gpus {self.ngpus}   {gpu_ids}

        """)
        input_terachem += input_system

        # Setting Block
        if (self.l_spin_pol):
            functional = f"u{self.functional}"
        else:
            functional = self.functional

        input_setting = textwrap.dedent(f"""\
        method {functional}
        basis {self.basis_set}
        charge {molecule.charge}
        spinmult {self.spin_mult}

        """)
        input_terachem += input_setting

        # DFT Block
        input_dft = textwrap.dedent(f"""\
        convthre {self.scf_wf_tol}
        maxit {self.scf_max_iter}
        """)
        input_terachem += input_dft

        if (restart):
            if (self.l_spin_pol):
                input_guess = textwrap.dedent(f"""\
                guess ./ca0 ./cb0
                """)
            else:
                input_guess = textwrap.dedent(f"""\
                guess ./c0
                """)
            input_terachem += input_guess

        # TODO: pointcharges? in qmmm?

        # Write 'input.tcin' file
        file_name = "input.tcin"
        with open(file_name, "w") as f:
            f.write(input_terachem)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run DFT calculation and save the output files to qm_log directory

            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Run TeraChem method
        qm_command = os.path.join(self.qm_path, "terachem")
        # OpenMP setting
        os.environ["OMP_NUM_THREADS"] = "1"
        command = f"{qm_command} input.tcin > log"
        os.system(command)
        # Copy the output file to 'qm_log' directory
        tmp_dir = os.path.join(base_dir, "qm_log")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_QM(self, molecule):
        """ Read the output files to get BO information

            :param object molecule: Molecule object
        """
        # Read 'log' file
        file_name = "log"
        with open(file_name, "r") as f:
            log_out = f.read()

        # Energy
        tmp_e = 'FINAL ENERGY:\s+([-]\S+)'
        energy = re.findall(tmp_e, log_out)
        energy = np.array(energy, dtype=np.float64)
        molecule.states[0].energy = energy[0]

        # Force
        tmp_g = 'Gradient units are Hartree/Bohr\n[-]+\n\s+dE/dX\s+dE/dY\s+dE/dZ' + \
	          '\n\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)' * molecule.nat_qm
        grad = re.findall(tmp_g, log_out)
        grad = np.array(grad[0], dtype=np.float64)
        grad = grad.reshape(molecule.nat_qm, 3, order='C')
        molecule.states[0].force[0:molecule.nat_qm] = - np.copy(grad)


