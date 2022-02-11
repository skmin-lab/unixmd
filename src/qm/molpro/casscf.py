from __future__ import division
from qm.molpro.molpro import Molpro
from misc import call_name
import os, shutil, re, textwrap
import numpy as np

class CASSCF(Molpro):
    """ Class for CASSCF method of Molpro

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string guess: Initial guess for (SA-)CASSCF method
        :param string guess_file: File containing initial guesses for (SA-)CASSCF calculations
        :param integer scf_max_iter: Maximum number of HF iterations
        :param double scf_en_tol: Energy convergence threshold for HF iterations
        :param double scf_rho_tol: Density convergence threshold for HF iterations
        :param integer mcscf_max_iter: Maximum number of (SA-)CASSCF iterations
        :param double mcscf_en_tol: Energy convergence threshold for (SA-)CASSCF iterations
        :param double mcscf_grad_tol: Gradient convergence threshold for (SA-)CASSCF iterations
        :param double mcscf_step_tol: Step length convergence threshold for (SA-)CASSCF iterations
        :param integer active_elec: Number of electrons in active space
        :param integer active_orb: Number of orbitals in active space
        :param double cpscf_grad_tol: Gradient convergence threshold for CP-CASSCF equations
        :param string qm_path: Path for QM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of Molpro
    """
    def __init__(self, molecule, basis_set="sto-3g", memory="500m", \
        guess="hf", guess_file="./wf.wfu", scf_max_iter=20, scf_en_tol=1E-8, scf_rho_tol=1E-6, \
        mcscf_max_iter=20, mcscf_en_tol=1E-8, mcscf_grad_tol=1E-6, mcscf_step_tol=1E-2, \
        active_elec=2, active_orb=2, cpscf_grad_tol=1E-7, \
        qm_path="./", nthreads=1, version="2015.1"):
        # Initialize Molpro common variables
        super(CASSCF, self).__init__(basis_set, memory, qm_path, nthreads, version)

        # Initialize Molpro CASSCF variables
        # Set initial guess for CASSCF calculation
        self.guess = guess.lower()
        self.guess_file = guess_file
        if not (self.guess in ["hf", "read"]):
            error_message = "Invalid initial guess for CASSCF!"
            error_vars = f"guess = {self.guess}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # HF calculation for initial guess of CASSCF calculation
        self.scf_max_iter = scf_max_iter
        self.scf_en_tol = scf_en_tol
        self.scf_rho_tol = scf_rho_tol

        # CASSCF calculation
        self.mcscf_max_iter = mcscf_max_iter
        self.mcscf_en_tol = mcscf_en_tol
        self.mcscf_grad_tol = mcscf_grad_tol
        self.mcscf_step_tol = mcscf_step_tol
        self.active_elec = active_elec
        self.active_orb = active_orb
        self.cpscf_grad_tol = cpscf_grad_tol

        # Check the closed shell for systems
        if (not int(molecule.nelec) % 2 == 0):
            error_message = "Only closed shell configuration supported, check charge!"
            error_vars = f"Molecule.nelec = {int(molecule.nelec)}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Molpro do not support periodic setting with CASSCF method

        # CASSCF calculation do not provide parallel computation
        # If your system provide parallel casscf, then this part should be removed
        if (self.nthreads > 1):
            error_message = "Parallel CASSCF not implemented!"
            error_vars = f"nthreads = {self.nthreads}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Calculate number of frozen, closed and occ orbitals in CASSCF method
        # No positive frozen core orbitals in CASSCF
        self.frozen_orb = 0
        self.closed_orb = int((int(molecule.nelec) - self.active_elec) / 2)
        self.occ_orb = int(self.closed_orb + self.active_orb)

        # Set 'l_nacme' with respect to the computational method
        # CASSCF can produce NACs, so we do not need to get NACME from CIoverlap
        # CASSCF can compute the gradient of several states simultaneously,
        #        but self.re_calc is set to be true to reduce cost.
        molecule.l_nacme = False
        self.re_calc = True

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from CASSCF method

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer,list bo_list: List of BO states for BO calculation
            :param double dt: Time interval
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        if (not calc_force_only):
            self.copy_files(istep)
        super().get_data(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, istep, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_QM(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def copy_files(self, istep):
        """ Copy necessary scratch files in previous step

            :param integer istep: Current MD step
        """
        # Copy required files to read initial guess
        if (self.guess == "read" and istep >= 0):
            # After T = 0.0 s
            shutil.copy(os.path.join(self.scr_qm_dir, "./wfu/wf.wfu"), \
                os.path.join(self.scr_qm_dir, "../wf.wfu"))

    def get_input(self, molecule, istep, bo_list, calc_force_only):
        """ Generate Molpro input files: molpro.inp

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Make 'molpro.inp' file
        input_molpro = ""

        # Scratch Block
        if (not calc_force_only):

            if (self.guess == "read"):
                wfu_dir = os.path.join(self.scr_qm_dir, "wfu")
                os.makedirs(wfu_dir)
                if (istep == -1):
                    if (os.path.isfile(self.guess_file)):
                        # Copy guess file to currect directory
                        shutil.copy(self.guess_file, os.path.join(wfu_dir, "wf.wfu"))
                        restart = "restart,2\n"
                        hf = False
                    else:
                        restart = ""
                        hf = True
                elif (istep >= 0):
                    # Move previous file to currect directory
                    os.rename("../wf.wfu", os.path.join(wfu_dir, "wf.wfu"))
                    restart = "restart,2\n"
                    hf = False
            elif (self.guess == "hf"):
                restart = ""
                hf = True

        else:

            # Read 'wf.wfu' for surface hopping dynamics when hop occurs
            restart = "restart,2\n"
            hf = False

        input_scr = textwrap.dedent(f"""\
        file,1,int.int,delete
        file,2,wf.wfu,unknown
        {restart}
        """)
        input_molpro += input_scr

        # Geometry Block
        input_geom = textwrap.dedent(f"""\
        geometry={{angstrom
        include,geometry.xyz,echo
        }}

        """)
        input_molpro += input_geom

        # Control Block
        input_control = textwrap.dedent(f"""\
        symmetry nosym
        basis={self.basis_set}

        """)
        input_molpro += input_control

        # HF Block: calculate energy option
        if (hf):
            input_hf = textwrap.dedent(f"""\
            {{hf,maxit={self.scf_max_iter},energy={self.scf_en_tol},accu={self.scf_rho_tol}
            start,2100.2
            orbital,2100.2
            wf,{int(molecule.nelec)},1,0
            }}

            """)
            input_molpro += input_hf

        # CASSCF Block: calculate energy option
        input_casscf = textwrap.dedent(f"""\
        {{mcscf,maxit={self.mcscf_max_iter},energy={self.mcscf_en_tol},gradient={self.mcscf_grad_tol},step={self.mcscf_step_tol}
        frozen,{self.frozen_orb}
        closed,{self.closed_orb}
        occ,{self.occ_orb}
        start,2140.2
        orbital,2140.2
        wf,{int(molecule.nelec)},1,0
        state,{molecule.nst}
        """)
        input_molpro += input_casscf

        # CASSCF Block: calculate gradient option
        casscf_grad = "\n".join([f"cpmcscf,grad,{ist + 1:d}.1,spin=0,accu={self.cpscf_grad_tol},record={5100.3 + float(ist):6.1f}" for ist in bo_list])
        if (not calc_force_only and self.calc_coupling):
            input_casscf_grad = f"""{casscf_grad}"""
        else:
            input_casscf_grad = f"""{casscf_grad}\n}}\n"""
        input_molpro += input_casscf_grad

        # CASSCF Block: calculate NAC option
        if (not calc_force_only and self.calc_coupling):
            kst = 0
            tmp_ind = 5100.3 + float(molecule.nst) - 1
            casscf_nac = ""
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    kst += 1
                    casscf_nac += f"\ncpmcscf,nacm,{ist + 1:d}.1,{jst + 1:d}.1,spin=0,accu={self.cpscf_grad_tol},record={tmp_ind + float(kst):6.1f}"
            input_casscf_nac = f"""{casscf_nac}\n}}\n"""
            input_molpro += input_casscf_nac

        # CASSCF Block: print energy option
        casscf_print = "".join([f"en({ist + 1}) = ENERGY({ist + 1})\n" for ist in range(molecule.nst)])
        input_casscf_print = f"""{casscf_print}\n"""
        input_molpro += input_casscf_print

        # CASSCF Block: print force option
        casscf_force = "".join([f"\n{{force\nsamc,{5100.3 + float(ist):6.1f}\n}}\n" for ist in bo_list])
        if (not calc_force_only and self.calc_coupling):
            kst = 0
            tmp_ind = 5100.3 + float(molecule.nst) - 1
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    kst += 1
                    casscf_force += f"\n{{force\nsamc,{tmp_ind + float(kst):6.1f}\n}}\n"
        input_casscf_force = f"""{casscf_force}\n"""
        input_molpro += input_casscf_force

        # Write 'molpro.inp' file
        file_name = "molpro.inp"
        with open(file_name, "w") as f:
            f.write(input_molpro)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run CASSCF calculation and save the output files to qm_log directory

            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Run Molpro method
        qm_command = os.path.join(self.qm_path, "molpro")
        # OpenMP setting
        os.environ["OMP_NUM_THREADS"] = "1"
        command = f"{qm_command} -m {self.memory} -I int -W wfu --no-xml-output -d int -o log -g -s molpro.inp > tmp_log"
        os.system(command)
        os.remove("tmp_log")
        # Copy the output file to 'qm_log' directory
        tmp_dir = os.path.join(base_dir, "qm_log")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_QM(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Read 'log' file
        file_name = "log"
        with open(file_name, "r") as f:
            log_out = f.read()

        # Energy
        if (not calc_force_only):
            tmp_e = 'SETTING EN\(\d+\)\s+[=]\s+([-]\S+)\s+HARTREE'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy, dtype=np.float64)
            for ist in range(molecule.nst):
                molecule.states[ist].energy = energy[ist]

        # Force
        for ist in bo_list:
            tmp_f = f'SA-MC GRADIENT FOR STATE {ist + 1:d}.1\n\n' + \
                '\s+Atom\s+dE\/dx\s+dE\/dy\s+dE\/dz\n\n' + \
                '\s+\d+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat_qm
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0], dtype=np.float64)
            force = force.reshape(molecule.nat_qm, 3, order='C')
            molecule.states[ist].force = - np.copy(force)

        # NAC
        if (not calc_force_only and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    tmp_c = f'SA-MC NACME FOR STATES {ist + 1:d}.1 - {jst + 1:d}.1\n\n' + \
                        '\s+Atom\s+dE\/dx\s+dE\/dy\s+dE\/dz\n\n' + \
                        '\s+\d+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat_qm
                    nac = re.findall(tmp_c, log_out)
                    nac = np.array(nac[0], dtype=np.float64)
                    nac = nac.reshape(molecule.nat_qm, 3, order='C')
                    molecule.nac[ist, jst] = nac
                    molecule.nac[jst, ist] = - nac


