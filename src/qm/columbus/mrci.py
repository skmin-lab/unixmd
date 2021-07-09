from __future__ import division
from qm.columbus.columbus import Columbus
from misc import data, amu_to_au, call_name
import os, shutil, re
import numpy as np

class MRCI(Columbus):
    """ Class for MRCI method of Columbus program

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param integer memory: Allocatable memory in the calculations
        :param string guess: Initial guess for MRCI method
        :param string guess_file: Initial guess file
        :param integer scf_en_tol: Energy convergence for SCF iterations
        :param integer scf_max_iter: Maximum number of SCF iterations
        :param integer mcscf_en_tol: Energy convergence for (SA-)CASSCF iterations
        :param integer mcscf_max_iter: Maximum number of (SA-)CASSCF iterations
        :param integer mrci_en_tol: Energy convergence for MRCI iterations
        :param integer mrci_max_iter: Maximum number of MRCI iterations
        :param integer state_avg: Number of states to be averaged for (SA-)CASSCF and MRCI
        :param integer active_elec: Number of electrons in active space
        :param integer active_orb: Number of orbitals in active space
        :param integer frozen_core_orb: Number of frozen core orbitals in doubly occupied space
        :param integer frozen_virt_orb: Number of frozen virtual orbitals from the highest unoccupied space
        :param integer cpscf_grad_tol: Gradient tolerance for CP-MRCI equations
        :param integer cpscf_max_iter: Maximum number of iterations for CP-MRCI equations
        :param string qm_path: Path for QM binary
        :param string version: Version of Columbus program
    """
    def __init__(self, molecule, basis_set="6-31g*", memory=500, \
        guess="hf", guess_file="./mocoef", scf_en_tol=9, scf_max_iter=40, \
        mcscf_en_tol=8, mcscf_max_iter=100, mrci_en_tol=4, mrci_max_iter=None, \
        state_avg=None, active_elec=2, active_orb=2, frozen_core_orb=0, frozen_virt_orb=0, \
        cpscf_grad_tol=6, cpscf_max_iter=100, qm_path="./", version="7.0"):
        # Initialize Columbus common variables
        super(MRCI, self).__init__(molecule, basis_set, memory, qm_path, version)

        # Initialize Columbus MRCI variables
        # Set initial guess for MRCI calculation
        # read: Read MCSCF orbitals obtained from previous step -> MCSCF -> MRCI
        # hf: Start HF orbitals -> MCSCF -> MRCI
        self.guess = guess.lower()
        self.guess_file = guess_file
        if not (self.guess in ["hf", "read"]):
            error_message = "Invalid initial guess for MRCI!"
            error_vars = f"guess = {self.guess}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # HF calculation for initial guess of CASSCF or MRCI calculations
        self.scf_en_tol = scf_en_tol
        self.scf_max_iter = scf_max_iter

        # CASSCF/MRCI calculation
        # TODO : Currently, same active space is applied to the CASSCF and MRCI calculations
        self.active_elec = active_elec
        self.active_orb = active_orb

        self.state_avg = state_avg
        if (self.state_avg == None):
            # Set number of state-averaging to number of states
            self.state_avg = molecule.nst
        else:
            if (self.state_avg < molecule.nst):
                error_message = "Number of state-averaging must be equal or larger than number of states!"
                error_vars = f"state_avg = {self.state_avg}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # CASSCF calculation
        self.mcscf_en_tol = mcscf_en_tol
        self.mcscf_max_iter = mcscf_max_iter

        # MRCI calculation
        self.frozen_core_orb = frozen_core_orb
        self.frozen_virt_orb = frozen_virt_orb

        self.mrci_en_tol = mrci_en_tol
        self.mrci_max_iter = mrci_max_iter
        if (self.mrci_max_iter == None):
            # Set maximum number of MRCI iterations
            self.mrci_max_iter = 30 * self.state_avg

        self.cpscf_grad_tol = cpscf_grad_tol
        self.cpscf_max_iter = cpscf_max_iter

        # Calculate number of doubly occuplied, closed, and internal orbitals in HF, CASSCF, and MRCI method
        # Note that there is positive frozen core orbitals in MRCI
        #           and the number of auxiliary internal orbitals is set to zero
        self.docc_orb = int(int(molecule.nelec) / 2)
        self.closed_orb = int((int(molecule.nelec) - self.active_elec) / 2)
        self.internal_orb = self.closed_orb + self.active_orb + 0 - self.frozen_core_orb

        # Check the closed shell for systems
        if (not int(molecule.nelec) % 2 == 0):
            error_message = "Only closed shell configuration supported, check charge!"
            error_vars = f"Molecule.nelec = {int(molecule.nelec)}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set 'l_nacme' with respect to the computational method
        # MRCI can produce NACs, so we do not need to get NACME from CIoverlap
        # MRCI can compute the gradient of several states simultaneously,
        #      but self.re_calc is set to be true to reduce cost.
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
            # The MCSCF orbitals from previous step are used as initial guess
            shutil.copy(os.path.join(self.scr_qm_dir, "./MOCOEFS/mocoef_mc.sp"), \
                os.path.join(self.scr_qm_dir, "../mocoef"))

    def get_input(self, molecule, istep, bo_list, calc_force_only):
        """ Generate Columbus input files: geom, prepin, stdin, mcscfin, transmomin, etc

            :param object molecule: Molecule object
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Generate 'geom' file used in Columbus
        geom = ""
        for iat in range(molecule.nat_qm):
            atom_num = list(data.keys()).index(f"{molecule.symbols[iat]}")
            tmp_atom = f' {molecule.symbols[iat]:5s}{atom_num:7.2f}' \
                + "".join([f'{molecule.pos[iat, isp]:15.8f}' for isp in range(molecule.ndim)]) \
                + f'{molecule.mass[iat] / amu_to_au:15.8f}' + "\n"
            geom += tmp_atom

        file_name = "geom"
        with open(file_name, "w") as f:
            f.write(geom)

        # Scratch Block
        hf = False
        mcscf = True
        if (self.guess == "read"):
            if (istep == -1):
                if (os.path.isfile(self.guess_file)):
                    # Copy guess file to currect directory
                    shutil.copy(self.guess_file, os.path.join(self.scr_qm_dir, "mocoef"))
                    restart = 1
                else:
                    restart = 0
                    hf = True
            elif (istep >= 0):
                # Move previous file to currect directory
                os.rename("../mocoef", os.path.join(self.scr_qm_dir, "mocoef"))
                restart = 1
        elif (self.guess == "hf"):
            restart = 0
            hf = True

        if (calc_force_only):
            restart = 1
            hf = False
            mcscf = False
            shutil.copy(os.path.join(self.scr_qm_dir, "./MOCOEFS/mocoef_mc.sp"), \
                os.path.join(self.scr_qm_dir, "mocoef"))

        # Generate new prepinp script
        shutil.copy(os.path.join(self.qm_path, "prepinp"), "prepinp_copy")

        file_name = "prepinp_copy"
        with open(file_name, "r") as f:
            prepinp = f.read()
            prepinp = prepinp.replace("( keys %sumformula )", "(sort keys %sumformula )", 1)

        file_name = "prepinp_fix"
        with open(file_name, "w") as f:
            f.write(prepinp)
        os.chmod("prepinp_fix", 0o755)

        # Generate 'prepin' file used in prepinp script of Columbus
        if (calc_force_only):
            # Here, y means overwritting of 'inpcol' file
            prepin = "1\ny\nc1\ngeom\n\n"
        else:
            prepin = "1\nc1\ngeom\n\n"

        prepin += self.basis_nums
        prepin += "\ny\n\n"

        file_name = "prepin"
        with open(file_name, "w") as f:
            f.write(prepin)

        os.system("./prepinp_fix < prepin > prepout")

        # Generate 'stdin' file used in colinp script of Columbus
        # DALTON, SCF input setting in colinp script of Columbus
        if (hf):
            stdin = f"\ny\n1\nn\nno\n2\nyes\n{self.docc_orb}\nyes\nyes\n{self.scf_max_iter}\n{self.scf_en_tol}\nno\n1\n\n"
        else:
            # Here, n in DALTON means no change of basis set; this choice appears when 'control.run' or 'mocoef' file exists
            stdin = f"\ny\n1\nn\nno\nn\n"

        # MCSCF input setting in colinp script of Columbus
        if (mcscf):
            stdin += f"3\nn\n1\n1\n{int(molecule.nelec)}\n1\n1\n0\n0\n{self.closed_orb}\n{self.active_orb}\nn\n" + "\t" * 14 + "\n"

        # MRCI input setting in colinp script of Columbus
        if (not calc_force_only):
            stdin += f"4\n\n2\ny\nn\n1\n{int(molecule.nelec)}\n1\n{self.frozen_core_orb}\n{self.frozen_virt_orb}\n" + \
                f"{self.internal_orb}\n{self.internal_orb - self.active_orb}\n0\n2\ny\n\nn\n1\n" + "\t" * 17 + "\n"

        # Job control setting in colinp script of Columbus
        if (calc_force_only):
            # Start from 'mocoef' file
            # Here, y in job control setting means discard of already existing 'control.run' file
            stdin += "5\n1\ny\n1\n3\n5\n11\n1\nn\n3\nn\n8\n4\n7\n\n"
        else:
            if (restart == 1):
                # Start from MCSCF calculation
                stdin += "5\n1\n1\n3\n5\n11\n1\nn\n3\nn\n8\n4\n7\n\n"
            else:
                # Start from SCF calculation
                stdin += "5\n1\n1\n2\n3\n5\n11\n1\nn\n3\nn\n8\n4\n7\n\n"

        file_name = "stdin"
        with open(file_name, "w") as f:
            f.write(stdin)

        os.system(f"{self.qm_path}/colinp < stdin > stdout")

        if (not calc_force_only):

            # Manually modify input files
            # Modify 'mcscfin' files
            file_name = "mcscfin"
            with open(file_name, "r") as f:
                mcscfin = f.readlines()

            mcscf_length = len(mcscfin)
            target_line = mcscf_length - 3

            new_mcscf = ""
            for i in range(target_line):
                if ("niter" in mcscfin[i]):
                    new_mcscf += f"  niter={self.mcscf_max_iter},\n"
                elif ("tol(1)" in mcscfin[i]):
                    new_mcscf += f"  tol(1)=1.e-{self.mcscf_en_tol},\n"
                else:
                    new_mcscf += mcscfin[i]

            new_mcscf += f"  NAVST(1) = {self.state_avg},\n"
            for i in range(self.state_avg):
                new_mcscf += f"  WAVST(1,{i + 1})=1 ,\n"
            new_mcscf += " &end\n"

            os.rename("mcscfin", "mcscfin.old")

            file_name = "mcscfin"
            with open(file_name, "w") as f:
                f.write(new_mcscf)

            # Manually modify input files
            # Modify 'ciudgin' files
            file_name = "ciudgin"
            with open(file_name, "r") as f:
                ciudgin = f.readlines()

            ciudg_length = len(ciudgin)

            # For 'NVBKMN', 'NVCIMN', 'NVCIMX', 'NVRFMX', 'NVBKMX', these variables affect
            # the iteration process, so the default values are used.
            new_ciudg = ""
            for i in range(ciudg_length):
                if ("NROOT" in ciudgin[i]):
                    new_ciudg += f" NROOT = {self.state_avg}\n"
                elif ("NVBKMN" in ciudgin[i]):
                    new_ciudg += f" NVBKMN = {self.state_avg}\n"
                elif ("RTOLBK" in ciudgin[i]):
                    new_ciudg += " RTOLBK = "
                    new_ciudg += "1e-4," * self.state_avg
                    new_ciudg += "\n"
                elif ("NITER" in ciudgin[i]):
                    new_ciudg += f" NITER = {self.mrci_max_iter}\n"
                elif ("NVCIMN" in ciudgin[i]):
                    new_ciudg += f" NVCIMN = {self.state_avg + 2}\n"
                elif ("RTOLCI" in ciudgin[i]):
                    new_ciudg += " RTOLCI = "
                    new_ciudg += f"1e-{self.mrci_en_tol}," * self.state_avg
                    new_ciudg += "\n"
                elif ("NVCIMX" in ciudgin[i]):
                    new_ciudg += f" NVCIMX = {self.state_avg + 5}\n"
                elif ("NVRFMX" in ciudgin[i]):
                    new_ciudg += f" NVRFMX = {self.state_avg + 5}\n"
                elif ("NVBKMX" in ciudgin[i]):
                    new_ciudg += f" NVBKMX = {self.state_avg + 5}\n"
                else:
                    new_ciudg += ciudgin[i]

            os.rename("ciudgin", "ciudgin.old")

            file_name = "ciudgin"
            with open(file_name, "w") as f:
                    f.write(new_ciudg)

            # Modify 'cidenin' files
            file_name = "cidenin"
            with open(file_name, "r") as f:
                cidenin = f.readlines()

            ciden_length = len(cidenin)

            new_ciden = ""
            for i in range(ciden_length):
                if ("lroot" in cidenin[i]):
                    new_ciden += f"  lroot = {self.state_avg}\n"
                else:
                    new_ciden += cidenin[i]

            os.rename("cidenin", "cidenin.old")

            file_name = "cidenin"
            with open(file_name, "w") as f:
                    f.write(new_ciden)

            # Manually modify input files
            # Modify 'cigrdin' files
            file_name = "cigrdin"
            with open(file_name, "r") as f:
                cigrdin = f.readlines()

            new_cigrd = ""
            for line in cigrdin:
                if ("nmiter" in line):
                    new_cigrd += f" nmiter= {self.cpscf_max_iter}, print=0, fresdd=1,\n"
                elif ("rtol" in line):
                    new_cigrd += f" rtol=1e-{self.cpscf_grad_tol}, dtol=1e-6,\n"
                else:
                    new_cigrd += line

            os.rename("cigrdin", "cigrdin.old")

            file_name = "cigrdin"
            with open(file_name, "w") as f:
                f.write(new_cigrd)

            # Copy 'daltcomm' files
            shutil.copy("daltcomm", "daltcomm.new")

        # Write 'transmomin' files
        transmomin = "CI\n"
        # Gradient part
        for ist in bo_list:
            transmomin += f"1  {ist + 1}  1  {ist + 1}\n"

        # NAC part
        if (not calc_force_only):
            for i in range(molecule.nst):
                for j in range(i + 1, molecule.nst):
                    transmomin += f"1  {i + 1}  1  {j + 1}\n"

        file_name = "transmomin"
        with open(file_name, "w") as f:
            f.write(transmomin)

    def run_QM(self, base_dir, istep, bo_list):
        """ Run MRCI calculation and save the output files to qm_log directory

            :param string base_dir: Base directory
            :param integer istep: Current MD step
            :param integer,list bo_list: List of BO states for BO calculation
        """
        # Run Columbus method
        qm_command = os.path.join(self.qm_path, "runc")
        command = f"{qm_command} -m {self.memory} > runls"
        os.system(command)
        # Copy the output file to 'qm_log' directory
        tmp_dir = os.path.join(base_dir, "qm_log")
        if (os.path.exists(tmp_dir)):
            log_step = f"runls.{istep + 1}.{bo_list[0]}"
            shutil.copy("runls", os.path.join(tmp_dir, log_step))
        # Remove scratch 'WORK' directory
        tmp_dir = os.path.join(self.scr_qm_dir, "WORK")
        if (os.path.exists(tmp_dir)):
            shutil.rmtree(tmp_dir)

    def extract_QM(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Energy
        if (not calc_force_only):
            # Read 'ciudgsm.sp' file
            file_name = "LISTINGS/ciudgsm.sp"
            with open(file_name, "r") as f:
                log_out = f.read()

            tmp_e = ''
            for ist in range(self.state_avg):
                if (ist < molecule.nst):
                    tmp_e += ' mr-sdci [#]\s*\S+\s+\S+\s+([-]\S+)\s+[-]*\S+\s+[-]*\S+\s+[-]*\S+\s+[-]*\S+\s+\n'
                else:
                    tmp_e += ' mr-sdci [#]\s*\S+\s+\S+\s+[-]\S+\s+[-]*\S+\s+[-]*\S+\s+[-]*\S+\s+[-]*\S+\s+\n'
            tmp_e += '\n number of'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy[0], dtype=np.float64)
            for ist in range(molecule.nst):
                molecule.states[ist].energy = energy[ist]

        # Force
        for ist in bo_list:
            # Read 'cartgrd.drt1.state?.sp' file
            file_name = f"GRADIENTS/cartgrd.drt1.state{ist + 1}.sp"
            with open(file_name, "r") as f:
                log_out = f.read()
                log_out = log_out.replace("D", "E", molecule.nat_qm * molecule.ndim)

            tmp_f ='\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat_qm
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0], dtype=np.float64)
            force = force.reshape(molecule.nat_qm, 3, order='C')
            molecule.states[ist].force = - np.copy(force)

        # NAC
        if (not calc_force_only and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(ist + 1, molecule.nst):
                    # Read 'cartgrd.nad.drt1.state?.drt1.state?.sp' file
                    file_name = f"GRADIENTS/cartgrd.nad.drt1.state{ist + 1}.drt1.state{jst + 1}.sp"
                    with open(file_name, "r") as f:
                        log_out = f.read()

                    tmp_c =  '\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat_qm
                    nac = re.findall(tmp_c, log_out)
                    nac = np.array(nac[0], dtype=np.float64)
                    nac = nac.reshape(molecule.nat_qm, 3, order='C')
                    molecule.nac[ist, jst] = nac
                    molecule.nac[jst, ist] = - nac


