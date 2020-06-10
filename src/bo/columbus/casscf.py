from __future__ import division
from bo.columbus.columbus import Columbus
from bo.columbus.colbasis import *
from misc import data, au_to_A, A_to_au, amu_to_au
import os, shutil, re, textwrap
import numpy as np

class CASSCF(Columbus):
    """ Class for CASSCF method of Columbus program

        :param object molecule: molecule object
        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param integer active_elec: number of electrons in active space
        :param integer active_orb: number of orbitals in active space
        :param string qm_path: path for QM binary
        :param integer nthreads: number of threads in the calculations
        :param double version: version of Columbus program
    """
    def __init__(self, molecule, basis_set="6-31g*", memory="500", \
        active_elec=2, active_orb=2, qm_path="./", nthreads=1, version=7.0):
        # Initialize Columbus common variables
        super(CASSCF, self).__init__(molecule, basis_set, memory, qm_path, nthreads, version)

        # Initialize Columbus CASSCF variables
        # TODO : restart option? from mocoef file
#        self.max_iter = max_iter
#        self.scf_en_tol = scf_en_tol
#        self.scf_grad_tol = scf_grad_tol
#        self.scf_step_tol = scf_step_tol
        self.active_elec = active_elec
        self.active_orb = active_orb
#        self.cpscf_grad_tol = cpscf_grad_tol

        # CASSCF calculation do not support parallel computation
        if (self.nthreads > 1):
            raise ValueError ("Parallel CASSCF Not Implemented")

        # Calculate number of frozen, closed and occ orbitals in CASSCF method
        # Note that there is no positive frozen core orbitals in CASSCF
        self.frozen_orb = 0
        self.closed_orb = int((int(molecule.nelec) - self.active_elec) / 2)
        self.docc_orb = int(int(molecule.nelec) / 2)

        # Check the closed shell for systems
        if (not int(molecule.nelec) % 2 == 0):
            raise ValueError ("Only closed shell configuration Implemented")

        # Set 'l_nacme' with respect to the computational method
        # CASSCF can produce NACs, so we do not need to get NACME from CIoverlap
        # CASSCF can compute the gradient of several states simultaneously,
        #        but self.re_calc is set to be true to reduce cost.
        molecule.l_nacme = False
        self.re_calc = True

    def get_bo(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
        """ Extract energy, gradient and nonadiabatic couplings from CASSCF method

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, calc_force_only):
        """ Generate Columbus input files: geom, prepin, stdin, mcscfin, transmomin, etc

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Generate 'geom' file used in Columbus
        geom = ""
        for iat in range(molecule.nat):
            atom_num = list(data.keys()).index(f"{molecule.symbols[iat]}")
            tmp_atom = f' {molecule.symbols[iat]:5s}{atom_num:7.2f}' \
                + "".join([f'{molecule.pos[iat, isp]:15.8f}' for isp in range(molecule.nsp)]) \
                + f'{molecule.mass[iat] / amu_to_au:15.8f}' + "\n"
            geom += tmp_atom

        file_name = "geom"
        with open(file_name, "w") as f:
            f.write(geom)

        # Set basis sets information
        # TODO : move to columbus.py, this is common part
        if (self.basis_set == "cc-pvdz"):
            tmp_basis = "\n".join([f"{cc_pvdz[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "cc-pvtz"):
            tmp_basis = "\n".join([f"{cc_pvtz[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "cc-pvqz"):
            tmp_basis = "\n".join([f"{cc_pvqz[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "3-21g*"):
            tmp_basis = "\n".join([f"{t_21gs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "3-21+g*"):
            tmp_basis = "\n".join([f"{t_21pgs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-31g"):
            tmp_basis = "\n".join([f"{s_31g[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-31g*"):
            tmp_basis = "\n".join([f"{s_31gs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-31+g*"):
            tmp_basis = "\n".join([f"{s_31pgs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-311g*"):
            tmp_basis = "\n".join([f"{s_311gs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-311g+*"):
            tmp_basis = "\n".join([f"{s_311pgs[f'{itype}']}" for itype in self.atom_type])
        else:
            raise ValueError ("Basis set not yet implemented in Columbus input: add manually (colbasis.py)")

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
        prepin = "1\nc1\ngeom\n\n"
        prepin += tmp_basis
        prepin += "\ny\n\n"

        file_name = "prepin"
        with open(file_name, "w") as f:
            f.write(prepin)

        os.system("./prepinp_fix < prepin > prepout")

        # Generate 'stdin' file used in colinp script of Columbus
        # DALTON, SCF input setting in colinp script of Columbus
        stdin = f"\ny\n1\nn\nno\n2\nyes\n{self.docc_orb}\nyes\nno\n\n"
        # MCSCF input setting in colinp script of Columbus
        stdin += f"3\nn\n3\n1\n{int(molecule.nelec)}\n1\n1\n0\n0\n{self.closed_orb}\n{self.active_orb}\nn\n" + "\t" * 14 + "\n"
        # Job control setting in colinp script of Columbus
        if (calc_force_only):
            # Start from 'mocoef' file
            stdin += "5\n1\n1\n3\n11\n1\nn\n3\ny\n8\n4\n7\n\n"
        else:
            # Start from SCF calculation
            stdin += "5\n1\n1\n2\n3\n11\n1\nn\n3\ny\n8\n4\n7\n\n"

        file_name = "stdin"
        with open(file_name, "w") as f:
            f.write(stdin)

        os.system(f"{self.qm_path}/colinp < stdin > stdout")

        # Manually modify input files
        # Modify 'mcscfin' files
        file_name = "mcscfin"
        with open(file_name, "r") as f:
            mcscfin = f.readlines()

        mcscf_length = len(mcscfin)
        target_line = mcscf_length - 3

        new_mcscf = ""
        for i in range (target_line):
            new_mcscf += mcscfin[i]
        new_mcscf += f"  NAVST(1) = {molecule.nst},\n"
        for i in range (molecule.nst):
            new_mcscf += f"  WAVST(1,{i + 1})=1 ,\n"
        new_mcscf += " &end\n"

        os.rename("mcscfin", "mcscfin.old")

        file_name = "mcscfin"
        with open(file_name, "w") as f:
            f.write(new_mcscf)

        # Modify 'transmomin' files
        transmomin = "MCSCF\n"
        # Gradient part
        for ist in bo_list:
            transmomin += f"1  {ist + 1}  1  {ist + 1}\n"

        # NAC part
        if (not calc_force_only):
            for i in range (molecule.nst):
                for j in range (i):
                    transmomin += f"1  {i + 1}  1  {j + 1}\n"

        file_name = "transmomin"
        with open(file_name, "w") as f:
            f.write(transmomin)

        # Copy 'daltcomm' files
        shutil.copy("daltcomm", "daltcomm.new")

        # Copy 'mocoef' file to scratch directory
        if (calc_force_only):
            shutil.copy("MOCOEFS/mocoef_mc.sp", "mocoef")

    def run_QM(self, base_dir, istep, bo_list):
        """ Run CASSCF calculation and save the output files to QMlog directory

            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
        """
        # Run Columbus method
        qm_command = os.path.join(self.qm_path, "runc")
        command = f"{qm_command} -m {self.memory} > runls"
        os.system(command)
        # Copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("runls", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get BO information

            :param object molecule: molecule object
            :param integer,list bo_list: list of BO states for BO calculation
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        # Energy
        if (not calc_force_only):
            # Read 'mcscfsm.sp' file
            file_name = "LISTINGS/mcscfsm.sp"
            with open(file_name, "r") as f:
                log_out = f.read()

            for states in molecule.states:
                states.energy = 0.

            tmp_e = 'total\senergy[=]\s*([-]\S+)[,]'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy)
            energy = energy.astype(float)
            for ist in range(molecule.nst):
                molecule.states[ist].energy = energy[ist]

        # Force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        for ist in bo_list:
            # Read 'cartgrd.drt1.state?.sp' file
            file_name = f"GRADIENTS/cartgrd.drt1.state{ist + 1}.sp"
            with open(file_name, "r") as f:
                log_out = f.read()
                log_out = log_out.replace("D", "E", molecule.nat * molecule.nsp)

            tmp_f ='\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[bo_list[0]].force = - np.copy(force)

        # NAC
        if (not calc_force_only and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.nac[ist, jst, :, :] = 0.
                    elif (ist < jst):
                        # Read 'cartgrd.nad.drt1.state?.drt1.state?.sp' file
                        file_name = f"GRADIENTS/cartgrd.nad.drt1.state{jst + 1}.drt1.state{ist + 1}.sp"
                        with open(file_name, "r") as f:
                            log_out = f.read()
                            log_out = log_out.replace("D", "E", molecule.nat * molecule.nsp)

                        tmp_c =  '\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat
                        nac = re.findall(tmp_c, log_out)
                        nac = np.array(nac[0])
                        nac = nac.astype(float)
                        nac = nac.reshape(molecule.nat, 3, order='C')
                        molecule.nac[ist, jst] = np.copy(nac)
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]


