from __future__ import division
import os, shutil, re
from bo.molpro.molpro import Molpro
import numpy as np
import textwrap

class CASSCF(Molpro):
    """ Class for Molpro CASSCF method
                 bomd | sh | eh | nac | re_calc
        CASSCF :  o     o    o     T      T
    """
    def __init__(self, molecule, basis_set="sto-3g", memory="500m", \
        max_iter=20, scf_en_tol=1E-8, scf_grad_tol=1E-6, scf_step_tol=1E-2, \
        active_elec=2, active_orb=2, cpscf_grad_tol=1E-7, \
        qm_path="./", nthreads=1, version=2015.1):
        # Initialize Molpro common variables
        super().__init__(basis_set, memory, qm_path, nthreads, version)

        # Initialize Molpro CASSCF variables
        # note that Molpro do not need restart option since it automatically search
        # the MO files for restart if possible
        # in addition, Molpro do not provide periodic setting with CASSCF method
        self.max_iter = max_iter
        self.scf_en_tol = scf_en_tol
        self.scf_grad_tol = scf_grad_tol
        self.scf_step_tol = scf_step_tol
        self.active_elec = active_elec
        self.active_orb = active_orb
        self.cpscf_grad_tol = cpscf_grad_tol

        # casscf calculation do not provide parallel computation
        # TODO: if your system provide parallel casscf, then this part should be removed
        if (self.nthreads > 1):
            raise ValueError ("Parallel CASSCF Not Implemented")

        # calculate number of frozen, closed and occ orbitals in CASSCF method
        # no positive frozen core orbitals in CASSCF
        self.frozen_orb = 0
        self.closed_orb = int((int(molecule.nelec) - self.active_elec) / 2)
        self.occ_orb = int(self.closed_orb + self.active_orb)

        # set 'l_nacme' with respect to the computational method
        # CASSCF can produce NACs, so we do not need to get NACME from CIoverlap
        # CASSCF can compute the gradient of several states simultaneously,
        #        but self.re_calc is set to be true to reduce cost.
        molecule.l_nacme = False
        self.re_calc = True

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Get/Extract BO information from Molpro
        """
        super().get_bo(base_dir, calc_force_only)
        self.write_xyz(molecule)
        self.get_input(molecule, bo_list, calc_force_only)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule, bo_list, calc_force_only):
        """ Generate Molpro input files: molpro.inp
        """
        # make 'molpro.inp' file
        input_molpro = ""

        # Scratch Block
        input_scr = textwrap.dedent(f"""\
        file,1,int.int,delete
        file,2,wf.wfu,unknown
        restart,2

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

        # CASSCF Block: calculation energy option
        input_casscf = textwrap.dedent(f"""\
        {{mcscf,maxit={self.max_iter},energy={self.scf_en_tol},gradient={self.scf_grad_tol},step={self.scf_step_tol}
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

        # CASSCF Block: calculate nac option
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

        # write 'molpro.inp' file
        file_name = "molpro.inp"
        with open(file_name, "w") as f:
            f.write(input_molpro)

    def run_QM(self, base_dir, istep, bo_list):
        """ run Molpro calculation and save the output files
        """
        # run DFTB+ method
        qm_command = os.path.join(self.qm_path, "molpro")
        # openmp setting
        os.environ["OMP_NUM_THREADS"] = "1"
        command = f"{qm_command} -m {self.memory} -I int -W wfu --no-xml-output -d int -o log -g -s molpro.inp > tmp_log"
        os.system(command)
        os.remove("tmp_log")
        # copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ read the output files to get BO data
        """
        # read 'log' file
        file_name = "log"
        with open(file_name, "r") as f:
            log_out = f.read()

        # energy
        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            tmp_e = 'SETTING EN\(\d+\)\s+[=]\s+([-]\S+)\s+HARTREE'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy)
            energy = energy.astype(float)
            for ist in range(molecule.nst):
                molecule.states[ist].energy = energy[ist]

        # force
        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        for ist in bo_list:
            tmp_f = f'SA-MC GRADIENT FOR STATE {ist + 1:d}.1\n\n' + \
                '\s+Atom\s+dE\/dx\s+dE\/dy\s+dE\/dz\n\n' + \
                '\s+\d+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, 3, order='C')
            molecule.states[ist].force = - np.copy(force)

        # NAC
        if (not calc_force_only and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.nac[ist, jst, :, :] = 0.
                    elif (ist < jst):
                        tmp_c = f'SA-MC NACME FOR STATES {ist + 1:d}.1 - {jst + 1:d}.1\n\n' + \
                            '\s+Atom\s+dE\/dx\s+dE\/dy\s+dE\/dz\n\n' + \
                            '\s+\d+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\n' * molecule.nat
                        nac = re.findall(tmp_c, log_out)
                        nac = np.array(nac[0])
                        nac = nac.astype(float)
                        nac = nac.reshape(molecule.nat, 3, order='C')
                        molecule.nac[ist, jst] = np.copy(nac)
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]




