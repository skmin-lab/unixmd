from __future__ import division
from bo.model.model import Model
import os, shutil, re
import numpy as np

class Shin_Metiu(Model):
    """ Class for shin-metiu model calculation
        you can extract 'ENERGY.DAT', 'FORCE.DAT', and 'NAC.DAT' (or 'NACME.DAT')
        or directly set the values to the variables such as molecule.states.energy
                     bomd | sh | eh | nac | re_calc
        shin-metiu :  o     o    o     T      F
    """
    def __init__(self, molecule, qm_path="./"):
        # Initialize model common variables
        super().__init__(qm_path)

        # set 'l_nacme' with respect to the computational method
        # shin-metiu model can produce NACs, so we do not need to get NACME
        molecule.l_nacme = False

        # shin-metiu model can compute the gradient of several states simultaneously,
        self.re_calc = False

    def get_bo(self, molecule, base_dir, istep, bo_list, calc_force_only):
        """ Get/Extract BO information from shin-metiu model
        """
        super().get_bo(base_dir, calc_force_only)
        self.get_input(molecule)
        self.run_QM(base_dir, istep, bo_list)
        self.extract_BO(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule):
        """ Generate shin-metiu input files: metiu.in
        """
        # make 'metiu.in' file
        input_shin_metiu = ""
        input_shin_metiu += f"metiu\n{molecule.pos[0, 0]:15.8f}\n{molecule.nst}"

        # write 'metiu.in' file
        file_name = "metiu.in"
        with open(file_name, "w") as f:
            f.write(input_shin_metiu)

    def run_QM(self, base_dir, istep, bo_list):
        """ run shin-metiu calculation and save the output files
        """
        # run shin-metiu model
        qm_command = os.path.join(self.qm_path, "metiu.x")
        command = f"{qm_command} < metiu.in > metiu.log"
        os.system(command)
        # copy the output file to 'QMlog' directory
        tmp_dir = os.path.join(base_dir, "QMlog")
        if (os.path.exists(tmp_dir)):
            log_step = f"log.{istep + 1}.{bo_list[0]}"
            shutil.copy("metiu.log", os.path.join(tmp_dir, log_step))

    def extract_BO(self, molecule, bo_list, calc_force_only):
        """ read the output files to get BO data
        """
        # energy
        file_name = "ENERGY.DAT"
        with open(file_name, "r") as f:
            log_out = f.read()

        if (not calc_force_only):
            for states in molecule.states:
                states.energy = 0.

            tmp_e = '([-]*\S+)'
            energy = re.findall(tmp_e, log_out)
            energy = np.array(energy)
            energy = energy.astype(float)
            for ist in range(molecule.nst):
                molecule.states[ist].energy = energy[ist]

        # force
        file_name = "FORCE.DAT"
        with open(file_name, "r") as f:
            log_out = f.read()

        if (not calc_force_only):
            for states in molecule.states:
                states.force = np.zeros((molecule.nat, molecule.nsp))

        for ist in bo_list:
            tmp_f = f'{ist + 1:d}\n\s*([-]*\S+E[-]*\S+)'
            force = re.findall(tmp_f, log_out)
            force = np.array(force[0])
            force = force.astype(float)
            force = force.reshape(molecule.nat, molecule.nsp, order='C')
            molecule.states[ist].force = np.copy(force)

        # NAC
        file_name = "NAC.DAT"
        with open(file_name, "r") as f:
            log_out = f.read()

        if (not calc_force_only and self.calc_coupling):
            for ist in range(molecule.nst):
                for jst in range(molecule.nst):
                    if (ist == jst):
                        molecule.nac[ist, jst, :, :] = 0.
                    elif (ist < jst):
                        tmp_c = f'{ist + 1:d} {jst + 1:d}\n\s*([-]*\S+)'
                        nac = re.findall(tmp_c, log_out)
                        nac = np.array(nac[0])
                        nac = nac.astype(float)
                        nac = nac.reshape(molecule.nat, molecule.nsp, order='C')
                        molecule.nac[ist, jst] = np.copy(nac)
                    else:
                        molecule.nac[ist, jst] = - molecule.nac[jst, ist]



