from __future__ import division
from misc import au_to_A
import os, shutil

class QM_calculator(object):
    """ Class for quantum mechanics calculator such as QM, ML, etc
    """
    def __init__(self):
        # Save name of QM calculator and its method
        self.qm_prog = str(self.__class__).split('.')[1]
        self.qm_method = self.__class__.__name__

    def get_data(self, base_dir, calc_force_only):
        """ Make scratch directory and copy geometry file

            :param string base_dir: Base directory
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Make 'scr_qm' directory
        unixmd_dir = os.path.join(base_dir, "md")
        self.scr_qm_dir = os.path.join(unixmd_dir, "scr_qm")
        if (not calc_force_only):
            if (os.path.exists(self.scr_qm_dir)):
                shutil.rmtree(self.scr_qm_dir)
            os.makedirs(self.scr_qm_dir)
        # Move to the scratch directory
        os.chdir(self.scr_qm_dir)

    def write_xyz(self, molecule):
        """ Make current geometry file

            :param object molecule: Molecule object
        """
        file_name = "geometry.xyz"
        with open(file_name, "w") as ftj:
            ftj.write(f"{molecule.nat_qm}\n\n")
            for iat in range(molecule.nat_qm):
                ftj.write(f"{molecule.symbols[iat]:4}")
                ftj.write("".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A]) + "\n")

    def move_dir(self, base_dir):
        """ Move to the base directory

            :param string base_dir: Base directory
        """
        os.chdir(base_dir)



