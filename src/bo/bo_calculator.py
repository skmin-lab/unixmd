from __future__ import division
import os, shutil
from misc import au_to_A

class BO_calculator(object):
    """ Class for Born-Oppenheimer calculator such as QM, ML, etc
    """
    def __init__(self):
        pass

    def get_bo(self, base_dir, calc_force_only):
        """ Make scratch directory and copy geomtry file
        """
        # make 'scr_qm' directory
        unixmd_dir = os.path.join(base_dir, "md")
        self.scr_qm_dir = os.path.join(unixmd_dir, "scr_qm")
        if (not calc_force_only):
            if (os.path.exists(self.scr_qm_dir)):
                shutil.rmtree(self.scr_qm_dir)
            os.makedirs(self.scr_qm_dir)
        # move to the scratch directory
        os.chdir(self.scr_qm_dir)

    def write_xyz(self, molecule):
        """ make current geomtry file
        """
        file_name = "geometry.xyz"
        with open(file_name, "w") as ftj:
            ftj.write(f"{molecule.nat}\n\n")
            for iat in range(molecule.nat):
                ftj.write(f"{molecule.symbols[iat]:4}")
                ftj.write("".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A]) + "\n")

    def move_dir(self, base_dir):
        """ Move to the base directory
        """
        os.chdir(base_dir)



