from __future__ import division
from misc import au_to_A
import os, shutil

class MM_calculator(object):
    """ Class for molecular mechanics calculator such as MM, model, etc
    """
    def __init__(self):
        # Save name of MM calculator
        self.mm_prog = self.__class__.__name__

    def get_mm(self, base_dir):
        """ Make scratch directory and copy geometry file

            :param string base_dir: base directory
        """
        # Make 'scr_mm' directory
        unixmd_dir = os.path.join(base_dir, "md")
        self.scr_mm_dir = os.path.join(unixmd_dir, "scr_mm")
        if (os.path.exists(self.scr_mm_dir)):
            shutil.rmtree(self.scr_mm_dir)
        os.makedirs(self.scr_mm_dir)
        # Move to the scratch directory
        os.chdir(self.scr_mm_dir)

    def write_xyz(self, molecule):
        """ Make current geometry file

            :param object molecule: molecule object
        """
        # TODO : pass a choice for QMMM scheme as argument
        file_name = "geometry.xyz"
        with open(file_name, "w") as ftj:
            ftj.write(f"{molecule.nat}\n\n")
            for iat in range(molecule.nat):
                ftj.write(f"{molecule.symbols[iat]:4}")
                ftj.write("".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A]) + "\n")

    def move_dir(self, base_dir):
        """ Move to the base directory

            :param string base_dir: base directory
        """
        os.chdir(base_dir)


