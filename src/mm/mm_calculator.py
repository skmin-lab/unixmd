from __future__ import division
import os, shutil

class MM_calculator(object):
    """ Class for molecular mechanics calculator such as MM, model, etc
    """
    def __init__(self):
        # Save name of MM calculator
        self.mm_prog = self.__class__.__name__

    def get_data(self, base_dir, calc_force_only):
        """ Make scratch directory and copy geometry file

            :param string base_dir: Base directory
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Make 'scr_mm' directory
        unixmd_dir = os.path.join(base_dir, "md")
        self.scr_mm_dir = os.path.join(unixmd_dir, "scr_mm")
        if (not calc_force_only):
            if (os.path.exists(self.scr_mm_dir)):
                shutil.rmtree(self.scr_mm_dir)
            os.makedirs(self.scr_mm_dir)
        # Move to the scratch directory
        os.chdir(self.scr_mm_dir)

    def move_dir(self, base_dir):
        """ Move to the base directory

            :param string base_dir: Base directory
        """
        os.chdir(base_dir)


