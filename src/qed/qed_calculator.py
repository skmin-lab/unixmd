from __future__ import division
import os, shutil

class QED_calculator(object):
    """ Class for (cavity) quantum electrodynamics calculator
    """
    def __init__(self):
        # Save name of QED calculator
        self.qed_method = self.__class__.__name__

    def get_data(self, base_dir, calc_force_only):
        """ Make scratch directory

            :param string base_dir: Base directory
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        # Make 'scr_qed' directory
        unixmd_dir = os.path.join(base_dir, "md")
        self.scr_qed_dir = os.path.join(unixmd_dir, "scr_qed")
        if (not calc_force_only):
            if (os.path.exists(self.scr_qed_dir)):
                shutil.rmtree(self.scr_qed_dir)
            os.makedirs(self.scr_qed_dir)
        # Move to the scratch directory
        os.chdir(self.scr_qed_dir)

    def move_dir(self, base_dir):
        """ Move to the base directory

            :param string base_dir: Base directory
        """
        os.chdir(base_dir)


