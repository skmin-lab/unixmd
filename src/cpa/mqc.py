from __future__ import division
from misc import fs_to_au, au_to_A, call_name, typewriter
import textwrap, datetime
import numpy as np
import os, shutil


class MQC(object):
    """ Class for electronic propagator used in MQC dynamics with classical path approximation
    """
    def __init__(self):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

    def run_init(self):
        """ Initialize MQC dynamics
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self):
        """ Routine to print the initial information of dynamics
        """
        pass

    def touch_file(self):
        """ Routine to write PyUNIxMD output files
        """
        pass

    def write_md_output(self, unixmd_dir, istep):
        """ Write output files
        """
        pass
