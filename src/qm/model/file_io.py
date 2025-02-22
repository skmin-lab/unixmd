from __future__ import division
from qm.model.model import Model

class File_IO(Model):
    """ Class for reading the pre-calculated data

        :param object molecule: molecule object
    """
    def __init__(self, molecule):
        # Initialize model common variables
        super(File_IO, self).__init__(None)

        # Set 'l_nacme' with respect to the computational method
        # File_IO can provide NACMEs, so we do not need to get NACME
        molecule.l_nacme = False

        # File_IO can provide the gradient of several states simultaneously
        self.re_calc = False


