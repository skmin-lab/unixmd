from __future__ import division
from bo.bo_calculator import BO_calculator

class Molpro(BO_calculator):
    """ Class for common parts of Molpro program

        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param string qm_path: path for QM binary
        :param integer nthreads: number of threads in the calculations
        :param double version: version of Molpro program
    """
    def __init__(self, basis_set, memory, qm_path, nthreads, version):
        # Save name of BO calculator and its method
        super().__init__()

        # Initialize Molpro common variables
        self.basis_set = basis_set

        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version

        if (self.version != 2015.1):
            raise ValueError ("Other Versions Not Implemented")




