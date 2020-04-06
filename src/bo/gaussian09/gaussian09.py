from __future__ import division
from bo.bo_calculator import BO_calculator

class Gaussian09(BO_calculator):
    """ Class for common parts of Gaussian09 program

        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param integer nthreads: number of threads in the calculations
        :param string g09_root_path: path for Gaussian09 root
        :param string version: version of Gaussian09 program
    """
    def __init__(self, basis_set, memory, nthreads, g09_root_path, version):
        # Save name of BO calculator and its method
        super().__init__()

        # Initialize Gaussian09 common variables
        self.basis_set = basis_set

        self.memory = memory
        self.nthreads = nthreads
        self.g09_root_path = g09_root_path
        self.version = version

        # Print Gaussian09 Revision warning
        print("\n\n WARNING: The Gaussian09 implementation is based on the Revision A.02 version, not the latest one! \n\n")


