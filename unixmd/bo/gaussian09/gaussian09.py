from __future__ import division
from bo.bo_calculator import BO_calculator

class Gaussian09(BO_calculator):
    """ Class for Gaussian09 program
    """
    def __init__(self, basis_set, memory, nthreads, g09_root_path, version):
        # Initialize Gaussian09 common variables
        self.basis_set = basis_set

        self.memory = memory
        self.nthreads = nthreads
        self.g09_root_path = g09_root_path
        self.version = version

        # print Gaussian09 Revision warning
        print("\n\n WARNING: The Gaussian09 implementation is based on the Revision A.02 version, not the latest one! \n\n")


