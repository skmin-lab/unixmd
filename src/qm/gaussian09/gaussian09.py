from __future__ import division
from qm.qm_calculator import QM_calculator

class Gaussian09(QM_calculator):
    """ Class for common parts of Gaussian09

        :param string basis_set: Basis set information
        :param string memory: Allocatable memory
        :param integer nthreads: Number of threads in the calculations
        :param string g09_root_path: Path for Gaussian 09 root
        :param string version: Version of Gaussian09
    """
    def __init__(self, basis_set, memory, nthreads, g09_root_path, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize Gaussian09 common variables
        self.basis_set = basis_set

        self.memory = memory
        self.nthreads = nthreads
        self.g09_root_path = g09_root_path
        self.version = version

        # Print Gaussian09 Revision warning
        print("\n\n WARNING: The Gaussian09 implementation is based on the Revision A.02 version, not the latest one! \n\n")


