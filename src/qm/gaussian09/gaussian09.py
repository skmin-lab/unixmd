from __future__ import division
from qm.qm_calculator import QM_calculator
import os

class Gaussian09(QM_calculator):
    """ Class for common parts of Gaussian 09

        :param string basis_set: Basis set information
        :param string memory: Allocatable memory
        :param integer nthreads: Number of threads in the calculations
        :param string root_path: Path for Gaussian 09 root directory
        :param string version: Version of Gaussian 09
    """
    def __init__(self, basis_set, memory, nthreads, root_path, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize Gaussian09 common variables
        self.basis_set = basis_set

        self.memory = memory
        self.nthreads = nthreads

        self.root_path = root_path
        if (not os.path.isdir(self.root_path)):
            error_message = "Root directory for Gaussian 09 not found!"
            error_vars = f"root_path = {self.root_path}"
            raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.version = version

        # Print Gaussian09 Revision warning
        print ("\n\n WARNING: The Gaussian09 implementation is based on the Revision A.02 version, not the latest one! \n\n", flush=True)


