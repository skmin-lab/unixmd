from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name
import os

class QChem(QM_calculator):
    """ Class for common parts of Q-Chem

        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculation
        :param string root_path: Path for Q-Chem root directory
        :param integer nthreads: Number of threads in the calculation
        :param string version: Version of Q-Chem
    """
    def __init__(self, basis_set, memory, root_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize Q-Chem common variables
        self.basis_set = basis_set

        self.memory = memory
        self.root_path = root_path
        if (not os.path.isdir(self.root_path)):
            error_message = "Root directory for Q-Chem not found!"
            error_vars = f"root_path = {self.root_path}"
            raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.nthreads = nthreads
        self.version = version

        if (isinstance(self.version, str)):
            if (self.version not in ["5.2", "5.4"]):
                error_message = "Other versions not implemented!"
                error_vars = f"version = {self.version}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of version must be string!"
            error_vars = f"version = {self.version}"
            raise TypeError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

