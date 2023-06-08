from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name
import os

class GAMESS(QM_calculator):
    """ Class for common parts of GAMESS

        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string qm_path: Path for QM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of GAMESS, check $VERNO
    """
    def __init__(self, basis_set, memory, qm_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize GAMESS common variables
        self.basis_set = basis_set

        self.memory = memory
        self.qm_path = qm_path
        if (not os.path.isdir(self.qm_path)):
            error_message = "Directory for GAMESS binary not found!"
            error_vars = f"qm_path = {self.qm_path}"
            raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.nthreads = nthreads
        self.version = version

        if (isinstance(self.version, str)):
            if not (self.version in ["00"]):
                error_message = "Other versions not implemented!"
                error_vars = f"version = {self.version}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of version must be string!"
            error_vars = f"version = {self.version}"
            raise TypeError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")


