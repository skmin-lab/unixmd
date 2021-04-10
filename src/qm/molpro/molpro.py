from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name

class Molpro(QM_calculator):
    """ Class for common parts of Molpro

        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string qm_path: Path for QM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of Molpro
    """
    def __init__(self, basis_set, memory, qm_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize Molpro common variables
        self.basis_set = basis_set

        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version

        if (self.version != "2015.1"):
            error_message = "Other versions not implemented!"
            error_vars = f"version = {self.version}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")


