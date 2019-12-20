from __future__ import division
from bo.bo_calculator import BO_calculator


class columbus(BO_calculator):
    """ Class for Columbus program
    """
    def __init__(self, basis_set, memory, qm_path, nthreads, version):
        # Initialize Columbus common variables

        self.basis_set = basis_set
        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version

        if (self.version != 7.0):
            raise ValueError ("Other Versions Not Implemented")




