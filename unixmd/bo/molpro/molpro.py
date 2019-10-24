from __future__ import division
from bo.bo_calculator import BO_calculator

class Molpro(BO_calculator):
    """ Class for Molpro program
    """
    def __init__(self, basis_set, memory, qm_path, nthreads, version):
        # Initialize Molpro common variables
        self.basis_set = basis_set

        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version

        if (self.version != 2015.1):
            raise ValueError ("Other Versions Not Implemented")




