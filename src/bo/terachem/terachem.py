from __future__ import division
from bo.bo_calculator import BO_calculator

class TeraChem(BO_calculator):
    """ Class for TeraChem program
    """
    def __init__(self, functional, basis_set, qm_path, ngpus, \
        gpu_id, precision, version):
        # Initialize TeraChem common variables
        self.functional = functional
        self.basis_set = basis_set

        self.qm_path = qm_path
        self.ngpus = ngpus
        self.gpu_id = gpu_id

        self.precision = precision
        self.version = version

        if (not (self.version == 1.92 or self.version == 1.93)):
            raise ValueError ("Other Versions Not Implemented and Tested")



