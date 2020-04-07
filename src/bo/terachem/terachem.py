from __future__ import division
from bo.bo_calculator import BO_calculator

class TeraChem(BO_calculator):
    """ Class for common parts of TeraChem program

        :param string basis_set: basis set information
        :param string functional: functional in the calculations
        :param string precision: precision in the calculations
        :param string qm_path: path for QM binary
        :param integer ngpus: number of GPUs
        :param string gpu_id: ID of used GPUs
        :param double version: version of TeraChem program
    """
    def __init__(self, functional, basis_set, qm_path, ngpus, \
        gpu_id, precision, version):
        # Save name of BO calculator and its method
        super().__init__()

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



