from __future__ import division
from bo.bo_calculator import BO_calculator
from misc import call_name

class Turbomole(BO_calculator):
    """ Class for common parts of Turbomole program

        :param string functional: xc functional information
        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param string qm_path: path for QM binary
        :param integer nthreads: number of threads in the calculations
        :param double version: version of Turbomole program
    """
    def __init__(self, functional, basis_set, memory, qm_path, nthreads, version):
        # Save name of BO calculator and its method
        super().__init__()

        # Initialize Turbomole common variables
        self.functional = functional
        self.basis_set = basis_set

        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version

        if (self.version != 6.4):
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")


