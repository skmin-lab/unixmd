from __future__ import division
from bo.bo_calculator import BO_calculator

class Turbomole(BO_calculator):
    """ Class for TURBOMOLE program
    """
    def __init__(self, functional, basis_set, memory, qm_path, nthreads, version):
      #
        self.functional = functional
        self.basis_set = basis_set
        
        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version
        
        if (self.version != 6.4):
            raise ValueError ("TURBOMOMLE6.4 is ONLY available!")
      
