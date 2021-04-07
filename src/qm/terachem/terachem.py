from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name

class TeraChem(QM_calculator):
    """ Class for common parts of TeraChem

        :param string basis_set: Basis set information
        :param string functional: Exchange-correlation functional information
        :param string precision: Precision in the calculations
        :param string qm_path: Path for QM binary
        :param integer ngpus: Number of GPUs
        :param integer,list gpu_id: ID of used GPUs
        :param string version: Version of TeraChem
    """
    def __init__(self, functional, basis_set, qm_path, ngpus, \
        gpu_id, precision, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize TeraChem common variables
        self.functional = functional
        self.basis_set = basis_set

        self.qm_path = qm_path
        self.ngpus = ngpus
        self.gpu_id = gpu_id
        
        if (self.gpu_id == None):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) The data type of gpu_id should be given as list! {self.gpu_id}")

        if (len(self.gpu_id) != self.ngpus):
            raise ValueError (f"( {self.qm_method}.{call_name()} ) The length of gpu_id should be same to ngpus! ({self.gpu_id} != {self.ngpus})")

        self.precision = precision
        self.version = version

        if (not (self.version == "1.93")):
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")



