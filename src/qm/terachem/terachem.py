from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name

class TeraChem(QM_calculator):
    """ Class for common parts of TeraChem program

        :param string basis_set: Basis set information
        :param string functional: Functional in the calculations
        :param string precision: Precision in the calculations
        :param string qm_path: Path for QM binary
        :param integer ngpus: Number of GPUs
        :param string gpu_id: ID of used GPUs
        :param string version: Version of TeraChem program
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

        self.precision = precision
        self.version = version

        if (not (self.version == "1.99" or self.version == "1.93")):
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")



