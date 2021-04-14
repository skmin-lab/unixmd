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
        if (not os.path.isdir(self.qm_path)):
            error_message = "Directory for TeraChem binary not found!"
            error_vars = f"qm_path = {self.qm_path}"
            raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.ngpus = ngpus
        self.gpu_id = gpu_id
        
        if (self.gpu_id == None):
            error_message = "GPU ID must be set in running script!"
            error_vars = f"gpu_id = {self.gpu_id}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (isinstance(self.gpu_id, list)):
            if (len(self.gpu_id) != self.ngpus):
                error_message = "Number of elements for GPU ID must be equal to number of GPUs!"
                error_vars = f"len(gpu_id) = {len(self.gpu_id)}, ngpus = {self.ngpus}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of GPU ID must be list consisting of integer!"
            error_vars = f"gpu_id = {self.gpu_id}"
            raise TypeError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.precision = precision
        self.version = version

        if (isinstance(self.version, str)):
            if (self.version != "1.93"):
                error_message = "Other versions not implemented!"
                error_vars = f"version = {self.version}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of version must be string!"
            error_vars = f"version = {self.version}"
            raise TypeError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")


