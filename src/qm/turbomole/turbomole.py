from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name
import os

class Turbomole(QM_calculator):
    """ Class for common parts of Turbomole

        :param string functional: Exchange-correlation functional information
        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string root_path: Path for Turbomole root directory
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of Turbomole
    """
    def __init__(self, functional, basis_set, memory, root_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize Turbomole common variables
        self.functional = functional
        self.basis_set = basis_set

        self.memory = memory
        self.nthreads = nthreads
        self.version = version

        self.root_path = root_path
        if (not os.path.isdir(self.root_path)):
            error_message = "Root directory for Turbomole not found!"
            error_vars = f"root_path = {self.root_path}"
            raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            os.environ["TURBODIR"] = root_path
            self.scripts_path = os.path.join(self.root_path, "scripts/")
            if (self.nthreads == 1):
                self.qm_path = os.path.join(self.root_path, "bin/em64t-unknown-linux-gnu/")
            else:
                os.environ["PARA_ARCH"] = "SMP"
                os.environ["PARNODES"] = f"{self.nthreads}"
                self.qm_path = os.path.join(self.root_path, "bin/em64t-unknown-linux-gnu_smp/")

        if (isinstance(self.version, str)):
            if (self.version != "6.4"):
                error_message = "Other versions not implemented!"
                error_vars = f"version = {self.version}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of version must be string!"
            error_vars = f"version = {self.version}"
            raise TypeError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")


