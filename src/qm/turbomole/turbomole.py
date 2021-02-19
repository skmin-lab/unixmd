from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name
import os

class Turbomole(QM_calculator):
    """ Class for common parts of Turbomole program

        :param string functional: Exchange-correlation functional information
        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string qm_path: Path for QM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of Turbomole program
    """
    def __init__(self, functional, basis_set, memory, qm_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize Turbomole common variables
        self.functional = functional
        self.basis_set = basis_set

        self.memory = memory
        self.nthreads = nthreads
        self.version = version

        self.qm_path = qm_path
        os.environ["TURBODIR"] = qm_path
        self.qm_scripts_path = os.path.join(self.qm_path, "scripts/")
        if (self.nthreads == 1):
            self.qm_bin_path = os.path.join(self.qm_path, "bin/em64t-unknown-linux-gnu/")
        else:
            os.environ["PARA_ARCH"] = "SMP"
            os.environ["PARNODES"] = f"{self.nthreads}"
            self.qm_bin_path = os.path.join(self.qm_path, "bin/em64t-unknown-linux-gnu_smp/")

        if (self.version != "6.4"):
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")

