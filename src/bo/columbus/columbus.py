from __future__ import division
from bo.bo_calculator import BO_calculator
import os

class Columbus(BO_calculator):
    """ Class for common parts of Columbus program

        :param object molecule: molecule object
        :param string basis_set: basis set information
        :param string memory: allocatable memory in the calculations
        :param string qm_path: path for QM binary
        :param integer nthreads: number of threads in the calculations
        :param double version: version of Columbus program
    """
    def __init__(self, molecule, basis_set, memory, qm_path, nthreads, version):
        # Save name of BO calculator and its method
        super().__init__()

        # Initialize Columbus common variables
        self.basis_set = basis_set

        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version

        # qm_path should be saved in the environmental variable "COLUMBUS"
        os.environ["COLUMBUS"] = self.qm_path

        # Check the atomic species with sorted command
        self.atom_type = sorted(set(molecule.symbols))

        # Check basis sets
        for itype in self.atom_type:
            try:
                pass
            except:
                raise KeyError (f"Data not found in Columbus : Atom {itype} with basis {self.basis_set}")

        if (self.version != 7.0):
            raise ValueError ("Other Versions Not Implemented")


