from __future__ import division
from bo.bo_calculator import BO_calculator
import os

class Columbus(BO_calculator):
    """ Class for Columbus program
    """
    def __init__(self, molecule, basis_set, memory, qm_path, nthreads, version):
        # Initialize Columbus common variables
        self.basis_set = basis_set

        self.memory = memory
        self.qm_path = qm_path
        self.nthreads = nthreads
        self.version = version

        # qm_path should be saved in environmental variable "COLUMBUS"
        os.environ["COLUMBUS"] = self.qm_path

        # check the atomic species
        self.atom_type = sorted(set(molecule.symbols))

        # check basis sets
        for itype in self.atom_type:
            try:
                pass
            except:
                raise KeyError (f"Data not found in Columbus : Atom {itype} with basis {self.basis_set}")

        if (self.version != 7.0):
            raise ValueError ("Other Versions Not Implemented")




