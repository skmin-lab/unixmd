from __future__ import division
from bo.bo_calculator import BO_calculator
from bo.columbus.colbasis import *
from misc import call_name
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
                raise KeyError (f"( {self.qm_method}.{call_name()} ) Data not found: Atom {itype} with basis! {self.basis_set}")

        # Set basis sets information
        if (self.basis_set == "cc-pvdz"):
            self.basis_nums = "\n".join([f"{cc_pvdz[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "cc-pvtz"):
            self.basis_nums = "\n".join([f"{cc_pvtz[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "cc-pvqz"):
            self.basis_nums = "\n".join([f"{cc_pvqz[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "3-21g*"):
            self.basis_nums = "\n".join([f"{t_21gs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "3-21+g*"):
            self.basis_nums = "\n".join([f"{t_21pgs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-31g"):
            self.basis_nums = "\n".join([f"{s_31g[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-31g*"):
            self.basis_nums = "\n".join([f"{s_31gs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-31+g*"):
            self.basis_nums = "\n".join([f"{s_31pgs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-311g*"):
            self.basis_nums = "\n".join([f"{s_311gs[f'{itype}']}" for itype in self.atom_type])
        elif (self.basis_set == "6-311g+*"):
            self.basis_nums = "\n".join([f"{s_311pgs[f'{itype}']}" for itype in self.atom_type])
        else:
            raise ValueError (f"( {self.qm_method}.{call_name()} ) No basis set in colbasis.py! {self.basis_set}")

        if (self.version != 7.0):
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")


