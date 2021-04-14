from __future__ import division
from qm.qm_calculator import QM_calculator
from qm.columbus.colbasis import *
from misc import call_name
import os

class Columbus(QM_calculator):
    """ Class for common parts of Columbus

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param integer memory: Allocatable memory in the calculations
        :param string qm_path: Path for QM binary
        :param string version: Version of Columbus
    """
    def __init__(self, molecule, basis_set, memory, qm_path, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize Columbus common variables
        self.basis_set = basis_set.lower()

        self.memory = memory
        self.qm_path = qm_path
        if (not os.path.isdir(self.qm_path)):
            error_message = "Directory for Columbus binary not found!"
            error_vars = f"qm_path = {self.qm_path}"
            raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.version = version

        # qm_path should be saved in the environmental variable "COLUMBUS"
        os.environ["COLUMBUS"] = self.qm_path

        # Check the atomic species with sorted command
        self.atom_type = sorted(set(molecule.symbols))

        # Set basis sets information
        try:
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
                error_message = "Current basis set not implemented, check '$PYUNIXMDHOME/src/qm/columbus/colbasis.py'!"
                error_vars = f"basis_set = {self.basis_set}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        except KeyError:
            error_message = "Current atom type in the basis set not implemented, check '$PYUNIXMDHOME/src/qm/columbus/colbasis.py'!"
            error_vars = f"basis_set = {self.basis_set}"
            raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (isinstance(self.version, str)):
            if (self.version != "7.0"):
                error_message = "Other versions not implemented!"
                error_vars = f"version = {self.version}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of version must be string!"
            error_vars = f"version = {self.version}"
            raise TypeError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

