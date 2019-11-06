from __future__ import division
from bo.bo_calculator import BO_calculator

class DFTBplus(BO_calculator):
    """ Class for DFTBplus program
    """
    def __init__(self, molecule, sk_path, qm_path, nthreads, version):
        # Initialize DFTBplus common variables
        self.sk_path = sk_path
        self.qm_path = qm_path

        self.nthreads = nthreads
        self.version = version

        # check the atomic species
        self.atom_type = set(molecule.symbols)




