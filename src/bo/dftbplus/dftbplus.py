from __future__ import division
from bo.bo_calculator import BO_calculator

class DFTBplus(BO_calculator):
    """ Class for common parts of DFTB+ program

        :param object molecule: molecule object
        :param string sk_path: path for slater-koster files
        :param string qm_path: path for QM binary
        :param integer nthreads: number of threads in the calculations
        :param double version: version of DFTB+ program
    """
    def __init__(self, molecule, sk_path, qm_path, nthreads, version):
        # Initialize DFTB+ common variables
        self.sk_path = sk_path
        self.qm_path = qm_path

        self.nthreads = nthreads
        self.version = version

        # Check the atomic species
        self.atom_type = set(molecule.symbols)


