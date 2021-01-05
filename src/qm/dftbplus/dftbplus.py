from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name
import os

class DFTBplus(QM_calculator):
    """ Class for common parts of DFTB+ program

        :param object molecule: molecule object
        :param string sk_path: path for slater-koster files
        :param string install_path: path for DFTB+ install directory
        :param integer nthreads: number of threads in the calculations
        :param string version: version of DFTB+ program
    """
    def __init__(self, molecule, sk_path, install_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize DFTB+ common variables
        self.sk_path = sk_path
        self.install_path = install_path

        self.nthreads = nthreads
        self.version = version

        # Environmental variable setting for python scripts such as xyz2gen used in DFTB+
        if (self.version == "19.1" or self.version == "20.1"):
            self.qm_path = os.path.join(self.install_path, "bin")
            # Note that the python version can be changed according to the users setting
            lib_dir = os.path.join(self.install_path, "lib/python3.6/site-packages")
            if (not os.path.exists(lib_dir)):
                raise ValueError (f"( {self.qm_prog}.{call_name()} ) Please modify python version in interfacing script! {lib_dir}")
        else:
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")

        # Append following paths to PATH and PYTHONPATH variables
        os.environ["PATH"] += os.pathsep + os.path.join(self.qm_path)
        os.environ["PYTHONPATH"] += os.pathsep + os.path.join(lib_dir)

        # Check the atomic species
        self.atom_type = set(molecule.symbols[0:molecule.nat_qm])


