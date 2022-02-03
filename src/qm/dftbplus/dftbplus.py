from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name
import os

class DFTBplus(QM_calculator):
    """ Class for common parts of DFTB+

        :param object molecule: Molecule object
        :param string sk_path: Path for Slater-Koster files
        :param string install_path: Path for DFTB+ install directory
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of DFTB+
    """
    def __init__(self, molecule, sk_path, install_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize DFTB+ common variables
        self.sk_path = sk_path
        self.install_path = install_path
        if (not os.path.isdir(self.install_path)):
            error_message = "Install directory for DFTB+ not found!"
            error_vars = f"install_path = {self.install_path}"
            raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        self.nthreads = nthreads
        self.version = version

        # Environmental variable setting for Python scripts such as xyz2gen used in DFTB+
        if (isinstance(self.version, str)):
            if (self.version in ["19.1", "20.1", "21.1"]):
                self.qm_path = os.path.join(self.install_path, "bin")
                # Note that the Python version can be changed according to the users setting
                lib_dir = os.path.join(self.install_path, "lib/python3.6/site-packages")
                if (not os.path.exists(lib_dir)):
                    error_message = "Please use proper Python version number in '$PYUNIXMDHOME/src/qm/dftbplus/dftbplus.py'!"
                    error_vars = f"library directory = {lib_dir}"
                    raise FileNotFoundError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                error_message = "Other versions not implemented!"
                error_vars = f"version = {self.version}"
                raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of version must be string!"
            error_vars = f"version = {self.version}"
            raise TypeError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Append following paths to PATH and PYTHONPATH variables
        os.environ["PATH"] += os.pathsep + os.path.join(self.qm_path)
        os.environ["PYTHONPATH"] += os.pathsep + os.path.join(lib_dir)

        # Check the atomic species
        self.atom_type = set(molecule.symbols[0:molecule.nat_qm])


