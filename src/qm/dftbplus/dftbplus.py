from __future__ import division
from qm.qm_calculator import QM_calculator
from misc import call_name
import os

class DFTBplus(QM_calculator):
    """ Class for common parts of DFTB+ program

        :param object molecule: molecule object
        :param string sk_path: path for slater-koster files
        :param string qm_path: path for QM binary
        :param string script_path: path for DFTB+ python script (dptools)
        :param integer nthreads: number of threads in the calculations
        :param double version: version of DFTB+ program
    """
    def __init__(self, molecule, sk_path, qm_path, script_path, nthreads, version):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize DFTB+ common variables
        self.sk_path = sk_path
        self.qm_path = qm_path
        # Note that self.script_path should include path upto $DFTBPLUS/tools/dptools
        self.script_path = script_path

        self.nthreads = nthreads
        self.version = version

        # Environmental variable setting for python scripts such as xyz2gen used in DFTB+
        # Note that self.script_path should include upto install directory of DFTB+
        if (self.version == 19.1 or self.version == 20.1):
            script_dir = os.path.join(self.script_path, "bin")
            # Note that the python version can be changed according to the users setting
            lib_dir = os.path.join(self.script_path, "lib/python3.6/site-packages")
            if (not os.path.exists(lib_dir)):
                raise ValueError (f"( {self.qm_prog}.{call_name()} ) Please modify python version in interfacing script! {lib_dir}")
        else:
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")

        # Append following paths to PATH and PYTHONPATH variables
        os.environ["PATH"] += os.pathsep + os.path.join(script_dir)
        os.environ["PYTHONPATH"] += os.pathsep + os.path.join(lib_dir)

        # Check the atomic species
        self.atom_type = set(molecule.symbols[0:molecule.nat_qm])


