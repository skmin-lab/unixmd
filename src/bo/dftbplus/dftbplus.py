from __future__ import division
from bo.bo_calculator import BO_calculator
from misc import call_name

class DFTBplus(BO_calculator):
    """ Class for common parts of DFTB+ program

        :param object molecule: molecule object
        :param string sk_path: path for slater-koster files
        :param string qm_path: path for QM binary
        :param integer nthreads: number of threads in the calculations
        :param double version: version of DFTB+ program
    """
    def __init__(self, molecule, sk_path, qm_path, script_path, nthreads, version):
        # Save name of BO calculator and its method
        super().__init__()

        # Initialize DFTB+ common variables
        self.sk_path = sk_path
        self.qm_path = qm_path
        # Note that self.script_path should include path upto $DFTBPLUS/tools/dptools
        self.script_path = script_path

        self.nthreads = nthreads
        self.version = version

        # Environmental variable setting for python scripts such as xyz2gen used in DFTB+
        # Note that self.script_path should include upto $DFTBPLUS/tools/dptools
        if (self.version == 19.1):
            # version 19.1
            lib_dir = os.path.join(self.script_path, "build/lib")
            # The version number can be changed according to the users setting
            script_dir = os.path.join(self.script_path, "build/scripts-3.6")
            if (not os.path.exists(script_dir)):
                raise ValueError (f"( {self.qm_method}.{call_name()} ) Please modify version number to correct number that you used! {script_dir}")
        elif (self.version == 20.1):
            # version 20.1
            lib_dir = os.path.join(self.script_path, "src")
            script_dir = os.path.join(self.script_path, "bin")
        else:
            raise ValueError (f"( {self.qm_prog}.{call_name()} ) Other version not implemented! {self.version}")

        # Append following paths to PATH and PYTHONPATH variables
        os.environ["PATH"] += os.pathsep + os.path.join(script_dir)
        os.environ["PYTHONPATH"] += os.pathsep + os.path.join(lib_dir)

        # Check the atomic species
        self.atom_type = set(molecule.symbols)


