from __future__ import division
from qm.gamess.gamess import GAMESS
from misc import call_name

class SSR(GAMESS):
    """ Class for SSR method of GAMESS

        :param object molecule: Molecule object
        :param string basis_set: Basis set information
        :param string memory: Allocatable memory in the calculations
        :param string functional: Exchange-correlation functional information
        :param integer active_space: Active space for SSR calculation
        :param string qm_path: Path for QM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of GAMESS, check $VERNO
    """
    def __init__(self, molecule, basis_set="6-31g*", memory="50", functional="bhhlyp", \
        active_space=2, qm_path="./", nthreads=1, version="00"):
        # Initialize GAMESS common variables
        super(SSR, self).__init__(basis_set, memory, qm_path, nthreads, version)

        # Initialize GAMESS SSR variables
        self.functional = functional

        self.active_space = active_space
        if not (self.active_space in [2]):
            error_message = "Invalid active space for SSR!"
            error_vars = f"active_space = {self.active_space}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        if (self.active_space == 2):
            if (molecule.nst > 2):
                error_message = "SSR(3,2) calculation is not implemented in current version of pyUNIxMD!"
                error_vars = f"Molecule.nstates = {molecule.nst}"
                raise NotImplementedError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Check the closed shell for systems
        if (not int(molecule.nelec) % 2 == 0):
            error_message = "Only closed shell configuration supported, check charge!"
            error_vars = f"Molecule.nelec = {int(molecule.nelec)}"
            raise ValueError (f"( {self.qm_method}.{call_name()} ) {error_message} ( {error_vars} )")

        # Set 'l_nacme' with respect to the computational method
        # SSR can produce NACs, so we do not need to get NACME from CIoverlap
        # SSR can compute the gradient of several states simultaneously.
        molecule.l_nacme = False
        self.re_calc = False


