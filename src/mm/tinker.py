from __future__ import division
from mm.mm_calculator import MM_calculator
from misc import au_to_A, call_name

class Tinker(MM_calculator):
    """ Class for Tinker program

        :param object molecule: molecule object
    """
    def __init__(self, molecule, scheme=None):
        # Save name of MM calculator
        super().__init__()

        self.scheme = scheme

        if (not (self.scheme == "additivie" or self.scheme == "subtractive")):
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) Wrong QM/MM scheme given! {self.scheme}")

#    def get_mm(self, molecule, base_dir, istep, bo_list, dt, calc_force_only):
    def get_mm(self, molecule, base_dir):
        """ Extract energy and gradient from Tinker

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer istep: current MD step
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        super().get_mm(base_dir)
#        self.write_xyz(molecule)
#        self.get_input(molecule, bo_list)
#        self.run_QM(base_dir, istep, bo_list)
#        self.extract_BO(molecule, bo_list)
        self.move_dir(base_dir)

    def write_xyz(self, molecule):
        """ Make current geometry file with Tinker format

            :param object molecule: molecule object
        """
        # TODO : pass a choice for QMMM scheme as argument
        if (self.scheme == "additive"):
        elif (self.scheme == "subtractive"):

#        file_name = "tinker.xyz"
#        with open(file_name, "w") as ftj:
#            ftj.write(f"{molecule.nat}\n\n")
#            for iat in range(molecule.nat):
#                ftj.write(f"{molecule.symbols[iat]:4}")
#                ftj.write("".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A]) + "\n")


