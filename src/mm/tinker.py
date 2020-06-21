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

        if (not (self.scheme == "additive" or self.scheme == "subtractive")):
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) Wrong QM/MM scheme given! {self.scheme}")

    def get_mm(self, molecule, base_dir):
        """ Extract energy and gradient from Tinker

            :param object molecule: molecule object
            :param string base_dir: base directory
        """
        super().get_mm(base_dir)
        self.write_xyz(molecule)
#        self.get_input(molecule, bo_list)
#        self.run_QM(base_dir, istep, bo_list)
#        self.extract_BO(molecule, bo_list)
        self.move_dir(base_dir)

    def write_xyz(self, molecule):
        """ Make current geometry file with Tinker format

            :param object molecule: molecule object
        """
        if (self.scheme == "additive"):
            # Atoms in MM region
            file_name = "tinker.xyz.2"
            with open(file_name, "w") as ftj:
                ftj.write(f"{molecule.nat_mm}\n\n")
                for iat in range(molecule.nat_qm, molecule.nat):
                    ftj.write(f"{molecule.symbols[iat]:4}")
                    ftj.write("".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A]) + "\n")

        elif (self.scheme == "subtractive"):
            # Atoms in QM + MM region
            file_name = "tinker.xyz.12"
            with open(file_name, "w") as ftj:
                ftj.write(f"{molecule.nat}\n\n")
                for iat in range(molecule.nat):
                    ftj.write(f"{molecule.symbols[iat]:4}")
                    ftj.write("".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A]) + "\n")
            # Atoms in QM region
            file_name = "tinker.xyz.1"
            with open(file_name, "w") as ftj:
                ftj.write(f"{molecule.nat_qm}\n\n")
                for iat in range(molecule.nat_qm):
                    ftj.write(f"{molecule.symbols[iat]:4}")
                    ftj.write("".join([f"{i:15.8f}" for i in molecule.pos[iat] * au_to_A]) + "\n")



