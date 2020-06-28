from __future__ import division
from mm.mm_calculator import MM_calculator
from misc import au_to_A, call_name
import os, shutil

class Tinker(MM_calculator):
    """ Class for Tinker program

        :param object molecule: molecule object
    """
    def __init__(self, molecule, scheme=None, do_charge=False, do_vdw=False, periodic=False, \
        cell_par=[0., 0., 0., 0., 0., 0.], xyz_file="./tinker.xyz", key_file="./tinker.key"):
        # Save name of MM calculator
        super().__init__()

        self.scheme = scheme

        self.do_charge = do_charge
        self.do_vdw = do_vdw

        self.periodic = periodic
        self.cell_par = cell_par

        self.xyz_file = xyz_file
        self.key_file = key_file

        if (not (self.scheme == "additive" or self.scheme == "subtractive")):
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) Wrong QM/MM scheme given! {self.scheme}")

    def get_mm(self, molecule, base_dir):
        """ Extract energy and gradient from Tinker

            :param object molecule: molecule object
            :param string base_dir: base directory
        """
        super().get_mm(base_dir)
        self.get_input(molecule)
#        self.run_QM(base_dir, istep, bo_list)
#        self.extract_BO(molecule, bo_list)
        self.move_dir(base_dir)

    def get_input(self, molecule):
        """ Generate Tinker input files: tinker.xyz, tinker.key

            :param object molecule: molecule object
        """
        # Copy xyz file to currect directory
        if (os.path.isfile(self.xyz_file)):
            shutil.copy(self.xyz_file, os.path.join(self.scr_mm_dir, "tinker.xyz"))
        else:
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) Initial tinker.xyz file not given! {self.xyz_file}")

        # Copy key file to currect directory
        if (os.path.isfile(self.key_file)):
            shutil.copy(self.key_file, os.path.join(self.scr_mm_dir, "tinker.key"))
        else:
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) Initial tinker.key file not given! {self.key_file}")

        # Write xyz file using current positions
        if (self.scheme == "additive"):

            # Atoms in MM region
            input_xyz2 = ""

            input_xyz2 += f" {molecule.nat_mm}\n"
            # Information about periodicity, line_period is number of lines to be skipped in 'tinker.xyz' file
            line_period = 1
            if (self.periodic):
                line_period = 2
                input_xyz2 += " ".join([f"{ i:12.6f}" for i in self.cell_par]) + "\n"

            # Read 'tinker.xyz' file to obtain atom type and topology
            file_name = "tinker.xyz"
            with open(file_name, "r") as f_xyz:
                lines = f_xyz.readlines()
                iline = 1
                for line in lines:
                    # Skip first or second lines
                    if (iline > line_period + molecule.nat_qm):
                        ind = iline - line_period
                        input_geom = f" {ind - molecule.nat_qm:6d}"
                        input_geom += f" {molecule.symbols[ind - 1]:4}"
                        input_geom += "".join([f"{i:15.8f}" for i in molecule.pos[ind - 1] * au_to_A])
                        # Count index of column in coordinate lines
                        col = 0
                        field = line.split()
                        for element in field:
                            if (col > 4):
                                # Read right block; atom type and topology information
                                input_geom += f" {element:6s}"
                            col += 1
                        input_xyz2 += input_geom + "\n"
                    iline += 1

            # Make 'tinker.xyz.2' file
            file_name = "tinker.xyz.2"
            with open(file_name, "w") as f_xyz:
                f_xyz.write(input_xyz2)

        elif (self.scheme == "subtractive"):

            # Atoms in QM + MM region
            input_xyz12 = ""

            input_xyz12 += f" {molecule.nat}\n"
            # Information about periodicity, line_period is number of lines to be skipped in 'tinker.xyz' file
            line_period = 1
            if (self.periodic):
                line_period = 2
                input_xyz12 += " ".join([f"{ i:12.6f}" for i in self.cell_par]) + "\n"

            # Read 'tinker.xyz' file to obtain atom type and topology
            file_name = "tinker.xyz"
            with open(file_name, "r") as f_xyz:
                lines = f_xyz.readlines()
                iline = 1
                for line in lines:
                    # Skip first or second lines
                    if (iline > line_period):
                        ind = iline - line_period
                        input_geom = f" {ind:6d}"
                        input_geom += f" {molecule.symbols[ind - 1]:4}"
                        input_geom += "".join([f"{i:15.8f}" for i in molecule.pos[ind - 1] * au_to_A])
                        # Count index of column in coordinate lines
                        col = 0
                        field = line.split()
                        for element in field:
                            if (col > 4):
                                # Read right block; atom type and topology information
                                input_geom += f" {element:6s}"
                            col += 1
                        input_xyz12 += input_geom + "\n"
                    iline += 1

            # Make 'tinker.xyz.12' file
            file_name = "tinker.xyz.12"
            with open(file_name, "w") as f_xyz:
                f_xyz.write(input_xyz12)

            # Atoms in QM region
            input_xyz1 = ""

            input_xyz1 += f" {molecule.nat_qm}\n"
            # Information about periodicity, line_period is number of lines to be skipped in 'tinker.xyz' file
            line_period = 1
            if (self.periodic):
                line_period = 2
                input_xyz1 += " ".join([f"{ i:12.6f}" for i in self.cell_par]) + "\n"

            # Read 'tinker.xyz' file to obtain atom type and topology
            file_name = "tinker.xyz"
            with open(file_name, "r") as f_xyz:
                lines = f_xyz.readlines()
                iline = 1
                for line in lines:
                    # Skip first or second lines
                    if (iline in range(line_period + 1, molecule.nat_qm + line_period + 1)):
                        ind = iline - line_period
                        input_geom = f" {ind:6d}"
                        input_geom += f" {molecule.symbols[ind - 1]:4}"
                        input_geom += "".join([f"{i:15.8f}" for i in molecule.pos[ind - 1] * au_to_A])
                        # Count index of column in coordinate lines
                        col = 0
                        field = line.split()
                        for element in field:
                            if (col > 4):
                                # Read right block; atom type and topology information
                                input_geom += f" {element:6s}"
                            col += 1
                        input_xyz1 += input_geom + "\n"
                    iline += 1

            # Make 'tinker.xyz.1' file
            file_name = "tinker.xyz.1"
            with open(file_name, "w") as f_xyz:
                f_xyz.write(input_xyz1)

        # Set non-bonded interaction for the systems; charge term from 'tinker.key' file
        # This is mechanical embedding for charge-charge interaction
        # To deal with charge-charge interaction with electrostatic interaction,
        # set do_vdw=True in the initialization of theory object, not field object
        file_be = open('tinker.key', 'r')
        file_af = open('tmp.key', 'w')
        is_charge = False
        for line in file_be:
            if ("chargeterm" in line):
                is_charge = True
                line = ""
                if (not self.do_charge):
                    line = "chargeterm none\n"
            file_af.write(line)
        # If chargeterm keyword does not exist, add chargeterm keyword to last line
        if (not is_charge and not self.do_charge):
            line = "chargeterm none\n"
            file_af.write(line)
        file_be.close()
        file_af.close()
        os.rename('tmp.key', 'tinker.key')

        # Set non-bonded interaction for the systems; vdw term from 'tinker.key' file
        file_be = open('tinker.key', 'r')
        file_af = open('tmp.key', 'w')
        is_vdw = False
        for line in file_be:
            if ("vdwterm" in line):
                is_vdw = True
                line = ""
                if (not self.do_vdw):
                    line = "vdwterm none\n"
            file_af.write(line)
        # If vdwterm keyword does not exist, add vdwterm keyword to last line
        if (not is_vdw and not self.do_vdw):
            line = "vdwterm none\n"
            file_af.write(line)
        file_be.close()
        file_af.close()
        os.rename('tmp.key', 'tinker.key')


