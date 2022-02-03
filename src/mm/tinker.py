from __future__ import division
from mm.mm_calculator import MM_calculator
from misc import au_to_A, kcalmol_to_au, call_name
import os, shutil, re, textwrap
import numpy as np

class Tinker(MM_calculator):
    """ Class for Tinker

        :param object molecule: Molecule object
        :param string scheme: Type of QM/MM scheme
        :param string embedding: Charge embedding options
        :param string vdw: Van der Walls interactions
        :param boolean l_periodic: Use periodicity in the calculations
        :param double,list cell_par: Cell lattice parameters (lengths and angles)
        :param string xyz_file: Initial tinker.xyz file
        :param string key_file: Initial tinker.key file
        :param string mm_path: Path for MM binary
        :param integer nthreads: Number of threads in the calculations
        :param string version: Version of Tinker
    """
    def __init__(self, molecule, scheme=None, embedding=None, vdw=None, l_periodic=False, \
        cell_par=[0., 0., 0., 0., 0., 0.], xyz_file="./tinker.xyz", key_file="./tinker.key",
        mm_path="./", nthreads=1, version="8.7"):
        # Save name of MM calculator
        super().__init__()

        self.scheme = scheme
        if (self.scheme != None):
            self.scheme = self.scheme.lower()

        if not (self.scheme in ["additive", "subtractive"]):
            error_message = "Invalid scheme for QM/MM calculation!"
            error_vars = f"scheme = {self.scheme}"
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        self.embedding = embedding
        if (self.embedding != None):
            self.embedding = self.embedding.lower()

        if not (self.embedding in [None, "mechanical", "electrostatic"]):
            error_message = "Invalid charge embedding for QM/MM calculation!"
            error_vars = f"embedding = {self.embedding}"
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        self.vdw = vdw
        if (self.vdw != None):
            self.vdw = self.vdw.lower()

        if not (self.vdw in [None, "lennardjones"]):
            error_message = "Invalid van der Waals interaction for QM/MM calculation!"
            error_vars = f"vdw = {self.vdw}"
            raise ValueError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        self.l_periodic = l_periodic
        self.cell_par = cell_par

        self.xyz_file = xyz_file
        self.key_file = key_file

        self.mm_path = mm_path
        if (not os.path.isdir(self.mm_path)):
            error_message = "Directory for Tinker binary not found!"
            error_vars = f"mm_path = {self.mm_path}"
            raise FileNotFoundError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        self.nthreads = nthreads
        self.version = version

        if (isinstance(self.version, str)):
            if (self.version != "8.7"):
                error_message = "Other versions not implemented!"
                error_vars = f"version = {self.version}"
                raise ValueError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of version must be string!"
            error_vars = f"version = {self.version}"
            raise TypeError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        if (self.embedding == "electrostatic"):
            # Save current atom type for electrostatic embedding
            self.atom_type = np.zeros(molecule.nat, dtype=np.int32)

            # Read 'tinker.xyz' file to obtain atom type for QM part
            file_name = self.xyz_file
            with open(file_name, "r") as f_xyz:
                lines = f_xyz.readlines()
                # Check the number of lines; Read only non-periodic format for tinker.xyz file
                if (len(lines) != molecule.nat + 1):
                    error_message = "PyUNIxMD reads only non-periodic xyz file, the number of lines is not matched with non-periodic xyz file!"
                    error_vars = f"number of lines != {molecule.nat + 1}, xyz_file = {self.xyz_file}"
                    raise ValueError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

                iline = 1
                for line in lines:
                    # Skip first line
                    if (iline in range(2, molecule.nat + 2)):
                        ind = iline - 1
                        # Count index of column in coordinate lines
                        col = 0
                        field = line.split()
                        for element in field:
                            if (col == 5):
                                # Read atom type
                                self.atom_type[ind - 1] = int(element)
                            col += 1
                    iline += 1

            # TODO : add part to check necessary keyword such as parameters, etc

            # Read 'tinker.key' file to obtain charge for MM part
            file_name = self.key_file
            is_prm = False
            with open(file_name, "r") as f_key:
                lines = f_key.readlines()
                for line in lines:
                    if ("parameters" in line):
                        is_prm = True
                        field = line.split()
                        prm_file = field[1]
                if (not is_prm):
                    error_message = "To decide force field parameters, the keyword 'parameters' must exists in the key file!"
                    error_vars = f"key_file = {self.key_file}"
                    raise ValueError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

            # Save point charge values of current parameter file for electrostatic embedding
            point_charges = []

            file_name = prm_file
            with open(file_name, "r") as f_prm:
                lines = f_prm.readlines()
                for line in lines:
                    if ("charge" in line):
                        field = line.split()
                        if (len(field) == 3):
                            point_charges.append(field[2])

            # Assign point charge values for MM atoms
            for iat in range(molecule.nat_mm):
                ind = molecule.nat_qm + iat
                molecule.mm_charge[iat] = point_charges[self.atom_type[ind] - 1]

    def get_data(self, molecule, base_dir, bo_list, istep, calc_force_only):
        """ Extract energy and gradient from Tinker

            :param object molecule: Molecule object
            :param string base_dir: Base directory
            :param integer,list bo_list: List of BO states for BO calculation
            :param integer istep: Current MD step
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        super().get_data(base_dir, calc_force_only)
        if (not calc_force_only):
            self.get_input(molecule)
            self.run_MM(base_dir, istep)
        self.extract_MM(molecule, bo_list, calc_force_only)
        self.move_dir(base_dir)

    def get_input(self, molecule):
        """ Generate Tinker input files: tinker.xyz, tinker.key

            :param object molecule: Molecule object
        """
        # Copy xyz file to currect directory
        if (os.path.isfile(self.xyz_file)):
            shutil.copy(self.xyz_file, os.path.join(self.scr_mm_dir, "tinker.xyz"))
        else:
            error_message = "Initial tinker.xyz file not given, check file name!"
            error_vars = f"xyz_file = {self.xyz_file}"
            raise FileNotFoundError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        # Copy key file to currect directory
        if (os.path.isfile(self.key_file)):
            shutil.copy(self.key_file, os.path.join(self.scr_mm_dir, "tinker.key"))
        else:
            error_message = "Initial tinker.key file not given, check file name!"
            error_vars = f"key_file = {self.key_file}"
            raise FileNotFoundError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        # Write xyz file using current positions
        if (self.scheme == "additive"):

            # Atoms in MM region
            input_xyz2 = ""

            input_xyz2 += f" {molecule.nat_mm}\n"
            if (self.l_periodic):
                input_xyz2 += " ".join([f"{ i:12.6f}" for i in self.cell_par]) + "\n"

            # Read 'tinker.xyz' file to obtain atom type and topology
            file_name = "tinker.xyz"
            with open(file_name, "r") as f_xyz:
                lines = f_xyz.readlines()
                iline = 1
                for line in lines:
                    # Skip first line
                    if (iline > molecule.nat_qm + 1):
                        ind = iline - 1
                        input_geom = f" {ind - molecule.nat_qm:6d}"
                        input_geom += f" {molecule.symbols[ind - 1]:4}"
                        input_geom += "".join([f"{i:15.8f}" for i in molecule.pos[ind - 1] * au_to_A])
                        # Count index of column in coordinate lines
                        col = 0
                        field = line.split()
                        for element in field:
                            if (col == 5):
                                # Read right block; atom type information
                                input_geom += f" {element:6s}"
                            elif (col > 5):
                                # Read right block; topology information
                                input_geom += f" {int(element) - molecule.nat_qm:6d}"
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
            if (self.l_periodic):
                input_xyz12 += " ".join([f"{ i:12.6f}" for i in self.cell_par]) + "\n"

            # Read 'tinker.xyz' file to obtain atom type and topology
            file_name = "tinker.xyz"
            with open(file_name, "r") as f_xyz:
                lines = f_xyz.readlines()
                iline = 1
                for line in lines:
                    # Skip first line
                    if (iline > 1):
                        ind = iline - 1
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
            if (self.l_periodic):
                input_xyz1 += " ".join([f"{ i:12.6f}" for i in self.cell_par]) + "\n"

            # Read 'tinker.xyz' file to obtain atom type and topology
            file_name = "tinker.xyz"
            with open(file_name, "r") as f_xyz:
                lines = f_xyz.readlines()
                iline = 1
                for line in lines:
                    # Skip first line
                    if (iline in range(2, molecule.nat_qm + 2)):
                        ind = iline - 1
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
        file_be = open('tinker.key', 'r')
        file_af = open('tmp.key', 'w')
        is_charge = False
        for line in file_be:
            if ("chargeterm" in line):
                is_charge = True
                line = ""
                if (self.embedding == None):
                    line = "chargeterm none\n"
            file_af.write(line)
        # If chargeterm keyword does not exist, add chargeterm keyword to last line
        if (not is_charge and self.embedding == None):
            line = "chargeterm none\n"
            file_af.write(line)
        file_be.close()
        file_af.close()
        os.rename('tmp.key', 'tinker.key')

        # To avoid double counting, consider only charge-charge interactions between MM atoms
        if (self.embedding == "electrostatic"):
            # Save atom types for QM atoms
            qm_atom_type = set(self.atom_type[0:molecule.nat_qm])
            # Save QM atom types existing in 'tinker.key' file
            tmp_atom_type = []

            file_be = open('tinker.key', 'r')
            file_af = open('tmp.key', 'w')
            for line in file_be:
                if ("charge" in line):
                    if (line[0] != "#"):
                        field = line.split()
                        col = 0
                        for element in field:
                            if (element == "charge"):
                                ind_charge = col
                            col += 1
                        # Set charge to zero for QM atoms in electrostatic embedding
                        if (int(field[ind_charge + 1]) in qm_atom_type):
                            tmp_atom_type.append(int(field[ind_charge + 1]))
                            line = f"charge {int(field[ind_charge + 1])} 0.0\n"
                            file_af.write(line)
                    else:
                        # Write lines including comments
                        file_af.write(line)
                else:
                    # Write lines without 'charge' word
                    file_af.write(line)
            # Set charge to zero for remaining QM atoms
            for itype in qm_atom_type:
                if (not itype in tmp_atom_type):
                    line = f"charge {itype} 0.0\n"
                    file_af.write(line)
            file_af.write("\n")
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
                if (self.vdw == None):
                    line = "vdwterm none\n"
            file_af.write(line)
        # If vdwterm keyword does not exist, add vdwterm keyword to last line
        if (not is_vdw and self.vdw == None):
            line = "vdwterm none\n"
            file_af.write(line)
        file_be.close()
        file_af.close()
        os.rename('tmp.key', 'tinker.key')

        # Set periodicity for the systems
        file_be = open('tinker.key', 'r')
        file_af = open('tmp.key', 'w')
        for line in file_be:
            if ("a-axis" in line):
                line = ""
            elif ("b-axis" in line):
                line = ""
            elif ("c-axis" in line):
                line = ""
            elif ("alpha" in line):
                line = ""
            elif ("beta" in line):
                line = ""
            elif ("gamma" in line):
                line = ""
            file_af.write(line)

        if (self.l_periodic):
            # Add periodic keywords to last line when periodicity is used
            input_periodic = textwrap.dedent(f"""\
            a-axis {self.cell_par[0]}
            b-axis {self.cell_par[1]}
            c-axis {self.cell_par[2]}
            alpha {self.cell_par[3]}
            beta {self.cell_par[4]}
            gamma {self.cell_par[5]}

            """)
            file_af.write(input_periodic)

        file_be.close()
        file_af.close()
        os.rename('tmp.key', 'tinker.key')

        # Turn off unnecessary interactions in 'tinker.key' file
        file_name = "tinker.key"
        with open(file_name, "a") as f_xyz:
            input_interaction = textwrap.dedent(f"""\
            angangterm none
            chgdplterm none
            dipoleterm none
            extraterm none
            impropterm none
            imptorsterm none
            metalterm none
            mpoleterm none
            opbendterm none
            opdistterm none
            pitorsterm none
            polarizeterm none
            restrainterm none
            rxnfieldterm none
            solvateterm none
            strbndterm none
            strtorterm none
            tortorterm none
            ureyterm none

            """)
            f_xyz.write(input_interaction)

    def run_MM(self, base_dir, istep):
        """ Run Tinker calculation and save the output files to mm_log directory

            :param string base_dir: Base directory
            :param integer istep: Current MD step
        """
        # Set correct binary for Tinker calculation
        binary1 = os.path.join(self.mm_path, "testgrad")
        binary2 = os.path.join(self.mm_path, "testgrad.x")
        if (os.path.isfile(binary1)):
            mm_command = binary1
        elif (os.path.isfile(binary2)):
            mm_command = binary2
        else:
            error_message = "No testgrad (or testgrad.x) binary in the Tinker directory!"
            error_vars = f"mm_path = {self.mm_path}"
            raise FileNotFoundError (f"( {self.mm_prog}.{call_name()} ) {error_message} ( {error_vars} )")

        # OpenMP setting
        os.environ["OMP_NUM_THREADS"] = f"{self.nthreads}"
        # Run Tinker method
        if (self.scheme == "additive"):
            command = f"{mm_command} tinker.xyz.2 -k tinker.key y n n > tinker.out.2"
            os.system(command)
        elif (self.scheme == "subtractive"):
            command = f"{mm_command} tinker.xyz.12 -k tinker.key y n n > tinker.out.12"
            os.system(command)
            command = f"{mm_command} tinker.xyz.1 -k tinker.key y n n > tinker.out.1"
            os.system(command)
        # Copy the output file to 'mm_log' directory
        tmp_dir = os.path.join(base_dir, "mm_log")
        if (os.path.exists(tmp_dir)):
            if (self.scheme == "additive"):
                log_step = f"tinker.out.2.{istep + 1}"
                shutil.copy("tinker.out.2", os.path.join(tmp_dir, log_step))
            elif (self.scheme == "subtractive"):
                log_step = f"tinker.out.12.{istep + 1}"
                shutil.copy("tinker.out.12", os.path.join(tmp_dir, log_step))
                log_step = f"tinker.out.1.{istep + 1}"
                shutil.copy("tinker.out.1", os.path.join(tmp_dir, log_step))

    def extract_MM(self, molecule, bo_list, calc_force_only):
        """ Read the output files to get MM information

            :param object molecule: Molecule object
            :param integer,list bo_list: List of BO states for BO calculation
            :param boolean calc_force_only: Logical to decide whether calculate force only
        """
        if (self.scheme == "additive"):
            # Read 'tinker.out.2' file
            file_name = "tinker.out.2"
            with open(file_name, "r") as f:
                tinker_out2 = f.read()
        elif (self.scheme == "subtractive"):
            # Read 'tinker.out.12' file
            file_name = "tinker.out.12"
            with open(file_name, "r") as f:
                tinker_out12 = f.read()
            # Read 'tinker.out.1' file
            file_name = "tinker.out.1"
            with open(file_name, "r") as f:
                tinker_out1 = f.read()

        if (not calc_force_only):
            # Energy; initialize the energy at MM level
            mm_energy = 0.

            if (self.scheme == "additive"):

                tmp_e = 'Total Potential Energy :' + '\s+([-]*\S+)\s+' + 'Kcal/mole'
                energy = re.findall(tmp_e, tinker_out2)
                energy = np.array(energy, dtype=np.float64)
                mm_energy += energy[0]

            elif (self.scheme == "subtractive"):

                tmp_e = 'Total Potential Energy :' + '\s+([-]*\S+)\s+' + 'Kcal/mole'
                energy = re.findall(tmp_e, tinker_out12)
                energy = np.array(energy, dtype=np.float64)
                mm_energy += energy[0]

                tmp_e = 'Total Potential Energy :' + '\s+([-]*\S+)\s+' + 'Kcal/mole'
                energy = re.findall(tmp_e, tinker_out1)
                energy = np.array(energy, dtype=np.float64)
                mm_energy -= energy[0]

            # Add energy of MM part to total energy; kcal/mol to hartree
            for ist in range(molecule.nst):
                molecule.states[ist].energy += mm_energy * kcalmol_to_au

        # Force; initialize the force at MM level
        mm_force = np.zeros((molecule.nat, molecule.ndim))

        if (self.scheme == "additive"):

            tmp_g = 'Cartesian Gradient Breakdown over Individual Atoms :' + \
                '\n\n\s+Type\s+Atom\s+dE/dX\s+dE/dY\s+dE/dZ\s+Norm\n\n' + \
                '\s+Anlyt\s+\S+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\s+\S+\n' * molecule.nat_mm
            grad = re.findall(tmp_g, tinker_out2)
            grad = np.array(grad[0], dtype=np.float64)
            grad = grad.reshape(molecule.nat_mm, 3, order='C')
            # Tinker calculates gradient, not force
            for iat in range(molecule.nat_mm):
                mm_force[iat + molecule.nat_qm] -= grad[iat]

        elif (self.scheme == "subtractive"):

            tmp_g = 'Cartesian Gradient Breakdown over Individual Atoms :' + \
                '\n\n\s+Type\s+Atom\s+dE/dX\s+dE/dY\s+dE/dZ\s+Norm\n\n' + \
                '\s+Anlyt\s+\S+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\s+\S+\n' * molecule.nat
            grad = re.findall(tmp_g, tinker_out12)
            grad = np.array(grad[0], dtype=np.float64)
            grad = grad.reshape(molecule.nat, 3, order='C')
            # Tinker calculates gradient, not force
            mm_force -= np.copy(grad)

            tmp_g = 'Cartesian Gradient Breakdown over Individual Atoms :' + \
                '\n\n\s+Type\s+Atom\s+dE/dX\s+dE/dY\s+dE/dZ\s+Norm\n\n' + \
                '\s+Anlyt\s+\S+\s+([-]*\S+)\s+([-]*\S+)\s+([-]*\S+)\s+\S+\n' * molecule.nat_qm
            grad = re.findall(tmp_g, tinker_out1)
            grad = np.array(grad[0], dtype=np.float64)
            grad = grad.reshape(molecule.nat_qm, 3, order='C')
            # Tinker calculates gradient, not force
            for iat in range(molecule.nat_qm):
                mm_force[iat] += grad[iat]

        # Add force of MM part to total force; kcal/(mol*A) to hartree/bohr
        molecule.states[bo_list[0]].force += np.copy(mm_force) * kcalmol_to_au * au_to_A


