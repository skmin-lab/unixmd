from __future__ import division
from misc import data, eps, A_to_au, fs_to_au, call_name
import textwrap
import numpy as np

class State(object):
    """ Class for BO states

        :param integer ndim: Dimension of space
        :param integer nat: Number of atoms
    """
    def __init__(self, ndim, nat):
        # Initialize variables
        self.energy = 0.
        self.energy_old = 0.
        self.force = np.zeros((nat, ndim))
        self.coef = 0. + 0.j
        self.multiplicity = 1


class Molecule(object):
    """ Class for a molecule object including State objects

        :param string geometry: A string containing atomic positions and velocities
        :param integer ndim: Dimension of space
        :param integer nstates: Number of BO states
        :param boolean l_qmmm: Use the QM/MM scheme
        :param integer natoms_mm: Number of atoms in the MM region
        :param integer ndof: Degrees of freedom (if model is False, the molecular DoF is given.)
        :param string unit_pos: Unit of atomic positions
        :param string unit_vel: Unit of atomic velocities
        :param double charge: Total charge of the system
        :param boolean l_model: Is the system a model system?
    """
    def __init__(self, geometry, ndim=3, nstates=3, l_qmmm=False, natoms_mm=None, ndof=None, \
        unit_pos='angs', unit_vel='au', charge=0., l_model=False):
        # Save name of Molecule class
        self.mol_type = self.__class__.__name__

        # Initialize input values
        self.ndim = ndim
        self.nst = nstates
        self.l_model = l_model

        # Conversion unit
        self.unit_pos = unit_pos.lower()
        if not (self.unit_pos in ["angs", "au"]):
            error_message = "Invalid unit for position!"
            error_vars = f"unit_pos = {self.unit_pos}"
            raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.unit_vel = unit_vel.lower()
        if not (self.unit_vel in ["angs/ps", "angs/fs", "au"]):
            error_message = "Invalid unit for velocity!"
            error_vars = f"unit_vel = {self.unit_vel}"
            raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Initialize geometry
        self.pos = []
        self.vel = []
        self.mass = []
        self.symbols = []
        self.read_geometry(geometry)

        # Initialize QM/MM method
        self.l_qmmm = l_qmmm
        self.nat_mm = natoms_mm
        if (self.l_qmmm):
            if (self.nat_mm == None):
                error_message = "Number of atoms in MM region is essential for QMMM!"
                error_vars = f"natoms_mm = {self.nat_mm}"
                raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")
            self.nat_qm = self.nat - self.nat_mm
        else:
            if (self.nat_mm != None):
                error_message = "Number of atoms in MM region is not necessary!"
                error_vars = f"natoms_mm = {self.nat_mm}"
                raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")
            self.nat_qm = self.nat

        # Initialize system charge and number of electrons
        if (not self.l_model):
            self.charge = charge
            self.get_nr_electrons()
        else:
            self.charge = 0.
            self.nelec = 0

        # Initialize degrees of freedom
        if (self.l_model):
            if (ndof == None):
                self.ndof = self.nat * self.ndim
            else:
                self.ndof = ndof
        else:
            if (ndof == None):
                if (self.nat == 1):
                    error_message = "Too small number of atoms, check geometry! Or Check l_model and ndof!"
                    error_vars = f"nat = {self.nat}"
                    raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")
                elif (self.nat == 2):
                    # Diatomic molecules
                    self.ndof = 1
                else:
                    # Non-linear molecules
                    self.ndof = self.ndim * self.nat - self.ndim * (self.ndim + 1) / 2
            else:
                self.ndof = ndof

        # Initialize BO states
        self.states = []
        for ist in range(self.nst):
            self.states.append(State(self.ndim, self.nat))

        # Initialize couplings
        self.nacme = np.zeros((self.nst, self.nst))
        self.nacme_old = np.zeros((self.nst, self.nst))
        self.socme = np.zeros((self.nst, self.nst), dtype=np.complex128)
        self.socme_old = np.zeros((self.nst, self.nst), dtype=np.complex128)

        # Initialize other properties
        self.nac = np.zeros((self.nst, self.nst, self.nat_qm, self.ndim))
        self.nac_old = np.zeros((self.nst, self.nst, self.nat_qm, self.ndim))
        self.rho = np.zeros((self.nst, self.nst), dtype=np.complex128)

        self.ekin = 0.
        self.ekin_qm = 0.
        self.epot = 0.
        self.etot = 0.

        self.l_nacme = False

        # Initialize point charges for QM/MM calculations
        if (self.l_qmmm):
            self.mm_charge = np.zeros(self.nat_mm)

    def read_geometry(self, geometry):
        """ Routine to read the geometry in extended xyz format.\n
            Example:\n\n
            geometry = '''\n
                       2\n
                       Hydrogen\n
                       H 0.0 0.0 0.0 0.0 0.0 0.0\n
                       H 0.0 0.0 0.8 0.0 0.0 0.0\n
                       '''\n
            self.read_geometry(geometry)

            :param string geometry: Cartesian coordinates for position and initial velocity in the extended xyz format
        """
        f = geometry.split('\n')

        # Read the number of atoms
        l_read_nr_atoms = False
        count_line = 0
        for line_number, line in enumerate(f):
            llength = len(line.split())
            if (not l_read_nr_atoms and llength == 0):
                # Skip the blank lines
                continue
            elif (count_line == 0 and llength == 1):
                # Read the number of atoms
                l_read_nr_atoms = True
                self.nat = int(line.split()[0])
                count_line += 1
            elif (count_line == 1):
                # Skip the comment line
                count_line += 1
            else:
                # Read the positions and velocities
                if (len(line.split()) == 0):
                    break
                assert (len(line.split()) == (1 + 2 * self.ndim))
                self.symbols.append(line.split()[0])
                self.mass.append(data[line.split()[0]])
                self.pos.append(list(map(float, line.split()[1:(self.ndim + 1)])))
                self.vel.append(list(map(float, line.split()[(self.ndim + 1):])))
                count_line += 1
        assert (self.nat == count_line - 2)

        self.symbols = np.array(self.symbols)
        self.mass = np.array(self.mass)

        # Conversion unit
        if (self.unit_pos == 'au'):
            fac_pos = 1.
        elif (self.unit_pos == 'angs'):
            fac_pos = A_to_au

        self.pos = np.array(self.pos) * fac_pos

        if (self.unit_vel == 'au'):
            fac_vel = 1.
        elif (self.unit_vel == 'angs/ps'):
            fac_vel = A_to_au / (1000.0 * fs_to_au)
        elif (self.unit_vel == 'angs/fs'):
            fac_vel = A_to_au / fs_to_au

        self.vel = np.array(self.vel) * fac_vel

    def adjust_nac(self):
        """ Adjust phase of nonadiabatic couplings
        """
        for ist in range(self.nst):
            for jst in range(ist, self.nst):
                ovlp = 0.
                snac_old = 0.
                snac = 0.

                snac_old = np.sum(self.nac_old[ist, jst] ** 2)
                snac = np.sum(self.nac[ist, jst] ** 2)

                snac_old = np.sqrt(snac_old)
                snac = np.sqrt(snac)

                if (np.sqrt(snac * snac_old) < eps):
                    ovlp = 1.
                else:
                    dot_nac = 0.
                    dot_nac = np.sum(self.nac_old[ist, jst] * self.nac[ist, jst])
                    ovlp = dot_nac / snac / snac_old

                if (ovlp < 0.):
                    self.nac[ist, jst] = - self.nac[ist, jst]
                    self.nac[jst, ist] = - self.nac[jst, ist]

    def get_nacme(self):
        """ Get NACME from nonadiabatic couplings
        """
        for ist in range(self.nst):
            for jst in range(ist + 1, self.nst):
                self.nacme[ist, jst] = np.sum(self.nac[ist, jst] * self.vel[0:self.nat_qm])
                self.nacme[jst, ist] = - self.nacme[ist, jst]

    def update_kinetic(self):
        """ Get kinetic energy
        """
        self.ekin = np.sum(0.5 * self.mass * np.sum(self.vel ** 2, axis=1))

        if (self.l_qmmm):
            # Calculate the kinetic energy for QM atoms
            self.ekin_qm = np.sum(0.5 * self.mass[0:self.nat_qm] * np.sum(self.vel[0:self.nat_qm] ** 2, axis=1))
        else:
            self.ekin_qm = self.ekin

    def reset_bo(self, calc_coupling):
        """ Reset BO energies, forces and nonadiabatic couplings

            :param boolean calc_coupling: Check whether the dynamics includes coupling calculation
        """
        for states in self.states:
            states.energy = 0.
            states.force = np.zeros((self.nat, self.ndim))

        if (calc_coupling):
            if (self.l_nacme):
                self.nacme = np.zeros((self.nst, self.nst))
            else:
                self.nac = np.zeros((self.nst, self.nst, self.nat_qm, self.ndim))

    def backup_bo(self):
        """ Backup BO energies and nonadiabatic couplings
        """
        for states in self.states:
            states.energy_old = states.energy
        self.nac_old = np.copy(self.nac)
        self.nacme_old = np.copy(self.nacme)

    def get_nr_electrons(self):
        """ Get the number of electrons
        """
        sym_list = list(data.keys())
        self.nelec = 0.
        for iat in range(self.nat_qm):
            self.nelec += float(sym_list.index(self.symbols[iat]))

        self.nelec -= self.charge

    def get_coefficient(self, coef, istate):
        """ Get initial coefficient
            
            :param coef: Initial BO coefficient
            :type coef: double, list or complex, list
            :param integer istate: Initial adiabatic state
        """
        if (coef == None):
            if (istate == None): 
                error_message = "Either initial states or coefficients must be given!"
                error_vars = f"(MQC) istate = {istate}, (MQC) init_coef = {coef}"
                raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                if (istate >= self.nst):
                    error_message = "Index for initial state must be smaller than number of states! The index for ground state is zero"
                    error_vars = f"(MQC) istate = {istate}, nstates = {self.nst}"
                    raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")
                else:
                    self.states[istate].coef = 1. + 0.j
                    self.rho[istate, istate] = 1. + 0.j
        else:
            if (isinstance(coef, list)):
                if (len(coef) != self.nst):
                    error_message = "Number of initial coefficients must be equal to number of states!"
                    error_vars = f"(MQC) len(init_coef) = {len(coef)}, nstates = {self.nst}"
                    raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")
                else:
                    for ist in range(self.nst):
                        if (isinstance(coef[ist], float)):
                            self.states[ist].coef = coef[ist] + 0.j
                        elif (isinstance(coef[ist], complex)):
                            self.states[ist].coef = coef[ist]
                        else:
                            error_message = "Type of initial coefficients must be float or complex!"
                            error_vars = f"(MQC) init_coef[{ist}] = {coef[ist]}"
                            raise TypeError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")

                    norm = 0.
                    for ist in range(self.nst):
                        for jst in range(self.nst):
                            self.rho[ist, jst] = self.states[ist].coef.conjugate() * self.states[jst].coef
                        norm += self.rho.real[ist, ist]

                    if (abs(norm - 1.) >= eps):
                        error_message = "Norm for electronic wave function should be 1.0!"
                        error_vars = f"(MQC) init_coef = {coef}"
                        raise ValueError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                error_message = "Type of initial coefficients must be list!"
                error_vars = f"(MQC) init_coef = {coef}"
                raise TypeError (f"( {self.mol_type}.{call_name()} ) {error_message} ( {error_vars} )")

    def print_init(self, mm):
        """ Print initial information about molecule.py

            :param object mm: MM object containing MM calculation infomation
        """
        geom_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Initial Coordinate (au)":>45s}
        {"-" * 68}
        {"X":>16s}{"Y":>15s}{"Z":>15s}{"Mass":>16s}
        """)

        for nth, atoms in enumerate(self.symbols):
            geom_info += f"  {atoms:3s}"
            for isp in range(self.ndim):
                geom_info += f"{self.pos[nth, isp]:15.8f}"
            geom_info += f"{self.mass[nth]:15.5f}\n"
        print (geom_info, flush=True)

        vel_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Initial Velocity (au)":>44s}
        {"-" * 68}
        {"X":>16s}{"Y":>15s}{"Z":>15s}
        """)

        for nth, atoms in enumerate(self.symbols):
            vel_info += f"  {atoms:3s}"
            for isp in range(self.ndim):
                vel_info += f"{self.vel[nth, isp]:15.8f}"
            vel_info += f"\n"
        print (vel_info, flush=True)

        ### TODO: multiplicity
        molecule_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Molecule Information":>43s}
        {"-" * 68}
          Number of Atoms (QM)     = {self.nat_qm:>16d}
        """)
        if (self.l_qmmm and mm != None):
            molecule_info += f"  Number of Atoms (MM)     = {self.nat_mm:>16d}\n"
        molecule_info += textwrap.indent(textwrap.dedent(f"""\
          Degrees of Freedom       = {int(self.ndof):>16d}
          Charge                   = {int(self.charge):>16d}
          Number of Electrons      = {int(self.nelec):>16d}
          Number of States         = {self.nst:>16d}
        """), "  ")
        ### TODO: Model case
        print (molecule_info, flush=True)


