from __future__ import division
from misc import data, eps, A_to_au, call_name
import textwrap
import numpy as np

class State(object):
    """ Class for BO states

        :param integer nsp: dimension of space where the molecule is
        :param integer nat: number of atoms
    """
    def __init__(self, nsp, nat):
        # Initialize variables
        self.energy = 0.
        self.energy_old = 0.
        self.force = np.zeros((nat, nsp))
        self.coef = 0. + 0.j
        self.multiplicity = 1


class Molecule(object):
    """ Class for a molecule object including State objects

        :param string geometry: initial cartesian coordinates for position and velocities in the extended xyz format
        :param integer nsp: dimension of space where the molecule is
        :param integer nstates: number of BO states
        :param integer dof: degrees of freedom (if model is False, molecular dof is given)
        :param string unit_pos: unit of position (A = angstrom, au = atomic unit [bohr])
        :param string unit_vel: unit of velocity (au = atomic unit, A/ps = angstrom per ps, A/fs = angstromm per fs)
        :param double charge: total charge of the system
        :param boolean model: is the system a model system?
    """
    def __init__(self, geometry, nsp=3, nstates=3, dof=None, \
        unit_pos='A', unit_vel='au', charge=0., model=False):
        # Save name of Molecule class
        self.mol_type = self.__class__.__name__

        # Initialize input values
        self.nsp = nsp
        self.nst = nstates
        self.model = model

        # Initialize geometry
        self.pos = []
        self.vel = []
        self.mass = []
        self.symbols = []
        self.read_geometry(geometry, unit_pos, unit_vel)

        # Initialize system charge and number of electrons
        if (not model):
            self.charge = charge
            self.get_nr_electrons()
        else:
            self.charge = 0.
            self.nelec = 0

        # Initialize degrees of freedom
        if (model):
            if (dof == None):
                self.dof = self.nat * self.nsp
            else:
                self.dof = dof
        else:
            if (dof == None):
                if (self.nat == 1):
                    raise ValueError (f"( {self.mol_type}.{call_name()} ) Too small number of atoms! {self.nat}")
                elif (self.nat == 2):
                    # Diatomic molecules
                    self.dof = 1
                else:
                    # Non-linear molecules
                    self.dof = self.nsp * self.nat - self.nsp * (self.nsp + 1) / 2
            else:
                self.dof = dof

        # Initialize BO states
        self.states = []
        for ist in range(self.nst):
            self.states.append(State(self.nsp, self.nat))

        # Initialize couplings
        self.nacme = np.zeros((self.nst, self.nst))
        self.nacme_old = np.zeros((self.nst, self.nst))
        self.socme = np.zeros((self.nst, self.nst), dtype=np.complex_)
        self.socme_old = np.zeros((self.nst, self.nst), dtype=np.complex_)

        # Initialize other properties
        self.nac = np.zeros((self.nst, self.nst, self.nat, self.nsp))
        self.nac_old = np.zeros((self.nst, self.nst, self.nat, self.nsp))
        self.rho = np.zeros((self.nst, self.nst), dtype=np.complex_)

        self.ekin = 0.
        self.epot = 0.
        self.etot = 0.

        self.l_nacme = False

    def read_geometry(self, geometry, unit_pos, unit_vel):
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
            :param string unit_pos: unit of position (A = angstrom, au = atomic unit [bohr])
            :param string unit_vel: unit of velocity (au = atomic unit, A/ps = angstrom per ps, A/fs = angstromm per fs)
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
                assert (len(line.split()) == (1 + 2 * self.nsp))
                self.symbols.append(line.split()[0])
                self.mass.append(data[line.split()[0]])
                self.pos.append(list(map(float, line.split()[1:(self.nsp + 1)])))
                self.vel.append(list(map(float, line.split()[(self.nsp + 1):])))
                count_line += 1
        assert (self.nat == count_line - 2)

        self.symbols = np.array(self.symbols)
        self.mass = np.array(self.mass)

        # Conversion unit
        if (unit_pos == 'au'):
            fac_pos = 1.
        elif (unit_pos == 'A'):
            fac_pos = A_to_au
        else:
            raise ValueError (f"( {self.mol_type}.{call_name()} ) Invalid unit for position! {unit_pos}")

        self.pos = np.array(self.pos) * fac_pos
        if (unit_vel == 'au'):
            fac_vel = 1.
        elif (unit_vel == 'A/ps'):
            fac_vel = A_to_au / (1000.0 * fs_to_au)
        elif (unit_vel == 'A/fs'):
            fac_vel = A_to_au / fs_to_au
        else:
            raise ValueError (f"( {self.mol_type}.{call_name()} ) Invalid unit for velocity! {unit_vel}")

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

        self.nac_old = np.copy(self.nac)

    def get_nacme(self):
        """ Get NACME from nonadiabatic couplings
        """
        self.nacme = np.sum(np.multiply(self.nac, self.vel), axis = (2, 3))

    def update_kinetic(self):
        """ Get kinetic energy
        """
        self.ekin = np.sum(0.5 * self.mass * np.sum(self.vel ** 2, axis=1))

    def backup_bo(self):
        """ Backup BO energies and nonadiabatic couplings
        """
        for states in self.states:
            states.energy_old = states.energy
        self.nacme_old = np.copy(self.nacme)

    def get_nr_electrons(self):
        """ Get the number of electrons
        """
        sym_list = list(data.keys())
        self.nelec = 0.
        for iat in range(self.nat):
            self.nelec += float(sym_list.index(self.symbols[iat]))

        self.nelec -= self.charge

    def print_init(self):
        """ Print initial information about molecule.py
        """
        geom_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Initial Coordinate (au)":>45s}
        {"-" * 68}
        {"X":>16s}{"Y":>15s}{"Z":>15s}{"Mass":>16s}
        """)

        for nth, atoms in enumerate(self.symbols):
            geom_info += f"  {atoms:3s}"
            for isp in range(self.nsp):
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
            for isp in range(self.nsp):
                vel_info += f"{self.vel[nth, isp]:15.8f}"
            vel_info += f"\n"
        print (vel_info, flush=True)

        ### TODO: multiplicity
        molecule_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Molecule Information":>43s}
        {"-" * 68}
          Number of Atoms          = {self.nat:>16d}
          Degrees of Freedom       = {int(self.dof):>16d}
          Charge                   = {int(self.charge):>16d}
          Number of Electrons      = {int(self.nelec):>16d}
          Number of States         = {self.nst:>16d}
        """)
        ### TODO: Model case
        print (molecule_info, flush=True)


