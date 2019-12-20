from __future__ import division
import numpy as np
import textwrap
from misc import data, A_to_au, eps

class State(object):
    """ Class for BO states
    """
    def __init__(self, nsp, nat):
        self.energy = 0.
        self.energy_old = 0.
        self.force = np.zeros((nat, nsp)) 
        self.coef = 0. + 0.j
        self.multiplicity = 1


class Molecule(object):
    """ Class for a molecule object including State objects.
        nsp - dimension of dynamics
        nat - number of atoms
        states - list of state objects
        nstates - number of states
    """
    def __init__(self, geometry, nsp=3, nstates=3, dof=None, \
        unit_pos='A', unit_vel='au', charge=0., model=False):
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

        # Initialize system charge and nr. of electrons
        if (not model):
            self.charge = charge
            self.get_nr_electrons()
        else:
            self.charge = 0
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
                    raise ValueError ("Too small NATOMS {self.nat}")
                elif (self.nat == 2):
                    # diatomic molecules
                    self.dof = 1
                else:
                    # non-linear molecules
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
        """ Routine to read the geometry in extended xyz format.
            Example:

            geometry = '''
                       2
                       Hydrogen
                       H 0.0 0.0 0.0 0.0 0.0 0.0
                       H 0.0 0.0 0.8 0.0 0.0 0.0
                       '''
            self.read_geometry( geometry )

            geometry - docstring in extended xyz format
        """
        f = geometry.split('\n')
        # remove any empty lines
        f = filter(None,f)
        for line_number, line in enumerate(f):
            if line_number == 0:
                assert len(line.split()) == 1
                self.nat = int(line.split()[0])
            elif (line_number == 1):
                pass
            else:
                if len(line.split()) == 0: break
                assert len(line.split()) == (1 + 2 * self.nsp)
                self.symbols.append(line.split()[0])
                self.mass.append(data[line.split()[0]])  
                self.pos.append(list(map(float, line.split()[1:(self.nsp + 1)])))
                self.vel.append(list(map(float, line.split()[(self.nsp + 1):]))) 

        self.symbols = np.array(self.symbols)
        self.mass = np.array(self.mass)

        # unit conversion
        if (unit_pos == 'au'):
            fac_pos = 1.
        elif (unit_pos == 'A'):
            fac_pos = A_to_au
        else:
            raise ValueError ("Invalid unit of position")

        self.pos = np.array(self.pos) * fac_pos
        if (unit_vel == 'au'):
            fac_vel = 1.
        elif (unit_vel == 'A/ps'):
            fac_vel = A_to_au / (1000.0 * fs_to_au)
        elif (unit_vel == 'A/fs'):
            fac_vel = A_to_au / fs_to_au
        else:
            raise ValueError ("Invalid unit of velocity")

        self.vel = np.array(self.vel) * fac_vel

    def adjust_nac(self):
        """ Adjust phase of NACs
        """
        p_name = "ADJUST_NAC"
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

                if(ovlp < 0.):
                    #print(f"{p_name} : the sign of NAC changed {ist} {jst}")
                    self.nac[ist, jst] = - self.nac[ist, jst]
                    self.nac[jst, ist] = - self.nac[jst, ist]

        self.nac_old = np.copy(self.nac) 

    def get_nacme(self):
        """ get NACME from NACs
        """
        self.nacme = np.sum(np.multiply(self.nac, self.vel), axis = (2, 3))

    def update_kinetic(self):
        """ get kinetic energy
        """
        self.ekin = np.sum(0.5 * self.mass * np.sum(self.vel ** 2, axis=1))

    def backup_bo(self):
        """ backup BO energies and couplings
        """
        for states in self.states:
            states.energy_old = states.energy
        self.nacme_old = np.copy(self.nacme)

    def get_nr_electrons(self):
        """ get the number of electrons
        """
        sym_list = list(data.keys())
        self.nelec = 0.
        for iat in range(self.nat):
            self.nelec += float(sym_list.index(self.symbols[iat]))

        self.nelec -= self.charge

    def print_init(self):
        """ print initial information about molecule.py
        """
        geom_info = textwrap.dedent(f"""\
        {'-'*118}
        {'Initial coordinate and velocity':>69}
        {'-'*118}
        {'':>8}{'X':^18}{'Y':^12}{'Z':^18}{'Vel_X':^12}{'Vel_Y':^18}{'Vel_Z':^12}{'Mass(H)':^20}
        """)
        for nth, atoms in enumerate(self.symbols):
            geom_info += f"{atoms:^6s}"
            for isp in range(self.nsp):
                geom_info += f"{self.pos[nth][isp]:15.8f}"
            for isp in range(self.nsp):
                geom_info += f"{self.vel[nth][isp]:15.8f}"
            geom_info += f"{self.mass[nth]:18.8f}\n"
        geom_info += f"{'-'*118}\n"
        print(geom_info, flush=True)
        
        molecule_info = textwrap.dedent(f"""\
        Charge  = {int(self.charge)}
        Total number of atoms   = {self.nat}
        Total number of states  = {self.nst}

        Model   = {self.model}
        Degrees of freedom   = {self.dof}
        """)
        print(molecule_info, flush=True)
