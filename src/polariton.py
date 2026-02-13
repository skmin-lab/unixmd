from __future__ import division
from misc import data, eps, A_to_au, fs_to_au, eV_to_au, call_name
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
        self.force = np.zeros((nat, ndim))
        self.multiplicity = 1


class Polaritonic_State(object):
    """ Class for polaritonic states

        :param integer ndim: Dimension of space
        :param integer nat: Number of atoms
    """
    def __init__(self, ndim, nat):
        # Initialize variables
        self.energy = 0.
        self.energy_old = 0.
        self.force = np.zeros((nat, ndim))
        # Initialize the electronic coefficients (adiabatic and diabatic)
        self.coef_a = 0. + 0.j
        self.coef_d = 0. + 0.j


class Polariton(object):
    """ Class for a polariton object including State and Polaritonic_State objects

        :param string geometry: A string containing atomic positions and velocities
        :param integer ndim: Dimension of space
        :param integer nstates: Number of BO states
        :param boolean l_qmmm: Use the QM/MM scheme
        :param integer natoms_mm: Number of atoms in the MM region
        :param integer ndof: Degrees of freedom (if model is False, the molecular DoF is given.)
        :param integer nphotons: Number of quantized photons inside the cavity
        :param double photon_freq: Frequency of the photon inside the cavity
        :param double,list field_pol_vec: Field polarization vector
        :param string unit_pos: Unit of atomic positions
        :param string unit_vel: Unit of atomic velocities
        :param string unit_freq: Unit of photon frequency
        :param double charge: Total charge of the system
        :param boolean l_model: Is the system a model system?
    """
    def __init__(self, geometry, ndim=3, nstates=1, l_qmmm=False, natoms_mm=None, ndof=None, \
        nphotons=1, photon_freq=0.1, field_pol_vec=None, unit_pos='angs', unit_vel='au', \
        unit_freq='ev', charge=0., l_model=False):
        # Save name of Polariton class
        self.pol_type = self.__class__.__name__

        # Initialize input values
        self.ndim = ndim
        self.nst = nstates
        self.l_model = l_model

        # Conversion unit
        self.unit_pos = unit_pos.lower()
        if not (self.unit_pos in ["angs", "au"]):
            error_message = "Invalid unit for position!"
            error_vars = f"unit_pos = {self.unit_pos}"
            raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.unit_vel = unit_vel.lower()
        if not (self.unit_vel in ["angs/ps", "angs/fs", "au"]):
            error_message = "Invalid unit for velocity!"
            error_vars = f"unit_vel = {self.unit_vel}"
            raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")

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
                raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
            self.nat_qm = self.nat - self.nat_mm
        else:
            if (self.nat_mm != None):
                error_message = "Number of atoms in MM region is not necessary!"
                error_vars = f"natoms_mm = {self.nat_mm}"
                raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
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
                    raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
                elif (self.nat == 2):
                    # Diatomic molecules
                    self.ndof = 1
                else:
                    # Non-linear molecules
                    self.ndof = self.ndim * self.nat - self.ndim * (self.ndim + 1) / 2
            else:
                self.ndof = ndof

        # Initialize polariton input values
        self.nphotons = nphotons

        # Initialize field polarization vector in the cavity
        self.field_pol_vec = field_pol_vec
        if (self.field_pol_vec == None):
            error_message = "The field polarization vector in the cavity must be set in running script!"
            error_vars = f"field_pol_vec = {self.field_pol_vec}"
            raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")

        if (isinstance(self.field_pol_vec, list)):
            if (len(self.field_pol_vec) != self.ndim):
                error_message = "Number of elements for field polarization vector must be equal to dimension of the space!"
                error_vars = f"len(field_pol_vec) = {len(self.field_pol_vec)}"
                raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            error_message = "Type of field polarization vector must be list consisting of float!"
            error_vars = f"field_pol_vec = {self.field_pol_vec}"
            raise TypeError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Unit of photon frequency
        self.unit_freq = unit_freq.lower()
        if (self.unit_freq == "au"):
            self.photon_freq = photon_freq
        elif (self.unit_freq == "ev"):
            self.photon_freq = photon_freq * eV_to_au
        else:
            error_message = "Invalid unit for photon frequency!"
            error_vars = f"unit_freq = {unit_freq}"
            raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")

        # Initialize BO states
        self.states = []
        for ist in range(self.nst):
            self.states.append(State(self.ndim, self.nat))

        # Initialize polaritonic states
        self.pst = self.nst * (self.nphotons + 1)

        self.pol_states = []
        for ist in range(self.pst):
            self.pol_states.append(Polaritonic_State(self.ndim, self.nat))

        # Initialize couplings
        self.nacme = np.zeros((self.nst, self.nst))
        self.nacme_old = np.zeros((self.nst, self.nst))
        self.socme = np.zeros((self.nst, self.nst), dtype=np.complex128)
        self.socme_old = np.zeros((self.nst, self.nst), dtype=np.complex128)

        self.pnacme = np.zeros((self.pst, self.pst))

        # Initialize other properties
        self.nac = np.zeros((self.nst, self.nst, self.nat_qm, self.ndim))
        self.nac_old = np.zeros((self.nst, self.nst, self.nat_qm, self.ndim))

        self.pnac = np.zeros((self.pst, self.pst, self.nat_qm, self.ndim))
        self.rho_a = np.zeros((self.pst, self.pst), dtype=np.complex128)
        self.rho_d = np.zeros((self.pst, self.pst), dtype=np.complex128)

        self.tdp = np.zeros((self.nst, self.nst, self.ndim))
        self.tdp_old = np.zeros((self.nst, self.nst, self.ndim))
        self.tdp_grad = np.zeros((self.nst, self.nst, 3, self.nat_qm, self.ndim))

        self.ekin = 0.
        self.ekin_qm = 0.
        self.epot = 0.
        self.etot = 0.

        self.l_nacme = False
        self.l_pnacme = False

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
        """ Adjust phase of nonadiabatic couplings (vectorized)
        """
        # Compute norms for all state pairs: sum over atom and coordinate axes
        snac_old = np.sqrt(np.sum(self.nac_old ** 2, axis=(2, 3)))  # (nst, nst)
        snac = np.sqrt(np.sum(self.nac ** 2, axis=(2, 3)))  # (nst, nst)
        # Compute dot products for all state pairs
        dot_nac = np.sum(self.nac_old * self.nac, axis=(2, 3))  # (nst, nst)
        # Compute overlap, avoiding division by zero
        norm_prod = snac * snac_old
        ovlp = np.where(norm_prod < eps, 1., dot_nac / np.maximum(norm_prod, eps))
        # Create sign matrix: -1 where overlap < 0, +1 otherwise (upper triangle)
        sign_upper = np.where(np.triu(ovlp, k=0) < 0., -1., 1.)
        # Symmetrize the sign matrix
        sign_matrix = np.triu(sign_upper, k=1) + np.triu(sign_upper, k=1).T + np.diag(np.diag(sign_upper))
        # Apply sign flip
        self.nac *= sign_matrix[:, :, np.newaxis, np.newaxis]

    def adjust_tdp(self):
        """ Adjust phase of transition dipole moments (vectorized)
        """
        # Compute norms for all state pairs: sum over coordinate axis
        stdp_old = np.sqrt(np.sum(self.tdp_old ** 2, axis=2))  # (nst, nst)
        stdp = np.sqrt(np.sum(self.tdp ** 2, axis=2))  # (nst, nst)
        # Compute dot products for all state pairs
        dot_tdp = np.sum(self.tdp_old * self.tdp, axis=2)  # (nst, nst)
        # Compute overlap, avoiding division by zero
        norm_prod = stdp * stdp_old
        ovlp = np.where(norm_prod < eps, 1., dot_tdp / np.maximum(norm_prod, eps))
        # Create sign matrix: -1 where overlap < 0, +1 otherwise (upper triangle)
        sign_upper = np.where(np.triu(ovlp, k=0) < 0., -1., 1.)
        # Symmetrize the sign matrix
        sign_matrix = np.triu(sign_upper, k=1) + np.triu(sign_upper, k=1).T + np.diag(np.diag(sign_upper))
        # Apply sign flip to tdp
        self.tdp *= sign_matrix[:, :, np.newaxis]
        # Apply sign flip to tdp_grad
        self.tdp_grad *= sign_matrix[:, :, np.newaxis, np.newaxis, np.newaxis]

    def get_nacme(self):
        """ Get NACME from nonadiabatic couplings
        """
        # Vectorized: contract nac with velocity over atom and coordinate indices
        vel_qm = self.vel[0:self.nat_qm]
        nacme_full = np.tensordot(self.nac, vel_qm, axes=([2, 3], [0, 1]))
        # Enforce antisymmetry: nacme[i,j] = -nacme[j,i], diagonal = 0
        self.nacme = np.triu(nacme_full, k=1)
        self.nacme -= self.nacme.T

    def get_pnacme(self):
        """ Get pNACME from polaritonic nonadiabatic couplings
        """
        # Vectorized: contract pnac with velocity over atom and coordinate indices
        vel_qm = self.vel[0:self.nat_qm]
        pnacme_full = np.tensordot(self.pnac, vel_qm, axes=([2, 3], [0, 1]))
        # Enforce antisymmetry: pnacme[i,j] = -pnacme[j,i], diagonal = 0
        self.pnacme = np.triu(pnacme_full, k=1)
        self.pnacme -= self.pnacme.T

    def update_kinetic(self):
        """ Get kinetic energy
        """
        self.ekin = 0.5 * np.einsum('i,ij->', self.mass, self.vel ** 2)

        if (self.l_qmmm):
            # Calculate the kinetic energy for QM atoms
            self.ekin_qm = 0.5 * np.einsum('i,ij->', self.mass[0:self.nat_qm], self.vel[0:self.nat_qm] ** 2)
        else:
            self.ekin_qm = self.ekin

    def reset_bo(self, calc_coupling, calc_tdp, calc_tdp_grad):
        """ Reset BO energies, forces, nonadiabatic couplings and transition dipoles

            :param boolean calc_coupling: Check whether the dynamics includes coupling calculation
            :param boolean calc_tdp: Check whether the dynamics includes transition dipole calculation
            :param boolean calc_tdp_grad: Check whether the dynamics includes transition dipole gradient calculation
        """
        for states in self.states:
            states.energy = 0.
            states.force = np.zeros((self.nat, self.ndim))

        if (calc_coupling):
            self.nacme = np.zeros((self.nst, self.nst))
            if (not self.l_nacme):
                self.nac = np.zeros((self.nst, self.nst, self.nat_qm, self.ndim))

        if (calc_tdp):
            self.tdp = np.zeros((self.nst, self.nst, self.ndim))

        if (calc_tdp_grad):
            self.tdp_grad = np.zeros((self.nst, self.nst, 3, self.nat_qm, self.ndim))

    def reset_qed(self, calc_coupling):
        """ Reset polaritonic state energies, forces and nonadiabatic couplings

            :param boolean calc_coupling: Check whether the dynamics includes coupling calculation
        """
        for states in self.pol_states:
            states.energy = 0.
            states.force = np.zeros((self.nat, self.ndim))

        if (calc_coupling):
            self.pnacme = np.zeros((self.pst, self.pst))
            if (not self.l_pnacme):
                self.pnac = np.zeros((self.pst, self.pst, self.nat_qm, self.ndim))

    def backup_bo(self, calc_coupling, calc_tdp):
        """ Backup nonadiabatic couplings and transition dipoles

            :param boolean calc_coupling: Check whether the dynamics includes coupling calculation
            :param boolean calc_tdp: Check whether the dynamics includes transition dipole calculation
        """
        if (calc_coupling):
            self.nacme_old = np.copy(self.nacme)
            if (not self.l_nacme):
                self.nac_old = np.copy(self.nac)

        if (calc_tdp):
            self.tdp_old = np.copy(self.tdp)

    def backup_qed(self):
        """ Backup polaritonic state energies for propagation of auxiliary trajectories
        """
        for states in self.pol_states:
            states.energy_old = states.energy

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
            
            :param coef: Initial coefficient for polaritonic states
            :type coef: double, list or complex, list
            :param integer istate: Initial polaritonic state
        """
        if (coef == None):
            if (istate == None): 
                error_message = "Either initial state or coefficients must be given!"
                error_vars = f"(MQC) istate = {istate}, (MQC) init_coef = {coef}"
                raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                if (istate >= self.pst):
                    error_message = "Index for initial state must be smaller than number of states! The index for ground state is zero"
                    error_vars = f"(MQC) istate = {istate}, nstates = {self.nst}, nphotons = {self.nphotons}"
                    raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
                else:
                    self.pol_states[istate].coef_a = 1. + 0.j
                    self.rho_a[istate, istate] = 1. + 0.j
        else:
            if (isinstance(coef, list)):
                if (len(coef) != self.pst):
                    error_message = "Number of initial coefficients must be equal to number of states!"
                    error_vars = f"(MQC) len(init_coef) = {len(coef)}, nstates = {self.nst}, nphotons = {self.nphotons}"
                    raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
                else:
                    for ist in range(self.pst):
                        if (isinstance(coef[ist], float)):
                            self.pol_states[ist].coef_a = coef[ist] + 0.j
                        elif (isinstance(coef[ist], complex)):
                            self.pol_states[ist].coef_a = coef[ist]
                        else:
                            error_message = "Type of initial coefficients must be float or complex!"
                            error_vars = f"(MQC) init_coef[{ist}] = {coef[ist]}"
                            raise TypeError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")

                    norm = 0.
                    for ist in range(self.pst):
                        for jst in range(self.pst):
                            self.rho_a[ist, jst] = self.pol_states[ist].coef_a.conjugate() * self.pol_states[jst].coef_a
                        norm += self.rho_a.real[ist, ist]

                    if (abs(norm - 1.) >= eps):
                        error_message = "Norm for electronic wave function should be 1.0!"
                        error_vars = f"(MQC) init_coef = {coef}"
                        raise ValueError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")
            else:
                error_message = "Type of initial coefficients must be list!"
                error_vars = f"(MQC) init_coef = {coef}"
                raise TypeError (f"( {self.pol_type}.{call_name()} ) {error_message} ( {error_vars} )")

    def print_init(self, mm):
        """ Print initial information about polariton.py

            :param object mm: MM object containing MM calculation information
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
        polariton_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Polariton Information":>43s}
        {"-" * 68}
          Number of Atoms (QM)     = {self.nat_qm:>16d}
        """)
        if (self.l_qmmm and mm != None):
            polariton_info += f"  Number of Atoms (MM)     = {self.nat_mm:>16d}\n"
        polariton_info += textwrap.indent(textwrap.dedent(f"""\
          Degrees of Freedom       = {int(self.ndof):>16d}
          Charge                   = {int(self.charge):>16d}
          Number of Electrons      = {int(self.nelec):>16d}
          Number of States         = {self.nst:>16d}

          Number of Photons        = {self.nphotons:>16d}
          Number of Total States   = {self.pst:>16d}
          Photon Frequency (eV)    = {self.photon_freq / eV_to_au:>16.6f}
        """), "  ")
        # Print field polarization vector components
        field_info = " " * 8
        field_info += "  ".join([f'{value:8.4f}' for value in self.field_pol_vec]) + "\n"
        polariton_info += f"  Field Polarization Vec.  =\n"
        polariton_info += field_info
        ### TODO: Model case
        print (polariton_info, flush=True)


