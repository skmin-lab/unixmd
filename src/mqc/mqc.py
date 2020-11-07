from __future__ import division
from misc import fs_to_au, call_name
import textwrap
import numpy as np

class MQC(object):
    """ Class for nuclear/electronic propagator used in MQC dynamics

        :param object molecule: molecule object
        :param integer istate: initial adiabatic state
        :param double dt: time interval
        :param integer nsteps: nuclear step
        :param integer nesteps: electronic step
        :param string propagation: propagation scheme
        :param string solver: propagation solver
        :param boolean l_pop_print: logical to print BO population and coherence
        :param boolean l_adjnac: logical to adjust nonadiabatic coupling
        :param coefficient: initial BO coefficient
        :type coefficient: double, list or complex, list
        :param string unit_dt: unit of time step (fs = femtosecond, au = atomic unit)
    """
    def __init__(self, molecule, istate, dt, nsteps, nesteps, \
        propagation, solver, l_pop_print, l_adjnac, coefficient, unit_dt):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize Molecule class object
        self.mol = molecule

        # Initialize input values
        self.istate = istate
        self.nsteps = nsteps
        self.nesteps = nesteps

        # Decide unit of time step
        if (unit_dt == 'au'):
            self.dt = dt
        elif (unit_dt == 'fs'):
            self.dt = dt * fs_to_au
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid unit for time step! {unit_dt}")

        # Check number of state and initial state
        if (self.istate >= self.mol.nst): 
            raise ValueError (f"( {self.md_type}.{call_name()} ) Index for initial state must be smaller than number of states! {self.istate}")

        # None for BOMD case
        self.propagation = propagation
        if not (self.propagation in [None, "coefficient", "density"]): 
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'propagation'! {self.propagation}")
        self.solver = solver
        if not (self.solver in [None, "rk4"]): 
            raise ValueError (f"( {self.md_type}.{call_name()} ) Invalid 'solver'! {self.solver}")

        self.l_pop_print = l_pop_print

        self.l_adjnac = l_adjnac

        self.rforce = np.zeros((self.mol.nat, self.mol.nsp))

        # Initialize coefficients and densities
        self.mol.get_coefficient(coefficient, self.istate)

    def cl_update_position(self):
        """ Routine to update nuclear positions

        """
        self.calculate_force(self.mol)

        self.mol.vel += 0.5 * self.dt * self.rforce / np.column_stack([self.mol.mass] * self.mol.nsp)
        self.mol.pos += self.dt * self.mol.vel

    def cl_update_velocity(self):
        """ Routine to update nuclear velocities

        """
        self.calculate_force(self.mol)

        self.mol.vel += 0.5 * self.dt * self.rforce / np.column_stack([self.mol.mass] * self.mol.nsp)
        self.mol.update_kinetic()

#    def calculate_temperature(self):
#        """ Routine to calculate current temperature
#        """
#        pass
#        #self.temperature = self.mol.ekin * 2 / float(self.mol.dof) * au_to_K

    def calculate_force(self):
        """ Routine to calculate the forces
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self, qm, mm, thermostat, debug):
        """ Routine to print the initial information of dynamics

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
            :param object thermostat: thermostat type
            :param integer debug: verbosity level for standard output
        """
        # Print molecule information: coordinate, velocity
        self.mol.print_init(mm)

        # Print dynamics information
        dynamics_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Dynamics Information":>43s}
        {"-" * 68}
          QM Program               = {qm.qm_prog:>16s}
          QM Method                = {qm.qm_method:>16s}
        """)
        if (self.mol.qmmm and mm != None):
            dynamics_info += textwrap.indent(textwrap.dedent(f"""\
              MM Program               = {mm.mm_prog:>16s}
              QMMM Scheme              = {mm.scheme:>16s}
            """), "  ")
            # Print charge embedding in MM program
            if (mm.embedding != None):
                dynamics_info += f"  Charge Embedding         = {mm.embedding:>16s}\n"
            else:
                dynamics_info += f"  Charge Embedding         = {'No':>16s}\n"
            # Print vdw interaction in MM program
            if (mm.vdw != None):
                dynamics_info += f"  VDW Interaction          = {mm.vdw:>16s}\n"
            else:
                dynamics_info += f"  VDW Interaction          = {'No':>16s}\n"
                                 
        dynamics_info += textwrap.indent(textwrap.dedent(f"""\

          MQC Method               = {self.md_type:>16s}
          Time Interval (fs)       = {self.dt / fs_to_au:16.6f}
          Initial State (0:GS)     = {self.istate:>16d}
          Nuclear Step             = {self.nsteps:>16d}
        """), "  ")
        if (self.md_type != "BOMD"):
            dynamics_info += f"  Electronic Step          = {self.nesteps:>16d}\n"
            dynamics_info += f"  Propagation Scheme       = {self.propagation:>16s}\n"

        # Print surface hopping variables
        if (self.md_type == "SH" or self.md_type == "SHXF"):
            dynamics_info += f"\n  Velocity Rescale in Hop  = {self.vel_rescale:>16s}\n"

        # Print XF variables
        if (self.md_type == "SHXF" or self.md_type == "EhXF"):
            # Print density threshold used in decoherence term
            dynamics_info += f"\n  Density Threshold        = {self.threshold:>16.6f}"
            if (self.md_type == "SHXF" and self.one_dim):
                # Print reduced mass
                dynamics_info += f"\n  Reduced Mass             = {self.aux.mass[0]:16.6f}"
            # Print wsigma values
            if (isinstance(self.wsigma, float)):
                dynamics_info += f"\n  Sigma                    = {self.wsigma:16.3f}\n"
            elif (isinstance(self.wsigma, list)):
                dynamics_info += f"\n  Sigma (1:N)              =\n"
                nlines = int(self.mol.nat_qm / 6)
                if (self.mol.nat_qm % 6 != 0):
                    nlines += 1
                wsigma_info = ""
                for iline in range(nlines):
                    iline1 = iline * 6
                    iline2 = (iline + 1) * 6
                    if (iline2 > self.mol.nat_qm):
                        iline2 = self.mol.nat_qm
                    wsigma_info += f"  {iline1 + 1:>3d}:{iline2:<3d};"
                    wsigma_info += "".join([f'{sigma:7.3f}' for sigma in self.wsigma[iline1:iline2]])
                    wsigma_info += "\n"
                dynamics_info += wsigma_info

        print (dynamics_info, flush=True)

        # Print thermostat information
        if (thermostat != None):
            thermostat.print_init()
        else:
            thermostat_info = "  No Thermostat: Total energy is conserved!\n"
            print (thermostat_info, flush=True)

    def check_qmmm(self, qm, mm):
        """ Routine to check compatibility between QM and MM objects

            :param object qm: qm object containing on-the-fly calculation infomation
            :param object mm: mm object containing MM calculation infomation
        """
        # Now check MM object
        if (mm.mm_prog == "Tinker"):
            # Now check QM object
            if (qm.qm_prog == "dftbplus"):
                if (qm.qm_method == "SSR"):
                    do_qmmm = True
                else:
                    do_qmmm = False
            else:
                do_qmmm = False
        else:
            do_qmmm = False

        if (do_qmmm):
            if (qm.embedding != mm.embedding):
                raise ValueError (f"( {self.md_type}.{call_name()} ) Inconsistent charge embedding options between QM and MM objects! {qm.embedding} and {mm.embedding}")
        else:
            raise ValueError (f"( {self.md_type}.{call_name()} ) Not compatible objects in QMMM! {qm.qm_prog}.{qm.qm_method} and {mm.mm_prog}")


