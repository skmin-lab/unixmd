from __future__ import division
from misc import fs_to_au, au_to_A, call_name, typewriter
import textwrap, datetime
import numpy as np
import os, shutil


class MQC(object):
    """ Class for electronic propagator used in MQC dynamics with classical path approximation

        :param object molecule: Molecule object
        :param integer istate: Initial adiabatic state
        :param double dt: Time interval
        :param integer nsteps: Nuclear step
        :param integer nesteps: Electronic step
        :param string elec_object: Electronic equation of motions
        :param string propagator: Electronic propagator
        :param boolean l_print_dm: Logical to print BO population and coherence
        :param boolean l_adj_nac: Logical to adjust nonadiabatic coupling
        :param init_coef: Initial BO coefficient
        :type init_coef: Double, list or complex, list
        :param string unit_dt: Unit of time step (fs = femtosecond, au = atomic unit)
        :param integer out_freq: Frequency of printing output
        :param integer verbosity: Verbosity of output
    """
    def __init__(self, molecule, istate, dt, nsteps, nesteps, \
        elec_object, propagator, l_print_dm, l_adj_nac, init_coef, unit_dt, out_freq, verbosity):
        # Save name of MQC dynamics
        self.md_type = self.__class__.__name__

        # Initialize Molecule object
        self.mol = molecule

        # Initialize input values
        self.istate = istate
        self.nsteps = nsteps
        self.nesteps = nesteps

        # Initialize time step
        self.istep = -1
        self.fstep = -1

        # Decide unit of time step
        self.unit_dt = unit_dt.lower()
        if (self.unit_dt == 'au'):
            self.dt = dt
        elif (self.unit_dt == 'fs'):
            self.dt = dt * fs_to_au
        else:
            error_message = "Invalid unit for time step!"
            error_vars = f"unit_dt = {unit_dt}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.elec_object = elec_object
        if (self.elec_object != None):
            self.elec_object = self.elec_object.lower()

        if not (self.elec_object in ["coefficient", "density"]):
            error_message = "Invalid electronic object!"
            error_vars = f"elec_object = {self.elec_object}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.propagator = propagator
        if (self.propagator != None):
            self.propagator = self.propagator.lower()

        if not (self.propagator in ["rk4", "exponential"]):
            error_message = "Invalid electronic propagator!"
            error_vars = f"propagator = {self.propagator}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        if (self.propagator == "exponential" and self.elec_object != "coefficient"):
            error_message = "exponential propagator is incompatible with objects other than coefficient"
            error_vars = f"elec_object = {self.elec_object}, propagator = {self.propagator}"
            raise ValueError (f"( {self.md_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.l_print_dm = l_print_dm

        self.l_adj_nac = l_adj_nac

        self.out_freq = out_freq
        self.verbosity = verbosity

        # Initialize coefficients and densities
        self.mol.get_coefficient(init_coef, self.istate)

    def run_init(self):
        """ Initialize MQC dynamics
        """
        pass

    def update_potential(self):
        """ Routine to update the potential of molecules
        """
        pass

    def print_init(self):
        """ Routine to print the initial information of dynamics
        """
        pass

    def touch_file(self):
        """ Routine to write PyUNIxMD output files
        """
        pass

    def write_md_output(self, unixmd_dir, istep):
        """ Write output files
        """
        pass
