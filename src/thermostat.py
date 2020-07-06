from __future__ import division
from misc import eps, au_to_K, call_name
import textwrap
import numpy as np

# In a model example, a specific routine should be devised separately.
class thermo(object):
    """ Object class for a thermostat. The type of given thermostat is classified by its subclasses.

        :param double temperature: temperature (K) set in the NVT ensemble
    """
    def __init__(self, temperature):
        # Save name of thermostat type
        self.thermostat_type = self.__class__.__name__

        # Initialize input values
        self.temp = temperature


class rescale1(thermo):
    """ Rescale the velocities in a given period

        :param double temperature: the temperature (K) set in the NVT ensemble
        :param integer nrescale: period for rescaling step
    """
    def __init__(self, temperature=300.0, nrescale=20):
        # Initialize input values
        super().__init__(temperature)
        self.nrescale = nrescale
        self.istep = -1

    def run(self, molecule, md):
        """ Control the temperature

            :param object molecule: molecule object
            :param object md: MQC object, the MD theory
        """
        self.istep += 1
        if (not (self.istep + 1) % self.nrescale == 0):
            return

        ctemp = molecule.ekin * 2 / float(molecule.dof) * au_to_K
        if (ctemp < eps):
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) Too small current temperature! {ctemp}")

        alpha = np.sqrt(self.temp / ctemp)
        molecule.vel *= alpha

        # Rescale the auxiliary velocities for DISH-XF
        if (md.md_type == "SHXF"):
            md.aux.vel *= alpha
            md.aux.vel_old *= alpha

    def print_init(self):
        """ Print information about thermostat
        """
        thermostat_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Thermostat Information":>44s}
        {"-" * 68}
          Thermostat               = {"rescale ver. 1":>16s}
          Target Temperature (K)   = {self.temp:>16.3f}
          Rescale Step             = {self.nrescale:>16d}
        """)
        print (thermostat_info, flush=True)


class rescale2(thermo):
    """ Rescale the velocities when the temerature is out of a given range

        :param double temperature: the temperature (K) set in the NVT ensemble
        :param double dtemperature: the trigger temperature difference (K)
    """
    def __init__(self, temperature=300.0, dtemperature=100.0):
        # Initialize input values
        super().__init__(temperature)
        self.dtemp = dtemperature

    def run(self, molecule, md):
        """ Control the temperature

            :param object molecule: molecule object
            :param object md: MQC object, the MD theory
        """
        ctemp = molecule.ekin * 2 / float(molecule.dof) * au_to_K
        if (ctemp < eps):
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) Too small current temperature! {ctemp}")

        if (abs(self.temp - ctemp) > self.dtemp):
            alpha = np.sqrt(self.temp / ctemp)
            molecule.vel *= alpha

            # Rescale the auxiliary velocities for DISH-XF
            if (md.md_type == "SHXF"):
                md.aux.vel *= alpha
                md.aux.vel_old *= alpha

    def print_init(self):
        """ Print information about thermostat
        """
        thermostat_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Thermostat Information":>44s}
        {"-" * 68}
          Thermostat               = {"rescale ver. 2":>16s}
          Target Temperature (K)   = {self.temp:>16.3f}
          Temperature Range (K)    = {self.dtemp:>16.3f}
        """)
        print (thermostat_info, flush=True)


