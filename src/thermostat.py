from __future__ import division
from misc import au_to_K, eps
from mqc.shxf import SHXF
import textwrap
import numpy as np

# In a model example, a specific routine should be devised separately.
class thermo(object):
    """ Object class for a thermostat. The type of given thermostat is classified by its subclasses.

        :param double temperature: temperature (K) set in the NVT ensemble
    """
    def __init__(self, temperature):
        # Initialize input values
        self.temp = temperature


class none(thermo):
    """ No temperature control; the NVE case
    """
    def __init__(self):
        pass

    def run(self, molecule, md):
        """ Control the temperature (dummy in this case)

            :param object molecule: molecule object
            :param object md: MQC object, the MD theory
        """
        pass

    def print_init(self):
        """ Print information about thermostat
        """
        thermostat_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Thermostat Information":>44s}
        {"-" * 68}
          Thermostat               = {"no thermostat":>16s}
        """)
        print (thermostat_info, flush=True)

        # TODO : efficient method for printing thermostat information
#        print ("NVE: Total energy is conserved!\n", flush=True)
#        print (f"Temp = {self.temp}", flush=True)


class rescale1(thermo):
    """ Rescale the velocities in a given period

        :param double temperature: the temperature (K) set in the NVT ensemble
        :param integer nrescale: period for rescaling step
    """
    def __init__(self, temperature=300.0, nrescale=20):
        # Initialize input values
        self.nrescale = nrescale
        self.istep = -1
        super().__init__(temperature)

    def run(self, molecule, md):
        """ Control the temperature

            :param object molecule: molecule object
            :param object md: MQC object, the MD theory
        """
        # TODO : p_name?
        p_name = "RESCALE1"
        self.istep += 1
        if (not (self.istep + 1) % self.nrescale == 0):
            return

        ctemp = molecule.ekin * 2 / float(molecule.dof) * au_to_K
        if (ctemp < eps):
            raise ValueError(f"{p_name} Current temperature too small or zero {ctemp}")

        alpha = np.sqrt(self.temp / ctemp)
        molecule.vel *= alpha

        # TODO : import direct md object?
        if (isinstance(md, SHXF)):
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

        # TODO : efficient method for printing thermostat information
#        print ("NVT: rescale type 1", flush=True)
#        print (f"Temperature is rescaled as \"{self.temp} K\" at each \"{self.nrescale} step\"!\n", flush=True)
#        print (f"nrescale = {self.nrescale}", flush=True)


class rescale2(thermo):
    """ Rescale the velocities when the temerature is out of a given range

        :param double temperature: the temperature (K) set in the NVT ensemble
        :param double dtemperature: the trigger temperature difference (K)
    """
    def __init__(self, temperature=300.0, dtemperature=100.0):
        # Initialize input values
        self.dtemp = dtemperature
        super().__init__(temperature)

    def run(self, molecule, md):
        """ Control the temperature

            :param object molecule: molecule object
            :param object md: MQC object, the MD theory
        """
        # TODO : p_name?
        p_name = "RESCALE2"

        ctemp = molecule.ekin * 2 / float(molecule.dof) * au_to_K
        if (ctemp < eps):
            raise ValueError(f"{p_name} Current temperature too small or zero {ctemp}")

        if (abs(self.temp - ctemp) > self.dtemp):
            alpha = np.sqrt(self.temp / ctemp)
            molecule.vel *= alpha

            # TODO : import direct md object?
            if (isinstance(md, SHXF)):
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

        # TODO : efficient method for printing thermostat information
#        print ("NVT: rescale type 2", flush=True)
#        print (f"Temperature is rescaled as \"{self.temp} K\" when temp. is ouf of range \"[{self.temp}-{self.dtemp},{self.temp}-{self.dtemp}]\"!\n", flush=True)
#        print (f"dT   = {self.dtemp}", flush=True)



