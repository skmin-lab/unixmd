from __future__ import division
from misc import au_to_K, eps
from mqc.shxf import SHXF
import numpy as np
import textwrap

# In a model example, a specific routine should be devised separately.
class thermo(object):
    """ check the thermostat (NVT)
    """
    def __init__(self, temperature):
        # Initialize input values
        self.temp = temperature

class none(thermo):
    """ do not rescale the velocities (i.e. NVE)
    """
    def __init__(self):
        pass

    def run(self, molecule, md):
        pass

    def print_init(self):
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
    """ rescale the velocities with every nrescale step
    """
    def __init__(self, temperature=300.0, nrescale=20):
        # Initialize input values
        self.nrescale = nrescale
        self.istep = -1
        super().__init__(temperature)

    def run(self, molecule, md):
        p_name = "RESCALE1"
        self.istep += 1
        if(not (self.istep + 1) % self.nrescale == 0):
            return

        ctemp = molecule.ekin * 2 / float(molecule.dof) * au_to_K
        if(ctemp < eps):
            raise ValueError(f"{p_name} Current temperature too small or zero {ctemp}")

        alpha = np.sqrt(self.temp / ctemp)
        molecule.vel *= alpha

        if(isinstance(md, SHXF)):
            md.aux.vel *= alpha
            md.aux.vel_old *= alpha

    def print_init(self):
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
    """ rescale the velocities when the temerature
        gets out of target temperature by dtemperature
    """
    def __init__(self, temperature=300.0, dtemperature=100.0):
        # Initialize input values
        self.dtemp = dtemperature
        super().__init__(temperature)

    def run(self, molecule, md):
        p_name = "RESCALE2"

        ctemp = molecule.ekin * 2 / float(molecule.dof) * au_to_K
        if(ctemp < eps):
            raise ValueError(f"{p_name} Current temperature too small or zero {ctemp}")

        if(abs(self.temp - ctemp) > self.dtemp):
            alpha = np.sqrt(self.temp / ctemp)
            molecule.vel *= alpha

            if(isinstance(md, SHXF)):
                md.aux.vel *= alpha
                md.aux.vel_old *= alpha

    def print_init(self):
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



