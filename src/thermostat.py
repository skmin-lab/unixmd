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

# TODO: unit of coup_param
class Berendsen(thermo):
    """ Rescale the velocities by Berendsen thermostat
        
        :param double temperature: the temperature (K) set in the NVT ensemble
        :param double coup_prm: the coupling parameter
    """
    def __init__(self, temperature=300.0, coup_prm=None):
        # Initialize 
        super().__init__(temperature)
        self.coup_prm = coup_prm
        
    def run(self, molecule, md):
        """
        """
        ctemp = molecule.ekin * 2 /float(molecule.dof) * au_to_K
        alpha = np.sqrt(1.0 + (md.dt * self.coup_prm) * (self.temp/ctemp - 1.0))
        
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
          Thermostat               = {"Berendsen":>16s}
          Target Temperature (K)   = {self.temp:>16.3f}
          Coupling parameter       = {self.coup_prm:>16.3f}
        """)
        print (thermostat_info, flush=True)

class NHC(thermo):
    """ Rescale the velocities by Nose-Hoover chain thermostat
        
        :param double temperature: the temperature (K) set in the NVT ensemble
        :param double coup_prm: the coupling parameter
        :param integer chainL: the number of particles in the thermostat chain
        :param integer order: the order in the thermostat chain
        :param integer nstep: the total propagation step
    """
    def __init__(self, temperature=300.0, coup_prm=None, chainL=3, order=3, nstep=1):
        # Initialize 
        super().__init__(temperature)
        self.coup_prm = coup_prm * 0.0000046
        self.chainL = chainL
        self.nstep = nstep

        # order can be only 3 or 5.
        if (order != 3 and order != 5):
            pass
#            raise ValueError (f"( {self.class_name}.{call_name()} ) Invalid order! {self.order}")
        else:
            self.order = order

    def run(self, molecule, md):
        """
        """
        # position: x, velocity: v, 
        x = np.zeros(self.chainL)
        v = np.zeros(self.chainL)
        g = np.zeros(self.chainL)
        q = np.zeros(self.chainL)

        w = np.zeros(self.order)
        if (len(w) == 3):
            w[0] = 1.0/(2.0 - 2.0**(1.0/3.0))
            w[1] = 1.0 - 2.0 * w[0]
            w[2] = w[0]
        else:
            w[0] = 1.0/(4.0 - 4.0**(1.0/3.0))
            w[1:4] = w[0]
            w[2] = 1.0 - 4.0 * w[0]
        
        wdti = np.zeros(self.order)
        wdti2 = np.zeros(self.order)
        wdti4 = np.zeros(self.order)
        wdti8 = np.zeros(self.order)

        wdti = w*md.dt/self.nstep
        wdti2 = wdti/2.0
        wdti4 = wdti/4.0
        wdti8 = wdti/8.0

        #TODO: restart
        x[:] = 1.0

        # particles in the chain
        npart = self.chainL 
        npart1 = npart - 1

        mass_nhc = np.zeros(npart)
  
        # unit is atomic unit
        ctemp = molecule.ekin * 2 /float(molecule.dof) 

        # mass of extended variables 1: q_1 = d.o.f*k*T/w_p**2
        q[0] = molecule.dof * ctemp / self.coup_prm**2
        for ipart in range(1,npart):
        # mass of extended variables i: q_i = k*T/w_p**2, i>1
            q[ipart] = ctemp / self.coup_prm**2

        alpha = 1.0
        g[0] = (2.0 * molecule.ekin - molecule.dof * ctemp) / q[0]
        for istep in range(self.nstep):
            for iorder in range(self.order):
                # update the thermostat velocities
                v[-1] += g[-1] * wdti4[iorder]
                for ipart in range(npart1):
                    aa = np.exp(-wdti8[iorder] * v[npart1-ipart])
                    v[npart1-ipart-1] = v[npart1-ipart-1] * aa**2 \
                                      + wdti4[iorder] * g[npart1-ipart-1] * aa
                # update the particle velocities
                aa = np.exp(-wdti2[iorder] * v[0])
                alpha *= aa

                # update the thermostat forces 
                g[0] = (alpha**2 * molecule.ekin * 2.0 - ctemp * molecule.dof)/q[0]

                # update thermostat positions
                for ipart in range(npart):
                    x[ipart] += v[ipart] * wdti2[iorder]

                # update thermostat velocities
                for ipart in range(npart1):
                    aa = np.exp(-wdti8[iorder] * v[ipart+1])
                    v[ipart] = v[ipart] * aa*2 + wdti4[iorder] * g[ipart] * aa
                    g[ipart+1] = (q[ipart] * v[ipart]**2 - ctemp)/q[ipart+1]

                v[-1] += g[-1] * wdti4[iorder]

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
          Thermostat                 = {"Nose-Hoover chain":>16s}
          Target Temperature (K)     = {self.temp:>16.3f}
          Coupling Parameter (cm^-1) = {self.coup_prm:>16.3f}
          Order                      = {self.order:>16.3f}
          Integrator Steps           = {self.nstep:>16.3f}
        """)
        print (thermostat_info, flush=True)
