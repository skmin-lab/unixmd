from __future__ import division
from misc import eps, au_to_K, call_name, fs_to_au, cm_to_au
import textwrap
import numpy as np

# In a model example, a specific routine should be devised separately.
class Thermostat(object):
    """ Object class for a thermostat. The type of given thermostat is classified by its subclasses.

        :param double temperature: Target temperature (K) of the thermostat
    """
    def __init__(self, temperature):
        # Save name of thermostat type
        self.thermostat_type = self.__class__.__name__

        # Initialize input values
        self.temp = temperature


class Rescale1(Thermostat):
    """ Rescale the velocities in a given period

        :param double temperature: Target temperature (K) of the thermostat
        :param integer nrescale: The number of MD steps between rescalings
    """
    def __init__(self, temperature=300., nrescale=20):
        # Initialize input values
        super().__init__(temperature)
        self.nrescale = nrescale
        self.istep = -1

    def run(self, md):
        """ Control the temperature

            :param object md: MQC object, the MD theory
        """
        self.istep += 1
        if (not (self.istep + 1) % self.nrescale == 0):
            return

        ctemp = md.mol.ekin * 2 / float(md.mol.ndof) * au_to_K
        if (ctemp < eps):
            error_message = "Too small current temperature!"
            error_vars = f"current temperature {ctemp} < {eps}"
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) {error_message} ( {error_vars} )")

        alpha = np.sqrt(self.temp / ctemp)
        md.mol.vel *= alpha

        # Rescale the auxiliary velocities
        if (md.md_type in ["SHXF"]):
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


class Rescale2(Thermostat):
    """ Rescale the velocities when the temerature is out of a given range

        :param double temperature: Target temperature (K) of the thermostat
        :param double dtemperature: Threshold temperature difference (K)
    """
    def __init__(self, temperature=300., dtemperature=100.):
        # Initialize input values
        super().__init__(temperature)
        self.dtemp = dtemperature

    def run(self, md):
        """ Control the temperature

            :param object md: MQC object, the MD theory
        """
        ctemp = md.mol.ekin * 2 / float(md.mol.ndof) * au_to_K
        if (ctemp < eps):
            error_message = "Too small current temperature!"
            error_vars = f"current temperature {ctemp} < {eps}"
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) {error_message} ( {error_vars} )")

        if (abs(self.temp - ctemp) > self.dtemp):
            alpha = np.sqrt(self.temp / ctemp)
            md.mol.vel *= alpha

            # Rescale the auxiliary velocities
            if (md.md_type in ["SHXF"]):
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


class Berendsen(Thermostat):
    """ Rescale the velocities by Berendsen thermostat
        
        :param double temperature: Target temperature (K) of the thermostat
        :param double coupling_parameter: The characteristic tims to damp temperature toward target temperature.
        :param double coupling_strength: Dimensionless coupling strength for the thermostat
    """
    def __init__(self, temperature=300., coupling_parameter=None, coupling_strength=None):
        # Initialize input values
        super().__init__(temperature)

        self.coup_prm = coupling_parameter
        self.coup_str = coupling_strength

        if (self.coup_prm == None and self.coup_str == None):
            error_message = "Either coupling paramter or coupling strength must be set!"
            error_vars = f"coupling_parameter = {self.coup_prm}, coupling_strength = {self.coup_str}"
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) {error_message} ( {error_vars} )")
        elif (self.coup_prm != None and self.coup_str != None):
            error_message = "Only coupling paramter or coupling strength can be set!"
            error_vars = f"coupling_parameter = {self.coup_prm}, coupling_strength = {self.coup_str}"
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) {error_message} ( {error_vars} )")

    def run(self, md):
        """ Control the temperature

            :param object md: MQC object, the MD theory
        """
        if (self.coup_prm != None):
            self.coup_str = md.dt / (self.coup_prm * fs_to_au)

        ctemp = md.mol.ekin * 2 / float(md.mol.ndof) * au_to_K
        alpha = np.sqrt(1. + self.coup_str * (self.temp / ctemp - 1.))

        md.mol.vel *= alpha

        # Rescale the auxiliary velocities
        if (md.md_type in ["SHXF"]):
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
        """)

        if (self.coup_prm != None):
            thermostat_info += textwrap.indent(textwrap.dedent(f"""\
              Coupling parameter (fs)  = {self.coup_prm:>16.3f}
            """), "  ")

        if (self.coup_str != None):
            thermostat_info += textwrap.indent(textwrap.dedent(f"""\
              Coupling strength        = {self.coup_str:>16.3f}
            """), "  ")

        print (thermostat_info, flush=True)


class NHC(Thermostat):
    """ Rescale the velocities by Nose-Hoover chain thermostat
        
        :param double temperature: Target temperature (K) of the thermostat
        :param double coupling_strength: Coupling strength
        :param double time_scale: Coupling time scale
        :param integer chain_length: The number of particles in the NHC chain
        :param integer order: The order of evolution operator
        :param integer nsteps: NHC propagation step
    """
    def __init__(self, temperature=300., coupling_strength=None, time_scale=None, chain_length=3, order=3, nsteps=1):
        # Initialize input values
        super().__init__(temperature)

        self.coup_str = coupling_strength
        self.time_scale = time_scale

        if (self.coup_str == None and self.time_scale == None):
            error_message = "Either coupling strength or time scale must be set!"
            error_vars = f"coupling_strength = {self.coup_str}, time_scale = {self.time_scale}"
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) {error_message} ( {error_vars} )")
        elif (self.coup_str != None and self.time_scale != None):
            error_message = "Only coupling strength or time scale can be set!"
            error_vars = f"coupling_strength = {self.coup_str}, time_scale = {self.time_scale}"
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.chain_length = chain_length
        self.nsteps = nsteps

        # position: x, velocity: v, gradient: g and mass: q of extended particles
        self.x = np.zeros(self.chain_length)
        self.v = np.zeros(self.chain_length)
        self.g = np.zeros(self.chain_length)
        self.q = np.zeros(self.chain_length)

        #TODO: restart
        self.x = np.ones(self.chain_length)

        # order can be only 3 or 5.
        self.order = order
        if (self.order == 3):
            self.w = np.zeros(self.order)

            self.w[0] = 1. / (2. - 2. ** (1. / 3.))
            self.w[1] = 1. - 2. * self.w[0]
            self.w[2] = self.w[0]
        elif (self.order == 5):
            self.w = np.zeros(self.order)

            self.w[0] = 1. / (4. - 4. ** (1. / 3.))
            self.w[1:4] = self.w[0]
            self.w[2] = 1. - 4. * self.w[0]
        else:
            error_message = "Invalid order for NHC thermostat!"
            error_vars = f"order = {self.order}"
            raise ValueError (f"( {self.thermostat_type}.{call_name()} ) {error_message} ( {error_vars} )")

    def run(self, md):
        """ Control the temperature

            :param object md: MQC object, the MD theory
        """
        wdti = np.zeros(self.order)
        wdti2 = np.zeros(self.order)
        wdti4 = np.zeros(self.order)
        wdti8 = np.zeros(self.order)

        wdti = self.w * md.dt / self.nsteps
        wdti2 = wdti / 2.
        wdti4 = wdti / 4.
        wdti8 = wdti / 8.

        # index for last particle in list
        npart1 = self.chain_length - 1

        # target temperature: unit is atomic unit
        ttemp = self.temp / au_to_K

        if (self.coup_str != None):
            coup_prm = self.coup_str * cm_to_au
        if (self.time_scale != None):
            coup_prm = 1 / (self.time_scale * fs_to_au)

        # mass of extended variables 1: q_1 = d.o.f*k*T/w_p**2
        self.q[0] = md.mol.ndof * ttemp / coup_prm ** 2
        
        # mass of extended variables i: q_i = k*T/w_p**2, i>1
        for ipart in range(1, self.chain_length):
            self.q[ipart] = ttemp / coup_prm ** 2

        alpha = 1.
        akin = 2 * md.mol.ekin

        # update the forces
        self.g[0] = (akin - md.mol.ndof * ttemp) / self.q[0]

        # start the multiple time step procedure
        aa = 0.
        for istep in range(self.nsteps):
            for iorder in range(self.order):
                # update the thermostat velocities
                self.v[npart1] += self.g[npart1] * wdti4[iorder]
                for ipart in range(npart1, 0, -1):
                    aa = np.exp(-wdti8[iorder] * self.v[ipart])
                    self.v[ipart - 1] = self.v[ipart - 1] * aa ** 2 \
                        + wdti4[iorder] * self.g[ipart - 1] * aa

                # update the particle velocities
                aa = np.exp(-wdti2[iorder] * self.v[0])
                alpha *= aa

                # update the thermostat forces 
                self.g[0] = (alpha ** 2 * akin - ttemp * md.mol.ndof) / self.q[0]

                # update thermostat positions
                for ipart in range(self.chain_length):
                    self.x[ipart] += self.v[ipart] * wdti2[iorder]

                # update thermostat velocities
                for ipart in range(npart1):
                    aa = np.exp(-wdti8[iorder] * self.v[ipart + 1])
                    self.v[ipart] = self.v[ipart] * aa ** 2 + wdti4[iorder] * self.g[ipart] * aa
                    self.g[ipart + 1] = (self.q[ipart] * self.v[ipart] ** 2 - ttemp) / self.q[ipart + 1]
                self.v[npart1] += self.g[npart1] * wdti4[iorder]

        md.mol.vel *= alpha

        # Rescale the auxiliary velocities
        if (md.md_type in ["SHXF"]):
            md.aux.vel *= alpha
            md.aux.vel_old *= alpha

    def print_init(self):
        """ Print information about thermostat
        """
        thermostat_info = textwrap.dedent(f"""\
        {"-" * 68}
        {"Thermostat Information":>44s}
        {"-" * 68}
          Thermostat                 = {"Nose-Hoover chain":>18s}
          Target Temperature (K)     = {self.temp:>18.3f}
        """)

        if (self.coup_str != None):
            thermostat_info += textwrap.indent(textwrap.dedent(f"""\
              Coupling Strength (cm^-1)  = {self.coup_str:>18.3f}
            """), "  ")

        if (self.time_scale != None):
            thermostat_info += textwrap.indent(textwrap.dedent(f"""\
              Time Scale (fs)            = {self.time_scale:>18.3f}
            """), "  ")

        thermostat_info += textwrap.indent(textwrap.dedent(f"""\
          Chain Length               = {self.chain_length:>18d}
          Order                      = {self.order:>18d}
          Integrator Steps           = {self.nsteps:>18d}
        """), "  ")

        print (thermostat_info, flush=True)
