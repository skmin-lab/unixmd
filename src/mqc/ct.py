from __future__ import division
import numpy as np
from mqc.mqc import MQC

"""
   coupled-trajectory mixed quantum-classical dynamics
"""
class CT(MQC):

    def __init__(self, molecule, istate=1, dt=0.5, nsteps=1000, nesteps=10000, \
        propagation="density", thermostat={"thermo":'none', "nrescale":0}, \
        l_adjnac=True, l_deco=False ):

        # Initialize input values
        self.l_deco = l_deco

        super().__init__(molecule, istate, dt, nsteps, nesteps, \
            propagation, thermostat, l_adjnac)

        pass     

    """
       Routine related to the SH works
    """
    def step(self):
        pass
 
    """
       Routine to calculate the forces
    """
    def calculate_force(self, states):
        super().calculate_force()
        pass




