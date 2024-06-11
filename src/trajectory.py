from __future__ import division
import numpy as np

class Trajectory (object):
    """ Class to save BOMD trajectory data for CPA dynamics

        :param object molecule: Molecule object
        :param integer nsteps: Total step of nuclear propagation
    """
    def __init__(self, molecule, nsteps):
        self.pos = np.zeros((nsteps + 1, molecule.nat, molecule.ndim))
        self.vel = np.zeros((nsteps + 1, molecule.nat, molecule.ndim))
        self.energy = np.zeros((nsteps + 1, molecule.nst))
        self.force = np.zeros((nsteps + 1, molecule.nat, molecule.ndim))
        self.nacme = np.zeros((nsteps + 1, molecule.nst, molecule.nst))
        
