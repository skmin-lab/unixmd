from __future__ import division
from qm.model.model import Model
from trajectory import Trajectory
import numpy as np

class FileIO(Model):
    """ Class for fetching pre-calculated sampling data used in CPA dynamics

        :param object molecule: molecule object
    """
    def __init__(self, molecule):
        # Initialize model common variables
        super(FileIO, self).__init__(None)

        # NACME is pre-caculated during BOMD sampling, so we do not need to get NACME
        molecule.l_nacme = False

        # Gradient of several states are pre-calculated during BOMD sampling
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, force and nonadiabatic coupling matrix elements from pre-caculated sampling data
            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        for ist in range(molecule.nst):
            molecule.states[ist].energy = np.copy(trajectory.energy[istep][ist])
        molecule.force = np.copy(trajectory.force[istep])
        molecule.nacme = np.copy(trajectory.nacme[istep])
