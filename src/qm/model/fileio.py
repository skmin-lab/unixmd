from __future__ import division
from qm.model.model import Model
import numpy as np

class File_IO(Model):
    """ Class for fetching pre-calculated sampling data used in CPA dynamics

        :param object molecule: molecule object
        :param object trajectory: trajectory object
    """
    def __init__(self, molecule, trajectory):
        # Initialize model common variables
        super(File_IO, self).__init__(None)

        # Initialize trajectory object
        self.trajectory = trajectory

        # NACME is pre-caculated during BOMD sampling, so we do not need to get NACME
        molecule.l_nacme = False

        # Gradient of several states are pre-calculated during BOMD sampling
        self.re_calc = False

    def get_data(self, molecule, trajectory, base_dir, bo_list, dt, istep, calc_force_only):
        """ Extract energy, force and nonadiabatic coupling matrix elements from pre-calculated sampling data
            :param object molecule: molecule object
            :param object trajectory: trajectory object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
        """
        for ist in range(molecule.nst):
            molecule.states[ist].energy = np.copy(trajectory.energy[istep][ist])
        molecule.states[0].force = np.copy(trajectory.force[istep])
        molecule.nacme = np.copy(trajectory.nacme[istep])

