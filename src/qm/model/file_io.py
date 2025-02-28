from __future__ import division
from qm.model.model import Model
import numpy as np

class File_IO(Model):
    """ Class for reading the pre-calculated data

        :param object molecule: molecule object
    """
    def __init__(self, molecule):
        # Initialize model common variables
        super(File_IO, self).__init__(None)

        # Set 'l_nacme' with respect to the computational method
        # File_IO can provide NACMEs, so we do not need to get NACVs
        molecule.l_nacme = True

        # File_IO can provide the gradient of several states simultaneously
        self.re_calc = False

    def get_data(self, molecule, base_dir, bo_list, dt, istep, calc_force_only, traj=None):
        """ Extract energy, gradient and nonadiabatic couplings from the pre-calculated data

            :param object molecule: molecule object
            :param string base_dir: base directory
            :param integer,list bo_list: list of BO states for BO calculation
            :param double dt: time interval
            :param integer istep: current MD step
            :param boolean calc_force_only: logical to decide whether calculate force only
            :param object traj: Trajectory object containing the calculator and trajectory
        """
        # Assign the calculator information from trajectory to molecule object
        for ist in range(molecule.nst):
            molecule.states[ist].energy = traj.energy[istep][ist]
        molecule.states[0].force = np.copy(traj.force[istep])
        molecule.nacme = np.copy(traj.nacme[istep])


