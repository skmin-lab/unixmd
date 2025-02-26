from __future__ import division
from misc import call_name
import os
import pickle

class Trajectory(object):
    """ Class for a trajectory object including the calculator and trajectory

        :param integer samp_index: Initial index set from sampling data for the dynamics with CPA
        :param string samp_bin_dir: Name of directory where the calculator and trajectory are saved
    """
    def __init__(self, samp_index=None, samp_bin_dir="./"):
        # Save name of Trajectory class
        self.traj_type = self.__class__.__name__

        # Initialize input variables
        self.samp_index = samp_index

        if (self.samp_index == None):
            error_message = "samp_index must be given in running script!"
            error_vars = f"samp_index = {self.samp_index}"
            raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            if (not isinstance(self.samp_index, int)):
                error_message = "samp_index must be integer!"
                error_vars = f"samp_index = {self.samp_index}"
                raise TypeError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")
                if (self.samp_index < 0):
                    error_message = "samp_index must not be smaller than 0"
                    error_vars = f"samp_index = {self.samp_index}"
                    raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.samp_bin_dir = samp_bin_dir

        # Initialize the calculator and trajectory variables
        self.energy = []
        self.force = []
        self.nacme = []

        self.pos = []
        self.vel = []

    def read_QM_from_file(self, nsteps):
        """ Routine to read energy, force, and NACME from binary files

            :param integer nsteps: Total step of nuclear propagation
        """
        for istep in range(self.samp_index + 1, self.samp_index + nsteps + 1):
            QM_path = os.path.join(self.samp_bin_dir, f"QM.{istep}.bin")
            with open(QM_path, "rb") as f:
                data = pickle.load(f)

            self.energy.append(data["energy"])
            self.force.append(data["force"])
            self.nacme.append(data["nacme"])

        QM_path = os.path.join(self.samp_bin_dir, f"QM.{self.samp_index}.bin")
        with open(QM_path, "rb") as f:
            data = pickle.load(f)

        self.energy.append(data["energy"])
        self.force.append(data["force"])
        self.nacme.append(data["nacme"])

    def read_RV_from_file(self, nsteps):
        """ Routine to read atomic position and velocity from binary files

            :param integer nsteps: Total step of nuclear propagation
        """
        for istep in range(self.samp_index + 1, self.samp_index + nsteps + 1):
            RV_path = os.path.join(self.samp_bin_dir, f"RV.{istep}.bin")
            with open(RV_path, "rb") as f:
                data = pickle.load(f)

            self.pos.append(data["pos"])
            self.vel.append(data["vel"])

        RV_path = os.path.join(self.samp_bin_dir, f"RV.{self.samp_index}.bin")
        with open(RV_path, "rb") as f:
            data = pickle.load(f)

        self.pos.append(data["pos"])
        self.vel.append(data["vel"])


