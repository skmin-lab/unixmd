from __future__ import division
from misc import call_name
import os
import pickle

class Trajectory(object):
    """ Class for a trajectory object including the calculator and trajectory

    """
    def __init__(self, index=None, traj_bin_dir="./"):
        # Save name of Trajectory class
        self.traj_type = self.__class__.__name__

        # Initialize input variables
        self.index = index

        if (self.index == None):
            error_message = "index must be given in running script!"
            error_vars = f"index = {self.index}"
            raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            if (not isinstance(self.index, int)):
                error_message = "index must be integer!"
                error_vars = f"index = {self.index}"
                raise TypeError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")
                if (self.index < 0):
                    error_message = "index must not be smaller than 0"
                    error_vars = f"index = {self.index}"
                    raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.traj_bin_dir = traj_bin_dir

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
        for istep in range(self.index + 1, self.index + nsteps + 1):
            QM_path = os.path.join(self.traj_bin_dir, f"QM.{istep}.bin")
            with open(QM_path, "rb") as f:
                data = pickle.load(f)

            self.energy.append(data["energy"])
            self.force.append(data["force"])
            self.nacme.append(data["nacme"])

        QM_path = os.path.join(self.traj_bin_dir, f"QM.{self.index}.bin")
        with open(QM_path, "rb") as f:
            data = pickle.load(f)

        self.energy.append(data["energy"])
        self.force.append(data["force"])
        self.nacme.append(data["nacme"])

    def read_RV_from_file(self, nsteps):
        """ Routine to read atomic position and velocity from binary files

            :param integer nsteps: Total step of nuclear propagation
        """
        for istep in range(self.index + 1, self.index + nsteps + 1):
            RV_path = os.path.join(self.traj_bin_dir, f"RV.{istep}.bin")
            with open(RV_path, "rb") as f:
                data = pickle.load(f)

            self.pos.append(data["pos"])
            self.vel.append(data["vel"])

        RV_path = os.path.join(self.traj_bin_dir, f"RV.{self.index}.bin")
        with open(RV_path, "rb") as f:
            data = pickle.load(f)

        self.pos.append(data["pos"])
        self.vel.append(data["vel"])


