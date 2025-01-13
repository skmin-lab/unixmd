from __future__ import division
import numpy as np
import os
from misc import call_name
import pickle

class Trajectory(object):
    """ Class to save BOMD trajectory data for CPA

        :param integer cpa_index: Initial index for running MQC dynamics with sampling data
        :param string samp_dir: Path of sampling data directory
    """
    def __init__(self, cpa_index=None, samp_dir="./"):
        # Save name of Trajectory class
        self.traj_type = self.__class__.__name__

        # Initialize trajectory variables
        self.cpa_index = cpa_index

        if (self.cpa_index == None):
            error_message = "cpa_index isn't defined!"
            error_vars = f"cpa_index = {self.cpa_index}"
            raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            if (not isinstance(self.cpa_index, int)):
                error_message = "cpa_index should be integer!"
                error_vars = f"cpa_index = {self.cpa_index}"
                raise TypeError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")
                if (self.cpa_index < 0):
                    error_message = "cpa_index shouldn't be smaller than 0"
                    error_vars = f"cpa_index = {self.cpa_index}"
                    raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.samp_dir = samp_dir

        self.pos = []
        self.vel = []
        self.energy = []
        self.force = []
        self.nacme = []

    def read_RV_from_file(self, nsteps):
        """ Routine to save precomputed atomic position, velocities for CPA dynamics
            
            :param integer nsteps: Total step of nuclear propagation
        """
        for istep in range(self.cpa_index + 1, self.cpa_index + nsteps + 1):
            RV_path = os.path.join(self.samp_dir, f"RV.{istep}.bin")

            with open(RV_path, "rb") as f:
                data = pickle.load(f)

            self.pos.append(data["pos"])
            self.vel.append(data["vel"])

        RV_path = os.path.join(self.samp_dir, f"RV.{self.cpa_index}.bin")
        with open(RV_path, "rb") as f:
            data = pickle.load(f)

        self.pos.append(data["pos"])
        self.vel.append(data["vel"])

    def read_QM_from_file(self, nsteps):
        """ Routine to save precomputed energy, force, NACME for CPA dynamics
            
            :param integer nsteps: Total step of nuclear propagation
        """
        for istep in range(self.cpa_index + 1, self.cpa_index + nsteps + 1):
            QM_path = os.path.join(self.samp_dir, f"QM.{istep}.bin")

            with open(QM_path, "rb") as f:
                data = pickle.load(f)

            self.energy.append(data["energy"])
            self.force.append(data["force"])
            self.nacme.append(data["nacme"])

        QM_path = os.path.join(self.samp_dir, f"QM.{self.cpa_index}.bin")
        with open(QM_path, "rb") as f:
            data = pickle.load(f)

        self.energy.append(data["energy"])
        self.force.append(data["force"])
        self.nacme.append(data["nacme"])


