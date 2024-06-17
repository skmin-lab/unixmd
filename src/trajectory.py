from __future__ import division
import numpy as np
import os
import pickle

class Trajectory (object):
    """ Class to save BOMD trajectory data for CPA dynamics

        :param object molecule: Molecule object
        :param integer nsteps: Total step of nuclear propagation
        :param integer cpa_index: Initial index for running MQC dynamics with sampled data
        :param string samp_dir: Path of sampling data folder
    """
    def __init__(self, molecule, nsteps, cpa_index, samp_dir):

        #Initialize trajectory variables
        self.cpa_index = cpa_index
        
        self.samp_dir = samp_dir
        
        self.pos = np.zeros((nsteps + 1, molecule.nat, molecule.ndim))
        self.vel = np.zeros((nsteps + 1, molecule.nat, molecule.ndim))
        self.energy = np.zeros((nsteps + 1, molecule.nst))
        self.force = np.zeros((nsteps + 1, molecule.nat, molecule.ndim))
        self.nacme = np.zeros((nsteps + 1, molecule.nst, molecule.nst))

    def read_QM_from_file(self, nsteps):
        """ Routine to save precomputed atomic position, velocities for CPA dynamics
            
            :param integer nsteps: Total step of nuclear propagation
        """
        for istep in range(self.cpa_index, self.cpa_index + nsteps):
            RV_path = os.path.join(self.samp_dir, f"RP.{istep}.bin")
            with open(RV_path, "rb") as f:
                Data = pickle.load(f)

            save_step = istep - self.cpa_index
            self.pos[save_step] = Data["pos"]
            self.vel[save_step] = Data["vel"]

        RV_path = os.path.join(self.samp_dir, f"RP.{self.cpa_index-1}.bin")
        with open(RV_path, "rb") as f:
            Data = pickle.load(f)

        self.pos[nsteps] = Data["pos"]
        self.vel[nsteps] = Data["vel"]
        

    def read_RV_from_file(self):
        """ Routine to save precomputed energy, force, NACME for CPA dynamics
            
            :param integer nsteps: Total step of nuclear propagation
        """
        for istep in range(self.cpa_index, self.cpa_index + nsteps):
            QM_path = os.path.join(self.samp_dir, f"QM.{istep}.bin")
            with open(QM_path, "rb") as f:
                Data = pickle.load(f)

            save_step = istep - self.cpa_index
            self.energy[save_step] = Data["energy"]
            self.force[save_step] = Data["force"]
            self.nacme[save_step] = Data["nacme"]

        QM_path = os.path.join(self.samp_dir, f"QM.{self.cpa_index-1}.bin")
        with open(QM_path, "rb") as f:
            Data = pickle.load(f)

        self.energy[nsteps] = Data["energy"]
        self.force[nsteps] = Data["force"]
        self.nacme[nsteps] = Data["nacme"]

