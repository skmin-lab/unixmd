from __future__ import division

class Trajectory(object)
    """ Class to save BOMD trajectory data for CPA dynamics

        :param object molecule: Molecule object
        :param integer cpa_index: Initial index for running MQC dynamics with sampling data
        :param string samp_dir: Path of sampling data directory
    """
    def __init__(self, molecule, cpa_index=None, samp_dir="./"):
        # Save name of Trajectory class
        self.traj_type = self.__class__.__name__

        #Initialize trajectory variables
        self.cpa_index = cpa_index

        if (cpa_index == None):
            error_message = "cpa_index isn't defined!"
            error_vars = f"cpa_index = {self.cpa_index}"
            raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")
        else:
            if (cpa_index < 0):
                error_message = "cpa_index shouldn't be smaller than 0"
                error_vars = f"cpa_index = {self.cpa_index}"
                raise ValueError (f"( {self.traj_type}.{call_name()} ) {error_message} ( {error_vars} )")

        self.samp_dir = samp_dir
    
        self.pos = []
        self.vel = []
        self.energy = []
        self.force = [] 
        self.nacme = []

