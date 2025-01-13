from __future__ import division

class Trajectory(object):
    """ Class to save BOMD trajectory data for CPA
    """
    def __init__(self):
        # Save name of Trajectory class
        self.traj_type = self.__class__.__name__
