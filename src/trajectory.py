from __future__ import division

class Trajectory(object):
    """ Class for a trajectory object including the information obtained from BOMD

    """
    def __init__(self):
        # Save name of Trajectory class
        self.traj_type = self.__class__.__name__


