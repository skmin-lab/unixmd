from __future__ import division
from bo.bo_calculator import BO_calculator

class Model(BO_calculator):
    """ Class for common parts of model calculation

        :param string qm_path: path for model calculation program
    """
    def __init__(self, qm_path):
        # Save name of BO calculator and its method
        super().__init__()

        # Initialize model common variables
        self.qm_path = qm_path


