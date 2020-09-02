from __future__ import division
from qm.qm_calculator import QM_calculator

class Model(QM_calculator):
    """ Class for common parts of model calculation

        :param string qm_path: path for model calculation program
    """
    def __init__(self, qm_path):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize model common variables
        self.qm_path = qm_path


