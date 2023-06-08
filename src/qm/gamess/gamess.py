from __future__ import division
from qm.qm_calculator import QM_calculator

class GAMESS(QM_calculator):
    """ Class for common parts of GAMESS

    """
    def __init__(self):
        # Save name of QM calculator and its method
        super().__init__()

        # Initialize GAMESS common variables


