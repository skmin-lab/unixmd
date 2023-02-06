from __future__ import division

class QED_calculator(object):
    """ Class for (cavity) quantum electrodynamics calculator
    """
    def __init__(self):
        # Save name of QED calculator
        self.qed_method = self.__class__.__name__


