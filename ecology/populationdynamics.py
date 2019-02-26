
"""
With python, we can create a model of population dynamics using a few different approaches.

In this example, we'll look at hte Leslie and Lefkobitch matrix as they apply to the field.
"""
class LMatrix():
    """
    Allow for any number of age classes or stages.
    """
    def __init__(self, stAges):
        """
        Initialize the variabels we'll use.
        """
        import numpy as np
        import numpy.matlib as M
        from numpy.matlib import rand, zeros, ones, empty, eye
        import interval

        self.stAges = stAges # number of age/stage classes
        self.m = zeros((self.stAges, self.stAges)) # our matrix
        self.step = 0 # begin at 0
        self.popvec = None
        self.survival = None
        self.recurrence = None
        self.fecundity = None

    def LM_AddFecundity(self, fvector):
        """
        Set fecundity values for an LMatrix
        Set first row of the matrix equal to values, and find mismatches
        between the length of the vector and the width of the matrix.
        Leave those positions unchanged.
        """
