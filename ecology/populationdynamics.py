
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
        if (fvector.shape[0] == self.stAges):
            self.m[0] = fvector
            self.fecundity = fvector
        else:
            print("Mismatch in size: %s vs. %s" %(self.stAges-1, fvector.shape[0]))

    def LM_AddSurvival(self, survival):
        """
        Add the values for survival that shift population members from one age/stage
        to the next. The values come in as the "survival vector" which is a numpy array.
        """
        if survival.shape[0] == (self.stAges -1)):
            for j in range(1, self.stAges):
                self.m[j, j-1] = survival[j-1]
                self.survival = survival
        else:
            print("Mismatch in size: %svs %s" %(self.stAges -1, survival.shape[0]))

    def LM_AddRecurrence(self, recur):
        """
        Add the values for survival of organisms remaining in the same stage.
        This is for stage-structured population models only. The input as the vector recur, and
        its values replace those in the m matrix along the main diagonal from [1, 1] to [N-1, N-1].
        """
        if (recur.shape[0] == (self.stAges -1)):
            for i in range(1, self.stAges):
                self.m[i, i] = recur[o-1]
                self.recurrence = recur
        else:
            print("Mismatch in size %s vs. %s" % (self.stAges - 1, recur.shape[0]))

    def LM_SetOneRelation(self, fromState, toState, value):
        """
        Set a relation that doens't fall on the survivla diagonal or the recurrence diagonal.
        This is useful for more complex stage-structured population modeling in which organisms
        from one stage may graduate to multiple other stages at defined rates.
        """
        i = interval.Interval.between(0, self.stAges-1)
        if (fromState in i) and (toState in i):
            print("current matrix " + str(self.m))
            self.m[toState, fromState] = value
            print("with interval " + str(self.m))
