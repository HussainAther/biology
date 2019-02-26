
"""
With python, we can create a model of population dynamics using a few different approaches.

In this example, we'll look at hte Leslie and Lefkobitch matrix as they apply to the field.
"""
class LMatrix():
    """
    Allow for any number of age classes or stages.
    """
    def __init__(self, stAges):
        import numpy as np
        import numpy.matlib as M
        from numpy.matlib import rand, zeros, ones, empty, eye
        import interval
