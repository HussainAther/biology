import multiprocessing
import numpy as np
import scipy.stats as st
import numba
import matplotlib.pyplot as plt
import seaborn as sns

"""
Gillespie algorithm simulation on a simple production of a protein.
"""

%load_ext line_profiler

rc = {"lines.linewidth": 2, 
      "axes.labelsize": 18, 
      "axes.titlesize": 18, 
      "axes.facecolor": "DFDFE5"}
sns.set_context("notebook", rc=rc)
sns.set_style("darkgrid", rc=rc)
