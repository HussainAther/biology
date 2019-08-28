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

# Each step
simple_update = np.array([[1, 0],
                          [-1, 0],
                          [0, 1],
                          [0, -1]], dtype=np.int)

def simple_propensity(params, population):
    """
    Return an array of propensities given a set of parameters
    and an array of populations.
    """
    # Unpack parameters
    beta_m, beta_p, gamma = params
    # Unpack population
    m, p = population
    return np.array([beta_m, 
                     m, 
                     beta_p * m, 
                     gamma * p])

def sample_discrete_scipy(probs):
    """
    Randomly sample an index with probability given by probs.
    """
    return st.rv_discrete(values=(range(len(probs)), probs)).rvs()

def sample_discrete(probs):
    """
    Randomly sample an index with probability given by probs.
    """
    # Generate random number
    q = np.random.rand()
    # Find index
    i = 0
    p_sum = 0.0
    while p_sum < q:
        p_sum += probs[i]
        i += 1
    return i - 1

# Make dummy probs
probs = np.array([0.1, 0.3, 0.4, 0.05, 0.15])

print("Result from scipy.stats:")
%timeit sample_discrete_scipy(probs)

print("\nResult from hand-coded method:")
%timeit sample_discrete(probs)
