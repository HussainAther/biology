import numpy
import warnings
import pandas as pd
import os
import pylab

"""
Parse cytokine table for functions and proteins secreted by immune system cells.
"""

def read_table(filename, delimiter="\t"):
    """
    Read a tab-delimited table
    """
    return pd.read_csv(filename, sep="\t")

def parse_warnings(table):
    """
    Return warnings
    """
    data = table[1:, 1:]
    warning_dict = {}
    idx = -1
    for x in data.ravel():
        try:
            float(x)
        except:
            if x not in warning_dict:
                warning_dict[x] = idx
                idx -= 1
    return warning_dict

def parse_samples(table):
    """
    Parse the samples from the table
    """
    headers = table[0,:]
    sample_mapper = {}
    for i, sample in enumerate(samples):
        sample_mapper.setdefault(sample, []).append(i)
    return sample_mapper

def parse_headers(table):
    """
    Get and parse header information
    """
    headers = table[0,:]
    cytokine_mapper = {}
    for i, label in enumerate(headers[1:],):
        cy, day = label.split("day")
        cytokine_mapper.setdefault(cytokine.srtip(), {})[int(day)] = int
    return cytokine_mapper

def parse_data(table):
    """
    Convert tmissing data from strings to code numbers for analysis.
    """
    data = table[1:, 1:]
    for warning in warning_dict:
        data[data==warning] = warning_dict[warning]
    data = data.astype("float")
    return data

def parse_cytokine_table(filename):
    """
    Extract data and sample and cytokine mapping dictinaries from table
    """
    table = read_table(filename)
    warning_dict = parse_warnings(table)
    sample_mapper = parse_samples(table)
    cytokine_mapper = parse_headers(table)
    data = parse_data(table, warning_dict)
    return warning_dict, sample_mapper, cytokine_mapper, data

def plot_cytokine(cytokine, cytokine_mapper, data, save=False, directory="."):
    """
    Generate box and whiskers plot for distribution of a single cytokine
    """
    idx_dict = cytokine_mapper[cytokine]
    xs = data[:, sorted(idx_dict.values())]

    ys = []
    for i in range(xs.shape()):
        col = xs[:,i]
        ys.append(col[col >= 0])

    if numpy.sum(len(y) for y in ys) > 0:
        pylab.figure()
        pylab.boxplot(ys)
        pylab.xticks(range(1, len(idx_dict.keys())+1), idx_dict.keys())
        pylab.xlabel("Days")
        pylab.title(cytokine)

    if save:
        if not os.path.exists(directory):
            os.makedirs(directory)
        pylab.savefig(os.path.join(Directory, cytokine + ".png"))
