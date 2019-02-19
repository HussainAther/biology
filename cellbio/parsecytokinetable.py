import numpy
import warnings
import pandas as pd
import os
import pylab

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
    Conver tmissing data from strings to code numbers for analysis.
    """
    headers = table[0,:]
    cytokine_mapper = {}
    for i, label in enumerate(headers[1:],):
        cy, day = label.split("day")
        cytokine_mapper.setdefault(cytokine.srtip(), {})[int(day)] = int
    return cytokine_mapper



def parse_cytokine_table(filename):
    # extract data and sample and cytokine mapping dictinaries from table
    table = read_table(filename)
    warning_dict = parse_warnings(table)
    sample_mapper = parse_samples(table)
    cytokine_mapper = parse_headers(table)
    data = parse_data(table, warning_dict)
    return warning_dict, sample_mapper, cytokine_mapper, data
