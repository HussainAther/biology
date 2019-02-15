#!/usr/bin/env python

# Draw markers in chromosome data from MySQL database

# standard library
import os
import sys
import re

# local stuff
import MySQLdb
from Bio.Graphics import BasicChromosome

# reportlab
from reportlab.lib import colors

def sortmarkers(crms,end):
    """
    Sort markers into chromosomes
    """
    i = 0
    crms_o = [[] for r in range(len(end))]
    crms_fo = [[] for r in range(len(end))]
    for crm in crms:
        for marker in crm:
            # add the marker start position at each chromosome.
            crms_fo[i].append(marker[1])
        crms_fo[i].sort() # Sort the marker positions.
        i += 1
    i=0
    for order in crms_fo:
        # Using the marker order set in crms_fo, fill crms_o
        # with all the marker information
        for pos in order:
            for mark in crms[i]:
                try:
                    if pos == mark[1]:
                        crms_o[i].append(mark)
                except:
                    pass
        i += 1
    return crms_o
