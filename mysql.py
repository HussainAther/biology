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

def getchromo(crms_o,end):
    """
    From an ordered list of markers,
    generate chromosomes.
    """
    chromo = [[] for r in range(len(end))]
    i = 0
    for crm_o in crms_o:
        j = 0
        if len(crm_o) > 1:
            for mark in crm_o:
                if mark==crm_o[0]: #first marker
                    chromo[i].append(("", None, mark[1]))
                    chromo[i].append((mark[0],colors.red,
                                    mark[2]-mark[1]))
                    ant = mark[2]
                elif mark==crm_o[-1]: #last marker
                    chromo[i].append(("", None, mark[1]-ant))
                    chromo[i].append((mark[0],colors.red,
                                    mark[2]-mark[1]))
                    chromo[i].append(("", None, end[i]-mark[2]))
                else:
                    chromo[i].append(("", None, mark[1]-ant))
                    chromo[i].append((mark[0],colors.red,
                                    ant=mark[2]))
        elif len(crm_o)==1: # For chromosomes with one marker
            chromo[i].append(("", None, crm_o[0][1]))
            chromo[i].append((crm_o[0][0],colors.red,
                  crm_o[0][2]-crm_o[0][1]))
            chromo[i].append(("", None, end[i]-crm_o[0][2]))
        else:
            # For chromosomes without markers
            # Add 3% of each chromosome
            chromo[i].append(("", None, int(0.03*end[i])))
            chromo[i].append(("", None, end[i]))
            chromo[i].append(("", none, int(0.03*end[i])))
        i += 1
        j += 1
    return chromo
