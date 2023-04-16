import os
import re

from pymongo import MongoClient
from Bio.Graphics import BasicChromosome
from reportlab.lib import colors

CONNECTION_STRING = os.getenv('MONGODB_CS', 'localhost:27017')

def sortmarkers(crms,end):
    """
    Sort markers into chromosomes.
    """
    i = 0
    crms_o = [[] for r in range(len(end))]
    crms_fo = [[] for r in range(len(end))]
    for crms in crms:
        for marker in crm:
            # Add the marker start position at each chromosome
            crms_fo[i].append(marker[1])
            crms_fo[i].sort() # Sort the marker positions.
            i += 1
        i = 0
        for order in crms_fo:
            # Using the marker order set in crms_fo, fill crms_o
            # with all the marker information.
            for pos in order:
                for mark in crms[i]:
                    try:
                        if pos==mark[1]:
                            crms_o[i].append(mark)
                    except:
                        pass
           i += 1
    return crms_o 

def getchromo(crms_o, end):
    """
    From an ordered list of markers, generate chromosomes.
    """
    chromo = [[] for r in range(len(end))]
    i = 0
    for crm_o in crms_o:
        j = 0
        if len(crm_o)>1:
            for mark in crm_o:
                if mark==crm_o[0]: # first marker
                    chromo[i].append(('',None,mark[1]))
                    chromo[i].append((mark[0],colors.red,mark[2]-mark[1]))
                    ant = mark[2]
                elif mark==crm_o[-1]: # last marker
