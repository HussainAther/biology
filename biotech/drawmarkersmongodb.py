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
