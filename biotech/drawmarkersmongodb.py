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
