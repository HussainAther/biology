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
                    chromo[i].append(('',None,mark[1]-ant))
                    chromo[i].append((mark[0],colors.red,mark[2]-mark[1]))
                    chromo[i].append(('',None,end[i]-mark[2]))
                else:
                    chromo[i].append(('',None,mark[1]-ant))
                    chromo[i].append((mark[0],colors.red,mark[2]-mark[1]))
                    ant=mark[2]
        elif len(crm_o)==1: # For chromosomes with one marker
            chromo[i].append(('',None,crm_o[0][1]))
            chromo[i].append((crm_o[0][0],colors.red,crm_o[0][2]-crm_o[0][1]))
            chromo[i].append(('',None,end[i]-crm_o[0][2]))
        else:
            # For chromosomes without markers
            # Add 3% of each chromosome.
            chromo[i].append(('',None,int(0.03*end[i])))
            chromo[i].append(('',None,end[i]))
            chromo[i].append(('',None,int(0.03*end[i])))
        i += 1
        j += 1
    return chromo

def addends(chromo):
    """
    Adds a 3% of blank region at both ends for better76graphic output.
    """
    size = 0
    for x in chromo:
        size += x[2]
    # get 3% of size of each chromosome:
    endsize = int(float(size)*.03)
    # add this size to both ends in chromo:
    chromo.insert(0,('', None, endsize))
    chromo.append(('', None, endsize))
    return chromo

def load_chrom(chr_name):
    """ 
    Generate a chromosome with information
    """
    cur_chromosome = BasicChromosome.Chromosome(chr_name[0])
    chr_segment_info = chr_name[1]
   
    for seg_info_num in range(len(chr_segment_info)):
        label, color, scale = chr_segment_info[seg_info_num]
        # make the top and bottom telomeres
        if seg_info_num == 0:
            cur_segment = BasicChromosome.TelomereSegment()
        elif seg_info_num == len(chr_segment_info) - 1:
            cur_segment = BasicChromosome.TelomereSegment(1)
        ## otherwise, they are just regular segments
        else:
            cur_segment = BasicChromosome.ChromosomeSegment()
        if label != "":
            cur_segment.label = label
            cur_segment.label_size = 12
        if color is not None:
            cur_segment.fill_color = color
        cur_segment.scale = scale
        cur_chromosome.add(cur_segment)
    cur_chromosome.scale_num = max(END) + (max(END)*.04)
    return cur_chromosome

def dblookup(atgids):
    """ 
    Code to retrieve all marker data from name using MongoDB.
    """
    client = MongoClient(CONNECTION_STRING)
    db = client.pr4
    collection = db.markers_map4
    markers = []
    for marker in atgids:
        mrk = collection.find_one(}"marker_id": marker})
        if mrk:
            markers.append((marker, (mrk['chromosome'],
                            mrk['start'], mrk['end'])))
        else:
            print('Marker {0} is not in the DB'.format(marker))
    return markers 
