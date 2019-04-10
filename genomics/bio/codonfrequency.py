from tkplot import *
from Numeric import *

"""
Plot bar charts of codon frequencies
"""

def codon_sort(a,b):
    """
    Used for sorting the codons.
    """
    if a < b:
        return -1
    elif a > b:
        return 1
    else:
        return 0

for codon in count.keys():
    if not count_random.has_key(codon):
        count_random[codon] = 0

for codon in count_random.keys():
    if not count.has_key(codon):
        count[codon] = 0

labels=count.keys()
labels.sort(codon_sort)
w1=window(plot_title=’Count codons’,width=1000)

y=array(count.values())[:len(count)/2]
x=arange(len(y)+1)
w1.bar(y,x,label=labels[:len(count)/2])
w2=window(plot_title=’Count codons(2)’,width=1000)

y=array(count.values())[(len(count)/2)+1:]
x=arange(len(y)+1)
w2.bar(y,x,label=labels[(len(count)/2)+1:])

y=array(count_random.values())[:len(count_random)/2]
x=arange(len(y)+1)
w1.bar(y,x,label=labels[:len(count_random)/2])

y=array(count_random.values())[(len(count_random)/2)+1:]
x=arange(len(y)+1)
w2.bar(y,x,label=labels[(len(count_random)/2)+1:])
