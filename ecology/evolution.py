from string import letters
from random import choice, random

"""
Evolution (evolution) algorithm on a string.
It's like Pokemon. Except not.
"""

target  = list("METHINKS IT IS LIKE A WEASEL")
charset = letters + " "
parent  = [choice(charset) for _ in range(len(target))]
minmutaterate  = .09
C = range(100)

perfectfitness = float(len(target))
 
def fitness(trial):
    """
    Sum of matching chars by position across a trial for our target.
    """
    return sum(t==h for t,h in zip(trial, target)) 

def mutaterate():
    """
    Less mutation the closer the fit of the parent.
    """
    return 1-((perfectfitness - fitness(parent)) / perfectfitness * (1 - minmutaterate))
 
def mutate(parent, rate):
    """
    Mutate with parent and rate.
    """
    return [(ch if random() <= rate else choice(charset)) for ch in parent]
 
def what():
    """
    An appropriate answer to any question.
    """
    print("#%-4i, fitness: %4.1f%%, '%s'" % (iterations, fitness(parent)*100./perfectfitness, "".join(parent)))
 
def mate(a, b):
    """
    Mate a and b individuals.
    """
    place = 0
    if choice(xrange(10)) < 7:
        place = choice(xrange(len(target)))
    else:
        return a, b
    return a, b, a[:place] + b[place:], b[:place] + a[place:]
  
