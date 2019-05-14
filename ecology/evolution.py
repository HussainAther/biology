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
