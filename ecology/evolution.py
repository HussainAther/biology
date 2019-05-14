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
 
