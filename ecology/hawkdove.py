import time

from random import choice, randint

"""
Evolutionary game theory (EGT) Hawk-dove (hawk dove hawkdove) game.
"""

startdoves = 100
starthawks = 100
startpop = starthawks + startdoves

rounds = 20
startenergy = 100;

minfoodround = 20
maxfoodround = 70

repenergy = 250
energylossround = 2
energybluffcost = 10
energylossfight = 200
livingenergy = 20

statusactive = "active"
statusasleep = "asleep"

hawk = "hawk"
dove = "dove"

agents = []

