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

# Graph stuff
graph_hawk_points = []
graph_dove_points = []

# Profiling


class Agent:
 	id = 0
	agent_type = None
	status = statusactive
	energy = startenergy



def main():
	init()

	current_round = 1
	death_count = 0
	dead_hawks  = 0
	dead_doves  = 0
	breed_count = 0
	main_tic = time.clock()

	while current_round <= rounds and len(agents) > 2:
		tic = time.clock()
		awakenAgents()
		food = getFood()

		# This could be optimized further by creating a list every time
		# that only has active agents, so it isn't iterating over entire list every time
		while True:
			agent, nemesis = getRandomAgents()
			if agent is None or nemesis is None: break
			compete(agent, nemesis, food)

		# Energy cost of 'living'
		for agent in agents:
			agent.energy += energylossround

		round_dead_hawks, round_dead_doves = cull()
		round_hawk_babies, round_dove_babies = breed()
		death_count += (round_dead_hawks + round_dead_doves)
		breed_count += (round_hawk_babies + round_dove_babies)


		toc = time.clock()

		print("ROUND %d" % current_round)
		print("Food produced          : %d" % food)
		print("Population             : Hawks-> %d, Doves-> %d" % (getAgentCountByType(hawk), getAgentCountByType(dove)))
		print("Dead hawks             : %d" % round_dead_hawks)
		print("Dead doves             : %d" % round_dead_doves)
		print("Hawk babies            : %s" % round_hawk_babies)
		print("Dove babies            : %s" % round_dove_babies)
		print("Hawks                  : %s" % getPercByType(hawk))
		print("Doves                  : %s" % getPercByType(dove))
		print("----")
		print("Round Processing time  : %s" % getTimeFormatted(toc - tic))
		print("Elapsed time           : %s\n" % getTimeFormatted(time.clock() - main_tic))

		# Plot
		graph_hawk_points.append(getAgentCountByType(hawk))
		graph_dove_points.append(getAgentCountByType(dove))

		current_round += 1


	main_toc = time.clock()

	print("=============================================================")
	print("Total dead agents      : %d" % death_count)
	print("Total breeding agents  : %d" % breed_count)
	print("Total rounds completed : %d" % (current_round - 1))
	print("Total population size  : %s" % len(agents))
	print("Hawks                  : %s" % getPercByType(hawk))
	print("Doves                  : %s" % getPercByType(dove))
	print("Processing time        : %s" % getTimeFormatted(main_toc - main_tic))
	print("=============================================================")


def init():

	for x in xrange(0,startdoves):
		a = Agent()
		a.agent_type = dove
		agents.append(a)

	for x2 in xrange(0,starthawks):
		a2 = Agent()
		a2.agent_type = hawk
		agents.append(a2)


