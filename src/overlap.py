import copy

import numpy as np

from deap import tools

def check(a, b, dist, grads):

	weakx = a if a.fitness < b.fitness else b

    if weakx.radius > dist:
    	wps = dm.denorm(a, b, waypts)
    	wpfit = np.array(map(feval, wps))
    	
    	if all(wpfit >= b.fitness.wvalues[0]):
    		return True

    return True

def fight(a, b, abdist, grads):
	if a.fitness > b.fitness and a.radius > abdist:
		return a
	else:
		

def trialByCombat(individuals, k):
	chosen = []
    for i in xrange(k):
        aspirants = tool.selRandom(individuals, 2)
        chosen.append(max(aspirants, key=attrgetter("fitness")))
    return chosen