import array
import copy
import random

import numpy as np
from scipy.spatial.distance import euclidean, pdist, squareform

from itertools import permutations as permut

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

from deap import tools
from deap import creator
from deap import base

import pdb

import dm

import sys
sys.path.append("../benchmarks/niching-benchmark-cec2013/python/")
from framework import pycec2013

def remove_overlaps(feval, individuals, hives):

	idx = range(0, len(hives))
	idx.sort(key=lambda i: individuals[i].fitness)

	nm = np.sqrt(len(individuals[0]))
	
	covs = [h.sigma*nm for h in hives]

	D = squareform(pdist(individuals))

	uniques = set(idx[:])
	# g = dm.inline_geom_sampler(1)
	# print g
	for i, j in permut(idx, 2):
		if covs[i] > D[i,j]:
			if individuals[i].fitness > individuals[j].fitness:
				wps = dm.denorm(individuals[i], individuals[j], [0.5])
				wpfit = np.array(map(feval, wps))
				if all(wpfit >= individuals[j].fitness.wvalues[0]):
					D[j,...] = np.inf
					uniques.discard(j)

	return uniques

class Test(object):

	def __init__(self, sg):
		self.sigma = sg


creator.create("Fitness", base.Fitness, weights=(1.0,))
creator.create("Individual", array.array, typecode='d',
               fitness=creator.Fitness, cpos=-1, radius=0, hive=None)

def main():

	ind_init = creator.Individual
	bench = pycec2013.Benchmark(4)

	hives = [Test(0.8) for i in range(10)]

	individuals = np.random.uniform(0,6, (10,2))
	pop = map(ind_init, individuals)
	fits = map(bench.evaluate, pop)
	for ind, fit in zip(pop, fits):
		ind.fitness.values = (fit,)	

	covs = [h.sigma*np.sqrt(2.) for h in hives]
	
	fig, ax = plt.subplots(1)
	for x, r in zip(pop, covs):
		c = mpatches.Circle(x, radius=r, alpha=0.6, fc="b", ec="b", lw=1)
		ax.add_patch(c)
	plt.scatter(individuals[:, 0], individuals[:, 1], s=100)

	for idx in remove_overlaps(bench.evaluate, pop, hives):
		plt.scatter(individuals[idx, 0], individuals[idx, 1], s=100, c='r')

	plt.show()
	print "End"

main()