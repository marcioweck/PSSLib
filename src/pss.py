# TODO license
# TODO explain

from itertools import izip, combinations

import numpy as np

import array
import random
import copy

from scipy.spatial.distance import euclidean, pdist, squareform
import scipy.stats as stats

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

from deap import base
from deap import creator
from deap import tools

import dm
import lhs

"""
Utils
"""

def sampler(dist_engine, xbounds, dim, n, convert=False):
    rngs = [(0., 1.) for i in range(dim)]
    P = np.array(lhs.lhs(dist_engine, rngs, siz=n)).T
    P = (xbounds[1] - xbounds[0])*P + xbounds[0]  # rebound

    if convert:
        if dim == 1:
            samples = map(creator.Individual, [np.array((p, )) for p in P])
        else:
            samples = map(creator.Individual, P)
        return samples

    return P


def ensure_bounds(lbounds, ubounds):
    rand_ = random.random
    def decorator(func):
        def wrappper(*args, **kargs):
            pop = func(*args, **kargs)
            for ind in pop:
                for i in range(len(ind)):
                    if ind[i] > ubounds[i]:
                        ind[i] = ubounds[i]*1.05 if rand_() > 0.8 \
                                                    else ubounds[i]*0.9
                    elif ind[i] < lbounds[i]:
                        ind[i] = lbounds[i]*1.05 if rand_() > 0.8 \
                                                    else lbounds[i]*0.9
            return pop
        return wrappper
    return decorator


def nearest_better_tree(X, fitness, copy_X=False, marker=0.0):
    """X are edge weights of fully connected graph"""
    if copy_X:
        X = X.copy()

    if X.shape[0] != X.shape[1]:
        raise ValueError("X needs to be square matrix of edge weights")

    n_vertices = X.shape[0]
    spanning_edges = []

    # exclude self connections:
    diag_indices = np.arange(n_vertices)
    X[diag_indices, diag_indices] = np.inf

    for i in xrange(n_vertices):
        sortedline = np.argsort(X[i])
        xid = next((j for j in sortedline if fitness[j] > fitness[i]), -1)

        if xid < 0:
            spanning_edges.append(((i, i), marker))
            continue  # go to the next

        new_edge = ((i, xid), X[i, xid])
        spanning_edges.append(new_edge)

        # remove all edges inside current tree
        X[i, xid] = np.inf
        X[xid, i] = np.inf

    return spanning_edges


"""
PSS framework
"""

def assign_fitness(population, origFit, distFit):
    for ind, ofit, dfit in izip(population, origFit, distFit):
        ind.fitness.values = (ofit, dfit)


def determine_edges(population, fitness):
    D = squareform(pdist(population))
    # TODO falta remover os colegas
    # take care here, these fitnesses are raw numbers not controlled by deap
    # nearest_better_tree select using maximization
    return nearest_better_tree(D, fitness)


import pdb


def seeds_selection(feval, population, nmin, useDM=True):
    # dim = bench.ndim
#    pdb.set_trace()
    fitnesses = [x.fitness for x in population]
    edges = determine_edges(population, fitnesses)
    dists = np.array([d for _, d in edges])
    idxs = np.arange(0, len(dists))

    assign_fitness(population, fitnesses, dists)

    mu_dist = np.mean(dists)

    seeds = dict()
    # Select the main seeds
    mask = ((dists >= 2.*mu_dist) + (dists == 0))

    sidxs = idxs[mask]
    for i in sidxs:
        population[i].id = i
        seeds[i] = (population[i], 0)

    stree = dict()
    for (a, b), d in edges:
        if b in stree:
            stree[b][a] = d
        else:
            stree[b] = {a: d}

    for sid in seeds.keys():
        try:
            adjnodes = stree[sid]
            c = min(1.5*mu_dist, max(adjnodes, key=adjnodes.get))
            seeds[sid] = (population[sid], c)
        except:
            seeds[sid] = (population[sid], mu_dist)

    if useDM and len(seeds) < nmin:
        # Random search
        interior_pts = dm.inline_geom_sampler(5)

        sidxs = idxs[-mask]
        shuffled = random.sample(sidxs, len(sidxs))
        for l in shuffled:
            xy, d = edges[l]
            ind = population[xy[0]]
            dom = population[xy[1]]
            wps = dm.denorm(ind, dom, interior_pts)
            wpfit = map(feval, wps)
            #print "--------"
            #print ind.fitness
            #print dom.fitness
            #print wpfit
            #print "--------"
            if any(wpfit < ind.fitness.wvalues[0]):
                ind.id = xy[0]
                seeds[ind.id] = (ind, d*0.6)
                try:
                    seeds[dom.id] = (dom, d*0.6)
                except:
                    pass

            if len(seeds) >= nmin:
                break
    return seeds.values(), edges, population


def check_overlap(a, b, dist, waypts):
    c1 = a.radius
    c2 = b.radius

    weakx = a if a.fitness < b.fitness else b
    wps = dm.denorm(a, b, waypts)
    wpfit = map(feval, wps)

    status = 0
    if c1 > dist:
        if a.fitness < b.fitness:
            pass # Reduce a.hive.sigma
        else:
            if all(wpfit > b.fitness.wvalues[0]):
                return 2 # remove b
            
    if c2 > dist:
        if b.fitness < a.fitness:
            pass # Reduce b.hive.sigma
        else:
            if all(wpfit > a.fitness.wvalues[0]):
                return 1 # remove a

    return status
 


def remove_overlaps(individuals, D, hives, ntests=1):
    dim = len(individuals[0])
    snorm = np.sqrt(dim)

    cs = [h.sigma * snorm for h in hives]
    for i, ind in enumerate(individuals):
        ind.radius = cs[i]
        ind.hive = hives[i]
        ind.cpos = i

    in_grads = dm.inline_geom_sampler(ntests)

    A = np.zeros(D.shape)
    for i, j in combinations(range(len(individuals)), 2):
         s = check_overlap(individuals[i], individuals[j], (cs[i], cs[j])
                      D[i,j], internal_grads)
         A[i,j] = s

    new_hives = list()
    selgrp = tools.selTournament(individuals, 0.5*len(individuals), 2)

    for a, b in combinations(selgrp):
        s = check_overlap(a, b, D[a.cpos,b.cpos], in_grads)
        if s != 1:
            new_hives.append(a.hive)
        if s != 2:
            new_hives.append(b.hive)


def generate(ind_init, hives):
    swarms = [h.generate(ind_init) for h in hives]
    return swarms


def updateHive(hive, swarm):
    hive.update(swarm)
    return hive.shouldContinue()


def generateHive(hclass, (xstart, sigma)):
    return hclass(xstart, sigma)


import sys
sys.path.append("../benchmarks/niching-benchmark-cec2013/python/")
from framework import pycec2013

import cma

creator.create("Fitness", base.Fitness, weights=(1.0, 1.0))
creator.create("Individual", array.array, typecode='d',
               fitness=creator.Fitness, cpos=-1, radius=0, hive=None)
creator.create("Hive", cma.Strategy)


def main():

    benchmark = pycec2013.Benchmark(4)

    toolbox = base.Toolbox()
    toolbox.register("generate", generate, creator.Individual)
    toolbox.register("update", map, updateHive)
    toolbox.register("feval", map, benchmark.evaluate)
    toolbox.register("fdist", nearest_better_tree)
    toolbox.register("restart", seeds_selection, benchmark.evaluate,
                     creator.Individual)
    toolbox.register("hive", generateHive, cma.Strategy)


    dim = benchmark.ndim
    nmin = benchmark.ngoptima
    leftfes = benchmark.max_fes
    ngoptima = 0
    max_ngoptima = nmin

    dists = [stats.uniform for i in range(dim)]

    lbounds = tuple(benchmark.get_lbounds())
    ubounds = tuple(benchmark.get_ubounds())
    min_ = min(lbounds)
    max_ = max(ubounds)

    samples = sampler(dists, (min_, max_), dim, 20*dim)

    population = [creator.Individual(x) for x in samples]

    seeds, edges = seeds_selection(benchmark.evaluate,
                                    population, 10)

    hives = [toolbox.hive(seed) for seed in seeds]

    while leftfes > 0  and ngoptima < max_ngoptima:

        swarms = toolbox.generate(hives)
        for swarm in swarms:
            fits = toolbox.feval(swarm)
            for ind, fit in izip(swarm, fits):
                ind.fitness.values = (fit,0)

        xstarts = (hive.centroid for hive in hives)
        centroids = [creator.Individual(x) for x in xstarts]

        cfit = toolbox.feval(centroids)
        for ind, fit in izip(centroids, cfit):
            centroids.fitness.values = (cfit, 0)

        leftfes -= len(swarms)*len(swarms[0]) + len(centroids)

        ngoptima = benchmark.count_goptima(centroids, cfit, 1e-5)

        checks = toolbox.update(hives, swarms)

        nextgen = [hive for ok, hive in izip(checks, hives) if ok]

        nextgen = remove_overlaps()

        del hives[:]
        hives = nextgen

        if len(hives) < nmin:
            print "more"
            samples = sampler(dists, (min_, max_), dim, 10*(len(hives)-nmin))

            comp = np.concatenate((samples, centroids))

            seeds, _, _ = toolbox.restart(comp, nmin)

            hives.extend((toolbox.hive(s) for s in seeds
                          if s.id < len(samples)))




    plotseeds(benchmark.evaluate, min_, max_, samples=centroids)

def plotseeds(feval, min_, max_, samples=None, edges=None, seeds=None):
    fig, ax = plt.subplots()
    X = np.arange(min_, max_, 0.01)
    Y = np.arange(min_, max_, 0.01)
    X, Y = np.meshgrid(X, Y)
    PTS = np.hstack((X.reshape(X.size, 1), Y.reshape(Y.size, 1)))
    Z = np.array(map(feval, PTS), copy=True)
    Z = Z.reshape(X.shape)
    plt.contour(X, Y, Z, zdir='z', cmap=cm.jet, offset=np.min(Z))
    if edges is not None:
        mudist = np.mean([d for _, d in edges])
        for (x, y), d in edges:
            if d < 2*mudist and d > 0:
                plt.plot([samples[x, 0], samples[y, 0]],
                         [samples[x, 1], samples[y, 1]], 'k')

    if samples is not None:
        sarr = np.array(samples)
        plt.scatter(sarr[:, 0], sarr[:, 1], c='r', s=100)

    if seeds is not None:
        for x, r in seeds: 
            c = mpatches.Circle(x, r, alpha=0.6, fc="b", ec="b", lw=1)
            ax.add_patch(c)
    
    plt.show()

if __name__ == "__main__":
    main()
