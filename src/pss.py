# TODO license
# TODO explain

import array
import random
import copy

from itertools import izip, combinations, chain
from itertools import permutations as permut

import numpy as np

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
import nbcdm
from archive import Archive
from utils import solowPolaskyQC, averageFitnessQC, ensure_bounds

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
    else:
        if dim == 1:
            P = np.array([np.array((x,)) for x in P])

    return P


"""
PSS framework
"""

import pdb


def remove_overlaps(feval, individuals, hives):

    if len(individuals) == 0: return []

    idx = range(0, len(hives))
    idx.sort(key=lambda i: individuals[i].fitness)

    nm = np.sqrt(len(individuals[0]))

    covs = [h.sigma*nm for h in hives]

    D = squareform(pdist(individuals))

    uniques = set(idx[:])

    for i, j in permut(idx, 2):
        if covs[i] > D[i,j]:
            if individuals[i].fitness > individuals[j].fitness:
                wps = dm.denorm(individuals[i], individuals[j], [0.5])
                wpfit = np.array(map(feval, wps))
                if all(wpfit >= individuals[j].fitness.wvalues[0]):
                    D[j,...] = np.inf
                    uniques.discard(j)

    return (i for i in uniques)


def generate(ind_init, hives):
    swarms = [h.generate(ind_init) for h in hives]
    return swarms


def updateHive(hive, swarm):
    hive.update(swarm)
    return hive.shouldContinue()


def generateHive(hclass, (xstart, sigma)):
    return hclass(copy.deepcopy(xstart), sigma)


import sys
sys.path.append("../benchmarks/niching-benchmark-cec2013/python/")
from framework import pycec2013

import cma


creator.create("Fitness", base.Fitness, weights=(1.0, 1.0))
creator.create("OneFitness", base.Fitness, weights=(1.0,))
creator.create("Individual", array.array, typecode='d', fitness=creator.Fitness)
creator.create("Centroid", array.array, typecode='d', fitness=creator.OneFitness)
creator.create("Hive", cma.Strategy)


def main():

    benchmark = pycec2013.Benchmark(17)

    lbounds = tuple(benchmark.get_lbounds())
    ubounds = tuple(benchmark.get_ubounds())
    min_ = min(lbounds)
    max_ = max(ubounds)

    toolbox = base.Toolbox()
    toolbox.register("generate", generate, creator.Individual)
    toolbox.register("update", map, updateHive)
    toolbox.register("feval", map, benchmark.evaluate)
    # toolbox.register("fdist", nearest_better_tree)
    toolbox.register("hive", generateHive, cma.Strategy)
    toolbox.register("bounds", ensure_bounds(lbounds, ubounds))
    toolbox.decorate("generate", toolbox.bounds)

    dim = benchmark.ndim
    nmin = benchmark.ngoptima
    leftfes = benchmark.max_fes
    ngoptima = 0
    max_ngoptima = benchmark.ngoptima

    def similarity_func(a, b):
        if np.isnan(np.sum(a)) or np.isnan(np.sum(b)):
            pdb.set_trace()

        d = euclidean(a, b)
        return d < 0.06

    hof = Archive(max_ngoptima, similarity_func)

    distribs = [stats.uniform for i in range(dim)]

    samples = sampler(distribs, (min_, max_), dim, 20*dim)

    #samples = np.loadtxt("/home/weckwar/inputs.txt", delimiter=',', ndmin=2)

    seeds, _ = nbcdm.raw_data_seeds_sel(benchmark.evaluate, samples, 20, useDM=True, maskmode='NEA1')
    # xpoints = np.array([x for x,c,_ in seeds])
    # np.savetxt("/home/weckwar/inputs.txt", xpoints, delimiter=',')
    #plotseeds(benchmark.evaluate, min_, max_, dim, samples=xpoints)
    #return

    hives = list()
    population = list()
    norm = float(np.sqrt(dim))
    for (xstart, c, (f1, f2)) in seeds:
        ind = creator.Individual(xstart)
        ind.fitness.values = (f1, f2)
        population.append(ind)
        hives.append(toolbox.hive((ind, c/norm)))

    verbose = True

    logbook = tools.Logbook()
    logbook.header = "gen", "nswarm", "ngoptima", "muerror", "dispersion"

    generation = 0
    logbook.record(gen=generation, nswarm=len(hives), ngoptima=ngoptima,
                   muerror=averageFitnessQC(population),
                   dispersion=solowPolaskyQC(population, 1.0/dim))

    while leftfes > 0  and ngoptima < max_ngoptima:

        swarms = toolbox.generate(hives)

        blob = list(chain(*swarms))
        D = squareform(pdist(blob))
        fitnesses = toolbox.feval(blob)

        nelem = len(swarms[0])
        for i, swarm in enumerate(swarms):
            k = i*nelem
            nbidx = np.arange(k, k+nelem)
            for j, ind in enumerate(swarm):
                D[k+j,nbidx] = np.inf
                sortedline = np.argsort(D[k+j,:])
                bestidx = next((l for l in sortedline
                                if fitnesses[l] > fitnesses[k+j]), -1)

                ind.fitness.values = (fitnesses[k+j], D[k+j, bestidx])

        checks = toolbox.update(hives, swarms)

        nextgen = [hives[i] for i, ok in enumerate(checks) if ok]

        xstarts = [creator.Centroid(x.centroid) for x in nextgen]
        cfit = toolbox.feval(xstarts)
        for x, fit in izip(xstarts, cfit):
            x.fitness.values = (fit,)

        uniques = list(remove_overlaps(benchmark.evaluate, xstarts, nextgen))
        hives = [nextgen[i] for i in uniques]
        xstarts = [xstarts[i]  for i in uniques]

        hof.update(xstarts)
        hfit = [x.fitness.values[0] for x in hof]

        ngoptima = benchmark.count_goptima(hof, hfit, 1e-5)

        if len(hives) < 2:
            samples = sampler(distribs, (min_, max_), dim, 2.*dim)

            seeds, _ = nbcdm.raw_data_seeds_sel(benchmark.evaluate, samples, 10.)

            for (xstart, c, (f1, f2)) in seeds:
                ind = creator.Individual(xstart)
                ind.fitness.values = (f1, f2)
                hives.append(toolbox.hive((ind, 0.5*c/norm)))

            leftfes -= len(samples)

        leftfes -= len(swarms)*nelem + len(xstarts)

        generation += 1
        logbook.record(gen=generation, nswarm=len(hives), ngoptima=ngoptima,
                   muerror=0,#averageFitnessQC(xstarts),
                   dispersion=0)#solowPolaskyQC(xstarts, 1.0/dim))
        print logbook.stream

    print "Used FEs: {0}".format(benchmark.max_fes - leftfes)
    print ngoptima
    for ind in hof:
        print "x: {0} -> {1}".format(ind, ind.fitness.values[0])
    plotseeds(benchmark.evaluate, min_, max_, dim, samples=hof)


def plotseeds(feval, min_, max_, dim, samples=None, edges=None, seeds=None):
    fig, ax = plt.subplots()

    if dim == 2:
        X = np.arange(min_, max_, 0.05)
        Y = np.arange(min_, max_, 0.05)
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

    if dim == 1:
        X = np.arange(min_, max_, 0.01)
        Y = [feval(x) for x in X]
        plt.plot(X,Y,'r')

        if samples is not None:
            F = [s.fitness.values[0] for s in samples]
            plt.scatter(samples, F, c='k', s = 10)

        if seeds is not None:
            for x, r, _ in seeds:
                c = mpatches.Circle((x, feval(x)), r, alpha=0.6, fc="b", ec="b", lw=1)
                ax.add_patch(c)

    plt.show()

if __name__ == "__main__":
    main()
