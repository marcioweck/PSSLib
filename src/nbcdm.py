import random
import copy

import numpy as np
from scipy.spatial.distance import euclidean, pdist, squareform

import dm

import pdb

_idxs = None

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


def nbtree(population, fitnesses):
    D = squareform(pdist(population))
    # TODO falta remover os colegas

    return nearest_better_tree(D, fitnesses)


def dm_plugin(population, edges, seeds, mask):
    # Random search
    interior_pts = dm.inline_geom_sampler(5)

    sidxs = _idxs[-mask]
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
            yield xy, d
            # ind.id = xy[0]
            # seeds[xy[0]] = (ind, d*0.6)
            # if xy[1] in seeds:
            #    seeds[xy[1]] = (dom, d*0.6)
            # else
            #     pass 

        if len(seeds) >= nmin:
                break


def adj_tree_plugin(edges):
    stree = dict()
    for (a, b), d in edges:
        if b in stree:
            stree[b][a] = d
        else:
            stree[b] = {a: d}

    return stree


def raw_data_seeds_sel(feval, individuals, nmin, useDM=True):
    """
    raw_data_seeds_sel selects seeds (good starting points) from
    a numpy.array matrix NxM with N individuals with M dimensions

    TODO explain parameters
    """

    fitnesses = [feval(x) for x in individuals]
    edges = nbtree(individuals, fitnesses)

    stree = adj_tree_plugin(edges)

    dists = [d for _, d in edges]
    np.arange(0, len(dists))
    mu_dist = np.mean(dists)

    seeds = dict()
    # Select the main seeds
    mask = ((dists >= 2.*mu_dist) + (dists == 0))

    sidxs = _idxs[mask]
    for i in sidxs:
        try:
            adjnodes = stree[sid]
            c = min(1.5*mu_dist, max(adjnodes, key=adjnodes.get))
            seeds[sid] = (individuals[sid], c, (fitnesses[sid], dists[sid]))
        except:
            seeds[sid] = (individuals[sid], mu_dist, (fitnesses[sid], dists[sid]))

    if useDM and len(seeds) < nmin:
        for (i,j), d in dm_plugin(individuals, edges, seeds, mask):
            seeds[i] = (individuals[i], 0.6*d, (fitnesses[i], d))
            if j in seeds:
                seeds[j] = (individuals[j], 0.6*d, (fitnesses[j], dists[j]))
        
    return seeds.values(), edges