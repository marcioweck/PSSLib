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
    return nearest_better_tree(D, fitnesses)


def dm_plugin(feval, fitnesses, population, edges, seeds, mask, nmin):
    global _idxs
    # Random search
    interior_pts = dm.inline_geom_sampler(5)

    sidxs = _idxs[mask]
    shuffled = random.sample(sidxs, len(sidxs))
    for l in shuffled:
        xy, d = edges[l]
        ind = population[xy[0]]
        dom = population[xy[1]]
        wps = dm.denorm(ind, dom, interior_pts)
        wpfit = np.array(map(feval, wps))

        if any(wpfit < fitnesses[xy[1]]):
            yield xy, d

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


def beta_plugin(nsamples, ndim):
    D2 = ndim**2
    S = nsamples

    return ((-4.69e-4 * D2 + 0.0263*ndim + 3.66/ndim - 0.457) * np.log10(S)
            + 7.51e-4 * D2 - 0.0421*ndim - 2.26/ndim + 1.83)

def raw_data_seeds_sel(feval, individuals, nmin, useDM=True, maskmode='NEA1'):
    """
    raw_data_seeds_sel selects seeds (good starting points) from
    a numpy.array matrix NxM with N individuals with M dimensions

    TODO explain parameters
    """
    global _idxs
    _idxs = np.arange(0, len(individuals))

    fitnesses = [feval(np.array(x)) for x in individuals]
    edges = nbtree(individuals, fitnesses)

    stree = adj_tree_plugin(edges)

    dists = np.array([d for _, d in edges])
    np.arange(0, len(dists))
    mu_dist = np.mean(dists)

    seeds = dict()
    beta = beta_plugin(len(individuals), len(individuals[0]))
    # Select the main seeds
    mask = ((dists >= 2.*mu_dist) + (dists == 0))
    if maskmode == 'NEA2':
        for k, adj in stree.iteritems():
            # print "log:"
            # print 'k {0}'.format(k)
            # print 'nadjs {0}'.format(len(adj))
            # print 'calc {0}'.format(edges[k][1]/np.median(adj.values()))
            # print 'beta {0}'.format(beta)
            # print
            # print
            if len(adj) >= 3 and edges[k][1]/np.median(adj.values()) > beta:
                # pdb.set_trace()
                mask[k] |= 1

    sidxs = _idxs[mask]
    for sid in sidxs:
        try:
            adjnodes = stree[sid]
            c = min(1.8*mu_dist, max(adjnodes.values())) if len(adjnodes) > 0 else 1.8*mu_dist
            seeds[sid] = (individuals[sid], c, (fitnesses[sid], dists[sid]))
        except:
            seeds[sid] = (individuals[sid], mu_dist, (fitnesses[sid], dists[sid]))

    if useDM and len(seeds) < nmin:
        for (i,j), d in dm_plugin(feval, fitnesses, individuals, edges, seeds, -mask, nmin):
            seeds[i] = (individuals[i], 0.6*d, (fitnesses[i], d))
            if j in seeds:
                seeds[j] = (individuals[j], 0.6*d, (fitnesses[j], dists[j]))

    return seeds.values(), edges