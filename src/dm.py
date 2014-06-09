from itertools import combinations, imap
from operator import attrgetter

import numpy as np

import pdb


def inline_geom_sampler(n):
    """
    TODO docstring here
    """
    point = (0., 1.)
    stack = [point]
    rstack = list()
    while(len(rstack) < n):
        p = stack.pop(0)
        m = sum(p)/2.
        stack.append((p[0], m))
        stack.append((m, p[1]))
        rstack.append(m)

    return rstack


def denorm(x, y, gradations):
    """
    TODO docstring here
    """
    return [np.array(x) + (np.array(y) - x) * grad for grad in gradations]


def create_tests(individuals, ntests=1):

    tests = []
    extendt = tests.extend
    for a, b in combinations(individuals, 2):
        extendt(denorm(a, b, inline_geom_sampler(ntests)))

    return tests


def check_tests(individuals, tests, ntests=1):
    """
    Individus e testes precisam ser instancias da classe Individual (deap)
    Retornar quem esta dominado
    """
    fails = np.zeros(len(individuals))
    # dominated = set()
    failed_tests = []

    for i, (a, b) in enumerate(combinations(individuals, 2)):
        minimum = min(a, b, key=attrgetter("fitness"))
        for t in xrange(ntests):
            if tests[(i*ntests)+t].fitness < minimum.fitness:
                break  # has valley, for now everything is ok
            failed_tests.append((i*ntests)*t)
        # if failed in all tests, i.e., no valley was found, we consider
        #  this solution dominated
        if len(failed_tests) == ntests:
            fails[minimum.pid] += 1

    # dominated = [ind for i, ind in enumerate(individuals) if fails[i] > 0]
    # basins = [ind for i, ind in enumerate(individuals) if fails[i] == 0]

    # return basins, dominated, failed_tests
    return fails, failed_tests
