import array
import random
import sys
sys.path.append("..")

from math import log

from deap import base
from deap import cma
from deap import creator
from deap import tools
from deap.benchmarks.tools import bound
from framework import main

creator.create("Fitness", base.Fitness, weights=(1.0,))
creator.create("Individual", array.array, typecode='d', fitness=creator.Fitness)

def uniform(lbounds, ubounds):
    return [random.uniform(low, high) for low, high in zip(lbounds, ubounds)]

def checkBounds(lbounds, ubounds):
    def decorator(func):
        def wrappper(*args, **kargs):
            pop = func(*args, **kargs)
            for ind in pop:
                for i in range(len(ind)):
                    if ind[i] > ubounds[i]:
                        ind[i] = ubounds[i]
                    elif ind[i] < lbounds[i]:
                        ind[i] = lbounds[i]
            return pop
        return wrappper
    return decorator

def init(benchmark):
    lbounds = tuple(benchmark.get_lbounds())
    ubounds = tuple(benchmark.get_ubounds())
    min_ = min(lbounds)
    max_ = max(ubounds)

    centroid = uniform(lbounds, ubounds)
    sigma = 0.3 * (max_ - min_)
    lambda_ = int(4 + 3 * log(benchmark.ndim))
    strategy = cma.Strategy(centroid=centroid, sigma=sigma, lambda_=lambda_)

    toolbox = base.Toolbox()
    toolbox.register("generate", strategy.generate, creator.Individual)
    toolbox.register("update", strategy.update)
    toolbox.register("bounds", checkBounds(lbounds, ubounds))
    toolbox.decorate("generate", toolbox.bounds)
    return toolbox

if __name__ == "__main__":
	main(init)