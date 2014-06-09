import array
import random
import sys
sys.path.append("..")

from deap import base
from deap import creator
from deap import tools
from framework import main

creator.create("Fitness", base.Fitness, weights=(1.0,))
creator.create("Individual", array.array, typecode='d', fitness=creator.Fitness)

def uniform(lbounds, ubounds):
    return (random.uniform(low, high) for low, high in zip(lbounds, ubounds))

def update(pop):
    pass

def init(benchmark):
    toolbox = base.Toolbox()
    lbounds = tuple(benchmark.get_lbounds())
    ubounds = tuple(benchmark.get_ubounds())

    toolbox.register("uniform", uniform, lbounds, ubounds)
    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.uniform)
    toolbox.register("generate", tools.initRepeat, list, toolbox.individual, 100)
    toolbox.register("update", update)
    return toolbox

if __name__ == "__main__":
    main(init)
