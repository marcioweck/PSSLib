from __future__ import division

import argparse
import csv
import sys
sys.path.append("..")

import pycec2013
import numpy

from collections import namedtuple

from deap import tools

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

# from itertools import chain
import pdb
import time

Result = namedtuple("Result", ("index", "accuracy", "pr", "sr", "ave_fes", "std_fes", "ave_ngo", "std_ngo"))


def optimize(init, benchmark, accuracy):
    toolbox = init(benchmark)
    leftfes = toolbox.update_fes(benchmark.max_fes)
    ngoptima = 0
    max_ngoptima = benchmark.ngoptima

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)

    # pdb.set_trace()
    # databk = list()
    # i = 0
    while leftfes > 0 and ngoptima < max_ngoptima:
        pop = toolbox.generate()

        fitnesses = toolbox.map(benchmark.evaluate, pop)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = (fit, 0.0)

        stats.compile(pop)
        leftfes -= len(fitnesses)

        # databk.extend(pop)
        # databk.sort(key=lambda x:x.fitness)
        # databk = databk[:300]

        # try:
        #     toolbox.update(pop)
        # except:
        #     break
        toolbox.update(pop, benchmark.evaluate)
        leftfes -= toolbox.reborn()

        selection = toolbox.retrieve()
        fitnesses = [sol.fitness.values[0] for sol in selection]

        ngoptima = benchmark.count_goptima(selection, fitnesses, accuracy)
        # print ngoptima
        # if ngoptima > 4:
        #     print len(selection)
    info_fes = [0.0]*benchmark.ndim; info_fes[0] = benchmark.max_fes - leftfes
    info_ngo = [0.0]*benchmark.ndim; info_ngo[0] = ngoptima
    selection.append(info_fes)
    selection.append(info_ngo)

    lseeds = numpy.array([x for x in selection])

    # databk = numpy.array(databk)
    # numpy.savetxt('/home/weckwar/Documents/logsv2/databk.dat', databk, delimiter=',')
    # # np.savetxt('/home/weckwar/Documents/res{0}_{1}.out', lseeds, delimiter=',')
    # plt.scatter(lseeds[:,0], lseeds[:,1], c='b')
    # plt.show()
    # print "fes {0}".format(benchmark.max_fes - leftfes)
    # print "ngo {0}".format(ngoptima)
    # print "seeds {0}".format(len(lseeds))
    return leftfes, ngoptima, lseeds

def execute_run(init, benchmark, nbrun, accuracy):

    fes = [benchmark.max_fes] * nbrun
    total_optima = 0
    successful = 0
    n_goptima = [0] * nbrun

    # if accuracy == 1e-5:
    #     ngo_iter = numpy.zeros(benchmark.max_fes/(benchmark.ngoptima*30))

    for i in range(nbrun):
        leftfes, ngoptima, lseeds = optimize(init, benchmark, accuracy)
        fes[i] -= max(leftfes, 0)
        total_optima += ngoptima
        successful += ngoptima == benchmark.ngoptima
        n_goptima[i] += ngoptima

        numpy.savetxt("/home/weckwar/Documents/logsv3/f{0}_iter{1}_acc{2}.log".format(benchmark.index, i, accuracy), lseeds, delimiter=',')

    pr = total_optima / (nbrun * benchmark.ngoptima)
    sr = successful / nbrun
    ave_fes = numpy.mean(fes)
    std_fes = numpy.std(fes)
    ave_ngo = total_optima/nbrun
    std_ngo = numpy.std(n_goptima)
    return Result(benchmark.index, accuracy, pr, sr, ave_fes, std_fes, ave_ngo, std_ngo)

def execute_framework(init, indices, accuracies, nbrun, file_):
    writer = csv.writer(file_)
    writer.writerow(Result._fields)
    for index in indices:
        bench = pycec2013.Benchmark(index)
        for accuracy in accuracies:
            result = execute_run(init, bench, nbrun, accuracy)
            writer.writerow(result)

def main(init):
    args = parser.parse_args()
    indices = args.function
    accuracies = args.accuracy
    nbrun = args.N
    filename = args.file

    if filename != "":
        file_ = open(filename, 'w')
    else:
        file_ = sys.stdout

    execute_framework(init, indices, accuracies, nbrun, file_)


parser = argparse.ArgumentParser(description="CEC2013 Niching Competition Framework for DEAP.")
parser.add_argument('--function', type=int, nargs='+', default=pycec2013.functions.keys(),
                    help='Indices of benchmarks to run. Def.: [1, 20]')
parser.add_argument('--accuracy', type=float, nargs='+', default=(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                    help="Optimal fitness accuracy level. Def.: [1e-5, 1e-1]")
parser.add_argument('-N', type=int, default=50,
                    help="Number of runs / benchmark / accuracy. Def.: 50")
parser.add_argument("--file", type=str, default="",
                    help="Path and name of the results file. Def.: stdout")

