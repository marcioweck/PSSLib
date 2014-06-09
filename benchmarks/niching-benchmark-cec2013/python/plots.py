#!python
"""Script that plots the benchmark functions defined by CEC2013
multimodal niching competition.
"""
import numpy
import pycec2013

from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

SAMPLES = 750

if __name__ == "__main__":

	for index, name in pycec2013.functions.items():
		bench = pycec2013.Benchmark(index)
		if bench.ndim > 2:
			continue

		fig = plt.figure()

		lbounds = tuple(bench.get_lbounds())
		ubounds = tuple(bench.get_ubounds())
		x = numpy.linspace(lbounds[0], ubounds[0], num=SAMPLES)
		if bench.ndim == 1:
			x = x.reshape(SAMPLES, 1)
			y = map(bench.evaluate, x)
			plt.plot(x, y)
		elif bench.ndim == 2:
			y = numpy.linspace(lbounds[1], ubounds[1], num=SAMPLES)
			x, y = numpy.meshgrid(x, y)
			pts = numpy.hstack((x.reshape(x.size, 1), y.reshape(y.size, 1)))
			z = numpy.array(map(bench.evaluate, pts), copy=False)
			z = z.reshape(x.shape)
			ax = fig.gca(projection='3d')
			ax.plot_surface(x, y, z, cmap=cm.jet)
			ax.contour(x, y, z, zdir='z', cmap=cm.jet, offset=numpy.min(z))

		plt.show()
