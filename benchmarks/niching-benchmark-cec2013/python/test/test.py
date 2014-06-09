import pycec2013

for i in range(1, 21):
	bench = pycec2013.Benchmark(i)
	print pycec2013.functions[i] 
	print "\t", "ndim:", bench.ndim
	print "\t", "rho:", bench.rho
	print "\t", "# global optima:", bench.ngoptima
	print "\t", "lower bounsd:", list(bench.get_lbounds())
	print "\t", "upper bounds:", list(bench.get_ubounds())
	print "\t", "max fes:", bench.max_fes
	print "\t", "global optima fitness:", bench.goptima_fitness
