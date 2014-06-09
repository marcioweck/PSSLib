# distutils: language = c++
# distutils: sources = ../c++/cec2013.cpp ../c++/cfunction.cpp
# distutils: include_dirs = ../c++

cimport numpy as np
import numpy
from libcpp.vector cimport vector

cdef extern from "cec2013.h":
    cdef cppclass CEC2013:
        CEC2013(int) except +
        long double evaluate(double*)
        
        double get_rho()
        double get_lbound(int)
        double get_ubound(int)
        double get_fitness_goptima()

        int count_goptima(vector[vector[double]],
                          vector[double],
                          double)

        int get_maxfes()
        int get_dimension()
        int get_no_goptima()

functions = {1 : "Five-Uneven-Peak Trap",
             2 : "Equal Maxima",
             3 : "Uneven Decreasing Maxima",
             4 : "Himmelblau",
             5 : "Six-Hump Camel Back",
             6 : "Shubert (2D)",
             8 : "Shubert (3D)",
             7 : "Vincent (2D)",
             9 : "Vincent (3D)",
             10 : "Modified Rastrigin - All Global Optima",
             11 : "Composition Function 1",
             12 : "Composition Function 2",
             13 : "Composition Function 3 (2D)",
             14 : "Composition Function 3 (3D)",
             16 : "Composition Function 3 (5D)",
             18 : "Composition Function 3 (10D)",
             15 : "Composition Function 4 (3D)",
             17 : "Composition Function 4 (5D)",
             19 : "Composition Function 4 (10D)",
             20 : "Composition Function 4 (20D)"
            }

functions_map = dict(zip(functions.values(), functions.keys()))

cdef class Benchmark:
    """Class that encapsulates the set of 20 multimodal benchmarks 
    proposed for the 2013 CEC Niching competition. The initialization
    takes an index between 1 and 20 corresponding to a specific function 
    and in some cases a specific number of dimensions. The name of the
    functions can be retrieved via the dictionnary *functions* at the root
    of the pycec2013 module, and the function index can be found based on
    the function name using *functions_map*.

    This class is a Cython wrapper around the C++ class CEC2013.
    """
    cdef CEC2013 *thisptr
    cdef public int index
    cdef public int max_fes
    cdef public int ndim
    cdef public int ngoptima
    cdef public double rho
    cdef public double goptima_fitness

    def __cinit__(self, index):
        self.index = index
        self.thisptr = new CEC2013(index)
        self.ndim = self.thisptr.get_dimension()
        self.rho = self.thisptr.get_rho()
        self.ngoptima = self.thisptr.get_no_goptima()
        self.max_fes = self.thisptr.get_maxfes()
        self.goptima_fitness = self.thisptr.get_fitness_goptima()
    
    def __dealloc__(self):
        del self.thisptr

    def __reduce__(self):
        return Benchmark, (self.index,)

    def evaluate(self, ind):
        """Evaluate the fitness of the individual *ind*. The individual
        must respect the buffer interface, i.e.: either an array.array
        or a numpy.ndarray."""
        cdef np.ndarray nparray = numpy.frombuffer(ind)
        return self.thisptr.evaluate(<double*>np.PyArray_DATA(nparray))

    def get_lbounds(self):
        """Return an iterator on the lower bounds for every dimension"""
        return (self.thisptr.get_lbound(dim) for dim in range(self.ndim))

    def get_ubounds(self):
        """Return an iterator on the upper bounds for every dimension."""
        return (self.thisptr.get_ubound(dim) for dim in range(self.ndim))

    def count_goptima(self, pop, fits, double epsilon):
        """Return the number of unique global optima in a population *pop*,
        with the fitnesses *fits*, given a level of accuracy *epsilon*.
        """
        cdef vector[vector[double]] pop_vec
        cdef vector[double] fits_vec
        cdef int i, j
        cdef int ndim = len(pop[0])

        pop_vec.resize(len(pop))
        fits_vec.resize(len(pop))
        for i, (ind, fit) in enumerate(zip(pop, fits)):
            pop_vec[i].resize(ndim)
            fits_vec[i] = fit
            for j, value in enumerate(ind):
                pop_vec[i][j] = value

        return self.thisptr.count_goptima(pop_vec, fits, epsilon)                

