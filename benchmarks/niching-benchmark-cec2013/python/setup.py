from distutils.core import setup
from Cython.Build import cythonize

import numpy

ext = cythonize("pycec2013.pyx")
ext[0].include_dirs.append(numpy.get_include())
setup(ext_modules=ext)