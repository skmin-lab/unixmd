from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

sourcefile = ["./unixmd/mqc/el_prop/el_propagator.pyx", "./unixmd/mqc/el_prop/rk4.c"]

extensions = [Extension("el_propagator", sources=sourcefile, include_dirs=[np.get_include()])]

setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules = extensions
)
