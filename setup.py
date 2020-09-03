from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

sourcefile1 = ["./src/mqc/el_prop/el_propagator.pyx", "./src/mqc/el_prop/rk4.c"]
sourcefile2 = ["./src/qm/cioverlap/cioverlap.pyx", "./src/qm/cioverlap/tdnac.c"]

extensions = [
    Extension("el_propagator", sources=sourcefile1, include_dirs=[np.get_include()]),
    Extension("cioverlap", sources=sourcefile2, include_dirs=[np.get_include()])
]

setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules = extensions
)
