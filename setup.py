from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

sourcefile1 = ["./src/mqc/el_prop/el_propagator.pyx", "./src/mqc/el_prop/rk4.c"]
sourcefile2 = ["./src/mqc/el_prop/el_propagator_xf.pyx", "./src/mqc/el_prop/rk4_xf.c"]
sourcefile3 = ["./src/mqc/el_prop/el_propagator_ct.pyx", "./src/mqc/el_prop/rk4_ct.c"]
sourcefile4 = ["./src/qm/cioverlap/cioverlap.pyx", "./src/qm/cioverlap/tdnac.c"]

extensions = [
    Extension("el_propagator", sources=sourcefile1, include_dirs=[np.get_include()]),
    Extension("el_propagator_xf", sources=sourcefile2, include_dirs=[np.get_include()]),
    Extension("el_propagator_ct", sources=sourcefile3, include_dirs=[np.get_include()]),
    Extension("cioverlap", sources=sourcefile4, include_dirs=[np.get_include()])
]

setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules = extensions
)
