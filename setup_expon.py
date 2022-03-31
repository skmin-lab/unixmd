from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np
import os

lib_type = os.environ.get('lib_type')
complier = os.environ.get('CC')

if lib_type == None:
    print("Library is not chosen. If you want to use exponetial propagator, you have to select library.")

if lib_type == "lapack":
    sourcefile1 = ["./src/mqc/el_prop/el_propagator.pyx", "./src/mqc/el_prop/rk4.c", "./src/mqc/el_prop/exponential.c"]
    if complier == "icc":
        lib_name = ["liblapack_ifort_zheev","ifcore","m"]
    else:
        lib_name = ["liblapack_gcc_zheev","gfortran","m"]
    lib_dir = ["./lib/"]
elif lib_type == "mkl":
    sourcefile1 = ["./src/mqc/el_prop/el_propagator.pyx", "./src/mqc/el_prop/rk4.c", "./src/mqc/el_prop/exponential.c"]
    lib_name = ["mkl_intel_lp64","mkl_sequential", "mkl_core",":libmkl_avx512.so.1"]
    lib_dir = ["${MKLROOT}/lib/intel64/"]
else:
    sourcefile1 = ["./src/mqc/el_prop/el_propagator.pyx", "./src/mqc/el_prop/rk4.c"]
    lib_name = []
    lib_dir = []

sourcefile2 = ["./src/mqc/el_prop/el_propagator_xf.pyx", "./src/mqc/el_prop/rk4_xf.c"]
sourcefile3 = ["./src/mqc/el_prop/el_propagator_ct.pyx", "./src/mqc/el_prop/rk4_ct.c"]
sourcefile4 = ["./src/qm/cioverlap/cioverlap.pyx", "./src/qm/cioverlap/tdnac.c"]

extensions = [
    Extension("el_propagator", sources=sourcefile1, include_dirs=[np.get_include()],libraries=lib_name,library_dirs=lib_dir),
    Extension("el_propagator_xf", sources=sourcefile2, include_dirs=[np.get_include()]),
    Extension("el_propagator_ct", sources=sourcefile3, include_dirs=[np.get_include()]),
    Extension("cioverlap", sources=sourcefile4, include_dirs=[np.get_include()])
]

setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules = extensions
)
