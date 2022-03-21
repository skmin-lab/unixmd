from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
import os

sourcefile1 = ["./src/mqc/el_prop/el_propagator.pyx", "./src/mqc/el_prop/rk4.c", "./src/mqc/el_prop/temp_expon.c"]
sourcefile2 = ["./src/mqc/el_prop/el_propagator_xf.pyx", "./src/mqc/el_prop/rk4_xf.c"]
sourcefile3 = ["./src/mqc/el_prop/el_propagator_ct.pyx", "./src/mqc/el_prop/rk4_ct.c"]
sourcefile4 = ["./src/qm/cioverlap/cioverlap.pyx", "./src/qm/cioverlap/tdnac.c"]

lib_type = "lapack"

if lib_type == "lapack":
    extensions = [
        Extension("el_propagator", sources=sourcefile1, include_dirs=[np.get_include()],libraries=["lapack","gfortran","m"],
            library_dirs=[".src/lib/"]),
        Extension("el_propagator_xf", sources=sourcefile2, include_dirs=[np.get_include()]),
        Extension("el_propagator_ct", sources=sourcefile3, include_dirs=[np.get_include()]),
        Extension("cioverlap", sources=sourcefile4, include_dirs=[np.get_include()])
    ]

elif lib_type == "mkl":
    extensions = [
        Extension("el_propagator", sources=sourcefile1, include_dirs=[np.get_include()],
            libraries=["mkl_intel_lp64","mkl_sequential", "mkl_core",":libmkl_avx512.so.1"] , 
            library_dirs=["${MKLROOT}/lib/intel64"]),
        Extension("el_propagator_xf", sources=sourcefile2, include_dirs=[np.get_include()]),
        Extension("el_propagator_ct", sources=sourcefile3, include_dirs=[np.get_include()]),
        Extension("cioverlap", sources=sourcefile4, include_dirs=[np.get_include()])
    ]
else:
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
