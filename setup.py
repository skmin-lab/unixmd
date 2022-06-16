from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import sys
import numpy as np

# Selects which math libraries to be used; Available options: None, lapack
math_lib = None
# Directories including the math libraries
math_lib_dir = None

sourcefile1 = ["./src/mqc/el_prop/el_propagator.pyx", "./src/mqc/el_prop/rk4.c"]
sourcefile2 = ["./src/mqc/el_prop/el_propagator_xf.pyx", "./src/mqc/el_prop/rk4_xf.c"]
sourcefile3 = ["./src/mqc/el_prop/el_propagator_ct.pyx", "./src/mqc/el_prop/rk4_ct.c"]
sourcefile4 = ["./src/qm/cioverlap/cioverlap.pyx", "./src/qm/cioverlap/tdnac.c"]

# External libraries to be linked
libs = []
# Directories for linking external libraries
lib_dirs = []
# Extra flags for compilation
extra_flags = []

if (math_lib == None):
    print ("No math libraries are used!", flush=True)
elif (math_lib == "lapack"):
    libs += ["lapack", "refblas"]
    lib_dirs += [math_lib_dir]
    extra_flags += ["-D HAVE_LAPACK"]
else:
    print (f"Wrong math_lib:{math_lib} is given!", flush=True)
    sys.exit(1)

extensions = [
    Extension("el_propagator", sources=sourcefile1, include_dirs=[np.get_include()]),
    Extension("el_propagator_xf", sources=sourcefile2, include_dirs=[np.get_include()]),
    Extension("el_propagator_ct", sources=sourcefile3, include_dirs=[np.get_include()]),
    Extension("cioverlap", sources=sourcefile4, include_dirs=[np.get_include()], \
        libraries=libs, library_dirs=lib_dirs, extra_compile_args=extra_flags)
]

setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules = extensions
)
