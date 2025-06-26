from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

# Selects the type of math libraries to be used; Available options: lapack, mkl
math_lib_type = "mkl"
#math_lib_type = "lapack"
# Directories including the math libraries
math_lib_dir = "${MKLROOT}/lib/intel64/"
#math_lib_dir = "/my_disk/my_name/lapack/"

sourcefile1 = ["./src/lib/mqc/el_propagator.pyx", "./src/lib/mqc/rk4.c", "./src/lib/mqc/exponential.c"]
sourcefile2 = ["./src/lib/mqc/el_propagator_xf.pyx", "./src/lib/mqc/rk4_xf.c","./src/lib/mqc/exponential_xf.c"]
sourcefile3 = ["./src/lib/mqc/el_propagator_ct.pyx", "./src/lib/mqc/rk4_ct.c"]
sourcefile4 = ["./src/lib/cioverlap/cioverlap.pyx", "./src/lib/cioverlap/tdnac.c"]

sourcefile1_qed = ["./src/lib/mqc_qed/el_propagator.pyx", "./src/lib/mqc_qed/rk4.c", "./src/lib/mqc_qed/exponential.c"]
sourcefile2_qed = ["./src/lib/mqc_qed/el_propagator_xf.pyx", "./src/lib/mqc_qed/rk4_xf.c", "./src/lib/mqc_qed/exponential_xf.c"]

# External libraries to be linked
libs = []
# Directories for linking external libraries
lib_dirs = []
# Extra flags for compilation
extra_flags = []

if (math_lib_type == "lapack"):
    libs += ["lapack", "refblas", "gfortran"]
    lib_dirs += [math_lib_dir]
    extra_flags += ["-D HAVE_LAPACK"]
elif (math_lib_type == "mkl"):
    libs += ["mkl_intel_lp64", "mkl_sequential", "mkl_core", ":libmkl_avx512.so.1"]
    lib_dirs += [math_lib_dir]
    extra_flags += ["-D HAVE_MKL"]
else:
    error_message = "Invalid type of math libraries is given!"
    error_vars = f"math_lib_type = {math_lib_type}"
    raise ValueError (f"( setup.py ) {error_message} ( {error_vars} )")

extensions = [
    # Electronic propagation in MQC dynamics
    Extension("libmqc", sources=sourcefile1,  include_dirs=[np.get_include()], \
        libraries=libs, library_dirs=lib_dirs),
    Extension("libmqcxf", sources=sourcefile2, include_dirs=[np.get_include()], \
        libraries=libs, library_dirs=lib_dirs),
    Extension("libctmqc", sources=sourcefile3, include_dirs=[np.get_include()]),
    Extension("libcioverlap", sources=sourcefile4, include_dirs=[np.get_include()], \
        libraries=libs, library_dirs=lib_dirs, extra_compile_args=extra_flags),
    # Electronic propagation in MQC_QED dynamics
    Extension("libmqc_qed", sources=sourcefile1_qed, include_dirs=[np.get_include()]),
    Extension("libmqcxf_qed", sources=sourcefile2_qed, include_dirs=[np.get_include()])
]

setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules = extensions
)
