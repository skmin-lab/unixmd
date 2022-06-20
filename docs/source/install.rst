==============================
Installation
==============================

PyUNIxMD is a Python based program with a little C code for time-consuming parts
(e.g., electronic propagation or TDNACs) interfaced via Cython, therefore compilation is needed.

Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Python 3.6 (or later)

-  Numpy

-  Cython https://cython.org

If you don't have Numpy or Cython, you can install them using :code:`pip` command.

.. code-block:: bash

   $ pip install --upgrade numpy Cython

Optional external libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  LAPACK/BLAS libraries

Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can obtain PyUNIxMD code by putting following command.

.. code-block:: bash

   $ git clone https://github.com/skmin-lab/unixmd.git

You can compile the C routines by typing the following
command in the top-level directory of the program which contains setup.py file.

.. code-block:: bash

   $ cd unixmd/
   $ python3 setup.py build_ext -b ./src/build

If you want to link the external math libraries for TDNACs, you may modify **math_lib_type** and **math_lib_dir**
defined in setup.py file as follows:

.. code-block:: python
   :linenos:

   from distutils.core import setup
   from distutils.extension import Extension
   from Cython.Distutils import build_ext

   import numpy as np

   # Selects which type of math libraries to be used; Available options: None, lapack
   math_lib_type = "lapack"

   # Directories including the math libraries
   math_lib_dir = "/my_disk/my_name/lapack/"

   ...

After successful compilation, you will need to add the source directory (:code:`$PYUNIXMDHOME/src`) to your Python path,
where :code:`$PYUNIXMDHOME` is an enviroment variable for the top-level directory.
