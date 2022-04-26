==========================
Installation
==========================

PyUNIxMD is a Python based program with a little C code for time-consuming
electronic propagation parts interfaced via Cython, therefore compilation is needed.

Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Python 3.6 (or later)

-  Numpy

-  Cython https://cython.org

If you don't have Numpy or Cython, you can install them using :code:`pip` command.

.. code-block:: bash

   $ pip install --upgrade numpy Cython

Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can obtain PyUNIxMD code by putting following command.

.. code-block:: bash

   $ git clone https://github.com/skmin-lab/unixmd.git

You can compile the C routines by typing the following
command in the top-level directory of the program which contains setup.py file.

.. code-block:: bash

   $ cd unixmd/
   $ python3 setup.py build_ext -b ./src/build

You will need to add the source directory(:code:`$PYUNIXMDHOME/src`) to your Python path, where :code:`$PYUNIXMDHOME` is an enviroment variable for the top-level directory.

+----------------------------+--------------------------------------------------+----------------+
| Option                     | Work                                             | Default        |
+============================+==================================================+================+
| **lib_type**               | library type                                     | *None*         |
| *(string)*                 |                                                  |                |
+----------------------------+--------------------------------------------------+----------------+
| **complier_name**          | complier name                                    | *GCC*          |
| *(string)*                 |                                                  |                |
+----------------------------+--------------------------------------------------+----------------+

If you want to use library for exponential propagator, you need to add :code:`lib_type` option.

.. code-block:: bash

   $ python3 setup.py build_ext -b ./src/build --lib_type=lapack

You can select lapack or mkl library.

In our program, lapack 3.10 version is provided.

If you want to change compiler such as Intel C Compiler, you need to add :code:`complier_name` option and :code:`LDSHARED="icc -shared" CC=icc` option before python3.

.. code-block:: bash

   $ LDSHARED="icc -shared" CC=icc python3 setup.py build_ext -b ./src/build --complier_name=icc