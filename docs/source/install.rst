==========================
Installation
==========================

PyUNIxMD is a Python based program with a little C code for time-consuming
electronic propagation parts interfaced via Cython, therefore compilation is needed.

Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^

Python 3.6 (or newer)

Numpy

Cython https://cython.org

If you don't have Numpy or Cython, you can install them using pip command.

.. code-block:: bash

   $ pip install --upgrade numpy Cython

Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can compile the C routines by typing the following
command in the root directory of the program which contains setup.py file.

.. code-block:: bash

   $ python3 setup.py build_ext -b ./src/build

You need to add the source directory(:code:`$PYUNIXMDHOME/src`) to your Python path, where :code:`$PYUNIXMDHOME` is the top-level directory where this file belongs.
