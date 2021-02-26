==========================
Build
==========================

UNI-xMD is python based program with a little C code for time-consuming
electronic propagation part interfaced via Cython, therefore compilation is needed.

Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^

Python 3.6 (or newer)

Numpy

Cython https://cython.org

If you don't have numpy or Cython, you can install them using pip command.

.. code-block:: bash

   $ pip install --upgrade numpy Cython

Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can compile electronic propagation routine by typing following
command in root directory of the program which contains setup.py file.

.. code-block:: bash

   $ python3 setup.py build_ext -b ./src/build


