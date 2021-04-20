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
