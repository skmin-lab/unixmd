==========================
Build
==========================

UNI-xMD is python based program with a little C code for time-consuming electronic propagation part interfaced via Cython,
therefore compiliation is needed.


Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^
Python 3.5 (or newer)

Numpy

Cython https://cython.org


If you don't have numpy or Cython, you can install them using pip command.

.. code-block:: bash
   
   $ pip install --upgade numpy Cython

Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can compile electronic propagation routine by typing following command in root directory of the program which contains setup.py file.

.. code-block:: bash

   $ python3 setup.py build_ext -b ./src/mqc/el_prop

================================
Code Flow and Program Structure
================================
Simple intro..?

Work flow
^^^^^^^^^^^^^^^^^^^^^^^^^^
Diagram 

Program structure
^^^^^^^^^^^^^^^^^^^^^^^^^^
Inheritance,...

==========================
Quick Start
==========================
| This is quick start.
| program is controlled by running script.
| Goto directory 
| $ cd [UNIXMDHOME]/quick_start/

Define system
^^^^^^^^^^^^^^^^^^^^^^^^^^
mol = Molecule(~~~)

Define MD method
^^^^^^^^^^^^^^^^^^^^^^^^^^
import and make object
md = EH(~~~)

Run MD
^^^^^^^^^^^^^^^^^^^^^^^^^^
each md module has 'run' method which actually run~~~
md.run(~~~)

Check output
^^^^^^^^^^^^^^^^^^^^^^^^^^
files~~~~~,simple explanation

