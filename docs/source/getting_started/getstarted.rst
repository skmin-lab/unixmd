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

Here, we provide executable running script file, which contains:

.. code-block:: python

   from molecule import Molecule
   import bo
   import mqc
   from thermostat import *
   from misc import data

   geom = """
   NUMBER_OF_ATOMS
   TITLE
   SYMBOL  COORDINATES  VELOCITIES
   """

   mol = Molecule(geometry=geom, nstates=NSTATES)

   qm = bo.PROG_NAME.QM_METHOD(ARGUMENTS)

   md = mqc.MD_METHOD(ARGUMETNS)

   bathT = THERMOSTAT(ARGUMENTS)

   md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir=INPUT_DIR)

If you execute this script, you can get output file.
UNI-xMD provides some output files such as **MDENERGY**, **MOVIE.xyz**, **FINAL.xyz** when all types of ``md``, which can be possible in UNI-xMD, are performed. 
In the case of ``Eh``, ``SH`` and ``SHXF``, **BOCOH**, **BOPOP** and **NACME** are addtionally generated. Escpecially, ``SH`` and ``SHXF`` generate more output file **SHPROB** and **SHSTATE**.

- MDENERGY  : energy which contains kinetic energy, potential energy of each adiabatic state and total energy 

.. code-block:: bash

    Here is data

- MOVIE.xyz : geometries at each step 

.. code-block:: bash

    Here is data

- FINAL.xyz : geometry at final step

.. code-block:: bash

    Here is data

- BOCOH     : off-diagonal term of adiabatic density matrix

.. code-block:: bash

    Here is data

- BOPOP     : adiabatic population

.. code-block:: bash

    Here is data

- NACME     : non-adiabatic coupling matrix

.. code-block:: bash

    Here is data

- SHPROB    : hopping probability between all of adiabatic states

.. code-block:: bash

    Here is data

- SHSTATE   : running state 
.. code-block:: bash

    Here is data

