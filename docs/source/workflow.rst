==========================
Workflow
==========================
Here, we explain how to run MD calculations with PyUNIxMD.

You will make a running script for the MD calculation you want to perform. In your running script, you will create PyUNIxMD objects successively.
A typical template of the running script is the following:

.. code-block:: python
   :linenos:

   from molecule import Molecule
   import qm, mqc
   from thermostat import *
   from misc import data

   geom = """
   <number of atoms>
   <comment>
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   ...
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   """

   mol = Molecule(geometry=geom, nstates=NSTATES)

   qm = qm.QM_PROG.QM_METHOD(ARGUMENTS)

   md = mqc.MDTYPE(ARGUMENTS)

   bathT = THERMOSTAT(ARGUMENTS)

   md.run(molecule=mol, theory=qm, thermostat=bathT)

**Line 1-4** import the PyUNIxMD packages for the below jobs.

**Line 6-12** set a target system you are interested in.
You need to prepare a string as an argument to specify initial geometry and velocities in extended XYZ format.
NSTATES means the number of adiabatic states considered in the dynamics calculations.
See :ref:`PyUNIxMD Objects <Objects Molecule>` for the list of parameters.

.. note:: The ``mol`` object must be created first because it will be used for making other objects.

**Line 14** determines an electronic structure calculation program and its method to obtain QM information such as energies, forces, and nonadiabatic coupling vectors.
QM_PROG is the directory name where the QM interface package is. QM_METHOD is the name of the Python class specifying one of QM methods provided with that interface package. See :ref:`PyUNIxMD Objects <Objects QM_calculator>` for the list.

**Line 16** determines a dynamics method you want to use. MDTYPE is the name of Python class specifying MQC methods. See :ref:`PyUNIxMD Objects <Objects MQC>` for the list.

**Line 18** sets a thermostat. THERMOSTAT is one of the options. See :ref:`PyUNIxMD Objects <Objects Thermostat>` for the list. 

**Line 20** runs the dynamics calculation. 

Finally, you will execute your running script.

.. code-block:: bash

   $ python3 running_script.py

Running MD calculations with PyUNIxMD, you will obtains the following output file according to the MD method.

+-----------+------+----+----+------+
|           | BOMD | Eh | SH | SHXF |
+===========+======+====+====+======+
| MDENERGY  | o    | o  | o  | o    |
+-----------+------+----+----+------+
| MOVIE.xyz | o    | o  | o  | o    |
+-----------+------+----+----+------+
| FINAL.xyz | o    | o  | o  | o    |
+-----------+------+----+----+------+
| BOCOH *   | x    | o  | o  | o    |
+-----------+------+----+----+------+
| BOPOP *   | x    | o  | o  | o    |
+-----------+------+----+----+------+
| NACME     | x    | o  | o  | o    |
+-----------+------+----+----+------+
| SHPROB    | x    | x  | o  | o    |
+-----------+------+----+----+------+
| SHSTATE   | x    | x  | o  | o    |
+-----------+------+----+----+------+
| DOTPOPD   | x    | x  | x  | o    |
+-----------+------+----+----+------+

.. note:: If you put **propagation** = *"density"* when setting the MD method, PyUNIxMD provides 'BOCOH' and 'BOPOP'.
   However, if you put **propagation** = *"coefficient"* when setting the MD method, PyUNIxMD provides 'BOCOEF' rather than 'BOCOH' and 'BOPOP'.

- MDENERGY

This file shows MD energies and energies of adiabatic states

.. code-block:: bash

   <MD step>   <kinetic energy>   <potential energy>   <total MD energy>   <adiabatic energy 1>   <adiabatic energy 2> ... <adiabatic energy last>
   <MD step>   <kinetic energy>   <potential energy>   <total MD energy>   <adiabatic energy 1>   <adiabatic energy 2> ... <adiabatic energy last>
   ...

- MOVIE.xyz

This file contains positions and velocities at each time step.
For the ease of visualization, those information are written successively in the extended XYZ format.

.. code-block:: bash

   <number of atoms>
   Step:     0
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   ...
   <number of atoms>
   Step:     1
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   ...

- FINAL.xyz

This file contains the final position and velocity of an MD calculation.

.. code-block:: bash

   <number of atoms>
   Step:    <last MD step>
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>
   ...
   <symbol> <X> <Y> <Z> <V_X> <V_Y> <V_Z>

- BOPOP

This file shows the adiabatic populations (diagonal elements of the density matrix) at each time step.

.. code-block:: bash

   <MD step> <population of state 1> <population of state 2> ... <population of last state> 
   <MD step> <population of state 1> <population of state 2> ... <population of last state> 
   ... 

- BOCOH 

This file shows off-diagonal elements of the density matrix at each time step. Only the upper triangular portions are given because of hermiticity. The real and imaginary part of each element are written alternately.

.. code-block:: bash

   <MD step> <Re. element 1, 2> <Im. element 1, 2> <Re. element 1, 3> <Im. element 1, 3> ... <Re. element last-1, last> <Im. element last-1, last> 
   <MD step> <Re. element 1, 2> <Im. element 1, 2> <Re. element 1, 3> <Im. element 1, 3> ... <Re. element last-1, last> <Im. element last-1, last> 
   ... 

- NACME

This file shows nonadiabatic coupling matrix elements at each time step. Only the upper triangular portions are given because of antihermiticity.

.. code-block:: bash

   <MD step> <element 1, 2> <element 1, 3> ... <element last-1, last> 
   <MD step> <element 1, 2> <element 1, 3> ... <element last-1, last> 
   ... 

- SHPROB

This file shows hopping probabilities from the running state of each time step.

.. code-block:: bash

   <MD step> <P(running -> 1)> <P(running -> 2)> ... <P(running -> last)>
   <MD step> <P(running -> 1)> <P(running -> 2)> ... <P(running -> last)>
   ... 

- SHSTATE

This file shows the running state at each time step.

.. code-block:: bash

   <MD step> <running>
   <MD step> <running>
   ... 

- DOTPOPD

This file shows time-derivative populations by decoherence at each time step.

.. code-block:: bash

   <MD step> <td population of state 1> <td population of state 2> ... <td population of last state> 
   <MD step> <td population of state 1> <td population of state 2> ... <td population of last state> 
   ... 

For a quick test of PyUNIxMD, see :ref:`Quick Start<Quick Start>` . Also, you can refer to scripts and log files in 'examples/' directory for practical calculations.

