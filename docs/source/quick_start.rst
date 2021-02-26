==========================
Quick Start
==========================

Here, we provide executable running script file, which contains:

.. code-block:: python

   from molecule import Molecule
   import qm, mqc
   from thermostat import *
   from misc import data

   geom = """
   NUMBER_OF_ATOMS
   TITLE
   SYMBOL  COORDINATES  VELOCITIES
   """

   mol = Molecule(geometry=geom, nstates=NSTATES)

   qm = qm.QM_PROG.QM_METHOD(ARGUMENTS)

   md = mqc.MDTYPE(ARGUMETNS)

   bathT = THERMOSTAT(ARGUMENTS)

   md.run(molecule=mol, theory=qm, thermostat=bathT, input_dir=INPUT_DIR)

If you execute this script, you can get output files listed in table:

+-----------+------+--------+----+
|           | BOMD | SH(XF) | Eh |
+===========+======+========+====+
| MDENERGY  | o    | o      | o  |
+-----------+------+--------+----+
| MOVIE.xyz | o    | o      | o  |
+-----------+------+--------+----+
| FINAL.xyz | o    | o      | o  |
+-----------+------+--------+----+
| BOCOH *   | x    | o      | o  |
+-----------+------+--------+----+
| BOPOP *   | x    | o      | o  |
+-----------+------+--------+----+
| NACME     | x    | o      | o  |
+-----------+------+--------+----+
| SHPROB    | x    | o      | x  |
+-----------+------+--------+----+
| SHSTATE   | x    | o      | x  |
+-----------+------+--------+----+

.. note:: \* If you set propagation="density", UNI-xMD provides **BOCOH** and **BOPOP**.
   However, if you set propagation="coefficient", UNI-xMD provides **BOCOEF** rather than **BOCOH** and **BOPOP**.

- MDENERGY : energy which contains kinetic energy, potential energy of each adiabatic state and total energy

.. code-block:: bash

   Here is data

- MOVIE.xyz : geometries at each step

.. code-block:: bash

   Here is data

- FINAL.xyz : geometry at final step

.. code-block:: bash

   Here is data

- BOCOH : off-diagonal term of adiabatic density matrix

.. code-block:: bash

   Here is data

- BOPOP : adiabatic population

.. code-block:: bash

   Here is data

- NACME : nonadiabatic coupling matrix element

.. code-block:: bash

   Here is data

- SHPROB : hopping probability between the adiabatic states

.. code-block:: bash

   Here is data

- SHSTATE : running state

.. code-block:: bash

   Here is data

