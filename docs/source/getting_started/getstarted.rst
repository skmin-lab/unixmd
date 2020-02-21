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
Code overview
================================

Program structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

The overall code structure is displayed in next figure.

.. image:: ./unixmd_structure.png
   :width: 400pt

UNI-xMD is mainly based on object-oriented programming, which is structured by 
several classes closely connected with each others. 
Central modules of UNI-xMD can be divided into followings.

Molecule: Describes overall molecule objects, including state objects.

MQC: Class for dynamics propagation. Contains subclasses for each methods
(ex)Ehrenfest, FSSH, etc...)

BO: Class for calculating dynamics properties using external software. 
Includes subclasses for each softwares (ex) g09, molpro, etc...)

Thermostat: Class for velocity rescaling in given condition.

Subclasses of MQC and BO classes are organized in inheritance structure.
This helps to simplify codes by inherit common arguments to different subclasses.

For the detailed information of each modules, check each sections in below.

Working process
^^^^^^^^^^^^^^^^^^^^^^^^^^

The general work flow for whole methods are controlled by run.py. To be specific, 
the user must state basic arguments for run method(md.run) in run.py. i.e, arguments 
for molecule, bo, mqc and thermostat objects are require to run UNI-xMD. For detailed 
information about each objects will be included in each section. 

UNI-xMD also provides template code for run.py. In addition, there are sample codes 
for each BO and MQC methods in example directory.

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

