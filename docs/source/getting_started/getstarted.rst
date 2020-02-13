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
Program Structure and Code Flow
================================

In this section, we want to describe the structure of UNI-xMD and how its code flows. 


Program structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

UNI-xMD is mainly based on object-oriented programming, which is structured by 
several classes closely connected with each others. 
Central modules of UNI-xMD can be devided into followings.

Molecule : Describes overall molecule objects, including state objects

MQC : Class for dynamics propagation. Contains subclasses for each methods 
(ex)Ehrenfest, FSSH, etc...)

BO : Class for calculating dynamics properties using external software. 
Includes subclasses for each softwares (ex) g09, molpro, etc...)

FileIO : Class for handling output files.

Subclasses of MQC and BO classes are orgazined in inheritance structure.
This helps to simplify codes by inherit common arguments to different subclasses.

For the detailed information of each modules, check each sections in below.

Work flow
^^^^^^^^^^^^^^^^^^^^^^^^^^

The overall code structure and flowchart are displayed in next figure.

.. image:: ./unixmd_structure.png
    :width: 400pt

The general code flow for whole methods are as follows.

1. Basic molecular informations in input file are stored in molecule class.

2. Initialize the dynamics steps  

3. Calculate dynamics properties from molecular properties using external software in BO module.

4. Propagation through time using calculated features in MQC module.

5. Repeat former steps until dynamics end step

In here, only general scheme is described. About detailed dynamic process for 
individual BO or MQC method will be stated in keyword section.

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

