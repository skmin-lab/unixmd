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

To run UNI-xMD, UNI-xMD requires creating some objects in running script.
A straightfoward way to perform a UNI-xMD calculation is as follows:

**1. First of all, you should add the corresponding directory to your python path.**

.. code-block:: bash
  
    $ export PYTHONPATH=$UNIXMDHOME/src:$PYTHONPATH
 
**2. You should import some library in running script.**

.. code-block:: python

    from molecule import Molecule
    import bo
    import mqc
    from thermostat import *
    from misc import data

**3. To run UNI-xMD, you should create several objects, which are** ``molecule`` **,** ``bo`` **,** ``md`` **and** ``thermostat`` **, in your running script. The important thing is that molecule object is created in the first place.**

- Define molecular infomation

.. code-block:: python

    geom = """
    NUMBER_OF_ATOMS
    TITLE
    SYMBOL  COORDINATES  VELOCITIES
    """
    
    mol = Molecule(geometry=geom, ARGUMENTS)

.. note:: molecule object should be already created before creating another objects such as ``bo``, ``md`` and ``thermostat``.

- Determine electronic structure calculation program and method to get energy, force and non-adiabatic coupling vector

.. code-block:: python
   
    qm = bo.QM_prog.QM_method(molecule=mol, ARGUMENTS)

QM_prog and QM_method is electronic structure calculation program and theory, respectively. They are listed in ???.

- Determine method of MD

.. code-block:: python
   
    md = mqc.XXX(molecule=mol, ARGUMENTS)

XXX can be replaced by BOMD, SH, Eh, SHXF which means Born-Opphenhimer molecular dynamics, surface hopping, Ehrenfest dynamics and decoherence induced surface hopping based on exact factorization, respectively.

- Choose thermostat type, there are three types of thermostat, that () states in detail

.. code-block:: python
   
    bathT = thermostat_type(temperature=300.0, ARGUMENTS)

thermostat_type is listed in ???.

- Put your objects into md method

.. code-block:: python
   
    md.run(molecule=mol, theory=qm, thermostat=bathT, ARGUMENTS)

**4. Execute your running script**

.. code-block:: bash

   $ python3 running_script_name

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

