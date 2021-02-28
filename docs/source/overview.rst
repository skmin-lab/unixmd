===========================
PyUNIxMD Overview
===========================

Features
---------------------------
| This is features.
| UNI-xMD program includes adiabatic, Ehrenfest, surface hopping and especially DISH-XF dynamics method.
| This package includes several scripts that interface the dynamics modules with some quantum chemistry calculation program such as Molpro, DFTB+, etc..

Authors
---------------------------
This is authors.


Acknowledgement
---------------------------
This is acknowledgement.


Code Overview
---------------------------

Program Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

The overall code structure is displayed in the next figure.

.. image:: ./unixmd_structure.png
   :width: 400pt

PyUNIxMD is mainly based on object-oriented programming, which is structured by
several classes closely connected with each other.
Central modules of UNI-xMD are divided into the following.

- Molecule: a class for describing overall molecule objects, including state objects.

- MQC: a class for dynamics propagation. Contains subclasses for each nonadiabatic dynamics such as Ehrenfest, surface hopping, etc.

- QM: a class for QM calculation of dynamics properties using external softwares such as Molpro, Gaussian 09, DFTB+, etc.

- MM: a class for MM calculation using external software such as Tinker.

- Thermostat: a class for controlling the temperature of a physical system.

Subclasses of MQC, QM and MM classes are organized in inheritance structure.
This helps to simplify the codes by inheriting common arguments to the subclasses.

For the detailed information of each module, check each section in below.

Working Process
^^^^^^^^^^^^^^^^^^^^^^^^^^

To run PyUNIxMD, it requires creating some objects in your running script.
A straightfoward way to perform a UNI-xMD calculation is as follows:

**1. First of all, you should add the corresponding directory to your python path.**

.. code-block:: bash

   $ export PYTHONPATH=$UNIXMDHOME/src:$PYTHONPATH
 
**2. You should import some libraries in the running script.**

.. code-block:: python

   from molecule import Molecule
   import qm, mqc
   from thermostat import *
   from misc import data

**3. To run PyUNIxMD, you should create several objects in your running script. The important
thing is that the object inherited from** ``Molecule`` **class must be created in the first place.**

- Define molecular infomation.

.. code-block:: python

   geom = """
   NUMBER_OF_ATOMS
   TITLE
   SYMBOL  COORDINATES  VELOCITIES
   """

   mol = Molecule(geometry=geom, ARGUMENTS)

.. note:: ``mol`` object must be created before creating another objects which describe QM, MQC and thermostat.

- Determine an electronic structure calculation program and a method to get the energy, force and the nonadiabatic coupling vector.

.. code-block:: python

   qm = qm.QM_prog.QM_method(molecule=mol, ARGUMENTS)

**QM_prog** and **QM_method** stand for an electronic structure calculation program and a theory, respectively. They are listed in ???.

- Determine a method for dynamics propagation.

.. code-block:: python

   md = mqc.MDTYPE(molecule=mol, ARGUMENTS)

**MDTYPE** can be replaced by BOMD, SH, Eh or SHXF which mean Born-Opphenhimer molecular dynamics, surface hopping,
Ehrenfest dynamics and decoherence induced surface hopping based on exact factorization, respectively.

- Choose a thermostat type. Currently, there are three types for a thermostat.

.. code-block:: python

   bathT = THERMOSTAT(temperature=300.0, ARGUMENTS)

**THERMOSTAT** is listed in ???.

- Put your objects into ``run`` method of ``md`` object.

.. code-block:: python

   md.run(molecule=mol, theory=qm, thermostat=bathT, ARGUMENTS)

**4. Execute your running script**

.. code-block:: bash

   $ python3 running_script.py


