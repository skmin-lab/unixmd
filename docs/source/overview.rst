===========================
PyUNIxMD Overview
===========================

Features
---------------------------
The features of PyUNIxMD are as follows.

- Conventional (non)adiabatic dynamics

  - Born-Oppenheimer molecular dynamics (BOMD)

  -  Ehrenfest dynamics

  -  Fewest switches surface hopping (FSSH) dynamics with ad hoc decoherence corrections

- Decoherence based on exact factorization

  -  Decoherence induced surface hopping based on exact factorization (DISH-XF)
  -  Coupled-trajectory mixed quantum/classical (CTMQC) method

- Numerical calculation of nonadiabatic couplings
- Accessible interface to external QM programs and built-in model Hamiltonians

  -  COLUMBUS: SA-CASSCF
  -  MOLPRO: SA-CASSCF
  -  Gaussian 09: TDDFT
  -  Q-Chem: TDDFT
  -  TURBOMOLE: TDDFT
  -  TeraChem: SI-SA-REKS (SSR)
  -  DFTB+: TDDFTB, LC-DFTB/SSR
  -  Model Hamiltonians: Tully, Shin-Metiu

- QM/MM functionalities
- Utility scripts in Python

Authors
---------------------------
The current version of PyUNIxMD has been developed by Seung Kyu Min, In Seong Lee, Jong-Kwon Ha, Daeho Han, Kicheol Kim, Tae In Kim, Sung Wook Moon in the Theorectical/Computational Chemistry Group for Excited State Phenomena of Ulsan National Institute of Science and Technology (UNIST). 

..
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

PyUNIxMD is an object-oriented program structured by
several modules closely connected with each other:

- Molecule: a class for molecules. A molecule object contains information of the electronic states as well as the geometry.

- MQC: a class for dynamics propagation for each nonadiabatic dynamics method such as Ehrenfest, surface hopping, etc.

- QM: a class for QM calculations of dynamics properties using external softwares such as Molpro, Gaussian 09, DFTB+, etc.

- MM: a class for MM calculations using external software such as Tinker.

- Thermostat: a class for controlling the temperature of a physical system.

PyUNIxMD takes advantage of the inheritance feature to simplify the codes by sharing the common parameters and methods.

For the detailed information of each module, check each section in below.

Working Process
^^^^^^^^^^^^^^^^^^^^^^^^^^

A straightfoward instruction to perform a PyUNIxMD calculation is as follows:

**1. First of all, you should add the corresponding directory to your python path.**

.. code-block:: bash

   $ export PYTHONPATH=$UNIXMDHOME/src:$PYTHONPATH
 
**2. You should import the PyUNIxMD packages in your running script.**

.. code-block:: python

   from molecule import Molecule
   import qm, mqc
   from thermostat import *
   from misc import data

**3. To run PyUNIxMD, you should create several PyUNIxMD objects in your running script. The important
thing is that you need to make a ** ``Molecule`` ** object first.**

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

**QM_prog** and **QM_method** stands for an electronic structure calculation program and a theory, respectively. For the list of the options, see the QM_calculator section of Chapter 4.

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


